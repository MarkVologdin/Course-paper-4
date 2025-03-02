#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NUM_EXPERIMENTS 4
#define DATA_POINTS 1000 // Количество точек данных в одном эксперименте
#define TIME_START 0.0
#define TIME_STOP 1000.0
#define M_PI 3.14159265358979323846

#define RadianToDegree  (180./M_PI)                         // Radians  --> Degree
#define D2R(Degrees)    ((Degrees)/RadianToDegree)          // Degree  --> Radians
#define RAD_TO_ARCSEC(RAD)  (RAD*206265)                    // Radians --> Arcseconds

// Блок определения места проведения эксперимента 

#define PHI 55 // широта москвы в градусах 
#define LAMBDA 37 // долгота москвы в градусах 
#define HEIGHT 30 // высота москвы в метрах

// Блок определения составляющих вектора возмущений сила тяжести \delta(g) = {g_E, g_N, g_UP}

#define g_E  4e-4 // не очень понимаю какие начальные значения положить
#define g_N  5e-4 // делаем из расчета 10 угловых секунд = 4.848*10^(-5) рад , а \delta(g) = g_0*(отклонение в радианах)
#define g_UP 0.0
#define DELTA_KSI 90

//Блок определения констант для белого шума 
#define MU 0
#define SIGMA 0.0001

// Определяем константы для формулы Гельмерта (GRS80/WGS84)
#define G_0 9.78030  // Ускорение свободного падения на экваторе (м/с^2)
#define K 0.005302 // Безразмерный параметр
#define E_2 0.000007 // Квадрат эксцентриситета Земли
#define A 6378137.0 // Большая полуось Земли (м)

double gravity_helmert(double phi, double height) {
    /*
     * Вычисляет нормальную силу тяжести по формуле Гельмерта.
     *
     * :param phi: Широта в градусах
     * :param height: Высота над уровнем моря в метрах
     * :return: Нормальная сила тяжести (м/с^2)
     */
    
    double lat_rad = D2R(phi); // Перевод широты в радианы

    // Гравитация на уровне геоида (уровня моря) по формуле Гельмерта
    double g_phi = G_0 * (1 + K * pow(sin(lat_rad), 2)) - E_2 * pow(sin(2 * lat_rad), 2);

    // Коррекция на высоту h (аппроксимация)
    double g_h = g_phi - (3.086 * 1e-6) * height;

    return g_h;
}

double sigmoid(double t, double start_heading, double stop_heading, double t_min, double t_max) {
    // Центрируем сигмоиду в середине диапазона времени
    double t0 = (t_min + t_max) / 2.0;
    double k = 0.03;  // Коэффициент крутизны (чем больше, тем резче переход)

    // Вычисляем значение сигмоиды
    double sigmoid_value = 1.0 / (1.0 + exp(-k * (t - t0)));

    // Масштабируем значение в диапазон [start_heading, stop_heading]
    return start_heading + (stop_heading - start_heading) * sigmoid_value;
}

// Структура для хранения одного считывания данных
typedef struct {
    double timestamp;  // Временная метка
    double phi;   // Широта (λ)
    double lambda;  // Долгота (φ)
    double height;   // Высота (h)
    double heading;    // Курс (ψ)
    double gamma;       // Крен (θ)
    double theta;      // Тангаж (γ)
    double f_z1;       // Сила f_z1
    double f_z2;       // Сила f_z2
    double f_z3;       // Сила f_z3
} DataPoint;

void reproject_vector(double heading, double gamma, double theta, double F_x[3], double F_z[3])
{

    double A_zx[3][3] = {
        {cos(heading)*cos(gamma)+sin(heading)*sin(theta)*sin(gamma), -sin(heading)*cos(gamma)+cos(heading)*sin(theta)*sin(gamma), -cos(theta)*sin(gamma)},
        {sin(heading)*cos(theta), cos(heading)*cos(theta), sin(theta)},
        {cos(heading)*sin(gamma)-sin(heading)*sin(theta)*cos(gamma), -sin(heading)*sin(gamma)-cos(heading)*sin(theta)*cos(gamma), cos(theta)*cos(gamma)},
    };

    //Применяем A к вектору F_x → F_z
    for (int i = 0; i < 3; i++) 
        {
        F_z[i] = 0;
        for (int j = 0; j < 3; j++)
            {
                F_z[i] += A_zx[i][j] * F_x[j];
            }
        }
    // for (int i = 0; i < 3; i++)
    // {
    //     printf("F_z[%d] = %lf\n", i, F_z[i]);
    // }
    

};


// Функция для инициализации данных эксперимента (для генерации сигмоидой добавить - double start_heading_value)
void initialize_data(DataPoint *data, int size, double stop_heading_value) {
    double g_0 = gravity_helmert(PHI, HEIGHT);

    double F_x[3] = {g_E, g_N, g_UP- g_0};
    double F_z[3] = {0, 0, 0};

    for (int i = 0; i < size; i++) {
        if(i == 0){
            data[i].timestamp = 0.0;
        }
        else{
            data[i].timestamp = data[i-1].timestamp + 0.01 ;
        };
        data[i].phi = PHI;
        data[i].lambda = LAMBDA;
        data[i].height = HEIGHT; 
        //проблема в масштабировании сигмоиды я засовывал по 0.01 в timestamp и получался очень маленький шаг, умножим на 100 мы получили шаг по секунде при характерном размере от 0 до 1000 
        // data[i].heading = D2R(sigmoid(data[i].timestamp * 100, start_heading_value, stop_heading_value, TIME_START, TIME_STOP));
        data[i].heading = D2R(stop_heading_value);
        //printf("heading = %lf\n", data[i].heading);
        data[i].gamma = D2R(0);
        data[i].theta = D2R(0);

        //происходит перепроектировка через A_zx(\psi,\theta,\gamma)
        reproject_vector(data[i].heading, data[i].gamma, data[i].theta, F_x, F_z);
        data[i].f_z1 = F_z[0];
        data[i].f_z2 = F_z[1];
        data[i].f_z3 = F_z[2];
        
    }
}

// Функция для генерации белого шума (нормального распределения => mu=0, sigma = 1)
double gaussian_noise(double mu, double sigma) {
    double U1 = (double)rand() / RAND_MAX;
    double U2 = (double)rand() / RAND_MAX;
    
    double Z = sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2); // Генерируем нормальное число
    return mu + sigma * Z;
}


void generate_white_noise(DataPoint *data, int size) {
    for (int i = 0; i < size; i++) {
        data[i].f_z1 = data[i].f_z1 + gaussian_noise(MU, SIGMA);
        data[i].f_z2 = data[i].f_z2 + gaussian_noise(MU, SIGMA);
        data[i].f_z3 = data[i].f_z3 + gaussian_noise(MU, SIGMA);
    }
}

// Функция для записи данных в файл
void save_data_to_file(const char *filename, DataPoint *data, int size) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Ошибка при открытии файла");
        return;
    }
    fprintf(file, "Time     Phi     Lambda      Height      Heading     Gamma   Theta   F_z1    F_z2    F_z3\n");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                data[i].timestamp, data[i].phi, data[i].lambda, data[i].height,
                data[i].heading, data[i].gamma, data[i].theta, 
                data[i].f_z1, data[i].f_z2, data[i].f_z3);
    }
    
    fclose(file);
}
// Функция для чтения данных из файла
int read_data(const char* filename, DataPoint* data) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Ошибка при открытии файла");
        return -1;
    }

    // Пропускаем заголовок
    char header[256];
    fgets(header, sizeof(header), file);

    int count = 0;
    while (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                  &data[count].timestamp, &data[count].phi, &data[count].lambda,
                  &data[count].height, &data[count].heading, &data[count].gamma,
                  &data[count].theta, &data[count].f_z1, &data[count].f_z2, &data[count].f_z3) == 10) {
        count++;
        if (count >= 10000) break;  // Ограничение на количество записей
    }

    fclose(file);
    return count;
}

//  Функция для осреднения последних трёх значений сил
void average_last_three_forces(DataPoint* data, int count, double* avg_fz1, double* avg_fz2, double* avg_fz3) {
    if (count < 3) {
        printf("Недостаточно данных для осреднения.\n");
        *avg_fz1 = *avg_fz2 = *avg_fz3 = 0.0;
        return;
    }

    *avg_fz1 = *avg_fz2 = *avg_fz3 = 0.0;

    for (int i = 0; i < count; i++) {
        *avg_fz1 += data[i].f_z1;
        *avg_fz2 += data[i].f_z2;
        *avg_fz3 += data[i].f_z3;
    }

    *avg_fz1 /= count;
    *avg_fz2 /= count;
    *avg_fz3 /= count;
}
int invert_2x2_matrix(double matrix[2][2], double inverse[2][2]) {
    // Вычисляем определитель
    double det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    if (det == 0) {
        printf("Ошибка: матрица необратима (det = 0).\n");
        return 0;  // Необратимая матрица
    }

    // Вычисляем обратную матрицу
    inverse[0][0] =  matrix[1][1] / det;
    inverse[0][1] = -matrix[0][1] / det;
    inverse[1][0] = -matrix[1][0] / det;
    inverse[1][1] =  matrix[0][0] / det;

    return 1;  // Успех
}
int invert_3x3_matrix(double matrix[3][3], double inverse[3][3]) {
    // Вычисляем определитель
    double det = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
                 matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
                 matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

    if (det == 0) {
        printf("Ошибка: матрица необратима (det = 0).\n");
        return 0;  // Необратимая матрица
    }

    // Вычисляем матрицу алгебраических дополнений и транспонируем (аджугатную матрицу)
    inverse[0][0] =  (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det;
    inverse[0][1] = -(matrix[0][1] * matrix[2][2] - matrix[0][2] * matrix[2][1]) / det;
    inverse[0][2] =  (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det;

    inverse[1][0] = -(matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) / det;
    inverse[1][1] =  (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det;
    inverse[1][2] = -(matrix[0][0] * matrix[1][2] - matrix[0][2] * matrix[1][0]) / det;

    inverse[2][0] =  (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det;
    inverse[2][1] = -(matrix[0][0] * matrix[2][1] - matrix[0][1] * matrix[2][0]) / det;
    inverse[2][2] =  (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det;

    return 1;  // Успех
}

//  Основная функция для обработки данных
void process_data(const char* file1, const char* file2, const char* file3) {
    DataPoint data1[DATA_POINTS], data2[DATA_POINTS];
    FILE *file = fopen(file3, "w");


    if (!file) {
        perror("Ошибка при открытии файла");
        return;
    }

    // Чтение данных из первого файла
    int count1 = read_data(file1, data1);
    if (count1 < 3) return;

    // Чтение данных из второго файла
    int count2 = read_data(file2, data2);
    if (count2 < 3) return;

    // Осреднённые значения сил
    double avg1_fz1, avg1_fz2, avg1_fz3;
    double avg2_fz1, avg2_fz2, avg2_fz3;

    // Осреднение для каждого файла
    average_last_three_forces(data1, count1, &avg1_fz1, &avg1_fz2, &avg1_fz3);
    average_last_three_forces(data2, count2, &avg2_fz1, &avg2_fz2, &avg2_fz3);

    //  Вычисление разниц между соответствующими силами (второй - первый)
    double diff_fz1 = avg2_fz1 - avg1_fz1;
    double diff_fz2 = avg2_fz2 - avg1_fz2;
    double diff_fz3 = avg2_fz3 - avg1_fz3;
    // получили первую дельту

    //  Вывод результатов
    printf("---------------------------------------------\n");
    printf("Осреднённые значения сил для %s: F_z1=%.3f, F_z2=%.3f, F_z3=%.3f\n", file1, avg1_fz1, avg1_fz2, avg1_fz3);
    printf("Осреднённые значения сил для %s: F_z1=%.3f, F_z2=%.3f, F_z3=%.3f\n", file2, avg2_fz1, avg2_fz2, avg2_fz3);
    printf("Разница между ними: F_z1=%.5f, F_z2=%.5f, F_z3=%.5f\n", diff_fz1, diff_fz2, diff_fz3);

    //  Матрица A(\ksi)
    double matrix_A_s_si[2][2] = { {cos(D2R(DELTA_KSI))-1, -sin(D2R(DELTA_KSI)) },
                            {sin(D2R(DELTA_KSI)), cos(D2R(DELTA_KSI))-1}};

                
    //  Матрица A^(-1)(\ksi)
    double inverse_matrix_A_s_si[2][2];

    invert_2x2_matrix(matrix_A_s_si, inverse_matrix_A_s_si);

    // Выведем матрицы в файл 
    fprintf(file, "Матрица A_s_si:\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            fprintf(file, "%.5f ", matrix_A_s_si[i][j]);
        }
        fprintf(file, "\n");
    }

    fprintf(file, "------------------------\n");
    fprintf(file, "Матрица A_s_si^(-1):\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            fprintf(file, "%.5f ", inverse_matrix_A_s_si[i][j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "------------------------\n");

    //  Вектор g_x_0
    double F_x_0[3] = {0, 0, G_0};
    double F_s_0[3] = {0, 0, 0};
    double F_z_0[3] = {0, 0, 0};

    //  Перепроектировка g_x_0 в g_s_0 помощью A_s_x(\psi,\gamma,\theta)
    reproject_vector(data1[0].heading, 0, 0, F_x_0, F_s_0);
    F_z_0[0] = F_s_0[0]*matrix_A_s_si[0][0] + F_s_0[1]*matrix_A_s_si[0][1];
    F_z_0[1] = F_s_0[0]*matrix_A_s_si[1][0] + F_s_0[1]*matrix_A_s_si[1][1];

    // получили \delta(\delta F_z)
    double Diff_f_z[2] = {0,0};
    Diff_f_z[0] = diff_fz1 - F_z_0[0];
    Diff_f_z[1] = diff_fz2 - F_z_0[1];

    // умножаем на обратную матрицу A^-1(\ksi)
    // invert_2x2_matrix(matrix_A_s_si, inverse_matrix_A_s_si);
    double g_E_1 = inverse_matrix_A_s_si[0][0] * Diff_f_z[0] + inverse_matrix_A_s_si[0][1] * Diff_f_z[1];
    double g_N_2 = inverse_matrix_A_s_si[1][0] * Diff_f_z[0] + inverse_matrix_A_s_si[1][1] * Diff_f_z[1];

    // нужна матрица A_s_x^0(\ksi, \gamma, \theta)
    double matrix_A_s_x[3][3] = { {cos(data1[0].heading), -sin(data1[0].heading), 0 },
                            {sin(data1[0].heading), cos(data1[0].heading), 0},
                            {0, 0, 1}};

    double inverse_matrix_A_s_x[3][3]; // обращаем матрицу A_s_x^0
    invert_3x3_matrix(matrix_A_s_x, inverse_matrix_A_s_x);

    // Выведем матрицы в файл A_s_x
    fprintf(file, "Матрица A_s_x:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            fprintf(file, "%.5f ", matrix_A_s_x[i][j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "------------------------\n");
    fprintf(file, "Матрица A_s_x^(-1):\n");

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            fprintf(file, "%.5f ", inverse_matrix_A_s_x[i][j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "------------------------\n");


    double g_E_exp = inverse_matrix_A_s_x[0][0] * g_E_1 + inverse_matrix_A_s_x[0][1] * g_N_2;
    double g_N_exp = inverse_matrix_A_s_x[1][0] * g_E_1 + inverse_matrix_A_s_x[1][1] * g_N_2;

    printf("Экспериментальные значения сил: g_E_exp=%.5f, g_N_exp=%.5f g_E=%.5f g_N=%.5f\n", g_E_exp, g_N_exp, g_E, g_N);
    printf("В угловых секундах: g_E_exp=%.3f, g_E=%.3f g_N_exp=%.3f g_N=%.3f\n", RAD_TO_ARCSEC(g_E_exp/G_0), RAD_TO_ARCSEC(g_E/G_0) , RAD_TO_ARCSEC(g_N_exp/G_0), RAD_TO_ARCSEC(g_N/G_0));
    printf("Ошибка в угловых секундах: D_g_E=%.3f, D_g_N=%.3f\n", RAD_TO_ARCSEC(g_E_exp/G_0) - RAD_TO_ARCSEC(g_E/G_0) , RAD_TO_ARCSEC(g_N_exp/G_0) - RAD_TO_ARCSEC(g_N/G_0));
}

int main() {
    srand(time(NULL)); // Инициализация генератора случайных чисел
    
    char filename[20];

    double F_x[3] = {0, 0, -10};
    double F_z[3] = {0, 0, 0};

    reproject_vector(0, 0, 0, F_x, F_z);

    // Динамическое выделение памяти для данных
    DataPoint* data = (DataPoint*)malloc(DATA_POINTS * sizeof(DataPoint));
    if (!data) {
        printf("Ошибка выделения памяти для данных.\n");
        return -1;
    }

    for (int exp = 0; exp < NUM_EXPERIMENTS; exp++) {
        initialize_data(data, DATA_POINTS, (exp + 1) * DELTA_KSI);
        generate_white_noise(data, DATA_POINTS);
        
        snprintf(filename, sizeof(filename), "experiment_%d.txt", exp + 1);
        save_data_to_file(filename, data, DATA_POINTS);
        
        printf("Данные эксперимента %d записаны в %s\n", exp + 1, filename);
    }

    printf("Данные для эксперимента сгенерированы с Mu = %d, Sigma = %lf\n", MU, SIGMA);
    printf("-------------------------------------------------------------\n");

    process_data("experiment_1.txt", "experiment_2.txt", "experiment_1_set.txt");
    process_data("experiment_2.txt", "experiment_3.txt", "experiment_2_set.txt");
    process_data("experiment_3.txt", "experiment_4.txt", "experiment_3_set.txt");
    // process_data("experiment_4.txt", "experiment_5.txt", "experiment_4_set.txt");
    // process_data("experiment_5.txt", "experiment_6.txt", "experiment_5_set.txt");
    // process_data("experiment_6.txt", "experiment_7.txt", "experiment_6_set.txt");

    // Освобождение динамической памяти
    free(data);
    
    return 0;
}