#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_EXPERIMENTS 7
#define DATA_POINTS 100 // Количество точек данных в одном эксперименте
#define TIME_START 0
#define TIME_STOP 100
#define M_PI 3.14159265358979323846

#define RadianToDegree  (180./M_PI)                         // Radians  --> Degree
#define D2R(Degrees)    ((Degrees)/RadianToDegree)          // Degree  --> Radians

// Блок определения места проведения эксперимента 

#define PHI 55 // широта москвы в градусах 
#define LAMBDA 37 // долгота москвы в градусах 
#define HEIGHT 30 // высота москвы в метрах

// Блок определения составляющих вектора возмущений сила тяжести \delta(g) = {g_E, g_N, g_UP}

#define g_E  4e-4 // не очень понимаю какие начальные значения положить
#define g_N  5e-4 // делаем из расчета 10 угловых секунд = 4.848*10^(-5) рад , а \delta(g) = g_0*(отклонение в радианах)
#define g_UP 0.0

//Блок определения констант для белого шума 
#define MU 0
#define SIGMA 1

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


// Функция для инициализации данных эксперимента
void initialize_data(DataPoint *data, int size, double start_heading_value, double stop_heading_value) {
    double g_0 = gravity_helmert(PHI, HEIGHT);

    double f_x1 = g_E;
    double f_x2 = g_N;
    double f_x3 = g_UP - g_0; 

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
        data[i].heading = D2R(sigmoid(data[i].timestamp, start_heading_value, stop_heading_value, TIME_START, TIME_STOP));
        data[i].gamma = D2R(0);
        data[i].theta = D2R(0);

        //происходит перепроектировка через A_zx(\psi,\theta,\gamma)
        data[i].f_z1 = f_x1*(sin(data[i].heading)*cos(data[i].theta)) + f_x2*(cos(data[i].heading)*cos(data[i].theta)) + f_x3*sin(data[i].theta);
        data[i].f_z2 = f_x1*(cos(data[i].heading)*sin(data[i].gamma)- sin(data[i].heading)*sin(data[i].theta)*cos(data[i].gamma)) + f_x2*(-sin(data[i].heading)*sin(data[i].gamma)- cos(data[i].heading)*sin(data[i].theta)*cos(data[i].gamma)) + f_x3*cos(data[i].theta)*cos(data[i].gamma);
        data[i].f_z3 = f_x1*(cos(data[i].heading)*cos(data[i].gamma)+ sin(data[i].heading)*sin(data[i].theta)*sin(data[i].gamma)) + f_x2*(-sin(data[i].heading)*cos(data[i].gamma)+ cos(data[i].heading)*sin(data[i].theta)*sin(data[i].gamma)) - f_x3*cos(data[i].theta)*sin(data[i].gamma);
        
    }
}

double sigmoid(double t, double start_heading, double stop_heading, double t_min, double t_max) {
    // Центрируем сигмоиду в середине диапазона времени
    double t0 = (t_min + t_max) / 2.0;
    double k = 0.1;  // Коэффициент крутизны (чем больше, тем резче переход)

    // Вычисляем значение сигмоиды
    double sigmoid_value = 1.0 / (1.0 + exp(-k * (t - t0)));

    // Масштабируем значение в диапазон [start_heading, stop_heading]
    return start_heading + (stop_heading - start_heading) * sigmoid_value;
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
    fprintf(file, "Time Phi Lambda Height Heading Gamma Theta F_z1 F_z2 F_z3\n");
    for (int i = 0; i < size; i++) {
        fprintf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                data[i].timestamp, data[i].phi, data[i].lambda, data[i].height,
                data[i].heading, data[i].gamma, data[i].theta, 
                data[i].f_z1, data[i].f_z2, data[i].f_z3);
    }
    
    fclose(file);
}

int main() {
    srand(time(NULL)); // Инициализация генератора случайных чисел
    
    char filename[20];
    DataPoint data[DATA_POINTS];
    
    for (int exp = 0; exp < NUM_EXPERIMENTS; exp++) {
        initialize_data(data, DATA_POINTS, exp*45.0, (exp+1)*45);
        generate_white_noise(data, DATA_POINTS);
        
        snprintf(filename, sizeof(filename), "experiment_%d.txt", exp + 1);
        save_data_to_file(filename, data, DATA_POINTS);
        
        printf("Данные эксперимента %d записаны в %s\n", exp + 1, filename);
    }
    
    return 0;
}
