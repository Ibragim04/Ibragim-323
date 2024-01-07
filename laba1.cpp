/*#include <iostream>
#include <fstream>
#include <cmath>

# define PI 3.141592653589793238463
# define epsilon 1e-5   

using namespace std;


// Функция для вычисления эксцентрической аномалии методом итераций
double Iteration(double M, double e) {
    
    double E_start = M, E_next = e * sin(E_start) + M; // начальное приближение E(0) = M
    while ((abs(E_next - E_start) < epsilon)) { // условие
        E_start = E_next; // Обновляем значение E
        E_next = e * sin(E_start) + M; // Вычисляем новое значение E
    }
    return E_next; // найденное приближение E
}

double HalfIteration(double M, double e) {
    
    double f = M - 2.0;
    double u = M + 1.0;
    double c = (f + u) / 2.0;

        do { //точность вычислений
                if ((f - e * sin(f) - M) * (c - e * sin(c) - M) < 0) {
                    u = c;
                }
                else {
                    f = c;
                }

            // Обновляем значение точки c
            c = (f + u) / 2.0;

        } while (fabs(u - f) < 2 * epsilon || fabs(c - e * sin(c) - M < epsilon));
        return c; // Возвращаем найденное значение
}

double GoldSection(double M, double e) {
    
    const double phi = (sqrt(5) + 1) / 2; // Золотое сечение

    double f = M - 3.0;
    double u = M + 3.0;
    double c, fc;

    do {
        // Вычисляем новую точку по формуле золотого сечения
        c = (f + (u - f) / phi);

        // Вычисляем значение функции в точке c
        fc = (c - e * sin(c) - M);

        // Проверяем знак изменения функции и обновляем границы интервала
        if ((f - e * sin(f) - M) * fc < 0)
            u = c;
        else
            f = c;
        
    } while (fabs(u - f) < 2 * epsilon || fabs(c - e * sin(c) - M < epsilon));
    return c; // Возвращаем найденное значение
}

double Newton(double M, double e) {

    double E_start = M;
    double E_next = E_start - ((E_start - e * sin(E_start) - M) / (1 - e * cos(E_start)));

    while (fabs(E_next - E_start) < epsilon) {
        E_start = E_next;
        E_next = E_start - ((E_start - e * sin(E_start) - M) / (1 - e * cos(E_start)));
    }
    return E_next;
}


int main() {
    double T = 43200.0;  // Период орбиты в с 
    double ra = 12739;   // Радиус апоцентра в км
    double rp = 2639;    // Радиус перицентра в км
    

    // Вычисляем параметры орбиты
    double a = (ra + rp) / 2;              // Большая полуось орбиты
    double e = (ra - rp) / (2 * a);      // Эксцентриситет орбиты
    double n = (2 * PI) / T;                 // Средняя угловая скорость

    ofstream data1;
    data1.open("Data_Itration.txt");

    if (!data1) {
        cout << "Error opening file.";
        return 0;
    }
 
    // Записываем заголовки столбцов в файл
    data1 << "t, c\t";
    data1 << "M(t), рад\t";
    data1 << "E(t), рад\t";
    data1 << "Theta(t), рад\t" << endl;

    // Вычисляем и записываем значения в файл
    for (int t = 0; t <= T; t++) {
        double M = n * t;                                                               // Средняя аномалия
        double E = Iteration(M, e);                                                     // Эксцентрическая аномалия - метод итераций
        double true_anomaly = atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2;           // Истинная аномалия
        if (true_anomaly < 0) { true_anomaly += 2 * PI; }                               // изменение истинной аномалии на интервал [0, 2π), т.к. не мб <0
        
        data1 << t << "\t" << M << "\t" << E << "\t" << true_anomaly << endl;
    }

    data1.close();
    cout << "The data was successfully written to the file.";

    ofstream data2;
    data2.open("Data_HalfIteration.txt");

    if (!data2) {
        cout << "Error opening file.";
        return 0;
    }

    data2 << "t, c\t";
    data2 << "M(t), рад\t";
    data2 << "E(t), рад\t";
    data2 << "Theta(t), рад\t" << endl;

    // Вычисляем и записываем значения в файл
    for (int t = 0; t <= T; t++) {
        double M = n * t;                                                               // Средняя аномалия
        double E = HalfIteration(M, e);                                                 // Эксцентрическая аномалия - метод деления на 2
        double true_anomaly = atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2;           // Истинная аномалия
        if (true_anomaly < 0) { true_anomaly += 2 * PI; }                               // изменение истинной аномалии на интервал [0, 2π), т.к. не мб <0

        data2 << t << "\t" << M << "\t" << E << "\t" << true_anomaly << endl;
    }

    data2.close();
    cout << "The data was successfully written to the file.";

 /*   ofstream data3;
    data3.open("Data_GoldSection.txt");

    if (!data3) {
        cout << "Error opening file.";
        return 0;
    }

    data3 << "t, c\t";
    data3 << "M(t), рад\t";
    data3 << "E(t), рад\t";
    data3 << "Theta(t), рад\t" << endl;

    // Вычисляем и записываем значения в файл
    for (int t = 0; t <= T; t++) {
        double M = n * t;                                                               // Средняя аномалия
        double E = GoldSection(M, e);                                                   // Эксцентрическая аномалия - метод золотого сечения
        double true_anomaly = atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2;           // Истинная аномалия
        if (true_anomaly < 0) { true_anomaly += 2 * PI; }                               // изменение истинной аномалии на интервал [0, 2π), т.к. не мб <0

        data3 << t << "\t" << M << "\t" << E << "\t" << true_anomaly << endl;
    }

    data3.close();
    cout << "The data was successfully written to the file.";

    ofstream data4;
    data3.open("Data_Newton.txt");

    if (!data4) {
        cout << "Error opening file.";
        return 0;
    }

    data4 << "t, c\t";
    data4 << "M(t), рад\t";
    data4 << "E(t), рад\t";
    data4 << "Theta(t), рад\t" << endl;

    // Вычисляем и записываем значения в файл
    for (int t = 0; t <= T; t++) {
        double M = n * t;                                                               // Средняя аномалия
        double E = Newton(M, e);                                                        // Эксцентрическая аномалия - метод Ньютона
        double true_anomaly = atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)) * 2;           // Истинная аномалия
        if (true_anomaly < 0) { true_anomaly += 2 * PI; }                               // изменение истинной аномалии на интервал [0, 2π), т.к. не мб <0

        data3 << t << "\t" << M << "\t" << E << "\t" << true_anomaly << endl;
    }

    data4.close();
    cout << "The data was successfully written to the file.";

    return 0;
}
*/