#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <SFML/Graphics.hpp>


const double m1 = 1.0, m2 = 1.0;
const double l1 = 1.0, l2 = 1.0;
const double m3 = 1.0, m4 = 1.0;
const double l3 = 1.0, l4 = 1.0;
const double g = 9.81;
const int windowWidth = 800;
const int windowHeight = 600;
const double scale = 100;


struct State {
    double theta1;
    double theta2;
    double omega1;
    double omega2;
};


void derivatives(const State& state, State& dstate) {
    double delta = state.theta2 - state.theta1;

    dstate.theta1 = state.omega1;
    dstate.theta2 = state.omega2;

    dstate.omega1 = m2 * l1 * state.omega1 * state.omega1 * sin(delta) * cos(delta) +
                    m2 * g * sin(state.theta2) * cos(delta) +
                    m2 * l2 * state.omega2 * state.omega2 * sin(delta) -
                    (m1 + m2) * g * sin(state.theta1);
    dstate.omega1 /= ((m1 + m2) * l1 - m2 * l1 * cos(delta) * cos(delta));

    dstate.omega2 = -m2 * l2 * state.omega2 * state.omega2 * sin(delta) * cos(delta) +
                    (m1 + m2) * (g * sin(state.theta1) * cos(delta) -
                    l1 * state.omega1 * state.omega1 * sin(delta) -
                    g * sin(state.theta2));
    dstate.omega2 /= ((m1 + m2) * l2 - m2 * l2 * cos(delta) * cos(delta));
}


void rk4(State& state, double dt) {
    State k1, k2, k3, k4, temp;

    derivatives(state, k1);
    temp = {state.theta1 + 0.5 * dt * k1.theta1,
            state.theta2 + 0.5 * dt * k1.theta2,
            state.omega1 + 0.5 * dt * k1.omega1,
            state.omega2 + 0.5 * dt * k1.omega2};

    derivatives(temp, k2);
    temp = {state.theta1 + 0.5 * dt * k2.theta1,
            state.theta2 + 0.5 * dt * k2.theta2,
            state.omega1 + 0.5 * dt * k2.omega1,
            state.omega2 + 0.5 * dt * k2.omega2};

    derivatives(temp, k3);
    temp = {state.theta1 + dt * k3.theta1,
            state.theta2 + dt * k3.theta2,
            state.omega1 + dt * k3.omega1,
            state.omega2 + dt * k3.omega2};

    derivatives(temp, k4);

    state.theta1 += dt * (k1.theta1 + 2 * k2.theta1 + 2 * k3.theta1 + k4.theta1) / 6;
    state.theta2 += dt * (k1.theta2 + 2 * k2.theta2 + 2 * k3.theta2 + k4.theta2) / 6;
    state.omega1 += dt * (k1.omega1 + 2 * k2.omega1 + 2 * k3.omega1 + k4.omega1) / 6;
    state.omega2 += dt * (k1.omega2 + 2 * k2.omega2 + 2 * k3.omega2 + k4.omega2) / 6;
}


void dopri5(State& state, double dt) {
    const double a2 = 1.0 / 5.0, a3 = 3.0 / 10.0, a4 = 4.0 / 5.0, a5 = 8.0 / 9.0, a6 = 1.0, a7 = 1.0;
    const double b21 = 1.0 / 5.0;
    const double b31 = 3.0 / 40.0, b32 = 9.0 / 40.0;
    const double b41 = 44.0 / 45.0, b42 = -56.0 / 15.0, b43 = 32.0 / 9.0;
    const double b51 = 19372.0 / 6561.0, b52 = -25360.0 / 2187.0, b53 = 64448.0 / 6561.0, b54 = -212.0 / 729.0;
    const double b61 = 9017.0 / 3168.0, b62 = -355.0 / 33.0, b63 = 46732.0 / 5247.0, b64 = 49.0 / 176.0, b65 = -5103.0 / 18656.0;
    const double b71 = 35.0 / 384.0, b72 = 0.0, b73 = 500.0 / 1113.0, b74 = 125.0 / 192.0, b75 = -2187.0 / 6784.0, b76 = 11.0 / 84.0;

    State k1, k2, k3, k4, k5, k6, k7, temp;

    derivatives(state, k1);

    temp = {state.theta1 + dt * b21 * k1.theta1,
            state.theta2 + dt * b21 * k1.theta2,
            state.omega1 + dt * b21 * k1.omega1,
            state.omega2 + dt * b21 * k1.omega2};
    derivatives(temp, k2);

    temp = {state.theta1 + dt * (b31 * k1.theta1 + b32 * k2.theta1),
            state.theta2 + dt * (b31 * k1.theta2 + b32 * k2.theta2),
            state.omega1 + dt * (b31 * k1.omega1 + b32 * k2.omega1),
            state.omega2 + dt * (b31 * k1.omega2 + b32 * k2.omega2)};
    derivatives(temp, k3);

    temp = {state.theta1 + dt * (b41 * k1.theta1 + b42 * k2.theta1 + b43 * k3.theta1),
            state.theta2 + dt * (b41 * k1.theta2 + b42 * k2.theta2 + b43 * k3.theta2),
            state.omega1 + dt * (b41 * k1.omega1 + b42 * k2.omega1 + b43 * k3.omega1),
            state.omega2 + dt * (b41 * k1.omega2 + b42 * k2.omega2 + b43 * k3.omega2)};
    derivatives(temp, k4);

    temp = {state.theta1 + dt * (b51 * k1.theta1 + b52 * k2.theta1 + b53 * k3.theta1 + b54 * k4.theta1),
            state.theta2 + dt * (b51 * k1.theta2 + b52 * k2.theta2 + b53 * k3.theta2 + b54 * k4.theta2),
            state.omega1 + dt * (b51 * k1.omega1 + b52 * k2.omega1 + b53 * k3.omega1 + b54 * k4.omega1),
            state.omega2 + dt * (b51 * k1.omega2 + b52 * k2.omega2 + b53 * k3.omega2 + b54 * k4.omega2)};
    derivatives(temp, k5);

    temp = {state.theta1 + dt * (b61 * k1.theta1 + b62 * k2.theta1 + b63 * k3.theta1 + b64 * k4.theta1 + b65 * k5.theta1),
            state.theta2 + dt * (b61 * k1.theta2 + b62 * k2.theta2 + b63 * k3.theta2 + b64 * k4.theta2 + b65 * k5.theta2),
            state.omega1 + dt * (b61 * k1.omega1 + b62 * k2.omega1 + b63 * k3.omega1 + b64 * k4.omega1 + b65 * k5.omega1),
            state.omega2 + dt * (b61 * k1.omega2 + b62 * k2.omega2 + b63 * k3.omega2 + b64 * k4.omega2 + b65 * k5.omega2)};
    derivatives(temp, k6);

    temp = {state.theta1 + dt * (b71 * k1.theta1 + b72 * k2.theta1 + b73 * k3.theta1 + b74 * k4.theta1 + b75 * k5.theta1 + b76 * k6.theta1),
            state.theta2 + dt * (b71 * k1.theta2 + b72 * k2.theta2 + b73 * k3.theta2 + b74 * k4.theta2 + b75 * k5.theta2 + b76 * k6.theta2),
            state.omega1 + dt * (b71 * k1.omega1 + b72 * k2.omega1 + b73 * k3.omega1 + b74 * k4.omega1 + b75 * k5.omega1 + b76 * k6.omega1),
            state.omega2 + dt * (b71 * k1.omega2 + b72 * k2.omega2 + b73 * k3.omega2 + b74 * k4.omega2 + b75 * k5.omega2 + b76 * k6.omega2)};
    derivatives(temp, k7);

    state.theta1 += dt * (b71 * k1.theta1 + b72 * k2.theta1 + b73 * k3.theta1 + b74 * k4.theta1 + b75 * k5.theta1 + b76 * k6.theta1);
    state.theta2 += dt * (b71 * k1.theta2 + b72 * k2.theta2 + b73 * k3.theta2 + b74 * k4.theta2 + b75 * k5.theta2 + b76 * k6.theta2);
    state.omega1 += dt * (b71 * k1.omega1 + b72 * k2.omega1 + b73 * k3.omega1 + b74 * k4.omega1 + b75 * k5.omega1 + b76 * k6.omega1);
    state.omega2 += dt * (b71 * k1.omega2 + b72 * k2.omega2 + b73 * k3.omega2 + b74 * k4.omega2 + b75 * k5.omega2 + b76 * k6.omega2);
}


void writeErrorToFile1(const State& state_rk4, const State& state_dopri5, double time, std::ofstream& file) {
    double error_theta1 = fabs(state_rk4.theta1 - state_dopri5.theta1);
    double error_theta2 = fabs(state_rk4.theta2 - state_dopri5.theta2);
    double error_omega1 = fabs(state_rk4.omega1 - state_dopri5.omega1);
    double error_omega2 = fabs(state_rk4.omega2 - state_dopri5.omega2);

    file << time << " " << error_theta1 << " " << error_theta2 << " " << error_omega1 << " " << error_omega2 << std::endl;
}


void writeErrorToFile2(const State& state_rk4, const State& state_dopri5, double time, std::ofstream& file) {
    double error_theta1 = fabs(state_rk4.theta1 - state_dopri5.theta1) / fabs(state_dopri5.theta1);
    double error_theta2 = fabs(state_rk4.theta2 - state_dopri5.theta2) / fabs(state_dopri5.theta2);
    double error_omega1 = fabs(state_rk4.omega1 - state_dopri5.omega1) / fabs(state_dopri5.omega1);
    double error_omega2 = fabs(state_rk4.omega2 - state_dopri5.omega2) / fabs(state_dopri5.omega2);

    file << time << " " << error_theta1 << " " << error_theta2 << " " << error_omega1 << " " << error_omega2 << std::endl;
}


void runVisualization(const State& initialState1, const State& initialState2, const std::string& method) {
    State state1 = initialState1;
    State state2 = initialState2;
    double dt = 0.01;

    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "Double Pendulum");
    window.setFramerateLimit(60);

    sf::VertexArray pendulum1(sf::LinesStrip, 3);
    pendulum1[0].position = sf::Vector2f(windowWidth / 2, windowHeight / 2);
    pendulum1[0].color = sf::Color::Blue;
    pendulum1[1].color = sf::Color::Blue;
    pendulum1[2].color = sf::Color::Blue;

    sf::VertexArray pendulum2(sf::LinesStrip, 3);
    pendulum2[0].position = sf::Vector2f(windowWidth / 2, windowHeight / 2);
    pendulum2[0].color = sf::Color::Red;
    pendulum2[1].color = sf::Color::Red;
    pendulum2[2].color = sf::Color::Red;

    sf::CircleShape mass1(10.0f);
    mass1.setFillColor(sf::Color::Blue);
    mass1.setOrigin(10.0f, 10.0f);

    sf::CircleShape mass2(10.0f);
    mass2.setFillColor(sf::Color::Blue);
    mass2.setOrigin(10.0f, 10.0f);

    sf::CircleShape mass3(10.0f);
    mass3.setFillColor(sf::Color::Red);
    mass3.setOrigin(10.0f, 10.0f);

    sf::CircleShape mass4(10.0f);
    mass4.setFillColor(sf::Color::Red);
    mass4.setOrigin(10.0f, 10.0f);

    sf::VertexArray trace1(sf::LinesStrip, 0);
    sf::VertexArray trace2(sf::LinesStrip, 0);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        double x1 = l1 * sin(state1.theta1) * scale + windowWidth / 2;
        double y1 = l1 * cos(state1.theta1) * scale + windowHeight / 2;
        double x2 = x1 + l2 * sin(state1.theta2) * scale;
        double y2 = y1 + l2 * cos(state1.theta2) * scale;

        pendulum1[1].position = sf::Vector2f(x1, y1);
        pendulum1[2].position = sf::Vector2f(x2, y2);

        mass1.setPosition(x1, y1);
        mass2.setPosition(x2, y2);

        trace1.append(sf::Vertex(sf::Vector2f(x2, y2), sf::Color(255, 165, 0)));

        double x3 = l3 * sin(state2.theta1) * scale + windowWidth / 2;
        double y3 = l3 * cos(state2.theta1) * scale + windowHeight / 2;
        double x4 = x3 + l4 * sin(state2.theta2) * scale;
        double y4 = y3 + l4 * cos(state2.theta2) * scale;

        pendulum2[1].position = sf::Vector2f(x3, y3);
        pendulum2[2].position = sf::Vector2f(x4, y4);

        mass3.setPosition(x3, y3);
        mass4.setPosition(x4, y4);

        trace2.append(sf::Vertex(sf::Vector2f(x4, y4), sf::Color(0, 255, 0)));

        window.clear(sf::Color::White);
        window.draw(pendulum1);
        window.draw(pendulum2);
        window.draw(mass1);
        window.draw(mass2);
        window.draw(mass3);
        window.draw(mass4);
        window.draw(trace1);
        window.draw(trace2);
        window.display();

        if (method == "rk4") {
            rk4(state1, dt);
            rk4(state2, dt);
        } else if (method == "d5") {
            dopri5(state1, dt);
            dopri5(state2, dt);
        }
    }
}


void calculateTotalError(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Не удалось открыть файл " << filename << " для чтения.\n";
        return;
    }

    double total_error_theta1 = 0.0, total_error_theta2 = 0.0, total_error_omega1 = 0.0, total_error_omega2 = 0.0;
    double time, error_theta1, error_theta2, error_omega1, error_omega2;

    while (file >> time >> error_theta1 >> error_theta2 >> error_omega1 >> error_omega2) {
        total_error_theta1 += error_theta1;
        total_error_theta2 += error_theta2;
        total_error_omega1 += error_omega1;
        total_error_omega2 += error_omega2;
    }

    file.close();

    std::cout << "Суммарная ошибка для файла " << filename << ":\n";
    std::cout << "Ошибка по theta1: " << total_error_theta1 << "\n";
    std::cout << "Ошибка по theta2: " << total_error_theta2 << "\n";
    std::cout << "Ошибка по omega1: " << total_error_omega1 << "\n";
    std::cout << "Ошибка по omega2: " << total_error_omega2 << "\n";
}


void runErrorCalculation(const State& initialState) {
    State state_rk4 = initialState;
    State state_dopri5 = initialState;
    double dt = 0.01;
    double time = 0.0;

    std::ofstream errorFile1("error1.txt");
    if (!errorFile1.is_open()) {
        std::cerr << "Не удалось открыть файл для записи ошибки.\n";
        return;
    }
    std::ofstream errorFile2("error2.txt");
    if (!errorFile2.is_open()) {
        std::cerr << "Не удалось открыть файл для записи ошибки.\n";
        return;
    }

    while (time <= 10.0) {
        rk4(state_rk4, dt);
        dopri5(state_dopri5, dt);

        writeErrorToFile1(state_rk4, state_dopri5, time, errorFile1);
        writeErrorToFile2(state_rk4, state_dopri5, time, errorFile2);

        time += dt;
    }

    errorFile1.close();
    errorFile2.close();
    
    std::cout << "Файл error1.txt: абсолютные ошибки\n";
    std::cout << "Файл error2.txt: относительные ошибки\n";
    
    calculateTotalError("error1.txt");
    calculateTotalError("error2.txt");
}

int main() {
    State initialState1 = {M_PI / 2, M_PI / 2, 0.0, 0.0};
    State initialState2 = {M_PI / 2, M_PI / 2, 0.0, 0.001};

    std::string mode;
    std::cout << "Введите режим ('show' или 'error'): ";
    std::cin >> mode;

    if (mode == "show") {
        std::string method;
        std::cout << "Введите метод ('rk4' или 'd5'): ";
        std::cin >> method;

        if (method == "rk4" || method == "d5") {
            runVisualization(initialState1, initialState2, method);
        } else {
            std::cout << "Неизвестный метод. По умолчанию используется 'rk4'.\n";
            runVisualization(initialState1, initialState2, "rk4");
        }
    } else if (mode == "error") {
        runErrorCalculation(initialState1);
    } else {
        std::cout << "Неизвестный режим. Завершение программы.\n";
    }

    return 0;
}

