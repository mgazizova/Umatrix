import numpy as np
from sympy.physics.wigner import gaunt
from cmath import sqrt
from math import pow
from sympy import pi


class Y:
    """
    Класс описывает сферические функции
    k - коэффициент при сферической функции в волновой функции
    array - массив для хранения произведений сферических функций
    """

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            try:
                if key == 'k':
                    self.k = value
                if key == 'l':
                    self.l = value
                if key == 'm':
                    self.m = value
                if key == 'array':
                    self.array = value
            finally:
                pass

    def __mul__(self, other):       # переопределенная функция для умножения объектов с типом Y
        return Y(
            k=self.k * other.k,
            array=[Y(l=self.l, m=self.m, k=1), Y(l=other.l, m=other.m, k=1)]
        )


def psi1():     # волновая функция
    y_array = np.array([Y(l=2, m=0, k=1)])
    return y_array


def psi2():     # волновая функция
    y_array = np.array([
        Y(l=2, m=-1, k=complex(0, 1 / sqrt(6))),
        Y(l=2, m=1, k=complex(0, 1 / sqrt(6))),
        Y(l=2, m=-2, k=complex(0, 1 / sqrt(3))),
        Y(l=2, m=2, k=complex(0, -1 / sqrt(3)))
    ])
    return y_array


def psi3():     # волновая функция
    y_array = np.array([
        Y(l=2, m=-1, k=complex(-1 / sqrt(6), 0)),
        Y(l=2, m=1, k=complex(1 / sqrt(6), 0)),
        Y(l=2, m=-2, k=complex(1 / sqrt(3), 0)),
        Y(l=2, m=2, k=complex(1 / sqrt(3), 0))
    ])
    return y_array


def multiply_y(ratio1, ratio2):     # умножение соответствующих волновых функций
    mul1_array = psi[ratio1]
    mul2_array = psi[ratio2]
    result_array = []
    for element1 in mul1_array:
        for element2 in mul2_array:
            result_array.append(element1 * element2)
    return result_array


def gaunt_integral(y1, y2, y3, y1_dash, y2_dash, y3_dash):      # расчет интеграла по d omega, d omega'
    return pow(-1, y1.m) * gaunt(y1.l, y2.l, y3.l, y1.m, y2.m, y3.m) * gaunt(y1_dash.l, y2_dash.l, y3_dash.l, y1_dash.m, y2_dash.m, y3_dash.m)


def find_m1(m2, m3):    # правило вычисления m1
    return 0 - m2 - m3


def calculation(ratio_u_matrix):
    omega_array = multiply_y(ratio_u_matrix[0], ratio_u_matrix[2])          # зависит от омега
    omega_array_dash = multiply_y(ratio_u_matrix[1], ratio_u_matrix[3])     # зависит от омега штрих
    u = {}      # словарь для хранения U-матрицы
    for l in l_array:       # цикл по l
        result = 0     # значение U-матрицы при заданном l
        for element in omega_array:
            m1 = find_m1(element.array[0].m, element.array[1].m)        # находим значение m1 для каждого произведения сферических функций, зависящего от omega
            for element_dash in omega_array_dash:
                m1_dash = find_m1(element_dash.array[0].m, element_dash.array[1].m)     # находим значение m1 для каждого произведения сферических функций, зависящего от omega'
                if m1 == -m1_dash:      # условие для соответствия коэффициентов m1
                    gaunt_ratio = gaunt_integral(       # расчет интеграла
                        Y(l=l, m=m1),
                        element.array[0],
                        element.array[1],
                        Y(l=l, m=m1_dash),
                        element_dash.array[0],
                        element_dash.array[1]
                    )
                    result += gaunt_ratio * element.k * element_dash.k  # домножаем на коэффициенты при сферических функциях
        result = result * 4 * pi / (2 * l + 1)
        u[l] = result       # добавляем значение U-матрицы в словарь для заданного l
    return u


u_matrix = np.array([           # коэффициенты U-матрицы
    [1, 1, 1, 1],
    [2, 2, 2, 2],
    [3, 3, 3, 3],
    [1, 1, 2, 2],
    [1, 1, 3, 3],
    [2, 2, 3, 3],
    [1, 2, 1, 2],
    [1, 3, 1, 3],
    [2, 3, 2, 3],
    [1, 2, 2, 1],
    [1, 3, 3, 1],
    [2, 3, 3, 2]
])


l_array = np.array([0, 2, 4])
psi = {1: psi1(), 2: psi2(), 3: psi3()}     # словарь для выбора волновой функции


if __name__ == '__main__':
    for item in u_matrix:       # цикл по строкам из массива коэффициентов U-матрицы
        u = calculation(item)
        print(item, u)





