import math
import matplotlib.pyplot as plt



def main():
    print()
    step, accuracy, num_func, x0, y0, a, b = read_data()

    modified_Euler_method(step, accuracy, get_function(num_func), x0, y0, a, b,get_exectly_function(num_func, x0, y0),num_func)

    runge(step, accuracy, get_function(num_func), x0, y0, a, b, get_exectly_function(num_func, x0, y0), num_func)

    if num_func!=3:
        miln(step, accuracy, get_function(num_func), x0, y0, a, b, get_exectly_function(num_func, x0, y0))
    else:
        print("Точного решения методом Милна для этого уравнения не существует")

def get_function(num):
    try:
        if num == 1:
            return lambda x, y: 4 * x + y / 3
        elif num == 2:
            return lambda x, y: x ** 2 + 3 * y / 2
        elif num == 3:
            return lambda x, y: x + 2 * y ** 2
        elif num == 4:
            return lambda x, y: x ** 2 - 2 * y
        else:
            return None
    except:
        print("Числа слишком большие")

def read_data():
    print("Выберите интересующее уравнение:")
    print("1: 4x + y/3")
    print("2: x^2 + 3y/2")
    print("3: x + 2y^2")
    print("4: x^2 - 2y")

    j = 0
    while j != 1:
        j = 1
        num_equation_str = input()

        if num_equation_str.isnumeric() == False:
            print("Вы должны ввести целое число")
            print("Повторите ввод:")
            j = 0
        else:
            num_equation = int(num_equation_str)
            if num_equation != 1 and num_equation != 2 and num_equation != 3 and num_equation != 4:
                j = 0
                print("Вы должны ввести номер уравнения соответственно 1, 2 или 3")
                print("Повторите ввод:")

    j = 0
    while j != 1:
        j = 1
        print("Введите x0:")
        x0_str = str(input())
        try:
            x0 = float(x0_str)
        except ValueError:
            print("X0 должно быть числом")
            j = 0

        print("Введите y0:")
        y0_str = str(input())
        try:
            y0 = float(y0_str)
        except ValueError:
            print("Y0 должно быть числом")
            j = 0

    j = 0
    while j != 1:
        j = 1
        print("Введите через пробел интервал дифферинцирования:")
        str_interval = input()
        str_interval = str_interval.replace(",", ".")
        str_interval.strip()
        low_border_str = str_interval.split()[0]
        height_border_str = str_interval.split()[1]

        try:
            low_border = float(low_border_str)
        except ValueError:
            print("Нижняя граница интервала не является числом")
            j = 0

        try:
            height_border = float(height_border_str)
        except ValueError:
            print("Верхняя граница интервала не является числом")
            j = 0

        if j == 0:
            print("Повторите ввод:")

    j = 0
    while j != 1:
        j = 1
        print("Введите через пробел шаг и точность:")
        str_step_acc = str(input())
        step_acc = str_step_acc.replace(",", ".")
        step_acc.strip()
        step_str = step_acc.split()[0]
        acc_str = step_acc.split()[1]

        try:
            step = float(step_str)
        except ValueError:
            print("Введенный шаг не является числом")
            j = 0

        try:
            acc = float(acc_str)
        except ValueError:
            print("Введенная точка не является числом")
            j = 0

        if j == 0:
            print("Повторите ввод:")

    return step, acc, num_equation, x0, y0, low_border, height_border

def runge(h, acc, func, x0, y0, a, b, func_exect,num_func):
    print("Метод Рунге-Кутта 4- го порядка")

    n = int((b - a) / h)
    p = 2
    x_h=x0
    y_h=y0
    x_seg = x0
    y_seg = y0
    step=h
    k=1
    r = 1.0
    while r>acc:
        if num_func != 3:
            list_f = []
            list_f.append(func(x0, y0))

            list_f_e = []
            list_f_e.append(func_exect(x0, y0))

        print("Решение с шагом:", step)
        print("| i |  x  |  y  |  k1  |  k2  |  k3  |  k4  |  y(i+1)  |")

        x_h = x0
        y_h = y0
        for i in range(0, n * k + 1):
            y_h = find_y_runge(i, x_h, y_h, step, func, True)
            x_h += step
            if num_func != 3:
                list_f.append(func(x_h, y_h))
                list_f_e.append(func_exect(x_h, y_h))


        print("| ",n * k + 1," | ", x_h, " | ", y_h, " |")
        print()
        print()
        print()

        step = step / 2
        k = k * 2

        x_seg = x0
        y_seg = y0
        for i in range(0, n * k + 1):
            y_seg = find_y_runge(i, x_seg, y_seg, step, func, False)
            x_seg += step

            if i>1000:
                print("Превышен лимит итераций, точность не достигнута")
                return

        r = abs(y_h - y_seg) / (2 ** p - 1)
    print("Итоговый ответ методом Рунге равен:", y_h)
    print("\n\n")
    if num_func!=3:
        drow_graph(list_f, list_f_e, "Рунге")
    else:
        print('Для этого уравнения не существует точного решения')


def find_y_runge(i, x, y, h, func, write):

    k1=h*func(x,y)
    k2=h*func(x+h/2, y+k1/2)
    k3=h*func(x+h/2,y+k2/2)
    k4=h*func(x+h, y+k3)
    y2=y+1/6 * (k1+2*k2+2*k3+k4)

    if write==True:
        print("| ", i ," | ", x," | ",y," | ",k1," | ",k2," | ",k3," | ",k4," | ",y2, " |")


    return y2


def modified_Euler_method(h, acc, func, x0, y0, a, b, func_exect, num_func):
    print("Усовершенствованный метод Эйлера")
    n = int((b - a) / h)
    x = x0
    y = y0

    p = 2

    step = h
    k = 1
    x_h = x
    y_h = y
    x_seg = x
    y_seg = y
    r = 1.0
    while r > acc:
        if num_func!=3:
            list_f=[]
            list_f.append(func(x,y))

            list_f_e=[]
            list_f_e.append(func_exect(x,y))

        print("Решение с шагом:",step)
        print("| i |  x  |  y  |  x + h/2  |  f(x,y)  |  y + h/2 * f(x,y)  |  f(x + h/2; y + h/2 * f(x,y)  |")

        x_h = x
        y_h = y
        for i in range(0, n * k):
            y_h = find_y_euler(i, x_h, y_h, step, func, True)
            x_h += step
            if num_func!=3:
                list_f.append(func(x_h,y_h))
                list_f_e.append(func_exect(x_h,y_h))

        print("| ", x_h, " | ", y_h, " |")
        print()
        print()
        print()

        step = step / 2
        k=k*2

        x_seg = x
        y_seg = y
        for i in range(0, n * k):
            y_seg = find_y_euler(i, x_seg, y_seg, step, func, False)
            x_seg += step
            if i>1000:
                print("Превышен лимит итераций, точность не достигнута")
                return

        r = abs(y_h - y_seg) / (2 ** p - 1)

    # r = 0.0
    # while r <= acc:
    #     print("| i |  x  |  y  |  x + h/2  |  f(x,y)  |  y + h/2 * f(x,y)  |  f(x + h/2; y + h/2 * f(x,y)  |")
    #
    #     x_h = x
    #     y_h = y
    #     for i in range(0, n * k + 1):
    #         y_h = find_y_euler(x_h, y_h, step, func, True)
    #         x_h += step
    #
    #     print("| ", x_h, " | ", y_h, " |")
    #
    #     step = step / 2
    #     k = k * 2
    #
    #     x_seg = x
    #     y_seg = y
    #     for i in range(0, n * k + 1):
    #         y_seg = find_y_euler(x_seg, y_seg, step, func, False)
    #         x_seg += step
    #
    #     r = abs(y_h - y_seg) / (2 ** p - 1)

    print("Итоговый ответ методом Эйлера равен:", y_h)
    print("\n \n")
    if num_func!=3:
        drow_graph(list_f, list_f_e, "Эйлер")
    else:
        print('Для этого уравнения не существует точного решения')




    # for i in range(0, n + 1):
    #     # step=h
    #     segments = 2
    #     y_h = find_y_euler(x, y, h, func, True)
    #
    #     next = True
    #     while next:
    #         # x_h=x
    #         # y_h=y
    #         # for k in range(0,segments//2):
    #
    #         # y_h = find_y_euler(x, y, step, func, True)
    #         # x_h+=step
    #
    #         x_seg = x
    #         y_seg = y
    #         for j in range(0, segments):
    #             y_seg = find_y_euler(x_seg, y_seg, h / segments, func, False)
    #             x_seg += h / segments
    #
    #         r = abs(y_h - y_seg) / (2 ** p - 1)
    #         if r <= acc:
    #             next = False
    #             x += h
    #             y = y_seg
    #         else:
    #             segments *= 2
    #             # step=step/2


def find_y_euler(i, x, y, h, func, write):
    x_h2 = x + h / 2
    f_xy = func(x, y)
    y_h2_f_xy = y + h / 2 * f_xy
    f_xh2_yh2 = func(x_h2, y_h2_f_xy)
    delt_y = f_xh2_yh2 * h


    if write==True:
        print("| ", i ," | ", x," | ",y," | ",x_h2," | ",f_xy," | ",y_h2_f_xy," | ",f_xh2_yh2," | ",delt_y)


    return y + delt_y


def y_runge_miln(x, y, h, func):

    k1=h*func(x,y)
    k2=h*func(x+h/2, y+k1/2)
    k3=h*func(x+h/2,y+k2/2)
    k4=h*func(x+h, y+k3)
    y2=y+1/6 * (k1+2*k2+2*k3+k4)
    return y2


def get_exectly_function(num, x0, y0):
    if num==1:
        c=(y0+12*x0+36)/math.exp(x0/3)
        return lambda x, y: c*math.exp(x/3) - 12*x-36
    elif num==2:
        c=(y0 + (2*(x0**2)/3) + 8*x0/9 + 16/27)/math.exp(3*x0/2)
        return lambda x,y: c*math.exp(3*x/2 - (2*(x**2)/3) - 8*x/9 - 16/27)
    elif num==4:
        c = (y0-x0**2 /2 +x0/2 - 0.25)*math.exp(2*x0)
        return lambda x, y: (c/math.exp(2*x))+(x**2 / 2) - x/2 + 0.25

def miln(h, acc, func, x0, y0, a, b, func_exect):
    print("Метод Милна")
    step=h
    k=1
    r=1.0
    n = int((b - a) / h)
    while r>acc:
        print("Решение с шагом:",step)
        print("| i |  x  |  y  |  y точное  |  погрешность  |")
        x=x0
        y=y0
        list_x=[]
        list_y = []
        list_f = []
        list_f_e=[]
        list_error = []
        list_x.append(x)
        list_y.append(y)
        f_xy = func(x, y)
        f_e_xy = func_exect(x, y)
        list_f.append(f_xy)
        list_f_e.append(f_e_xy)
        error=abs(f_e_xy-y)
        list_error.append(error)


        print("| ",0," | ", x," | ",y," | ",f_e_xy," | ",error," |")
        for i in range(1,4):
            y=y_runge_miln(x,y,step,func)
            x=x+step

            f_xy=func(x, y)
            f_e_xy=func_exect(x, y)
            error = abs(f_e_xy - y)

            print("| ",i," | ", x," | ",y," | ",f_e_xy," | ",error," |")

            list_x.append(x)
            list_y.append(y)
            list_f.append(f_xy)
            list_f_e.append(f_e_xy)
            list_error.append(error)


        for i in range(4, n*k+1):
            x+=step
            check=False
            y_pred=list_y[i-4]+(4*step/3)*(2*list_f[i-3]-list_f[i-2]+2*list_f[i-1])
            while check != True:
                f_pred=func(x,y_pred)
                y_cor=list_y[i-2]+step/3 * (list_f[i-2] +4*list_f[i-1]+f_pred)
                if abs(y_cor-y_pred)>=acc:
                    y_pred=y_cor
                else:
                    check=True
                    y=y_cor

            if i>1000:
                print("Превышен лимит итераций, точность не достигнута")
                return

            f_xy = func(x, y)
            f_e_xy = func_exect(x, y)
            error = abs(f_e_xy - y)

            list_x.append(x)
            list_y.append(y)
            list_f.append(f_xy)
            list_f_e.append(f_e_xy)
            list_error.append(error)

            print("| ",i," | ", x," | ",y," | ",f_e_xy," | ",error," |")

        r=max(list_error)
        if r>acc:
            step=step/2
            k=k*2
            print("\n\n")
        else:

            print("Итоговый ответ методом Милна:", y)
            print("\n\n")
            drow_graph(list_f,list_f_e, "Милн")

def drow_graph(f, f_e, name):
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.title('График функций')
    plt.plot(f, label=name)
    plt.plot(f_e, label='Точное значение')

    plt.legend()
    plt.show()

main()





