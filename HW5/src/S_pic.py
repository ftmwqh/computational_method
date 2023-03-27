import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']
def pic(filename1,filename2,filename3):
    f1=open(filename1,'r')
    f2=open(filename2,'r')
    f3 = open(filename3, 'r')
    lines1 = f1.readlines()
    lines2 = f2.readlines()
    lines3 = f3.readlines()
    m = len(lines1)
    print('height=', m)
    n = len(lines1[1].split(' '))
    print('width=', n)
    x = np.zeros(m, dtype=float)
    y=np.zeros(m,dtype=float)
    a=np.zeros(m-1, dtype=float)
    b= np.zeros(m-1, dtype=float)
    c=np.zeros(m-1,dtype=float)
    d = np.zeros(m-1, dtype=float)
    a_c = np.zeros(m - 1, dtype=float)
    b_c = np.zeros(m - 1, dtype=float)
    c_c = np.zeros(m - 1, dtype=float)
    d_c = np.zeros(m - 1, dtype=float)

    i = 0
    for line in lines1:
        list = line.strip('\n').split('\t')#除去换行符
        x[i] = list[0]
        y[i]=list[1]
        i += 1
        # 把文件读入到一个矩阵

    i=0
    for line in lines2:
        list = line.strip('\n').split('\t')#除去换行符
        a[i] = list[0]
        b[i]=list[1]
        c[i]=list[2]
        d[i]=list[3]
        i += 1
        # 把文件读入到一个矩阵

    i = 0
    for line in lines3:
        list = line.strip('\n').split('\t')  # 除去换行符
        a_c[i] = list[0]
        b_c[i] = list[1]
        c_c[i] = list[2]
        d_c[i] = list[3]
        i += 1
        # 把文件读入到一个矩阵
    axis_x = np.arange(x[0], x[-1] , 0.1)
    axis_y_1 = np.arange(x[0], x[-1] , 0.1)
    axis_y_2 = np.arange(x[0], x[-1], 0.1)
    for j in range(0,len(axis_y_1)):
        for i in range(0,20):
            if x[i]<=axis_x[j] and axis_x[j]<=x[i+1]:
                axis_y_1[j]=a[i]*axis_x[j]**3+b[i]*axis_x[j]**2+c[i]*axis_x[j]+d[i]
    for j in range(0, len(axis_y_2)):
        for i in range(0, 20):
            if x[i] <= axis_x[j] and axis_x[j] <= x[i + 1]:
                axis_y_2[j] = a_c[i] * axis_x[j] ** 3 + b_c[i] * axis_x[j] ** 2 + c_c[i] * axis_x[j] + d_c[i]



    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)

    ax1.scatter(x, y, label='压铁点', color='r')
    ax1.plot(axis_x, axis_y_1, label='S(x)', linewidth=1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title("(a)原始数据S(x)",y=-0.2)
    ax1.legend(loc='upper left')
    y[9]=10
    ax2.scatter(x, y, label='压铁点', color='r')
    ax2.plot(axis_x, axis_y_2, label='S(x)', linewidth=1)
    ax2.plot(axis_x, axis_y_1, label='原S(x)', linewidth=1)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title("(b)改变后S(x)",y=-0.2)
    ax2.legend(loc='upper left')
    plt.show()


pic('point.txt','coefficient.txt','coefficient_10.txt')