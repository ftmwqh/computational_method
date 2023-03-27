import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']
def pic_FFT(filename1):
    f1=open(filename1,'r')

    lines1 = f1.readlines()

    m = len(lines1)
    print('height=', m)
    n = len(lines1[1].split(' '))
    print('width=', n)
    x = np.arange(0,m)
    f=np.zeros(m,dtype=float)
    g=np.zeros(m, dtype=float)
    i_f= np.zeros(m, dtype=float)

    i = 0
    for line in lines1:
        list = line.strip('\n').split('\t')#除去换行符
        f[i] = list[0]
        g[i]=list[1]
        i_f[i]=list[2]
        print(list)
        i += 1
        # 把文件读入到一个矩阵

    fig=plt.figure()
    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3 = fig.add_subplot(223)

    ax1.plot(x, f, label='f', linewidth=1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title("原f(x)离散点",y=-0.2)
    ax1.legend(loc='upper left')

    ax2.plot(x, g, label='g', linewidth=1)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title("FFT结果g",y=-0.2)
    ax2.legend(loc='upper left')

    ax3.plot(x, i_f, label='IFFT(FFT(f))', linewidth=1)
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_title("IFFT结果", y=-0.2)
    ax3.legend(loc='upper left')
    plt.show()

def pic_IFFT(filename1,filename2):
    f1=open(filename1,'r')
    f2 = open(filename2, 'r')

    lines1 = f1.readlines()
    lines2 = f2.readlines()

    m = len(lines1)
    m_1=len(lines2)
    print('height=', m)
    n = len(lines1[1].split(' '))
    print('width=', n)
    x = np.arange(0,m)
    x_1=np.arange(0,m_1)
    f=np.zeros(m,dtype=float)
    g=np.zeros(m, dtype=float)
    i_f= np.zeros(m, dtype=float)
    if_025 = np.zeros(m_1, dtype=float)

    i = 0
    for line in lines1:
        list = line.strip('\n').split('\t')#除去换行符
        f[i] = list[0]
        g[i]=list[1]
        i_f[i]=list[2]
        print(list)
        i += 1
        # 把文件读入到一个矩阵

    i = 0
    for line in lines2:
        list = line.strip('\n').split('\t')  # 除去换行符
        if_025[i] = list[0]
        print(list)
        i += 1

    fig=plt.figure()
    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax1.plot(x, f, label='f', linewidth=1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title("原f(x)离散点",y=-0.2)
    ax1.legend(loc='upper left')

    ax2.plot(x, g, label='g', linewidth=1)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title("FFT结果g",y=-0.2)
    ax2.legend(loc='upper left')

    ax3.plot(x, i_f, label='IFFT(FFT(f))', linewidth=1)
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_title("IFFT结果", y=-0.2)
    ax3.legend(loc='upper left')

    ax4.plot(x, if_025, label='IFFT(FFT(f))', linewidth=1)
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.set_title("前25%IFFT结果", y=-0.2)
    ax4.legend(loc='upper left')
    plt.show()

pic_FFT('f_1_4.txt')
pic_FFT('f_1_7.txt')
pic_IFFT('f_2.txt','if_2_025.txt')