import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']
def pic(filename):
    f=open(filename,'r')
    lines = f.readlines()
    m = len(lines)
    print('height=', m)
    n = len(lines[1].split(' '))
    print('width=', n)
    x = np.zeros(m, dtype=float)
    y=np.zeros(m,dtype=float)

    i = 0
    for line in lines:
        list = line.strip('\n').split('\t')
        print(list)
        x[i] = list[0]
        y[i]=list[1]
        i += 1
        # 把文件读入到一个矩阵

    plt.scatter(x[0:49],y[0:49],color='b')
    plt.scatter(x[50:99],y[50:99],color='r')
    plt.scatter(x[100:149], y[100:149],color='g')
    plt.show()

print(np.sin(np.pi/2))
pic("coordinate.txt")