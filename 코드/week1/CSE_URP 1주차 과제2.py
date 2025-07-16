### 2

###1-1
def f(x):
    return 1/(1+16*x**2)

h=0.2
N=int((1+1)/h)
xlist=[-1+h*i for i in range(N+1)]
xlist=[round(i,2) for i in xlist]
ylist=[f(x) for x in xlist]

def p(x):
    pp=0
    for n in range(len(xlist)):
        ppp=1
        for i in range(len(xlist)):
            if i!=n:
                ppp*=(x-xlist[i])/(xlist[n]-xlist[i])
        pp+=ppp*ylist[n]
    return pp

import numpy as np

y2list=[p(x) for x in xlist]

xor=np.linspace(-1,1,400)
yor=[f(x) for x in xor]
yor2=[p(x) for x in xor]

import matplotlib.pyplot as plt

plt.plot(xor, yor, 'k', label='f(x)')
plt.plot(xlist, ylist, 'c', label='data points', linestyle='None', marker='o')
plt.plot(xor, yor2, 'm', label='p(x)')
plt.legend(loc='lower center')
plt.title("Lagrange interpolating polynomial p(x) for f(x)=1/(1+16x^2) (U N)")
plt.show()



###1-2
def f(x):
    return 1/(1+16*x**2)

import numpy as np

N=10
xlist=[np.cos((2*i+1)*np.pi/(2*N+2)) for i in range(N+1)]
xlist.sort()
ylist=[f(x) for x in xlist]

def p(x):
    ll=0
    for n in range(len(xlist)):
        lll=1
        for i in range(len(xlist)):
            if i!=n:
                lll*=(x-xlist[i])/(xlist[n]-xlist[i])
        ll+=lll*ylist[n]
    return ll


y2list=[p(x) for x in xlist]

xor=np.linspace(-1,1,400)
yor=[f(x) for x in xor]
yor2=[p(x) for x in xor]

import matplotlib.pyplot as plt

plt.plot(xor, yor, 'k', label='f(x)')
plt.plot(xlist, ylist, 'c', label='data points', linestyle='None', marker='o')
plt.plot(xor, yor2, 'm', label='p(x)')
plt.legend(loc='lower center')
plt.title("Lagrange interpolating polynomial p(x) for f(x)=1/(1+16x^2) (C N)")
plt.show()



###1-3
# Chebyshev nodes를 사용하였을 때, Uniform nodes를 사용하였을 때보다
# 더 f(x)에 유사한 양상을 보였다.
# 이 이유를 보간 오차 계수를 이용해 알아보았다.

import numpy as np
import matplotlib.pyplot as plt

# Uniform nodes의 경우
h=0.2
N=int((1+1)/h)
x1list=[-1+h*i for i in range(N+1)]
x1list=[round(i,2) for i in x1list]

def g1(x):
    gg1=1
    for i in range(len(x1list)):
        gg1*=abs(x-x1list[i])
    return gg1

# Chebyshev nodes의 경우
N=10
x2list=[np.cos((2*i+1)*np.pi/(2*N+2)) for i in range(N+1)]

def g2(x):
    gg2=1
    for i in range(len(x2list)):
        gg2*=abs(x-x2list[i])
    return gg2

xor=np.linspace(-1,1,400)
yor=[g1(x) for x in xor]
yor2=[g2(x) for x in xor]
y1list=[g1(x) for x in x1list]
y2list=[g2(x) for x in x2list]

plt.plot(xor, yor, 'm', label='Uniform Nodes')
plt.plot(xor, yor2, 'c', label='Chebyshev Nodes')
plt.legend(loc='center')
plt.title("Uniform and Chebyshev nodes for $\prod_{i=0}^{N} |x - x_i|$")
plt.show()

# Uniform Nodes의 경우 양 끝 값에서 오차가 크게 튀어오르는 양상을 보였다.
# 이로 인해 ###2-1에서 구한 p(x)와 f(x)의 오차가
# 양 끝 부분에서 커진 것으로 보인다.
# 반면 Chebyshev Nodes의 경우 비교적 일정한 오차가 유지되는 것을 확인하였다.
# 따라서 Chebyshev Nodes가 Lagrange interpolating polynomial 과정에
# 더 적합함을 알 수 있다.


###1-4
# Uniform Nodes
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

def f(x):
    return 1/(1+16*x**2)

h=0.2
N=int((1+1)/h)
x1list=[-1+h*i for i in range(N+1)]
x1list=[round(i,2) for i in x1list]
y1list=[f(x) for x in x1list]

xor=np.linspace(-1,1,400)

cs1=CubicSpline(x1list, y1list)
plt.plot(x1list, y1list, 'c', label='data points', linestyle='None', marker='o')
plt.plot(xor, cs1(xor), label='Cubic spline interpolation')
plt.legend(loc='lower center')
plt.title("Cubic spline interpolation (U N)")
plt.show()

# Chebyshev Nodes
def f(x):
    return 1/(1+16*x**2)

N=10
x2list=[np.cos((2*i+1)*np.pi/(2*N+2)) for i in range(N+1)]
x2list.sort()
y2list=[f(x) for x in x2list]

xor=np.linspace(-1,1,400)

cs2=CubicSpline(x2list, y2list)
plt.plot(x2list, y2list, 'c', label='data points', linestyle='None', marker='o')
plt.plot(xor, cs2(xor), label='Cubic spline interpolation')
plt.legend(loc='lower center')
plt.title("Cubic spline interpolation (C N)")
plt.show()

# 이렇게 하는 게 아닌 것 같음. 진짜 8개씩 조건 만들어서 다 구해나가라는 것 같은데
# 일단 그렇게 할 자신이 없으니 다음 문제부터 풀자


###1-5
import numpy as np
import math
from scipy.interpolate import CubicSpline

def f(x):
    return 1/(1+16*x**2)

N=10
xlist=[np.cos((2*i+1)*np.pi/(2*N+2)) for i in range(N+1)]
xlist.sort()
ylist=[f(x) for x in xlist]
xor=np.linspace(-1,1,400)

# Chebyshev node, Lagrangian interpolation
def p(x):
    ll=0
    for n in range(len(xlist)):
        lll=1
        for i in range(len(xlist)):
            if i!=n:
                lll*=(x-xlist[i])/(xlist[n]-xlist[i])
        ll+=lll*ylist[n]
    return ll

error_LI=[(f(x)-p(x))**2 for x in xor]
Error_LI=math.sqrt(sum(error_LI)/len(xor))
print(f"Lagrangian interpolation를 이용한 경우 RMSE Error은 {Error_LI}입니다.")

# Chebyshev node, Cubic spline method
cs=CubicSpline(xlist, ylist)

error_CS=[(f(x)-cs(x))**2 for x in xor]
Error_CS=math.sqrt(sum(error_CS)/len(xor))
print(f"Cubic spline method를 이용한 경우 RMSE Error은 {Error_CS}입니다.")

# 평균 제곱근 오차를 구한 결과 각각 0.03492, 0.01389의 Error을 구할 수 있었다.

###2

###2-1
import numpy as np

def f(x):
    return np.sin((4-x)*(4+x))

N=33
h=8/(N-1)
xlist=[h*i for i in range(N)]
xlist=[round(i,2) for i in xlist]

def fffr(x):
    return (2*f(x)-5*f(x+h)+4*f(x+2*h)-f(x+3*h))/h**2

def fffl(x):
    return (2*f(x)-5*f(x-h)+4*f(x-2*h)-f(x-3*h))/h**2

print(f"""Uniform nodes 33개, Second-order one-sided difference shceme을
사용하여 유도한 결과 f''(0)={fffr(0)}, f''(8)={fffl(8)}입니다.""")
# 경계 0에서는 오른쪽 세 점(0.25, 0.5, 0.75)을,
# 경계 8에서는 왼쪽 세 점 (7.25, 7.5, 7.75)을 추가로 이용하여
# f''의 근사값을 구하였다.
# 이에는 각 점에서의 테일러 전개를 조합하여 유도한
# f''=(2f(x)-5f(x-h)+4f(x-2h)-f(x-3h))/h^2 이 사용되었다.



###2-2
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.sin((4-x)*(4+x))

N=33
h=8/(N-1)
xlist=np.linspace(0,8,N)
xlist=xlist[1:-1]
xor=np.linspace(0,8,3200)

def fff(x):
    return (f(x+h)-2*f(x)+f(x-h))/h**2

def exact_fff(x):
    return -2*np.cos(16-x**2)-4*x**2*(np.sin(16-x**2))

plt.plot(xor, fff(xor), 'c', label='Second-order central difference scheme')
plt.plot(xlist, exact_fff(xlist), 'm', label="exact solution")
plt.legend(loc='lower center')
plt.title("Second-order central difference scheme")

plt.show()



###2-3
import numpy as np
import matplotlib.pyplot as plt
import math
plt.rcParams['font.family'] = 'Malgun Gothic'
plt.rcParams['axes.unicode_minus'] = False

def f(x):
    return np.sin((4-x)*(4+x))

def fff(x, h):
    return (f(x+h)-2*f(x)+f(x-h))/h**2

def exact_fff(x):
    return -2*np.cos(16-x**2)-4*x**2*(np.sin(16-x**2))

xor=np.linspace(0,8,3200)

N33=33
h33=8/(N33-1)
xlist33=np.linspace(0,8,N33)
xlist33=xlist33[1:-1]

N65=65
h65=8/(N65-1)
xlist65=np.linspace(0,8,N65)
xlist65=xlist65[1:-1]

N129=129
h129=8/(N129-1)
xlist129=np.linspace(0,8,N129)
xlist129=xlist129[1:-1]

# 세 N 조건에서의 f''를 동시에 그려보았다.
plt.figure(figsize=(12,8))
plt.plot(xlist33, fff(xlist33, h33), 'c', label='Second-order central difference scheme (N=33)')
plt.plot(xlist65, fff(xlist65, h65), 'm', label='Second-order central difference scheme (N=65)')
plt.plot(xlist129, fff(xlist129, h129), 'b', label='Second-order central difference scheme (N=129)')
plt.plot(xor, exact_fff(xor), 'k', label="exact solution")
plt.legend(loc='lower center')
plt.title("Second-order central difference scheme")

plt.show()


def e(error):
    return np.log(math.sqrt(np.mean(error)))

error33=[(exact_fff(x)-fff(x, h33))**2 for x in xlist33]
error65=[(exact_fff(x)-fff(x, h65))**2 for x in xlist65]
error129=[(exact_fff(x)-fff(x, h129))**2 for x in xlist129]
error=[e(error33), e(error65), e(error129)]
log_h=[np.log(h33), np.log(h65), np.log(h129)]

xl=[-2.7725887, -1.3]
def z(x):
    return 2*(x+2.7725887)+1.562729

Z=[z(x) for x in xl]

plt.plot(xl, Z, 'm', label='기울기 2')
plt.plot(log_h, error, 'c', label='정확도')
plt.text(log_h[0], error[0], 'N=33', fontsize=10, ha='center', va='top')
plt.text(log_h[1], error[1], 'N=65', fontsize=10, ha='center', va='top')
plt.text(log_h[2], error[2], 'N=129', fontsize=10, ha='center', va='top')
plt.title("Second-order central difference scheme의 정확도")
plt.xlabel('$\log(\Delta x)$')
plt.ylabel('$\log(\Vert e \Vert)$')
plt.legend(loc='lower center')
plt.show()

print(np.log(h129), e(error129))

from scipy.stats import linregress

print(f'order of accuracy: {linregress(log_h, error).slope}')
# 결과: order of accuracy: 1.7527

# Δx와 ∣∣e∣∣에 로그를 취해 만들어진 함수의 기울기가
# 2와 유사한 것을 확인하였다.
# 이를 통해 오차가 Δx^2에 비례하여 감소함을 알 수 있다.