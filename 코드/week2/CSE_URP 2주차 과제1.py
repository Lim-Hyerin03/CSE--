###3
import numpy as np
import matplotlib.pyplot as plt
import math
plt.rcParams['font.family'] = 'Malgun Gothic'
plt.rcParams['axes.unicode_minus'] = False

###2-1
def f(x):
    return np.sin((4-x)*(4+x))

N=33
h=8/(N-1)
xlist=np.linspace(0,8,N)

def fffr(x):
    global h
    return (2*f(x)-5*f(x+h)+4*f(x+2*h)-f(x+3*h))/h**2

def fffl(x):
    global h
    return (2*f(x)-5*f(x-h)+4*f(x-2*h)-f(x-3*h))/h**2

print(f"""Uniform nodes 33개, Second-order one-sided difference shceme을
사용하여 유도한 결과 f''(0)={fffr(0)}, f''(8)={fffl(8)}입니다.""")
# 경계 0에서는 오른쪽 세 점(0.25, 0.5, 0.75)을,
# 경계 8에서는 왼쪽 세 점 (7.25, 7.5, 7.75)을 추가로 이용하여
# f''의 근사값을 구하였다.
# 이에는 각 점에서의 테일러 전개를 조합하여 유도한
# f''=(2f(x)-5f(x-h)+4f(x-2h)-f(x-3h))/h^2 이 사용되었다.



###2-2
N=33
h=8/(N-1)
xlist=np.linspace(0,8,N)
xor=np.linspace(0,8,3200)

def fff(x):
    global h
    return (f(x+h)-2*f(x)+f(x-h))/h**2

def exact_fff(x):
    return -2*np.cos(16-x**2)-4*x**2*(np.sin(16-x**2))

def endpoints(x):
    x[0]=fffr(0)
    x[-1]=fffl(8)

fff_xlist=[fff(x) for x in xlist]
endpoints(fff_xlist)

plt.plot(xor, exact_fff(xor), 'k', label="exact solution")
plt.plot(xlist, fff_xlist, 'c', label='Second-order central difference scheme')
plt.legend(loc='lower center')
plt.title("Second-order central difference scheme")

plt.show()



###2-3
N33=33
h=8/(N33-1)
xlist33=np.linspace(0,8,N33)
fff_N33=[fff(x) for x in xlist33]
endpoints(fff_N33)

N65=65
h=8/(N65-1)
xlist65=np.linspace(0,8,N65)
fff_N65=[fff(x) for x in xlist65]
endpoints(fff_N65)
print(f'N=65, x=8에서의 함숫값: {fff_N65[-1]}')
print(exact_fff(8))

N129=129
h=8/(N129-1)
xlist129=np.linspace(0,8,N129)
fff_N129=[fff(x) for x in xlist129]
endpoints(fff_N129)


# 세 N 조건에서의 f''를 동시에 그려보았다.
plt.figure(figsize=(12,8))
plt.plot(xlist33, fff_N33, 'c', label='Second-order central difference scheme (N=33)')
plt.plot(xlist65, fff_N65, 'm', label='Second-order central difference scheme (N=65)')
plt.plot(xlist129, fff_N129, 'b', label='Second-order central difference scheme (N=129)')
plt.plot(xor, exact_fff(xor), 'k', label="exact solution")
plt.legend(loc='lower center')
plt.title("Second-order central difference scheme")
plt.show()
# N=65인 경우 x=8인 지점에서 함수값이 302로 튀는 것을 확인하였다. exact solution의 해는 -195이다.
# 이 오차가 아래의 오차 검증에도 큰 영향을 미쳐 오차 계산 시 이 점을 제외하고 진행하였다.


y_f=[np.sin((4-x)*(4+x)) for x in xlist65]
plt.plot(xlist65, y_f, c='b', label='sin((4-x)*(4+x))')
plt.scatter(xlist65, y_f, c='m', label='sin((4-x)*(4+x))')
plt.title("N=65일 때 f(x) 위에서 선택되는 점들")
plt.show()
# f(x)에 N=65인 경우 골라지는 점들을 찍어 확인해 보아도 원인을 모르겠다. 대체 왜 8에서 튀는 걸까.

 
error33=[(exact_fff(xlist33[i])-fff_N33[i])**2 for i in range(33)]
# x=8일 때 값이 튀는 것을 보정하기 위하여 range(65)가 아닌 range(64)로 진행하였다.
error65=[(exact_fff(xlist65[i])-fff_N65[i])**2 for i in range(64)]
error129=[(exact_fff(xlist129[i])-fff_N129[i])**2 for i in range(129)]
error=[np.log10(math.sqrt(sum(error33)*8/(N33-1))), np.log10(math.sqrt(sum(error65)*8/(N65-1))), np.log10(math.sqrt(sum(error129)*8/(N129-1)))]
log_h=[np.log10(8/(N33-1)), np.log10(8/(N65-1)), np.log10(8/(N129-1))]



plt.plot(log_h, error, 'c', label='정확도')
plt.text(log_h[0], error[0], 'N=33', fontsize=10, ha='center', va='top')
plt.text(log_h[1], error[1], 'N=65', fontsize=10, ha='center', va='top')
plt.text(log_h[2], error[2], 'N=129', fontsize=10, ha='center', va='top')
plt.title("Second-order central difference scheme의 정확도")
plt.xlabel('$\log(\Delta x)$')
plt.ylabel('$\log(\Vert e \Vert)$')
plt.show()

from scipy.stats import linregress

print(f'order of accuracy: {linregress(log_h, error).slope}')
# order of accuracy: 1.6388

# Δx와 ∣∣e∣∣에 로그를 취해 만들어진 함수의 기울기가
# 2와 유사한 것을 확인하였다.
# 이를 통해 오차가 Δx^2에 비례하여 감소함을 알 수 있다.