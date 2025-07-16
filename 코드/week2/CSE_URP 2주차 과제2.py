###4

###1
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Malgun Gothic'
plt.rcParams['axes.unicode_minus'] = False

def f(t, x):
    return t+2*x*t

###1-1,2
def list_maker(h):
    t=np.linspace(0, 2, int(2/h)+1)
    x1=np.zeros(int(2/h)+1)
    x2=np.zeros(int(2/h)+1)
    for n in range(int(2/h)):
        k1=f(t[n], x1[n])
        k2=f(t[n]+h, x1[n]+k1*h)
        x1[n+1]=x1[n]+h*(k1+k2)/2
    for n in range(int(2/h)):
        k1=f(t[n], x2[n])
        k2=f(t[n]+h/2, x2[n]+k1*h/2)
        k3=f(t[n]+h/2, x2[n]+k2*h/2)
        k4=f(t[n]+h, x2[n]+k3*h)
        x2[n+1]=x2[n]+(k1+2*k2+2*k3+k4)*h/6
    return t, x1, x2

t, x1, x2=list_maker(0.01)

plt.plot(t, x1, 'c')
plt.xlabel('t', fontsize=15)
plt.ylabel('x', fontsize=15)
plt.title('x(t) using RK2 with h=0.01', fontsize=15)
plt.show()

plt.plot(t, x2, 'm')
plt.xlabel('t', fontsize=15)
plt.ylabel('x', fontsize=15)
plt.title('x(t) using RK4 with h=0.01', fontsize=15)
plt.show()




###1-3
# 1-1, 1-2에서 그린 두 그래프 RK2와 RK4 함께 나타내기
plt.plot(t, x1, 'c', label='RK2')
plt.plot(t, x2, 'm', label='RK4')
plt.legend(loc='center')
plt.xlabel('t', fontsize=15)
plt.ylabel('x', fontsize=15)
plt.title('x(t) using RK2, RK4 with h=0.01', fontsize=15)
plt.show()
# 육안으로는 차이를 찾을 수 없을 정도로 겹쳐 보인다.

# error 계산 (모든 t에 대해 평균 오차)
def exact_sol(t):
    return abs(np.exp(t**2)-1)/2

x=[exact_sol(t) for t in t]

def error1(h):
    t, x1, x2=list_maker(h)
    N=0
    for n in range(int(2/h)):
        N+=(exact_sol(t[n])-x1[n])**2
    return np.log(np.sqrt(N*h/(int(2/h)+1)))

def error2(h):
    t, x1, x2=list_maker(h)
    N=0
    for n in range(int(2/h)):
        N+=(exact_sol(t[n])-x2[n])**2
    return np.log(np.sqrt(N*h/(int(2/h)+1)))

Error_RK2=[error1(0.1), error1(0.05), error1(0.01)]
Error_RK4=[error2(0.1), error2(0.05), error2(0.01)]
log_h=[np.log(0.1), np.log(0.05), np.log(0.01)]

plt.xlabel('logh', fontsize=12)
plt.ylabel('log||e||', fontsize=12)
plt.plot(log_h, Error_RK2, 'c', label='RK2')
plt.plot(log_h, Error_RK4, 'm', label='RK4')
plt.title('RK2와 RK4의 order of accuracy', fontsize=15)
plt.legend(loc='lower center')
plt.figtext(0.11, 0.02,
            '''RK2의 order of accuracy:2.2950,   RK4의 order of accuracy:4.2575''',
            fontweight='bold', fontsize=12)
plt.subplots_adjust(bottom=0.15)
plt.show()

# error 그래프 기울기 계산
def order_of_accuracy(log_h, log_error):
    return linregress(log_h, log_error).slope

slope_RK2 = order_of_accuracy(log_h, Error_RK2)
print(f'RK2의 order of accuracy: {slope_RK2}')
slope_RK4 = order_of_accuracy(log_h, Error_RK4)
print(f'RK4의 order of accuracy: {slope_RK4}')



###1-4
#1-3과의 다른 점은 t가 t=2로 고정되어 있다는 것이다.
#t=0에서 시작하여 근사를 통해 최종적으로 t=2의 값을 구하는 것이 목표인 만큼,
#최종 t=2에서의 error를 구하여 h의 크기가 error에 갖는 의미를 찾아보는 것이다.

def RK4_e_in2(h):
    t, x1, x2=list_maker(h)
    return np.log(abs(exact_sol(2)-x2[-1]))

RK4_Error_in2=[RK4_e_in2(0.1), RK4_e_in2(0.05), RK4_e_in2(0.01)]

plt.xlabel('logh', fontsize=12)
plt.ylabel('log|e|', fontsize=12)
plt.plot(log_h, RK4_Error_in2, 'b', label='RK2')
plt.scatter(log_h, RK4_Error_in2, color='b', s=40)
plt.title('RK4, t=2에서의 order of accuracy', fontsize=15)
plt.legend(loc='lower center')
plt.figtext(0.25, 0.02,
            'RK4, t=2에서의 order of accuracy: 3.9008',
            fontweight='bold', fontsize=12)
plt.subplots_adjust(bottom=0.15)
plt.text(log_h[0], RK4_Error_in2[0], 'h=0.1   ', ha='right', va='center')
plt.text(log_h[1], RK4_Error_in2[1], '   h=0.05', ha='left', va='center')
plt.text(log_h[2], RK4_Error_in2[2], '   h=0.01', ha='left', va='center')
plt.show()

print(f'RK4, t=2에서의 order of accuracy: {linregress(log_h, RK4_Error_in2).slope}')
# RK4, t=2에서의 order of accuracy: 3.9008





###2
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

###2-1
#상수 정의 (v0, h, phi는 아직 주어지지 않았으므로 0을 넣었다.)
g=9.81
B=4.1e-4
w=1800*2*np.pi/60
theta=np.radians(1)
v0=0
h=0
phi=0

# 초기 조건 설정
x0=0
y0=0
z0=h
vx0=v0*np.cos(theta)
vy0=0
vz0=v0*np.sin(theta)
S=np.array([x0, y0, z0, vx0, vy0, vz0])

# F(V)
def F(V):
    return 0.0039+0.0058/(1+np.exp((V-35)/5))

# x, y, z, vx, vy, vz 한 번에 묶기
def deriv(S):
    x, y, z, vx, vy, vz=S
    V=np.sqrt(vx**2+vy**2+vz**2)
    dS_dt=np.zeros(6)
    dS_dt[0]=vx
    dS_dt[1]=vy
    dS_dt[2]=vz
    dS_dt[3]=-F(V)*vx+B*w*(vz*np.sin(phi)-vy*np.cos(phi))
    dS_dt[4]=-F(V)*vy+B*w*vx*np.sin(phi)
    dS_dt[5]=-g-F(V)*V*vz-B*w*vx*np.sin(phi)
    return dS_dt

# 시간 설정, 시간에 따른 S 보존
dt=0.001
every_S=[S.copy()]
every_t=[0]

#RK4 (조건: x(t)=18.39까지)
while S[0]<18.39:
    k1=deriv(S)
    k2=deriv(S+dt*k1/2)
    k3=deriv(S+dt*k2/2)
    k4=deriv(S+dt*k3)
    S+=(k1+2*k2+2*k3+k4)*dt/6
    every_t.append(every_t[-1] + dt)

# x, y, z 각자만의 리스트로 정리
every_S=np.array(every_S)
x=every_S[:,0]
y=every_S[:,1]
z=every_S[:,2]
vx=every_S[:,3]
vy=every_S[:,4]
vz=every_S[:,5]
# 위의 6개의 리스트가 t에 대한 6개의 방정식의 값이며,
# 각 리스트의 마지막 값이 x(t)=18.39일 때의 값이 된다.




###2-2
g=9.81
B=4.1e-4
w=1800*2*np.pi/60
theta=np.radians(1)

dt=0.0001
every_t=[0]

def F(V):
    return 0.0039+0.0058/(1+np.exp((V-35)/5))

def deriv(S):
    x, y, z, vx, vy, vz=S
    global phi
    V=np.sqrt(vx**2+vy**2+vz**2)
    dS_dt=np.zeros(6)
    dS_dt[0]=vx
    dS_dt[1]=vy
    dS_dt[2]=vz
    dS_dt[3]=-F(V)*vx+B*w*(vz*np.sin(phi)-vy*np.cos(phi))
    dS_dt[4]=-F(V)*vy+B*w*vx*np.cos(phi)
    dS_dt[5]=-g-F(V)*V*vz-B*w*vx*np.sin(phi)
    return dS_dt

def RK4(S):
    while S[0]<18.39:
        k1=deriv(S)
        k2=deriv(S+dt*k1/2)
        k3=deriv(S+dt*k2/2)
        k4=deriv(S+dt*k3)
        S+=(k1+2*k2+2*k3+k4)*dt/6
        every_S.append(S.copy())
        every_t.append(every_t[-1] + dt)
    return every_S

# fastball
phi=np.radians(225)
S=np.array([0, 0, 1.7, 40*np.cos(theta), 0, 40*np.sin(theta)])
every_S=[S.copy()]
RK4(S)
every_S_fastball=np.array(every_S)
x_fastball=every_S_fastball[:,0]
y_fastball=every_S_fastball[:,1]
z_fastball=every_S_fastball[:,2]
print("포수에게 가장 가까운 두 좌표 (fastball)")
print(f'({x_fastball[-2]}, {y_fastball[-2]}, {z_fastball[-2]})')   # :18.3868 -0.2311 1.2242
print(f'({x_fastball[-1]}, {y_fastball[-1]}, {z_fastball[-1]})')   # :18.3908 -0.2312 1.2239

def xyz_3d(x, y, z , title):
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, c='c', label="fastball")
    ax.set_xlabel('X', fontsize=15, labelpad=20)
    ax.set_ylabel('Y', fontsize=15, labelpad=20)
    ax.set_zlabel('Z', fontsize=15, labelpad=20)
    ax.set_title(title, fontsize=15)
    ax.legend()
    plt.show()

xyz_3d(x_fastball, y_fastball, z_fastball, "fastball's trajectory" )


# curveball
phi=np.radians(45)
S=np.array([0, 0, 1.7, 30*np.cos(theta), 0, 30*np.sin(theta)])
every_S=[S.copy()]
RK4(S)
every_S_curveball=np.array(every_S)
x_curveball=every_S_curveball[:,0]
y_curveball=every_S_curveball[:,1]
z_curveball=every_S_curveball[:,2]

print("포수에게 가장 가까운 두 좌표 (curveball)")
print(f'({x_curveball[-2]}, {y_curveball[-2]}, {z_curveball[-2]})')   # :18.3883 0.3091 -0.0629
print(f'({x_curveball[-1]}, {y_curveball[-1]}, {z_curveball[-1]})')   # :18.3914 0.3092 -0.0634

def xyz_3d(x, y, z , title):
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, c='b', label="curveball")
    ax.set_xlabel('X', fontsize=15, labelpad=20)
    ax.set_ylabel('Y', fontsize=15, labelpad=20)
    ax.set_zlabel('Z', fontsize=15, labelpad=20)
    ax.set_title(title, fontsize=15)
    ax.legend()
    plt.show()

xyz_3d(x_curveball, y_curveball, z_curveball, "curveball's trajectory" )


# slider
phi=np.radians(0)
S=np.array([0, 0, 1.7, 30*np.cos(theta), 0, 30*np.sin(theta)])
every_S=[S.copy()]
RK4(S)
every_S_slider=np.array(every_S)
x_slider=every_S_slider[:,0]
y_slider=every_S_slider[:,1]
z_slider=every_S_slider[:,2]
print("포수에게 가장 가까운 두 좌표 (slider)")
print(f'({x_slider[-2]}, {y_slider[-2]}, {z_slider[-2]})')   # :18.3884 0.3091 -0.0629
print(f'({x_slider[-1]}, {y_slider[-1]}, {z_slider[-1]})')   # :18.3914 0.3092 -0.0635

def xyz_3d(x, y, z , title):
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, c='g', label="slider")
    ax.set_xlabel('X', fontsize=15, labelpad=20)
    ax.set_ylabel('Y', fontsize=15, labelpad=20)
    ax.set_zlabel('Z', fontsize=15, labelpad=20)
    ax.set_title(title, fontsize=15)
    ax.legend()
    plt.show()

xyz_3d(x_slider, y_slider, z_slider, "slider's trajectory" )



#screwball
phi=np.radians(135)
S=np.array([0, 0, 1.7, 30*np.cos(theta), 0, 30*np.sin(theta)])
every_S=[S.copy()]
RK4(S)
every_S_screwball=np.array(every_S)
x_screwball=every_S_screwball[:,0]
y_screwball=every_S_screwball[:,1]
z_screwball=every_S_screwball[:,2]
print("포수에게 가장 가까운 두 좌표 (screwball)")
print(f'({x_screwball[-2]}, {y_screwball[-2]}, {z_screwball[-2]})')   # :18.3892 0.0 0.2354
print(f'({x_screwball[-1]}, {y_screwball[-1]}, {z_screwball[-1]})')   # :18.3921 0.0 0.2348

def xyz_3d(x, y, z , title):
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, c='m', label="screwball")
    ax.set_xlabel('X', fontsize=15, labelpad=20)
    ax.set_ylabel('Y', fontsize=15, labelpad=20)
    ax.set_zlabel('Z', fontsize=15, labelpad=20)
    ax.set_title(title, fontsize=15)
    ax.legend()
    plt.show()

xyz_3d(x_screwball, y_screwball, z_screwball, "screwball's trajectory" )


# 4개 공 한 번에 그리기
def xyz_3d_all(xyz, title):
    fig=plt.figure(figsize=(10, 7))
    ax=fig.add_subplot(111, projection='3d')
    color_set=['c', 'b', 'g', 'm']
    for i, (x, y, z, label) in enumerate(xyz):
        ax.plot(x, y, z, label=label, c=color_set[i % len(color_set)])
    ax.set_xlabel('X', fontsize=15, labelpad=20)
    ax.set_ylabel('Y', fontsize=15, labelpad=20)
    ax.set_zlabel('Z', fontsize=15, labelpad=20)
    ax.set_title(title, fontsize=15)
    ax.legend()
    plt.show()

xyz_3d_all([
(x_fastball, y_fastball, z_fastball, 'Fastball'),
    (x_curveball, y_curveball, z_curveball, 'Curveball'),
    (x_slider, y_slider, z_slider, 'Slider'),
    (x_screwball, y_screwball, z_screwball, 'Screwball')
], "trajectory of 4 balls")