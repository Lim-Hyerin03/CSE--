###1

###1-1
import random
i=0
slist=[]
while i<100:
    x=random.randint(-1000, 1000)
    slist.append(x)
    i+=1
slist.sort()
print(slist)



### 1-2
# matrix C 만들기
import numpy as np
A=np.array([[1,2,3,4],[2,0,1,1],[6,3,9,8],[4,2,8,5]])
B=np.array([[6,0,3,-4],[-8,0,2,5],[1,-5,0,4],[-1,-4,-6,7]])
C=A@B

#det(A)\=0, det(B)\=0 확인
print(np.linalg.det(A))
print(np.linalg.det(B))

#C-1 구하기
# A×B=C 이므로 C-1=B-1×A-1
A_1=np.linalg.inv(A)
B_1=np.linalg.inv(B)
C_1=B_1@A_1
print(C_1)

#identity matrix
print(C@C_1)

#오차 판정
if np.allclose(C@C_1, np.eye(4))==True:
    print("구한 C×C-1가 indentity matrix와 유사합니다.")



###2
def f(x):
    return x**5-9*x**4-x**3+17*x**2-8*x-8

def ff(x):
    return 5*x**4-36*x**3-3*x**2-34*x-8

###2-1
#[-10,-1] (찾은 근: -1.3875...)
if f(-10)*f(-1)<=0:
    a=-10
    b=-1
    n11=0
while b-a>(0.1)*0.00000001:
    m=(a+b)/2
    if f(a)*f(m)<=0:
        b=m
    else:
        a=m
    n11+=1
print(f"이분법을 통해 구한 [-10,-1]의 근은 {m}입니다.")

#[-1,0] (찾은 근: -0.5104...)
if f(-1)*f(0)<=0:
    a=-1
    b=0
    n12=0
while b-a>0.00000001:
    m=(a+b)/2
    if f(a)*f(m)<=0:
        b=m
    else:
        a=m
    n12+=1
print(f"이분법을 통해 구한 [-1,0]의 근은 {m}입니다.")

#[0,10] (찾은 근: 8.9107...)
if f(0)*f(10)<=0:
    a=0
    b=10
    n13=0
while b-a>0.00000001:
    m=(a+b)/2
    if f(a)*f(m)<=0:
        b=m
    else:
        a=m
    n13+=1
print(f"이분법을 통해 구한 [0,10]의 근은 {m}입니다.")


###2-2
# x0=-10 (찾은 근: -1.3875...)
x0=-10
n21=0
while abs(f(x0))>0.00000001:
    x0=x0-f(x0)/ff(x0)
    n21+=1
print(f"뉴턴법을 통해 구한 x0=-10에서의 근은{x0}입니다.")
#print(n21)

# x0=-0.1 (찾은 근: -1.3875...)
x0=-0.1
n22=0
while abs(f(x0))>0.00000001:
    x0=x0-f(x0)/ff(x0)
    n22+=1
print(f"뉴턴법을 통해 구한 x0=-0.1에서의 근은{x0}입니다.")
#print(n22)

# x0=-0.1인 경우 -10에서와 마찬가지로 -1.3875...의 근을 찾게 됨.

#x0=10 (찾은 근: 8.9107...)
x0=10
n23=0
while abs(f(x0))>0.00000001:
    x0=x0-f(x0)/ff(x0)
    n23+=1
print(f"뉴턴법을 통해 구한 x0=10에서의 근은{x0}입니다.")
#print(n23)


###2-3 (찾은 근: -1.3875...)
x0=0
n3=0
while abs(f(x0))>0.00000001:
    x0=x0-f(x0)/ff(x0)
    n3+=1
print(f"뉴턴법을 통해 구한 x0=0에서의 근은{x0}입니다.")
#print(n3)
# x0=-0.5104가 아닌 x0=-1.3875의 근을 찾음.
# 0이 -1.3875로 수렴하는 영역에 속함을 알 수 있음.



###2-4
### x1=-10, x2=-9.9 (찾은 근: -1.3875...)
xn2=-10
xn1=-9.9
n41=0
while abs(f(xn1))>=0.00000001:
    x=xn1-f(xn1)*(xn1-xn2)/(f(xn1)-f(xn2))
    xn2=xn1
    xn1=x
    n41+=1
print(f"시컨트법을 통해 구한 x1=-10, x2=-9.9에서의 근은{xn1}입니다.")
#print(n41)

### x1=-0.1, x2=-0.2 (찾은 근: -0.5104...)
xn2=-0.1
xn1=-0.2
n42=0
while abs(f(xn1))>=0.00000001:
    x=xn1-f(xn1)*(xn1-xn2)/(f(xn1)-f(xn2))
    xn2=xn1
    xn1=x
    n42+=1
print(f"시컨트법을 통해 구한 x1=-0.1, x2=-0.2에서의 근은{xn1}입니다.")
#print(n42)

### x1=10, x2=9.9 (찾은 근: 8.9107...)
xn2=10
xn1=9.9
n43=0
while abs(f(xn1))>=0.00000001:
    x=xn1-f(xn1)*(xn1-xn2)/(f(xn1)-f(xn2))
    xn2=xn1
    xn1=x
    n43+=1
print(f"시컨트법을 통해 구한 x1=10, x2=9.9에서의 근은{xn1}입니다.")
#print(n43)


###2-5
print(f"""이분법을 이용한 경우 각각 {n11}, {n12}, {n13}번
뉴턴법을 이용한 경우 각각 {n21}, {n22}, {n23}번
시컨트법을 이용한 경우 각각 {n41}, {n42}, {n43}번
의 반복 계산 후 해를 구하였습니다.""")