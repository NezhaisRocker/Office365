
# coding: utf-8

# In[4]:

from scipy.optimize import minimize
from math import *
from numpy import *
from scipy.special import spence


# In[5]:
#define angle for facet
theta = 11.2/180*pi
#Volume of the island, Slope of island
s = tan(array([0,1.,2,3,2,1,0,-1.,-2,-3,-2,-1,0])*theta)
#define corner points of island
def mya(t,Vol,L):
    a = zeros(12)
    a[1] = t[0]
    a[2] = a[1]+t[1]
    a[3] = a[2]+L[0]
    a[4] = a[3]+t[2]
    a[5] = a[4]+t[3]
    t8   = t[0]+t[3]-t[4]+s[2]/s[1]*(t[1]+t[2]-t[5]-t[6])+s[3]/s[1]*(L[0]-L[1])
    a[6] = a[5]+(Vol-(s[1]*t[0]*(t[0]*0.5+t[1]+L[0]+t[2]+t[3])+\
           s[2]*t[1]*(t[1]*0.5+t[2]+t[3]+L[0])+s[3]*L[0]*(L[0]*0.5+t[2]\
           +t[3])+s[4]*t[2]*(t[2]*0.5+t[3])+s[5]*t[3]**2*0.5)-\
           (s[1]*t8*(t8*0.5+t[6]+L[1]+t[5]+t[4])+\
           s[2]*t[6]*(t[6]*0.5+L[1]+t[5]+t[4])+s[3]*L[1]*(L[1]*0.5+t[5]\
           +t[4])+s[4]*t[5]*(t[5]*0.5+t[4])+s[5]*t[4]**2*0.5))\
           /(s[1]*t[0]+s[2]*t[1]+s[3]*L[0]+s[4]*t[2]+s[5]*t[3])
    a[7] = a[6]+t[4]
    a[8] = a[7]+t[5]
    a[9] = a[8]+L[1]
    a[10] = a[9]+t[6]
    a[11] = a[10]+t8
    return a
#print mya([0.2,0.2,0.2],0.15)
def myh(t,Vol,L):
    a = array(mya(t,Vol,L))    
    h = zeros(11)
    for i in range(1,11):
        h[i] = sum([s[k]*(a[k]-a[k-1]) for k in range(1,i+1)])
    return h
#print h([0.2,0.2,0.2],0.15)
def shape(t,Vol,L):
    a = array(mya(t,Vol,L))
    b = sorted(a)
    return array_equal(a,b)


# In[6]:

def eint1(a,b,c,d):
    if a==b or c==d:
        return 0
    else:
        if b<c or d<a:
            return -1.5*a*c-log(abs(c-a))*a**2*0.5+log(abs(c-a))*c*a\
            -log(abs(c-a))*c**2*0.5 +1.5*d*a+log(abs(d-a))*a**2*0.5\
            -log(abs(d-a))*d*a+log(abs(d-a))*d**2*0.5+1.5*b*c\
            +log(abs(-b+c))*b**2*0.5-log(abs(-b+c))*c*b+log(abs(-b+c))*c**2*0.5\
            -1.5*b*d-log(abs(-b+d))*b**2*0.5+log(abs(-b+d))*b*d\
            -log(abs(-b+d))*d**2*0.5
        elif b==c:
            return 1.5*d*a+log(abs(d-a))*a**2*0.5-log(abs(d-a))*d*a\
            +log(abs(d-a))*d**2*0.5-1.5*a*b+1.5*b**2-log(abs(b-a))*a**2*0.5\
            +log(abs(b-a))*a*b-log(abs(b-a))*b**2*0.5-1.5*b*d\
            -log(abs(-b+d))*b**2*0.5+log(abs(-b+d))*b*d-log(abs(-b+d))*d**2*0.5
        elif d==a:
            return -1.5*a*c+1.5*a**2-log(-c+a)*a**2*0.5+log(-c+a)*a*c\
            -log(-c+a)*c**2*0.5+1.5*b*c-1.5*a*b-log(b-a)*a**2*0.5+log(b-a)*a*b\
            -log(b-a)*b**2*0.5+log(b-c)*b**2*0.5-log(b-c)*b*c+log(b-c)*c**2*0.5
        elif a==c:
            return 3*a*b-1.5*b**2+log(b-a)*a**2-2.*log(b-a)*a*b+log(b-a)*b**2\
            -1.5*a**2


# In[7]:

def dilog(x):
    return spence(x)
def eint2(a,b,c,d):
    c,d = sorted([c,d])
    if a == b:
        return 0
    else:
        if b<c:
            return (-4*a + 4*b - 2*(a - c)*log(-a + c)*(-1 + log(-a + d)) + \
            2*a*log(-a + d)-2*d*log(-a+d)-c*log(-a+d)**2+d*log(-a + d)**2+\
            2*(b - c)*log(-b + c)*(-1 + log(-b + d)) - 2*b*log(-b + d)+\
            2*d*log(-b + d) + c*log(-b + d)**2 - d*log(-b + d)**2)/2.\
            +(-c + d)*dilog((a-c)/(a - d)) + (c - d)*dilog((b-c)/(b - d))
        elif b==c and c<d:
            return (-12*a + 12*b + b*pi**2 - d*pi**2\
            -6*(a - b)*log(-a + b)*(-1 + log(-a + d)) + 6*a*log(-a + d) -\
            6*d*log(-a + d) - 3*b*log(-a + d)**2 + 3*d*log(-a + d)**2 - \
            6*b*log(-b + d) + 6*d*log(-b + d) + 3*b*log(-b + d)**2 - \
            3*d*log(-b + d)**2 - 6*(b - d)*dilog((a-b)/(a - d)))/6.
        elif (b==c and c==d) or (a==c and c==d):
            return -(a - b)*(2 + (-2 + log(-a + b))*log(-a + b))
        elif a==c and b==d:
            return ((a - b)*(-12 + pi**2 - 6*(-2 + log(-a + b))*log(-a + b)))/6.
        elif a==c and b<d:
            return -d*log(d-a)+log(-b+d)*log(b-a)*b-log(-b + d)*log(b - a)*d\
            + a*log(b - a) + log(d - a)*a + dilog((-d + b)/(-d + a))*a\
            - dilog((-d + b)/(-d + a))*d - 2.*a + 2.*b -log(d - a)*a*log(b - a)\
            + log(d - a)*d*log(b - a) - log(-b + d)*b + d*log(-b + d)\
            -log(b - a)* b
        elif c<a and a==d:
            return (12*b + c*pi**2 - a*(12 + pi**2) +3*((-a + c)\
            *(-2 + log(a - c))*log(a - c) -  2*(a - b)*log(-a + b)\
            *(-1 + log(b - c)) + 2*(-b + c)*log(b - c) +\
            (a - c)*log(b - c)**2) + 6*(a - c)*dilog((b-a)/(b-c)))/6.
        elif c<a and b==d:
            return -log(b - a)*log(a - c)*a + log(b - a)*log(a - c)*b \
            -log(b - c)*b*log(a - c)+log(b - c)**2*b+log(b - c)*c*log(a - c)\
            - log(b - c)**2*c+a*log(b - a)-log(b - a)*b+log(a - c)*a \
            + dilog((b-a) / (b - c))*b - dilog((b-a) / (b - c))*c \
            - c*log(a - c)-2.*a+2.*b-log(b - c)*b-pi**2*b/6. \
            +pi ** 2 * c / 6. + c *log(b - c)
        elif c<a and b<d:
            return log(a-c)*(a-c)+log(d-a)*(-d + a)-2.*a+2.*b+log(d-c)*log(a-c)\
            * (-d + c) - dilog((-d + a) / (-d + c))*(c-d)+log(d-a) * \
            log(a-c)*(d-a) + (c - b) * log(b - c) + dilog((-d + b)/(-d + c))\
            *(-d + c)+log(d - b)*(d - b)+log(d - c)*(d - c)*log(b - c) + \
            log(d - b) * log(b - c) * (-d + b)
        elif d<a and c<a:
            return (-4*a + 4*b - 2*(a - c)*log(a - c)*(-1 + log(a - d)) +\
            2*(a - d)*log(a - d) + (-c + d)*log(a - d)**2 +\
            2*(b - c)*log(b - c)*(-1 + log(b - d)) + 2*(-b + d)*log(b - d) + \
            (c - d)*log(b - d)**2 - 2*(c - d)*(dilog((a-c)/(a - d)) \
            - dilog((b-c)/(b - d))))/2.


# In[8]:

def eint3(a,b,c,d):
    c,d=sorted([c,d])    
    if a==b:
        return 0
    else:
        if b<c and c!=d:
            return ((-a + b)*(a + b + 3*(c + d)) - b**2*log(-b + c) +(a - d)\
            *(a + 2*c + d)*log(-a + d) + (-b + d)*(b + 2*c + d)*log(-b + d) +\
            log(-a + c)*((a - c)*(a + c + 2*d) + 2*(-a**2 + d**2)*log(-a + d) +\
            2*(c - d)*(c + d)*log(-c + d)) +log(-b + c)*(c**2 - 2*b*d + 2*c*d +\
            2*(b - d)*(b + d)*log(-b + d) +2*(-c**2 + d**2)*log(-c + d)) -\
             2*(c - d)*(c + d)*(dilog((a-d)/(c - d)) -dilog((b-d)/(c - d))))/4.
        if b<c and c==d:
            return log(c-a)**2*c**2/2.+log(c-a)*a*c-1.5*log(c-a)*c**2-1.5*c*a\
            -log(c-a)**2*a**2/2.+log(c-a)*a**2/2.-a**2/4.-log(c-b)**2*c**2/2.\
            -log(c-b)*b*c+1.5*log(c-b)*c**2+1.5*c*b+log(c-b)**2*b**2/2.\
            -log(c-b)*b**2/2.+b**2/4.
        if b==c and c!=d:
            return (-((a - b)*(a + 4*b+3*d))+(a-d)*(a + 2*b + d)*log(-a + d)+\
            (-b + d)*(3*b + d)*log(-b + d) +log(-a + b)*((a - b)*(a + b + 2*d)\
             + 2*(-a**2 + d**2)*log(-a + d) +2*(b - d)*(b + d)*log(-b + d)) +\
             2*(-b**2 + d**2)*dilog((a-d)/(b - d)))/4.
        if b==c and c==d:
            return log(b-a)**2*b**2/2.+log(b-a)*a*b-1.5*b**2*log(b-a)-\
            1.5*b*a+1.75*b**2-log(b-a)**2*a**2/2.+log(b-a)*a**2/2.-a**2/4.
        if a==c and b==d:
            return ((a-b)*(a+b)*(-12+pi**2-6*(-2+log(-a+b))*log(-a+b)))/12.
        if a==c and b<d:
            return (-3*(a - b)*(4*a + b + 3*d) + (a - d)*(a + d)*pi**2 +\
            6*(-a**2 + d**2)*log(-a + d)**2-3*(b-d)*(2*a+b+d)*log(-b+d)-\
            3*(a - b)*log(-a + b)*(-a - b - 2*d + 2*(a + b)*log(-b + d)) +  \
            3*(a - d)*log(-a + d)*(3*a + d + 2*(a + d)*log(-b + d)) +    \
            6*(-a**2 + d**2)*dilog((a-b)/(a - d)))/12.
        if a==c and c==d:
            return -((a-b)*(7*a+b-2*(3*a+b)*log(-a+b)+2*(a+b)*log(-a+b)**2))/4.
        if a>c and a==d:
            return (-((a - b)*(4*a + b + 3*c)) + (a - c)*(3*a + c)*log(a-c)+ \
            (-b + c)*(2*a+b+c)*log(b-c) + log(-a + b)*((a - b)*(a + b + 2*c)\
            + 2*(-a**2 + c**2)*log(a - c) +2*(b - c)*(b + c)*log(b - c)) +\
            2*(a - c)*(a + c)*dilog((b-c)/(a - c)))/4. 
        if a>c and b==d:
            return (-((a - b)*(a + 4*b + 3*c))+(a-c)*(a+2*b+c)*log(a-c) +\
             (-b + c)*(3*b+c)*log(b-c)+log(-a + b)*((a - b)*(a + b + 2*c)\
             + 2*(-a**2 + c**2)*log(a - c) +2*(b - c)*(b + c)*log(b - c)) + \
             2*(-b**2 + c**2)*dilog((a-c)/(b - c)))/4.
        if a>c and b<d:
            return (-((a - b)*(a + b + 3*(c + d))) + a**2*log(-a + d) + \
            2*a*c*log(-a + d) - 2*c*d*log(-a + d) - d**2*log(-a + d) + \
            log(a - c)*((a - c)*(a + c + 2*d) + 2*(-a**2 + c**2)*log(-a+d))-\
            b**2*log(-b + d) - 2*b*c*log(-b + d) + 2*c*d*log(-b + d) + \
            d**2*log(-b + d) + (b - c)*log(b - c)*  \
            (-b - c - 2*d + 2*(b + c)*log(-b + d)) - \
            2*(c - d)*(c + d)*(log(-a + d) - log(-b + d))*log(-c + d) +   \
            2*(c - d)*(c + d)*(dilog((c-a)/(c - d))-dilog((c-b)/(c - d))))/4.
        if d<a and c!=d:
            return (-((a-b)*(a+b+3*(c+d)))-b**2*log(b-d)+(a-c)*(a+2*d+c)\
            *log(a-c)+(-b+c)*(b+2*d+c)*log(b-c)+log(a-d)*((a-d)*(a+d+2*c)\
            +2*(-a**2+c**2)*log(a-c)+2*(d-c)*(c+d)*log(d-c))+\
            log(b-d)*(d**2-2*b*c+2*c*d+2*(b-c)*(b+c)*log(b-c)\
            +2*(-d**2+c**2)*log(d-c))-2*(d-c)*(c+d)*(dilog((a-c)/(d-c)) \
            -dilog((b-c)/(d-c))))/4.
        if d<a and c==d:
            return (-((a-b)*(a+b+6*c))+2*(a-c)*(a+3*c)*log(a-c)+2*(-a**2+c**2)\
            *log(a-c)**2+2*(b-c)*log(b-c)*(-b-3*c+(b+c)*log(b-c)))/4.


# In[18]:

def Eterm1(t,Vol,L):
    a = array(mya(t,Vol,L))    
    return 2./pi*sum([sum([s[i+1]*s[j+1]*eint1(a[i],a[i+1],a[j],a[j+1])\
    for j in range(11)]) for i in range(11)])
#print Eterm1([0.2,0.2,0.2],0.15)


# In[19]:

def Eterm2(t,Vol,L):
    a = array(mya(t,Vol,L))
    h = array(myh(t,Vol,L))
    M2 = empty((11,12,12),dtype=float) 
    for k in range(11):
        for j in range(12):
            for i in range(12):
                M2[k][j][i]=eint2(a[k],a[k+1],a[i],a[j])
    M3 = empty((11,12,12),dtype=float) 
    for k in range(11):
        for j in range(12):
            for i in range(12):
                M3[k][j][i]=eint3(a[k],a[k+1],a[i],a[j])
    return 4./pi**2*sum([sum([sum([\
           s[i+1]*s[j+1]*(h[k]-s[k+1]*a[k])*\
           (M2[k][i+1][j+1]-M2[k][i][j+1]-M2[k][i+1][j]+M2[k][i][j])\
           +s[i+1]*s[j+1]*s[k+1]*\
           (M3[k][i+1][j+1]-M3[k][i][j+1]-M3[k][i+1][j]+M3[k][i][j])
           for i in range(11)]) for j in range(11)]) for k in range(11)])
#print Eterm2([0.2,0.2,0.2],0.15)


# In[21]:

def Eterm3(t,Vol,L):
    a = array(mya(t,Vol,L))
    h = array(myh(t,Vol,L))
    return -4.*sum([s[i+1]**2*h[i]*(a[i+1]-a[i])+0.5*s[i+1]**3*(a[i+1]-a[i])**2\
           for i in range(11)])
#print Eterm3([0.2,0.2,0.2],0.15)


# In[17]:

def Energy_elastic(t,Vol,L):
    if shape(t,Vol,L)==True:
        return Eterm1(t,Vol,L)+Eterm2(t,Vol,L)+Eterm3(t,Vol,L)
    else:
        return 100
#print Energy_elastic([0.2,0.2,0.2],0.15)


# In[16]:

def Energy_surface(t,Vol,L):
    a = array(mya(t,Vol,L))   
    if shape(t,Vol,L)==True:
        return sum([(sqrt(1+s[i+1]**2)-1)*(a[i+1]-a[i])\
        for i in range(11)])
    else:
        return 100
#print Energy_surface([0.2,0.2,0.2],0.15)


# In[22]:

def Energy_total(t,Vol,L):
    return Energy_elastic(t,Vol,L)+Energy_surface(t,Vol,L)
#print Energy_total([0.2,0.2,0.2],0.15)
#x1=linspace(0.05,0.4,120)
#t9= time()
#Ex = array([Energy_total([0.1,y,0.3]) for y in x1])
#plot(x1,Ex)
#t10=time()
#print t10-t9


# In[23]:

def MinEnergy(x0,Vol,L):
    res=minimize(Energy_total,x0,args=(Vol,L,),method='L-BFGS-B',\
    options={'ftol':1e-9,'gtol':1e-6})
    return res
#t0 = time()
#print MinEnergy([ 0.06764238,  0.12133122,  0.31559877,  0.35082444,  0.35080851,
#        0.31560051,  0.12133631],0.5,[0.1,0.105])
#t1 = time()
#print t1-t0

# In[24 ]:


















