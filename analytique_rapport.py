import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
##### without fins #####
fig, ax = plt.subplots()
def dTdt(t, T):
    Pabs = 5
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    return (Pabs-h*A*(T-Tamb))/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 1000, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "5 W")



##### 10W #####
def dTdt(t, T):
    Pabs = 10
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    return (Pabs-h*A*(T-Tamb))/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 1000, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "10 W")


##### 15W #####

def dTdt(t, T):
    Pabs = 15
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    return (Pabs-h*A*(T-Tamb))/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 1000, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "15 W")


##### 20W #####

def dTdt(t, T):
    Pabs = 20
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    return (Pabs-h*A*(T-Tamb))/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 1000, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "20 W")

plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
#plt.show()


##### with fins ######
##### 5W #####
fig, ax = plt.subplots()
def dTdt(t, T):
    Pabs = 5
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    #calcul de P fins
    N = 16
    w = 0.05
    L = 0.02
    epa = 0.002
    Lc = L + epa/2
    Af = 2*w*Lc
    Ab = w**2 - N*w*epa
    At = N*Af+Ab
    Ac = w*epa
    peri = 2*epa+2*w
    kalu = 177
    T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
    m = np.sqrt((h*peri)/(kalu*Ac))
    nuf = np.tanh(m*Lc)/(m*Lc)
    nut = 1-(((N*Af)/At)*(1-nuf))
    Pfins = (T_0-Tamb)*nut*h*At
    print("hello")
    print(Pfins)
    return (Pabs-h*A*(T-Tamb)-Pfins)/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 180, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "5 W")


##### 10W #####
def dTdt(t, T):
    Pabs = 10
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    #calcul de P fins
    N = 16
    w = 0.05
    L = 0.02
    epa = 0.002
    Lc = L + epa/2
    Af = 2*w*Lc
    Ab = w**2 - N*w*epa
    At = N*Af+Ab
    Ac = w*epa
    peri = 2*epa+2*w
    kalu = 177
    T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
    m = np.sqrt((h*peri)/(kalu*Ac))
    nuf = np.tanh(m*Lc)/(m*Lc)
    nut = 1-(((N*Af)/At)*(1-nuf))
    Pfins = (T_0-Tamb)*nut*h*At
    print("hello")
    print(Pfins)
    return (Pabs-h*A*(T-Tamb)-Pfins)/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 180, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "10 W")


##### 15W #####
def dTdt(t, T):
    Pabs = 15
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    #calcul de P fins
    N = 16
    w = 0.05
    L = 0.02
    epa = 0.002
    Lc = L + epa/2
    Af = 2*w*Lc
    Ab = w**2 - N*w*epa
    At = N*Af+Ab
    Ac = w*epa
    peri = 2*epa+2*w
    kalu = 177
    T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
    m = np.sqrt((h*peri)/(kalu*Ac))
    nuf = np.tanh(m*Lc)/(m*Lc)
    nut = 1-(((N*Af)/At)*(1-nuf))
    Pfins = (T_0-Tamb)*nut*h*At
    print("hello")
    print(Pfins)
    return (Pabs-h*A*(T-Tamb)-Pfins)/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 180, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "15 W")



##### 20W #####
def dTdt(t, T):
    Pabs = 20
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    #calcul de P fins
    N = 16
    w = 0.05
    L = 0.02
    epa = 0.002
    Lc = L + epa/2
    Af = 2*w*Lc
    Ab = w**2 - N*w*epa
    At = N*Af+Ab
    Ac = w*epa
    peri = 2*epa+2*w
    kalu = 177
    T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
    m = np.sqrt((h*peri)/(kalu*Ac))
    nuf = np.tanh(m*Lc)/(m*Lc)
    nut = 1-(((N*Af)/At)*(1-nuf))
    Pfins = (T_0-Tamb)*nut*h*At
    print("hello")
    print(Pfins)
    return (Pabs-h*A*(T-Tamb)-Pfins)/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 180, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "20 W")

plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
#plt.show()


##### test differents nbr de fins
##### 20W #####

##### 4 fins #####
fig, ax = plt.subplots()
def dTdt(t, T):
    Pabs = 20
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    #calcul de P fins
    N = 4
    w = 0.05
    L = 0.02
    epa = 0.002
    Lc = L + epa/2
    Af = 2*w*Lc
    Ab = w**2 - N*w*epa
    At = N*Af+Ab
    Ac = w*epa
    peri = 2*epa+2*w
    kalu = 177
    T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
    m = np.sqrt((h*peri)/(kalu*Ac))
    nuf = np.tanh(m*Lc)/(m*Lc)
    nut = 1-(((N*Af)/At)*(1-nuf))
    Pfins = (T_0-Tamb)*nut*h*At
    print("hello")
    print(Pfins)
    return (Pabs-h*A*(T-Tamb)-Pfins)/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 180, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "4 fins")

##### 8 fins #####
def dTdt(t, T):
    Pabs = 20
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    #calcul de P fins
    N = 8
    w = 0.05
    L = 0.02
    epa = 0.002
    Lc = L + epa/2
    Af = 2*w*Lc
    Ab = w**2 - N*w*epa
    At = N*Af+Ab
    Ac = w*epa
    peri = 2*epa+2*w
    kalu = 177
    T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
    m = np.sqrt((h*peri)/(kalu*Ac))
    nuf = np.tanh(m*Lc)/(m*Lc)
    nut = 1-(((N*Af)/At)*(1-nuf))
    Pfins = (T_0-Tamb)*nut*h*At
    print("hello")
    print(Pfins)
    return (Pabs-h*A*(T-Tamb)-Pfins)/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 180, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "8 fins")

##### 12fins ######
def dTdt(t, T):
    Pabs = 20
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    #calcul de P fins
    N = 12
    w = 0.05
    L = 0.02
    epa = 0.002
    Lc = L + epa/2
    Af = 2*w*Lc
    Ab = w**2 - N*w*epa
    At = N*Af+Ab
    Ac = w*epa
    peri = 2*epa+2*w
    kalu = 177
    T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
    m = np.sqrt((h*peri)/(kalu*Ac))
    nuf = np.tanh(m*Lc)/(m*Lc)
    nut = 1-(((N*Af)/At)*(1-nuf))
    Pfins = (T_0-Tamb)*nut*h*At
    print("hello")
    print(Pfins)
    return (Pabs-h*A*(T-Tamb)-Pfins)/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 180, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "12 fins")


##### 16 fins #####
def dTdt(t, T):
    Pabs = 20
    R = 0.025
    H = 0.00127
    Tamb = 295
    rho = 2770
    V = np.pi*R**2*H
    Cp = 875
    A = 2*np.pi*R*H + 2*np.pi*R**2
    #calcul de h
    #prop a 300K
    g = 9.8
    Ts = T
    beta = 1/Ts
    Tinf = Tamb
    D = 2*R
    nu = 15.89*10**-6
    Pr = 0.707
    k = 26.3*10**-3
    Gr = (g*beta*(Ts-Tinf)*D**3)/nu**2
    Ra = Gr*Pr
    Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h = (Nu*k)/D
    #calcul de P fins
    N = 16
    w = 0.05
    L = 0.02
    epa = 0.002
    Lc = L + epa/2
    Af = 2*w*Lc
    Ab = w**2 - N*w*epa
    At = N*Af+Ab
    Ac = w*epa
    peri = 2*epa+2*w
    kalu = 177
    T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
    m = np.sqrt((h*peri)/(kalu*Ac))
    nuf = np.tanh(m*Lc)/(m*Lc)
    nut = 1-(((N*Af)/At)*(1-nuf))
    Pfins = (T_0-Tamb)*nut*h*At
    print("hello")
    print(Pfins)
    return (Pabs-h*A*(T-Tamb)-Pfins)/(rho*V*Cp)

T0 = 295

t = np.linspace(0, 180, 100)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
print(sol_m1.T[0])

ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "16 fins")

plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
plt.show()