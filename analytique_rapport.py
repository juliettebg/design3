import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
##### without fins #####
fig, ax = plt.subplots()
for PABS in [5, 10, 15, 20]:
    def dTdt(t, T):
        Pabs = PABS
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

    ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "{} W".format(PABS))
plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
plt.show()


fig, ax = plt.subplots()
###### nbr of fins ######
for n in [4, 8, 12, 16]:
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
        N = n
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
        return (Pabs-h*A*(T-Tamb)-(np.pi/4)*Pfins)/(rho*V*Cp)

    T0 = 295

    t = np.linspace(0, 180, 100)
    sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
    print(sol_m1.T[0])

    ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "{} fins".format(n))

 
plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
plt.show()



fig, ax = plt.subplots()
###### fin length ######
for l in [0.01, 0.02, 0.03, 0.04, 0.05]:
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
        L = l
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
        return (Pabs-h*A*(T-Tamb)-(np.pi/4)*Pfins)/(rho*V*Cp)

    T0 = 295

    t = np.linspace(0, 180, 100)
    sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)

    ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "{} mm".format(l*1000))

 
plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
plt.show()


##### thickness aluminium plate #####
fig, ax = plt.subplots()
###### with fins ######
for thickness in [0.00076, 0.00127, 0.00159, 0.00317]:
    def dTdt(t, T):
        Pabs = 20
        R = 0.025
        H = thickness
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
        return (Pabs-h*A*(T-Tamb)-(np.pi/4)*Pfins)/(rho*V*Cp)

    T0 = 295

    t = np.linspace(0, 180, 100)
    sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
    print(sol_m1.T[0])

    ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "{} mm".format(thickness))

 
plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
plt.show()


##### types of aluminium #####
fig, ax = plt.subplots()
###### with fins ######
#rho, cp, k
#380 https://www.matweb.com/search/datasheet_print.aspx?matguid=7441a0886ba142cc82eb7af5d4eece6d
for alu in [("pure", 2702, 903, 237), ("Alloy 2024-T6", 2770, 875, 177), ("Alloy 195", 2790, 883, 168), ("Alloy A380", 2740, 963, 109), ("Alloy 6061-T6", 2700, 896, 170)]:
    def dTdt(t, T):
        Pabs = 20
        R = 0.025 
        H = thickness
        Tamb = 295
        rho = alu[1]
        V = np.pi*R**2*H
        Cp = alu[2]
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
        kalu = alu[3]
        T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
        m = np.sqrt((h*peri)/(kalu*Ac))
        nuf = np.tanh(m*Lc)/(m*Lc)
        nut = 1-(((N*Af)/At)*(1-nuf))
        Pfins = (T_0-Tamb)*nut*h*At
        print("hello")
        print(Pfins)
        return (Pabs-h*A*(T-Tamb)-(np.pi/4)*Pfins)/(rho*V*Cp)

    T0 = 295

    t = np.linspace(0, 180, 100)
    sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
    print(sol_m1.T[0])

    ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "{}".format(alu[0]))

 
plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
plt.show()


fig, ax = plt.subplots()
###### with fins ######
for PABS in [5, 10, 15, 20]:
    def dTdt(t, T):
        Pabs = PABS
        R = 0.025
        H = 0.00127
        Tamb = 295
        rho = 2740
        V = np.pi*R**2*H
        Cp = 963
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
        L = 0.03
        epa = 0.002
        Lc = L + epa/2
        Af = 2*w*Lc
        Ab = w**2 - N*w*epa
        At = N*Af+Ab
        Ac = w*epa
        peri = 2*epa+2*w
        kalu = 109
        T_0 = T - (Pabs*H)/(kalu*np.pi*R**2)
        m = np.sqrt((h*peri)/(kalu*Ac))
        nuf = np.tanh(m*Lc)/(m*Lc)
        nut = 1-(((N*Af)/At)*(1-nuf))
        Pfins = (T_0-Tamb)*nut*h*At
        print("hello")
        print(Pfins)
        return (Pabs-h*A*(T-Tamb)-(np.pi/4)*Pfins)/(rho*V*Cp)

    T0 = 295

    t = np.linspace(0, 180, 100)
    sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
    print(sol_m1.T[0])

    ax.plot(t, [x-273 for x in sol_m1.T[0]], label = "{} W".format(PABS))

 
plt.legend()
ax.set_xlabel("Temps [s]")
ax.set_ylabel("Temprature [Celcius]")
plt.show()

###### radial #####
K = 109
rho = 2740
Cp = 963
max_iter_time = 181
alpha = K/(rho*Cp)
delta_x = 1
delta_t = (delta_x ** 2)/(4 * alpha)
#delta_t =1 
gamma = (alpha * delta_t) / (delta_x ** 2)
size = 500
sizey = 500
u = np.empty((max_iter_time, size, sizey))
u_initial = 295-273

u.fill(u_initial)
t = np.arange(0, max_iter_time-1, (max_iter_time-1)/180)
print(t)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
#print(sol_m1)
#print(len(sol_m1))
#print(sol_m1[6051])
def calculate(u):
    for k in range(0, len(sol_m1)-1, 1):
        print(k)
        print(sol_m1[k])
        for i in range(np.shape(u)[1]):
            for j in range(np.shape(u)[2]):
                if (i-size//2)**2 + (j-size//2)**2 <= (25)**2:
                    u[k][i][j] = sol_m1[k][0]-273
                    #print("ici!")
                    #print(sol_m1[k][0]-273)
                    
        for i in range(1, size-1, delta_x):
            for j in range(1, sizey-1, delta_x):
                u[k + 1, i, j] = gamma * (u[k][i+1][j] + u[k][i-1][j] + u[k][i][j+1] + u[k][i][j-1] - 4*u[k][i][j]) + u[k][i][j]

    return u

def plotheatmap(u_k, k):
    # Clear the current plot figure
    plt.clf()

    plt.title(f"Temperature at t = {k} s")
    plt.xlabel("x [mm*10]")
    plt.ylabel("y [mm*10]")

    # This is to plot u_k (u at time-step k)
    plt.pcolormesh(u_k, cmap=plt.cm.jet, vmin=295-273, vmax=91)
    plt.colorbar()

    return plt

# Do the calculation here
u = calculate(u)


def animate(k):
    plotheatmap(u[k], k)
animate
for k in [60]:
    plotheatmap(u[k], k)
    plt.show()

##### axial #####
K = 109
rho = 2740
Cp = 963
max_iter_time = 181
alpha = K/(rho*Cp)
delta_x = 1
delta_t = (delta_x ** 2)/(4 * alpha)
#delta_t =1 
gamma = (alpha * delta_t) / (delta_x ** 2)
size = 500
sizey = 30
u = np.empty((max_iter_time, size, sizey))
u_initial = 295-273

u.fill(u_initial)
t = np.arange(0, max_iter_time-1, (max_iter_time-1)/180)
print(t)
sol_m1 = odeint(dTdt, y0 = T0, t=t, tfirst= True)
#print(sol_m1)
#print(len(sol_m1))
#print(sol_m1[6051])
def calculate(u):
    for k in range(0, len(sol_m1)-1, 1):
        print(k)
        print(sol_m1[k])
        u[k, size//2-30:size//2+20, :1] = sol_m1[k][0]-273
                    
        for i in range(1, size-1, delta_x):
            for j in range(1, sizey-1, delta_x):
                u[k + 1, i, j] = gamma * (u[k][i+1][j] + u[k][i-1][j] + u[k][i][j+1] + u[k][i][j-1] - 4*u[k][i][j]) + u[k][i][j]

    return u

def plotheatmap(u_k, k):
    # Clear the current plot figure
    plt.clf()

    plt.title(f"Temperature at t = {k} s")
    plt.xlabel("z [mm*10]")
    plt.ylabel("y [mm*10]")

    # This is to plot u_k (u at time-step k)
    plt.pcolormesh(u_k, cmap=plt.cm.jet, vmin=295-273, vmax=91)
    plt.colorbar()

    return plt

# Do the calculation here
u = calculate(u)


def animate(k):
    plotheatmap(u[k], k)
animate
for k in [60]:
    plotheatmap(u[k], k)
    plt.show()