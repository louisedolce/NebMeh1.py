import math as m
epsilon = 23.43929111 #наклон эклиптики к экватору
k = 0.01720209895  #постоянная гаусса
moo = 132510000000 #гравитационной параметр
file = open('nebo13.txt', 'r')
ok = file.read().split('\n')
omega = float(ok[0])
Bomega = float(ok[1])
naklon_i = float(ok[2])
e = float(ok[3])
a = float(ok[4])
M_0 = float(ok[5])
t_0 = float(ok[6])
T = float(ok[7])
t_i = float(ok[8])
X = float(ok[9])
Y = float(ok[10])
Z = float(ok[11])
file.close()

alpha_1 = m.sin(m.radians(Bomega))*m.sin(m.radians(omega))
alpha_2 = m.sin(m.radians(Bomega))*m.cos(m.radians(omega))
beta_1 = m.cos(m.radians(Bomega))*m.sin(m.radians(omega))
beta_2 = m.cos(m.radians(Bomega))*m.cos(m.radians(omega))
gamma_1 = m.sin(m.radians(naklon_i))*m.sin(m.radians(omega))
gamma_2 = m.sin(m.radians(naklon_i))*m.cos(m.radians(omega))
Px = beta_2 - alpha_1 * m.cos(m.radians(naklon_i))
Py = (alpha_2 + beta_1 * m.cos(m.radians(naklon_i)))*m.cos(m.radians(epsilon)) - gamma_1*m.sin(m.radians(epsilon))
Pz = (alpha_2 + beta_1 * m.cos(m.radians(naklon_i)))*m.sin(m.radians(epsilon)) + gamma_1*m.cos(m.radians(epsilon))
Qx = -beta_1 - alpha_2 * m.cos(m.radians(naklon_i))
Qy = (-alpha_1 + beta_2 * m.cos(m.radians(naklon_i)))*m.cos(m.radians(epsilon)) - gamma_2*m.sin(m.radians(epsilon))
Qz = (-alpha_1 + beta_2 * m.cos(m.radians(naklon_i)))*m.sin(m.radians(epsilon)) + gamma_2*m.cos(m.radians(epsilon))

Sum_P = m.pow(Px,2) + m.pow(Py,2) + m.pow(Pz,2) #проверяем контрольные соотношения
Sum_Q = m.pow(Qx,2) + m.pow(Qy,2) + m.pow(Qz,2)
Mult_PQ = Px*Qx+Py*Qy+Pz*Qz
if abs(Sum_Q-1)<0.00005 and abs(Sum_P-1)<0.00005 and abs(Mult_PQ)<0.00005:
    print("Проверка пройдена") #если не выходит за пределы заданной точности
else:
    print("Проверка P и Q не пройдена")

if T == 0:
    n = (m.pow(moo, 1/2)) / (m.pow(a, 3/2)) #среднее суточное движение рад/сутки
    n = round(n, 10)
    M_i = M_0 + n*(t_i - T) #средняя аномалия на моменты t_i
    M_i = round(M_i, 8)
else:
    n = k / (m.pow(a, 3 / 2))
    n = round(n, 10)
    M_i = M_0 + n*(t_i - T)
    M_i = round(M_i, 8)

eps = 0.00000001 #задаем некую точность для итерационного метода
E = M_i #для нулевого приближения
E_1 = 0
if e < 0.4:
    while abs(E-E_1) > eps:
        E_1 = E
        E = M_i + e*m.sin(E_1)
else:
    while abs(E - E_1) > eps:
        E_1 = E
        E = E_1 + (M_i+e*m.sin(E_1) - E_1)/(1 - e*m.cos(E_1))
E = round(E,8)

rsin=a*m.sqrt(1-m.pow(e,2))*m.sin(E) #находим после определение E
rcos=a*(m.cos(E)-e)
rsin,rcos=round(rsin,8),round(rcos,8)

x=Px*rcos+Qx*rsin #вычисляем гелиоцентрические экваториальные координаты тела
y=Py*rcos+Qy*rsin
z=Pz*rcos+Qz*rsin
x,y,z=round(x,8),round(y,8),round(z,8)

iks=x+X #получаем геоцентрические координаты небесного тела
igrek=y+Y
zet=z+Z
iks,igrek,zet=round(iks,8),round(igrek,8),round(zet,8)

Ro = m.sqrt(m.pow(iks,2)+m.pow(igrek,2)+m.pow(zet,2)) #расстояние с сферических координатах
Ro = round(Ro,8)

delta = m.asin(zet/Ro) #склонение
delta = m.degrees(delta)
deg_d = int(abs(delta)//1)
min_d = int(((abs(delta)-deg_d)*60)//1)
sec_d = round(((abs(delta)-deg_d)*60-min_d)*60,1)

sin_alfa = igrek/(Ro*m.cos(m.radians(delta))) #восхождение
cos_alfa = iks/(Ro*m.cos(m.radians(delta)))

if sin_alfa > 0 and cos_alfa > 0: #в какой четверти находится альфа
    alfa = m.degrees(m.asin(sin_alfa))
elif sin_alfa > 0 and cos_alfa < 0:
    alfa = 180 - m.degrees(m.asin(sin_alfa))
elif sin_alfa < 0 and cos_alfa < 0:
    alfa = 180 + abs(m.degrees(m.asin(sin_alfa)))
else:
    alfa = 360 - abs(m.degrees(m.asin(sin_alfa)))
hour_alfa = int((alfa/15)//1)
min_alfa = int(((alfa/15 - hour_alfa)*60)//1)
sec_alfa = round(((alfa/15 - hour_alfa)*60 - min_alfa)*60,2)

#ответ
rezultat1 = open('itog.txt', 'w')
rezultat1.write(str(Ro)+'а.е.'+'\n')
if delta>0:
    rezultat1.write(str(deg_d)+'°'+'\n') #склонение
    rezultat1.write(str(min_d)+"'"+'\n')
    rezultat1.write(str(sec_d)+'"'+'\n')
else:
    rezultat1.write('-'+str(deg_d)+'°'+'\n')
    rezultat1.write(str(min_d)+"'"+'\n')
    rezultat1.write(str(sec_d)+'"'+'\n')
rezultat1.write(str(hour_alfa)+'h'+'\n') #прямое восхождение
rezultat1.write(str(min_alfa)+'m'+'\n')
rezultat1.write(str(sec_alfa)+'s'+'\n')
rezultat1.close()