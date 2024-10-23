import numpy as np 
import math  
from odlib.odlib import *


k = 0.01720209895

def angularMomentum(position, velocity):
    h = np.cross(position, velocity)
    return h

def SemiMajAx(position, velocity):
    magPos = np.linalg.norm(position)
    magVel = np.linalg.norm(velocity)
    a = 1/((2/magPos)-(magVel**2))
    return a 

def eccentricity(h, a):
    # e = np.linalg.norm(np.cross(((-1)*Angular_Momentum), velocity)-(position/np.linalg.norm(position)))
    e = math.sqrt(1 - (((np.linalg.norm(h))**2)/a))
    return e 

def inclination(h):
    i = math.atan((math.sqrt((h[0])**2 + (h[1])**2))/(h[2]))
    incline = math.degrees(i) 
    return incline

def LongAscending(h, incline):
    MagH = np.linalg.norm(h)
    i = np.radians(incline)
    s = (h[0])/((MagH)*(math.sin(i)))
    l = (-h[1])/(MagH*(math.sin(i)))
    Big_omega = math.asin((h[0])/((MagH)*(math.sin(i)))) 
    if s <= 0 and l <= 0:  
        Big_omega = (-Big_omega) + math.pi
        # test = math.acos((h[0])/((MagH)*(math.sin(i))))
        # if test<= math.pi:
        #     Big_omega = ((-1)*Big_omega + math.pi) 
    elif s > 0 and l < 0:
        # test = math.acos((h[0])/((MagH)*(math.sin(i))))
        # if test>math.pi:
        Big_omega = math.pi - Big_omega 
    elif s < 0 and l > 0:
        Big_omega = 2*math.pi - Big_omega 
    else: 
        Big_omega = math.asin((h[0])/((MagH)*(math.sin(i))))
            

    Big_Omega = math.degrees(Big_omega)
    return Big_Omega 
    
def Angle_U(position, incline):
    i = math.radians(incline)
    u = math.asin((position[2])/((np.linalg.norm(position))*(math.sin(i))))
    return math.degrees(u)

def truAnomaly(a, e, h, position, velocity):
    p = (np.dot(position, velocity))/(np.linalg.norm(position))
    MagH = np.linalg.norm(h) 
    E = e
    s = (((a*(1-(E**2))))/(E*MagH))*p 
    Nu = math.asin(((a*(1-(E**2)))*(p))/(E*MagH))
    l = (1/E)*(((a*(1-(E**2)))/(np.linalg.norm(position))) - 1)
    if l <= 0 :
        if s < 0:   
            Nu = 2*math.pi - Nu
        # test = math.acos((h[0])/((MagH)*(math.sin(i))))
        # if test<= math.pi:
        #     Big_omega = ((-1)*Big_omega + math.pi) 
        else: 
            Nu = Nu 
        # test = math.acos((h[0])/((MagH)*(math.sin(i))))
        # if test>math.pi:
    else:
        if s < 0: 
            Nu = 2*math.pi + Nu 
        else:
            Nu = Nu
    
    return math.degrees(Nu)

def Agr_Perhelion(u, Nu):
    small_Omega = u - Nu
    if small_Omega<0:
        SmallO = 360 - abs(small_Omega) 
    else:
        SmallO = small_Omega
    return SmallO 

def EccentricAnomaly(e, position, a, Nu):
    # E_rads = math.acos((1/e)*((1-((np.linalg.norm(position))/(a)))))  
    # l = ((1/e)*((1-((np.linalg.norm(position))/(a)))))
    # if l < 0:
    #     if 0 < Nu < 180:
    #         E_rads = E_rads 
    #     else:
    #         E_rads = 2*math.pi - E_rads 
    # else: 
    #     if 180 < Nu < 360:
    #         E_rads = 2*math.pi - E_rads
    #     else: 
    #         E_rads = E_rads
    if (Nu <= 180):
        E = math.acos((1/e)*(1 - ((np.linalg.norm(position))/a))) # Eq 52
    else:
        E = 2*math.pi - math.acos((1/e)*(1 - ((np.linalg.norm(position))/a))) # Eq 52
    return E

gauss_conversion = (2*math.pi)/365.256893

def TimeOfPerihelion(a, M, J_o):
    k = 0.01720209895
    T = J_o - M*(math.sqrt((a**3)/(k**2)))
    return T 

def MeanAnomaly(J_o, T, a):
    k = 0.01720209894 
    # n = (k*(math.sqrt(1/(a**3))))
    P = math.sqrt((a**3)*365.2568983)
    n = ((2*math.pi)/P)
    return np.degrees(n*((J_o - T)))

def MeanAnomaly2(position, a, e, Nu):
    E_0 = EccentricAnomaly(e, position, a, Nu)
    return np.degrees(E_0) - e*np.degrees(math.sin(E_0))

def newton_raphson(E_o, tol, max_iter, e, M):
    E = E_o
    repeat = 0
    while repeat < max_iter:
        fE = E - (e * math.sin(E)) - math.radians(M)
        dfE = 1 - (e * math.cos(E))
        if dfE == 0:
            print("Derivative is zero. No solution found.")
        if abs(fE/dfE) < tol:
            return E
        else:
            E -= (fE / dfE)
        repeat += 1
    return E  

def getDate_obs1(fileName):
    data = np.array(open(fileName).readlines())
    date = np.array(data[0].split())
    years = date[0].astype(float)
    months = date[1].astype(float)
    days = date[2].astype(float)
    hours = date[3].astype(float)
    minutes = date[4].astype(float)
    seconds = date[5].astype(float)
    return np.array([years, months, days, hours, minutes, seconds])  

def getDate_obs2(fileName):
    data = np.array(open(fileName).readlines())
    date = np.array(data[4].split())
    years = date[0].astype(float)
    months = date[1].astype(float)
    days = date[2].astype(float)
    hours = date[3].astype(float)
    minutes = date[4].astype(float)
    seconds = date[5].astype(float)
    return np.array([years, months, days, hours, minutes, seconds]) 

def getDate_obs3(fileName):
    data = np.array(open(fileName).readlines())
    date = np.array(data[8].split())
    years = date[0]
    months = date[1]
    days = date[2]
    hours = date[3]
    minutes = date[4]
    seconds = date[5]
    return np.array([years, months, days, hours, minutes, seconds]) 

def RA_DEC1(fileName):

    data = np.array(open(fileName).readlines())
    c = np.array((data[2]).split())
    # d = np.array((c[0]).split(":"))
    e = ((c[0]).astype(float))
    f = c[1].astype(float)
    #to get DEC in decimal degrees 
    # c1 = np.array((data[2]).split())
    # d = np.array(((c1[1]).strip("-")).split(":"))
    # if float(d[0]) < 0:
    #     f = (float(d[0]) - (float(d[1]))/60 - (float(d[2]))/3600)
    # else:
    #     f = (float(d[0]) + (float(d[1]))/60 + (float(d[2]))/3600)
    return [e, f]

def RA_DEC2(fileName):
    data = np.array(open(fileName).readlines())
    c = np.array((data[6]).split())
    # d = np.array((c[0]).split(":"))
    e = c[0].astype(float)
    f = c[1].astype(float)
    # data = np.array(open(fileName).readlines())
    #     #conditional if RA and DEC already in degrees 
    
    # c = np.array((data[6]).split())
    # d = np.array((c[0]).split(":"))
    # e = ((float(d[0]))*15 + 15*(float(d[1])/60) + 15*(float(d[2])/3600))
        
    #     #to get DEC in decimal degrees    
    # c1 = np.array((data[6]).split())
    # d = np.array(((c1[1]).strip("+")).split(":"))
    # if float(d[0]) < 0:
    #     f = (float(d[0]) - (float(d[1]))/60 - (float(d[2]))/3600)
    # else:
    #     f = (float(d[0]) + (float(d[1]))/60 + (float(d[2]))/3600)
    return [e, f]

def RA_DEC3(fileName):

    data = np.array(open(fileName).readlines())
    c = np.array((data[10]).split())
    # d = np.array((c[0]).split(":"))
    e = c[0].astype(float)
    f = c[1].astype(float)
    # data = np.array(open(fileName).readlines())

    #     #conditional if RA and DEC already in degrees 
    
    # c = np.array((data[10]).split())
    # d = np.array((c[0]).split(":"))
    # e = ((float(d[0]))*15 + 15*(float(d[1])/60) + 15*(float(d[2])/3600))
        

    #     #to get DEC in decimal degrees 
            
    # c1 = np.array((data[10]).split())
    # d = np.array(((c1[1]).strip("+")).split(":"))
    # if float(d[0]) < 0:
    #     f = (float(d[0]) - (float(d[1]))/60 - (float(d[2]))/3600)
    # else:
    #     f = (float(d[0]) + (float(d[1]))/60 + (float(d[2]))/3600)
    return [e, f]

def ligthtimeCorrection1(J_1, J_2, J_3, rhomags):
    cAU = 173.144643267
    #k = 0.0172020989484
    t_1 = J_1 - (rhomags[0]/cAU)
    t_2 = J_2 - (rhomags[1]/cAU)
    t_3 = J_3 - (rhomags[2]/cAU)
    return [t_1, t_2, t_3]

def ligthtimeCorrection(J_1, J_2, J_3, rhomags):
    cAU = 173.144643267
    k = 0.0172020989484
    t_1 = J_1 - (rhomags[0]/cAU)
    t_2 = J_2 - (rhomags[1]/cAU)
    t_3 = J_3 - (rhomags[2]/cAU)
    tau1 = (t_1-t_2)*k
    tau3 = (t_3-t_2)*k
    tau0 = (t_3-t_1)*k 
    return [tau1, tau3, tau0]

def getSunVec1(fileName, i):
    data = np.array(open(fileName).readlines())
    vec = np.array(data[i].split())
    vec = vec.astype(float)
    return vec

def getSunVec2(fileName, i):
    data = np.array(open(fileName).readlines())
    vec = np.array(data[i].split())
    vec = vec.astype(float)
    return vec

def getSunVec3(fileName, i):
    data = np.array(open(fileName).readlines())
    np.array(data[i].split())
    position = position.astype(float)
    return position
     
# def JulianDate(DateNtime):
#     UT = (DateNtime[3].astype(float)) + ((DateNtime[4].astype(float)))/60 + ((DateNtime[5].astype(float)))/3600
#     abc = 7*((DateNtime[0].astype(float)) + ((DateNtime[1].astype(float)) + 9)//12)
#     #abc = (7(float(Years) + ((float(months) + 9)/12)))
#     J_o = 367*((DateNtime[0].astype(float))) - (abc/(4)) + ((275*((DateNtime[1].astype(float))))//9) + (DateNtime[2].astype(float)) + 1721013.5 
#     JD = (UT/24) + J_o
#     return JD

def JulianDate(DateNtime):
    UT = (DateNtime[3].astype(float)) + (DateNtime[4].astype(float)/60) + (DateNtime[5].astype(float)/3600)
    j0 = 367*(DateNtime[0].astype(float)) - int(7*(DateNtime[0].astype(float) + int((DateNtime[1].astype(float)+9)/12))/4) + int(275*(DateNtime[1].astype(float))/9) + DateNtime[2].astype(float) + 1721013.5
    jd = j0 + (UT/24)
    return jd

def TimeOfPerihelion(a, M, J_o):
    T = J_o - M*(math.sqrt((a**3)/(k**2)))
    return T 

# def getSUNvec(fileName):
#     data = np.array(open(fileName).readlines())
#     position = np.array(data[7].split())
#     R = position.astype(float)
#     return R 

def Generate_ra_dec(a, e, E, incline, small_Omega, Big_Omega, R):
    Epsilon = math.radians(23.44)

    x_e = a*(np.cos(E)) - a*e 
    y_e = a*(np.sqrt(1 - (e)**2))*(np.sin(E))
    r_e = np.array([[x_e], [y_e], [0]])
   

    rot_matrix1 = np.array([[math.cos(math.radians(small_Omega)), (-1)*(np.sin(np.radians(small_Omega))), 0], [np.sin(np.radians(small_Omega)), np.cos(np.radians(small_Omega)), 0], [0, 0, 1]])


    rot_matrix_2 = np.array([[1, 0, 0], [0, np.cos(np.radians(incline)), (-1)*(np.sin(np.radians(incline)))], [0, np.sin(np.radians(incline)), np.cos(np.radians(incline))]])


    rot_matrix_3 = np.array([[np.cos(np.radians(Big_Omega)), (-1)*(np.sin(np.radians(Big_Omega))), 0], [np.sin(np.radians(Big_Omega)), np.cos(np.radians(Big_Omega)), 0], [0, 0, 1]])

  

    rot_matrix_4 = np.array([[1, 0, 0], [0, np.cos(Epsilon), (-1)*(np.sin(Epsilon))], [0, np.sin(Epsilon), np.cos(Epsilon)]])

    r_1= np.dot(rot_matrix1, r_e)
    r_2 = np.dot(rot_matrix_2, r_1)
    x_1 = np.dot(rot_matrix_3, r_2)

    r_eq = np.dot(rot_matrix_4, np.dot(rot_matrix_3, np.dot(rot_matrix_2, np.dot(rot_matrix1, r_e))))

    R_eq = R
   
    r_eq = [r_eq[0][0], r_eq[1][0], r_eq[2][0]]
    
    ert_ast = np.add(R_eq, r_eq)
    ert_ast_unit = (ert_ast/(np.linalg.norm(ert_ast)))
    
    dec = np.arcsin(ert_ast_unit[2])

    s = ((ert_ast_unit[0])/(np.cos(dec)))
    ra = np.arccos((ert_ast_unit[0])/(np.cos(dec)))
    if s<=0:  
        test = math.asin(s)
        if test<= math.pi:
            ra = (ra + math.pi) 
    else: 
        test = math.asin(s)
        if test>math.pi :
           ra = 2*math.pi - ra  
    
    return math.degrees(dec), math.degrees(ra)

def error_percent(file, e, a, incline, Big_Omega, Nu, small_Omega, M, T):
    data = np.array(open(file).readlines())
    vec1 = np.array(data[4].split())
    vec2 = [e, a, incline, Big_Omega, Nu, small_Omega, M, T]
    vec1 = [value.replace(r"[a-zA-Z]",'') for value in vec1]
    difference = []
    for i in range(len(vec1)):
        difference.append(abs((float(vec1[i])-float(vec2[i]))))
    #difference = vec1 - vec2
    error = [] 
    for i in range(len(difference)):
        error.append((difference[i])/(float(vec1[i])))
    percent = (np.array(error))*100 
    return percent 

def error(calc, exp):#non-degrees
    return abs(calc-exp)/exp*100
def degError(calc, exp):
    return abs(calc-exp)/360*100

def getFG(taus, i, r_2):
    r2mag = np.linalg.norm(r_2)
    u = 1/((r2mag)**3)
    f = 1 - (1/2)*u*(taus[i]**2) 
    g = taus[i] - (1/6)*u*(taus[i]**3) 
    return f, g 

def get_r2dot(f1, f3, g1, g3, SunVec3, SunVec1, rho1, rho3, rho_hat_2, rho2mag, SunVec2):
    r1 = rho1 - SunVec1
    r2 =np.multiply(rho_hat_2, rho2mag) - SunVec2
    r3 = rho3 - SunVec3

    d1 = -(f3)/(f1*g3 - f3*g1)
    d3 = (f1)/(f1*g3 - f3*g1)
    r2dot = d1*r1 + d3*r3
    return r2dot

def ToEquitorial(vec):
    Epsilon = math.radians(23.44)
    rot_matrix_4 = np.array([[1, 0, 0], [0, np.cos(Epsilon), (-1)*(np.sin(Epsilon))], [0, np.sin(Epsilon), np.cos(Epsilon)]])
    rot_matrix = np.linalg.inv(rot_matrix_4)
    eq_vec = np.dot(rot_matrix, vec)
    return eq_vec

def OrbitalElements(position, velocity, corrected_jTimes, DateNtime):
        
        h = angularMomentum(position, velocity)
        a = SemiMajAx(position, velocity)
        e = eccentricity(h, a)
        # a = 3.423119774454217
        # e = 0.6770088849489315
        incline = inclination(h)
        Big_Omega = LongAscending(h, incline)
        u = Angle_U(position, incline)
        #Nu = truAnomaly(a, e, h, position, velocity)
        Nu = get_anom(position, velocity, a, h, e)
        small_Omega = Agr_Perhelion(u, Nu)
        # J_o = JulianDate(DateNtime)
        J_o = DateNtime
        T = corrected_jTimes[1] #of the second observation
        # M = (MeanAnomaly(J_o, T, a)) + 
        M = MeanAnomaly2(position, a, e, Nu)
        # if (Nu <= 180):
        #     E = math.acos((1/e)*(1 - (np.linalg.norm(position)/a))) # Eq 52
        # else:
        #     E = 2*pi - math.acos((1/e)*(1 - (np.linalg.norm(position)/a))) # Eq 52
        #     E = E * (180/pi)
        # M = (E - (e*math.sin(E*(pi/180)))*(180/pi))  # Eq 53

        # M = MeanAnomaly2(position, a, e, Nu)
        #print(math.degrees(MeanAnomaly(J_o, T, a)))
        #print('mean anomaly 2 '  + str(MeanAnomaly2(position, a, e, Nu)))
        #E = newton_raphson(0.0, 0.00000001, 10000, e, M)
        E = EccentricAnomaly(e, position, a, Nu)
        # R = getSUNvec(fileName)
        # DEC_RA = Generate_ra_dec(a, e, E, incline, small_Omega, Big_Omega, R)
        # print(f"e: {e}")
        # print(f"h: {h}")
        # print(f"semimajor axis: {a}")
        # print(f"inclination: {incline}")
        # print(f"Eccentric Anomaly: {E}")
        # print(f"ascending node longitude: {Big_Omega}")
        # print(f"true anomaly: {Nu}")
        # print(f"argument of perihelion: {small_Omega}")
        # print(f"mean anomaly: {M}") 
        # print("******************errors*************************")
        # print(f"semi major error = {error(a, 3.423119774454217E+00)}")
        # print(f"eccentricity error = {error(e, 6.770088849489315E-01)}")
        # print(f"inclination error = {error(incline, 5.616390117735835E+01)}")
        # print(f"logitude of ascending node error = {error(Big_Omega, 2.530523773999922E+02)}")
        # print(f"True anomaly error = {error(Nu, 7.243923573300094E+01)}")
        # print(f"argument of perihelion error = {error(small_Omega, 3.008338981468623E+02)}")
        # print(f"Mean Anomaly error = {error(M,1.303659573706276E+01)}")

        #print("time of perihelion", TimeOfPerihelion(a, M, J_o))

        return a, e, incline, Big_Omega, small_Omega, M

def SEL(taus,Sun2,rhohat2,Ds):
    roots = [0.,0.,0.] #for up to three real, positive roots
    rhos = []
    A1 = taus[1]/taus[2]
    B1 = (A1*((taus[2])**2 - (taus[1])**2))/6
    A3 = (-taus[0])/(taus[2])
    B3 = (A3*((taus[2])**2 - (taus[0])**2))/6
    A = (A1*Ds[1] - Ds[2] + A3*Ds[3])/(-Ds[0])
    B = (B1*Ds[1] + B3*Ds[3])/(-Ds[0])
    E = -2*(np.dot(rhohat2, Sun2))
    F = (np.linalg.norm(Sun2))**2
    a = -(A**2 + A*E + F)
    b = -(2*A*B + B*E)
    c = -B**2
    polyArray = [c, 0, 0, b, 0, 0, a, 0, 1]
    #print('coefficients', a, b, c)
    roots = np.polynomial.polynomial.polyroots(polyArray)
    #print(f"roots: {roots}")
    real_positive_roots = []

    for value in roots:
        if np.imag(value) == 0 and np.real(value) > 0.05:
            real_positive_roots.append(np.real(value))
    
    # rhos1 = [] 
    for i in range(len(real_positive_roots)):
        rhos.append(A + (B/(real_positive_roots[i]**3)))
   
    # print('real positive', real_positive_roots)
    # print('rho', rhos)

    r_2 = []
    rho2mag = []

    for i in range(len(rhos)):
        if rhos[i] > 0.05:
            rho2mag.append(rhos[i])
            r_2.append(real_positive_roots[i])
            # r_2 = r2_roots[i]
            # rho2mag = rhos[i]

    if len(r_2) == 1:
        return [r_2[0], rho2mag[0]] 
    else:
        indices = int(input("which rho value 0 indexed"))
        return [r_2[indices], rho2mag[indices]]      

def Final(fileName, rms_ra, rms_dec):
#getting the observation data and time from reading the input file
    DateNtime1 = getDate_obs1(fileName)
    DateNtime2 = getDate_obs2(fileName)
    DateNtime3 = getDate_obs3(fileName)
#converting to julian date 
    J_1 = JulianDate(DateNtime1)
    J_2 = JulianDate(DateNtime2)
    J_3 = JulianDate(DateNtime3)

    #print(f"Julian date {J_1, J_2, J_3}")
#calculating the tau values to input into the lagrange function, to get r2 roots and rho2 initial value 
    tau1 = k*((J_1 - J_2))
    tau3 = k*((J_3 - J_2))
    tau0 = k*((J_3 - J_1))

    taus = np.array([tau1, tau3, tau0])
    #print(f"taus {taus}")
#read file and got the ra and dec for each observation to calculate the unit rho vector for each observation times
    RADEC1 = (RA_DEC1(fileName))
    ras1 = np.random.normal(RADEC1[0], rms_ra, 200000)
    decs1 = np.random.normal(RADEC1[1], rms_dec, 200000)
    
    RADEC2 = (RA_DEC2(fileName))
    ras2 = np.random.normal(RADEC2[0], rms_ra, 200000)
    decs2 = np.random.normal(RADEC2[1], rms_dec, 200000)
    
    RADEC3 = (RA_DEC3(fileName))
    ras3 = np.random.normal(RADEC3[0], rms_ra, 200000)
    decs3 = np.random.normal(RADEC3[1], rms_dec, 200000)

    #print(RADEC1, RADEC2, RADEC3)
    semi_major_arr = []
    eccent_arr = []
    inclin_arr = []
    longit_asc_arr = []
    arg_perih_arr = []
    mean_anom_arr = []
    for i in range(200000):
        rho_hat_1 = [np.cos(np.radians((ras1[i])))*np.cos(np.radians((decs1[i]))), np.sin(np.radians((ras1[i])))*np.cos(np.radians((decs1[i]))), np.sin(np.radians(decs1[i]))]
        rho_hat_2 = [np.cos(np.radians((ras2[i])))*np.cos(np.radians((decs2[i]))), np.sin(np.radians((ras2[i])))*np.cos(np.radians((decs2[i]))), np.sin(np.radians(decs2[i]))]
        rho_hat_3 = [np.cos(np.radians((ras3[i])))*np.cos(np.radians((decs3[i]))), np.sin(np.radians((ras3[i])))*np.cos(np.radians((decs3[i]))), np.sin(np.radians(decs3[i]))]
        
        #print(np.linalg.norm(rho_hat_1), np.linalg.norm(rho_hat_2), np.linalg.norm(rho_hat_3))
    #sunvector read from file for further computation and calculating the D's and array to use in the lagrange function defined SEL
        SunVec1 = getSunVec1(fileName, 3)
        SunVec2 = getSunVec1(fileName, 7)
        SunVec3 = getSunVec1(fileName, 11)
        #print(SunVec1, SunVec2, SunVec3)

        D_o = np.dot(rho_hat_1, np.cross(rho_hat_2, rho_hat_3))
        D_21 = np.dot(np.cross(rho_hat_1, SunVec1), rho_hat_3)
        D_22 = np.dot(np.cross(rho_hat_1, SunVec2), rho_hat_3)
        D_23 = np.dot(np.cross(rho_hat_1, SunVec3), rho_hat_3)
        #print(f"Ds {D_o, D_21, D_22, D_23}")
        
        D_11 = np.dot(np.cross(SunVec1, rho_hat_2), rho_hat_3)
        D_12 = np.dot(np.cross(SunVec2, rho_hat_2), rho_hat_3)
        D_13 = np.dot(np.cross(SunVec3, rho_hat_2), rho_hat_3)
        #print(D_11, D_12, D_13)

        D_31 = np.dot(rho_hat_1, np.cross(rho_hat_2, SunVec1))
        D_32 = np.dot(rho_hat_1, np.cross(rho_hat_2, SunVec2))
        D_33 = np.dot(rho_hat_1, np.cross(rho_hat_2, SunVec3))
        #print(D_31, D_32, D_33)

        Ds = np.array([D_o, D_21, D_22, D_23])
    
        r2_rhos = (SEL(taus,SunVec2,rho_hat_2,Ds))
        #print(f" lagrange {r2_rhos}")
    
        r2mag = np.array(r2_rhos[0])
        rho2mag = np.array(r2_rhos[1])

        #print(f"rho2mag {rho2mag}")

    #definig initial f and g variables to be further used to define small d's and C's, then using the d arrays to find the r2 dot as a linear combination of r1 and r3
        f1_g1 = np.array(getFG(taus, 0, r2mag))
        f3_g3 = np.array(getFG(taus, 1, r2mag))
        f1 = f1_g1[0]
        g1 = f1_g1[1]
        f3 = f3_g3[0]
        g3 = f3_g3[1]
        #print(f"FG {f1, g1, f3, g3}")
    
        Cs = [(g3/(f1*g3 - g1*f3)), -1, (-g1/(f1*g3 - g1*f3))]
        #print(f"Cs {Cs}")

    # rho magnitudes used in calculating r1 and r3 vectors from the sun, rho and position triangle
        rho1mag = (Cs[0]*D_11 + Cs[1]*D_12 + Cs[2]*D_13)/(Cs[0]*D_o)
        rho3mag = (Cs[0]*D_31 + Cs[1]*D_32 + Cs[2]*D_33)/(Cs[2]*D_o)

        rhoMags = np.array([rho1mag, rho2mag, rho3mag])

        #print(f"rhomags {rhoMags}")
    
        rho1 = np.multiply(rho_hat_1, rho1mag)
        rho3 = np.multiply(rho_hat_3, rho3mag)

        #print(f"rho 1 rho 3{rho1, rho3}")

        r1 = rho1 - SunVec1
        r2 = np.multiply(rho_hat_2, rho2mag) - SunVec2
        # print("r2", r2)
        r3 = rho3 - SunVec3

        # print(r1, r2, r3)
        

        ds = [(-f3/(f1*g3 - f3*g1)), (f1/(f1*g3 - f3*g1))]

        # print(ds)

    #initial r2 dot value for the fourth order f&g iteration 
        # print('dsss', ds[0], ds[1])
        # print('r1', r1)
        # print('r3', r3)
        r2dot = np.add(np.dot(ds[0], r1),np.dot(ds[1], r3))

        # print('r2 dot', r2dot)
        # print(rho_hat_1, rho_hat_2, rho_hat_3)


        rho2mag_new = 0 
        iterations = 0 
        diff = 100
        # print("************************************************")
        # print('Main Iteration Loop')
        while diff >= 10E-9 :
            prev = rhoMags[1] 
            corrected_times = ligthtimeCorrection(J_1, J_2, J_3, rhoMags)
            # print('corrected_times', corrected_times)
            a = get_semi_major(np.linalg.norm(r2), np.linalg.norm(r2dot))
            e = get_eccent(np.cross(r2, r2dot), a)
            FG = fgfunc(corrected_times[0], corrected_times[1], r2, r2dot, np.linalg.norm(r2), a, e)
            #FG = fg(corrected_times[0], corrected_times[1], r2, r2dot)
            # print('Fg', FG)
            ds = [(-FG[1]/(FG[0]*FG[3] - FG[1]*FG[2])), (FG[0]/(FG[0]*FG[3] - FG[1]*FG[2]))]
            # print('ds', ds)
            Cs = [(FG[3]/(FG[0]*FG[3] - FG[2]*FG[1])), -1, (-FG[2]/(FG[0]*FG[3] - FG[2]*FG[1]))]
            # print('Cs', Cs)
            # print('r2', r2)
            # print('a', a)
            # print('c', Cs)
            rho1mag = (Cs[0]*D_11 + Cs[1]*D_12 + Cs[2]*D_13)/(Cs[0]*D_o)
            # print((Cs[0]*D_11 + Cs[1]*D_12 + Cs[2]*D_13))
            # print((Cs[0]*D_o))
            rho2mag = (Cs[0]*D_21 + Cs[1]*D_22 + Cs[2]*D_23)/(Cs[1]*D_o)
            rho3mag = (Cs[0]*D_31 + Cs[1]*D_32 + Cs[2]*D_33)/(Cs[2]*D_o)
        
            diff = abs(prev - rho2mag)
            # print(diff)
            rhoMags = np.array([rho1mag, rho2mag, rho3mag])
            # print('rhoMags', rhoMags)
            # print(f"julian dates {J_1, J_2, J_3}")
            
            rho1 = np.multiply(rho_hat_1, rho1mag)
            rho2 = np.multiply(rho_hat_2, rho2mag)
            rho3 = np.multiply(rho_hat_3, rho3mag)
            # print(5)
            r1 = rho1 - SunVec1
            # print('rho_hat_2', rho_hat_2)
            # print('rho2mag', rho2mag)
            r2 = rho2 - SunVec2
            r3 = rho3 - SunVec3
            # print('r values', r1, r2, r3)

            r2dot = np.add(np.multiply(ds[0], r1), np.multiply(ds[1], r3))
            # print('r2dot', r2dot)

            iterations += 1

    #final position and velocity vectors form the second observation in equitorial
        #print(r2, iterations, r2dot)
    #converting to ecliptic coordinates (rotation matrix)
        position = ToEquitorial(r2)
        velocity = ToEquitorial(r2dot)
        # print(np.linalg.norm(position))
        # print(np.linalg.norm([-2.953795826337923E-01, -4.469102330078713E-01, 2.398463566701330E-01]))
        # velocity = [-9.411194138477498E-03, -7.669174894923517E-03, 1.394054361217390E-02]

        #print(f"real position: {position}, velocity: {velocity}")

    #ALL orbital elements from the OD codes 1 to 4 
        # print("************************************************")
        # print("Orbital Elements")
    #updated time for the second observation from 00:00 to 7:00 and updated outputs 
        # time = np.array([2021.00, 7.00, 24.00, 7, 0, 0])

        corrected_jTimes = ligthtimeCorrection1(J_1, J_2, J_3, rhoMags)
        # J_o = JulianDate(time)

        
        # corrected_JTimes = ligthtimeCorrection1(J_1, J_o, J_3, rhoMags)
    
        a, e, i, big_Omega, small_Omega, M = OrbitalElements(position, velocity, corrected_jTimes, corrected_jTimes[1])
        semi_major_arr.append(a)
        eccent_arr.append(e)
        inclin_arr.append(i)
        longit_asc_arr.append(big_Omega)
        arg_perih_arr.append(small_Omega)
        mean_anom_arr.append(M)
    
    #print(f"real position: {position}, velocity: {velocity}")

    return position, velocity, semi_major_arr, eccent_arr, inclin_arr, longit_asc_arr, arg_perih_arr, mean_anom_arr
    #"D:/OD_lib/Horizons_data"