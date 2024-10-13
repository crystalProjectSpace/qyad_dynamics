import numpy as np
# барометрическое уравнение
def Ro(H):
    return 1.2 * np.exp(-H / 7500)
# тестовые АДХ 
def getADX(alpha):
    return {
        'Cxa': abs(0.135 * np.cos(alpha)) + abs(0.27 * np.sin(alpha)),
        'Cya': 0.02 * (alpha * 57.3)
    }
# ускорение свободного падения
g = 9.81

# на основе вектора входных параметров и программы управления сформировать законы управления
def setControls(state, ctrl, prog):
    t = state[0]      # время 
    th = state[2]     # угол наклона траектории к местному горизонту
    d_aoa = state[4]   # производная от угла атаки
    d_roll = state[5] # производная от угла крена
    aoa = state[6]    # угол атаки
    roll = state[7]   # угол крена

    K_AOA = prog.k_prop_aoa
    K_D_AOA = prog.k_diff_aoa
    K_ROLL = prog.k_prop_roll
    K_D_ROLL = prog.k_diff_roll
    K_THRUST = prog.k_prop_thrust
    
    delta_aoa = prog.get_alpha(t) - aoa 
    delta_roll = prog.get_psi(t) - roll
    delta_th = th - prog.get_th(t)
    
    signal_aoa = delta_aoa * K_AOA - K_D_ROLL * d_aoa
    signal_roll = delta_roll * K_AOA - K_D_ROLL * d_roll
    signal_thrust = max(0, 1 - delta_th * K_THRUST)

    _rev1 = 0.5 + signal_aoa + signal_roll
    _rev2 = 0.5 + signal_aoa - signal_roll
    _rev3 = 0.5 - signal_aoa + signal_roll
    _rev4 = 0.5 - signal_aoa - signal_roll

    if t < 51.0:
        ctrl['rev1'] = min(1, max(0, _rev1) * signal_thrust) 
        ctrl['rev2'] = min(1, max(0, _rev2) * signal_thrust)
        ctrl['rev3'] = min(1, max(0, _rev3) * signal_thrust)
        ctrl['rev4'] = min(1, max(0, _rev4) * signal_thrust)
    else:
        ctrl['rev1'] = 0.55
        ctrl['rev2'] = 0.55
        ctrl['rev3'] = 0.55
        ctrl['rev4'] = 0.55
# получить производные от кинематических функций
def derivatives(state, ctrl, prms):
    V = state[1]
    Th = state[2]
    Psi = state[3]
    wAlpha = state[4]
    wGamma = state[5]
    alpha = state[6]
    gamma = state[7]
    Y = state[9]
    R = prms['w'] * (ctrl['rev1'] + ctrl['rev2'] + ctrl['rev3'] + ctrl['rev4'])
    Q = Ro(Y) * V * V * 0.5
    adx = getADX(alpha)
    Cxa = adx['Cxa']
    Cya = adx['Cya']

    CA = np.cos(alpha)
    SA = np.sin(alpha)
    CG = np.cos(gamma)
    SG = np.sin(gamma)
    CTH = np.cos(Th)
    STH = np.sin(Th)
    CPSI = np.cos(Psi)
    SPSI = np.sin(Psi)
    XA = Cxa * Q * prms['S']
    YA = Cya * Q * prms['S']
    G = prms['m'] * g
    MV = prms['m'] * V
    
    AX = (-R * SA - XA - G * STH) / prms['m']
    dTH = (R * CA * CG + YA * CG - G * CTH) / MV
    dPSI = (R * CA * SG + YA * SG) / MV

    dAlpha = prms['w'] * (ctrl['rev1'] + ctrl['rev2'] - ctrl['rev3'] - ctrl['rev4']) * prms['lx'] / prms['Jz']
    dGamma = prms['w'] * (ctrl['rev1'] + ctrl['rev3'] - ctrl['rev2'] - ctrl['rev4']) * prms['ly'] / prms['Jx']

    return [
        0,
        AX,
        dTH,
        dPSI,
        dAlpha,
        dGamma,
        wAlpha,
        wGamma,
        V * CTH * CPSI,
        V * STH,
        V * CTH * SPSI,
        0,
        0,
        0,
        0,
    ]
# Уравнений движения квадкоптера по заданному вектору начальных состояний, параметрам массы/геометрии и закону управления
def integrate(state, prms, prog, dT, tMax):
    dT_05 = dT * 0.5
    currentControls = {
        'rev1': 0.4,
        'rev2': 0.4,
        'rev3': 0.6,
        'rev4': 0.6,
    }
    t = 0
    i = 0
    result = [state.copy()]

    while t < tMax:
        _state = result[i].copy() 
        setControls(_state, currentControls, prog)

        K1 = derivatives(state, currentControls, prms)
        state_K1 = arrayCombination(_state, K1, 1, dT)
        K2 = derivatives(state_K1, currentControls, prms)
        derivSumm = arraySumm(K1, K2)
        
        state2 = arrayCombination(_state, derivSumm, 1, dT_05)
        t += dT
        state2[0] = t
        state2[11] = currentControls['rev1']
        state2[12] = currentControls['rev2']
        state2[13] = currentControls['rev3']
        state2[14] = currentControls['rev4']
        result.append(state2)
        i += 1

    return result
# сформировать текстовый вывод из массива результатов интегрирования
def create_text_output(arr):
    result = ""
    for row in arr:
        t = row[0]
        v = row[1]
        th = row[2] * 57.3
        psi = row[3] * 57.3
        daoa = row[4] * 57.3
        droll = row[5] * 57.3
        aoa = row[6] * 57.3
        roll = row[7] * 57.3
        x = row[8]
        y = row[9]
        z = row[10]
        r1 = row[11]
        r2 = row[12]
        r3 = row[13]
        r4 = row[14]
        temp_str = f"{t:.2f}\t{v:.1f}\t{th:.2f}\t{psi:.2f}\t{daoa:.3f}\t{droll:.3f}\t{aoa:.2f}\t{roll:.2f}\t{x:.1f}\t{y:.1f}\t{z:.1f}\t{r1:.2f}\t{r2:.2f}\t{r3:.2f}\t{r4:.2f}" 
        result = result + '\n' + temp_str
    return result
# сумма двух массивов
def arraySumm(U, V):
    return [u + v for u, v in zip(U, V)]
# линейная комбинация двух массивов с заданными коэффициентами
def arrayCombination(U, V, kU, kV):
    return [u * kU + v * kV for u, v in zip(U, V)]
# класс функции с заданными программами управления
class ControlProg:
    k_prop_aoa = 0
    k_diff_aoa = 0
    k_prop_roll = 0
    k_diff_roll = 0
    k_prop_thrust = 0

    @staticmethod
    def get_alpha(t):
        if (t < 20):
            return -(12.5 / 57.3)
        if (t < 25):
            return -(0 / 57.3)
        if (25 < t < 30):
            return 25 / 57.3
        if (30 < t < 35):
            return 45 / 57.3   
        if (35 < t < 40):
            return 55 / 57.3
        if (40 < t < 45):
            return 75 / 57.3
        if (45 < t < 60):
            return 85 / 57.3        

    @staticmethod
    def get_psi(t):
        return 0

    @staticmethod
    def get_th(t):
        if t < 20:
            return 0
        if 20 < t < 25:
            return -20 / 57.3
        if 25 < t < 30:
            return -30 / 47.3
        if 30 < t < 35:
            return -40 / 47.3
        if 35 < t < 40:
            return -40 / 47.3  
        if 40 < t < 45:
            return -55 / 47.3
        if 45 < t < 70:
            return -75 / 47.3        

    def __init__(self, k_prop_aoa, k_diff_aoa, k_prop_roll, k_diff_roll, k_prop_thrust):
        self.k_prop_aoa = k_prop_aoa
        self.k_diff_aoa = k_diff_aoa
        self.k_prop_roll = k_prop_roll
        self.k_diff_roll = k_diff_roll
        self.k_prop_thrust = k_prop_thrust

def simulation_run():
    test_prog = ControlProg(2.75, 12, 2.75, 12, 5.85)

    test_state = [0, 5, 0, 0, 0, 0, 0, 0, 0, 500, 0, 0.5, 0.5, 0.5, 0.5]

    test_prms = {
        'm': 2.5,
        'S': 0.032,
        'Jz': 1.75,
        'Jx': 0.5,
        'w': 12.5,
        'lx': 0.2,
        'ly': 0.15
    }

    result = integrate(test_state, test_prms, test_prog, 0.01, 60)
    res_str = create_text_output(result)
    output = open('result.txt', 'wt')
    output.write(res_str)
    output.close()
    print('simulation_complete')
print('test')
simulation_run()