import numpy as np

def Ro(H):
    return 1.2 * np.exp(-H / 7500)

def getADX(alpha):
    return {
        'Cxa': abs(0.135 * np.cos(alpha)) + abs(0.27 * np.sin(alpha)),
        'Cya': 0.02 * (alpha * 57.3)
    }

KTH = 2.5
KPSI = 2.5
g = 9.81

def setControls(state, ctrl, prog):
    t = state[0]
    Th = state[2]
    dAoa = state[4]
    dGamma = state[5]
    aoa = state[6]    
    roll = state[7]
    
    deltaA = (prog.getAlpha(t) - aoa) - dAoa * 4
    deltaPsi = (prog.getPsi(t) - roll) - dGamma * 4
    deltaTh = Th - prog.getTh(t)

    _rev1 = 0.5 + deltaA * KTH + deltaPsi * KPSI
    _rev2 = 0.5 + deltaA * KTH - deltaPsi * KPSI
    _rev3 = 0.5 - deltaA * KTH + deltaPsi * KPSI
    _rev4 = 0.5 - deltaA * KTH - deltaPsi * KPSI
    
    kThrust = max(0, 1 - 5.85 * deltaTh)
    
    ctrl['rev1'] = min(1, max(0, _rev1) * kThrust) 
    ctrl['rev2'] = min(1, max(0, _rev2) * kThrust)
    ctrl['rev3'] = min(1, max(0, _rev3) * kThrust)
    ctrl['rev4'] = min(1, max(0, _rev4) * kThrust) 

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

def arraySumm(U, V):
    return [u + v for u, v in zip(U, V)]

def arrayCombination(U, V, kU, kV):
    return [u * kU + v * kV for u, v in zip(U, V)]

test_state = [0, 5, 0, 0, 0, 0, 0, 0, 0, 500, 0, 0.5, 0.5, 0.5, 0.5]

test_params = {
    'm': 2.5,
    'S': 0.032,
    'Jz': 1.75,
    'Jx': 0.5,
    'w': 12.5,
    'lx': 0.2,
    'ly': 0.15
}

class TestProg:
    @staticmethod
    def get_alpha(t):
        return -(0.5 / 57.3)

    @staticmethod
    def get_psi(t):
        if t < 20:
            return 0
        if 20 < t < 30:
            return 15 / 57.3
        if t > 30:
            return 0

    @staticmethod
    def get_th(t):
        if t < 30:
            return 0
        if 30 < t < 35:
            return -20 / 57.3
        if t > 35:
            return 0

test_prog = TestProg()