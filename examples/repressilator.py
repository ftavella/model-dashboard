import numpy as np

variables ={
    "B": 0.0,
    "Bn": 0.0,
    "C": 0.0,
    "Cn": 0.0,
    "T": 0.27,
    "Tn": 0.27,
    "N": 0.0,
}

parameters = {
    "ks": {'init': 1.5, 'min': 0.0, 'max': 4.0},
    "W": {'init': 0.6, 'min': 0.0, 'max': 3.0},
    "Wn": {'init': 0.9, 'min': 0.0, 'max': 3.0},
    "Ttot": {'init': 0.27, 'min': 0.0, 'max': 3.0},
    # Degradation
    "aA": {'init': 0.018, 'min': 0.0, 'max': 1.0},
    "bA": {'init': 0.036, 'min': 0.0, 'max': 1.0},
    "KA": {'init': 24.0, 'min': 0.0, 'max': 40.0},
    "nA": {'init': 21.0, 'min': 0.0, 'max': 40.0},
    # Cdc25
    "aT": {'init': 0.16, 'min': 0.0, 'max': 1.0},
    "bT": {'init': 0.7, 'min': 0.0, 'max': 2.0},
    "KT": {'init': 28.0, 'min': 0.0, 'max': 40.0},
    "nT": {'init': 6.8, 'min': 0.0, 'max': 10.0},
    # Wee1
    "aW": {'init': 0.08, 'min': 0.0, 'max': 1.0},
    "bW": {'init': 0.43, 'min': 0.0, 'max': 2.0},
    "KW": {'init': 35.0, 'min': 0.0, 'max': 40.0},
    "nW": {'init': 3.1, 'min': 0.0, 'max': 6.0},
    # Nucleus
    "aN": {'init': 0.0, 'min': 0.0, 'max': 1.0},
    "bN": {'init': 1.0, 'min': 0.0, 'max': 2.0},
    "KN": {'init': 1.0, 'min': 0.0, 'max': 2.0},
    "nN": {'init': 5.0, 'min': 0.0, 'max': 10.0},
    "alpha": {'init': 0.083, 'min': 0.0, 'max': 1.0},
    "beta": {'init': 16.0, 'min': 0.0, 'max': 30.0},
    "gamma": {'init': 10.44, 'min': 0.0, 'max': 30.0},
    # Import and Export
    "aI": {'init': 0.001, 'min': 0.0, 'max': 1.0},
    "KI": {'init': 1.0, 'min': 0.0, 'max': 3.0},
    "nI": {'init': 4.0, 'min': 0.0, 'max': 6.0},
    "aIC": {'init': 0.1, 'min': 0.0, 'max': 1.0},
    "bIC": {'init': 0.0333, 'min': 0.0, 'max': 1.0},
    "aIT": {'init': 1.2, 'min': 0.0, 'max': 3.0},
    "bIT": {'init': 0.03, 'min': 0.0, 'max': 1.0},
    "kiB": {'init': 1.468, 'min': 0.0, 'max': 2.0},
    "keB": {'init': 0.204, 'min': 0.0, 'max': 2.0},
    "kiC": {'init': 1.428, 'min': 0.0, 'max': 2.0},
    "keC": {'init': 0.168, 'min': 0.0, 'max': 2.0},
    "kiT": {'init': 1.092, 'min': 0.0, 'max': 2.0},
    "keT": {'init': 0.390, 'min': 0.0, 'max': 2.0},
}


def Hi(x,K,n):
    """Increasing Hill function"""
    a = np.power(abs(x),n)
    b = np.power(K,n)
    y = np.divide(a,a+b)
    return y


def Hd(x,K,n):
    """Decreasing Hill function"""
    a = np.power(abs(x),n)
    b = np.power(K,n)
    y = np.divide(b,a+b)
    return y


def equations(t, y, p):
    """ODE equations for the sink source model"""
    B, Bn, C, Cn, T, Tn, N = y

    # Interactions 
    H_T = (p.aT + p.bT * Hi(C, p.KT, p.nT))
    H_Tn = (p.aT + p.bT * Hi(Cn, p.KT, p.nT)) # C instead of Cn in Hi
    H_W = (p.aW + p.bW * Hd(C, p.KW, p.nW))
    H_Wn = (p.aW + p.bW * Hd(Cn, p.KW, p.nW)) # C instead of Cn in Hd
    H_A = (p.aA + p.bA * Hi(C, p.KA, p.nA))
    H_N = (p.aN + p.bN * Hi(N, p.KN, p.nN))

    # Import and Export functions
    I_C = (p.aI + (1.0 - p.aI) * Hd(N, p.KI, p.nI)) * (p.aIC + p.bIC * Cn)
    E_C = 1.0
    I_T = (p.aI + (1.0 - p.aI) * Hd(N, p.KI, p.nI)) * (p.aIT + p.bIT * C)
    E_T = 1.0

    # Total Cdk1:CyclinB in the cytoplasm
    dB = p.ks - H_A * B
    dB += - p.kiB * I_C * B + p.keB * E_C * Bn
    # Total Cdk1:CyclinB in the nucleus
    dBn = p.kiB * I_C * B - p.keB * E_C * Bn
    # Active Cdk1:CyclinB in the cytoplasm
    dC = p.ks - H_A * C
    dC += - p.kiC * I_C * C + p.keC * E_C * Cn
    dC += T * H_T * (B - C)
    dC += - p.W * H_W * C 
    # Active Cdk1:CyclinB in the nucleus
    dCn = p.kiC * I_C * C - p.keC * E_C * Cn
    dCn += Tn * H_Tn * (Bn - Cn)
    dCn += - p.Wn * H_Wn * Cn 
    # Active Cdc25 in the cytoplasm
    dT = - p.kiT * I_T * T + p.keT * E_T * Tn
    # Active Cdc25 in the nucleus
    dTn = p.kiT * I_T * T - p.keT * E_T * Tn
    # Nucleus
    dN = p.alpha * (C + Cn) + p.beta * H_N - p.gamma * N

    return np.array([dB, dBn, dC, dCn, dT, dTn, dN])
