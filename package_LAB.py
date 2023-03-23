def LEADLAG_RT(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "LL_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    "LL_RT" = Lead-lag real-time
    The function "LL_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    

    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1 / (1 + K)) * PV[-1] + ((Kp * K) / (1 + K)) * ((1 + Tlead / Ts) * MV[-1] - (Tlead / Ts) * MV[-2]))
            elif method == 'EFD':
                PV.append((1 - K) * PV[-1] + Kp * K * ((Tlead / Ts) * MV[-1] + (1 - Tlead / Ts) * MV[-2]))
            elif method == 'TRAP':  #the equation is not in the slides
                PV.append(0)      
            else:   #if typo -> EBD (default)
                PV.append((1 / (1 + K)) * PV[-1] + Kp * K * ((Tlead * MV[-1]) / Ts) + (1 - Tlead / Ts) * MV[-2])
    else:
        PV.append(Kp*MV[-1])

def PID_RT(SP,E,MV,MVP,MVI,MVD,man_mode,MVMan,MVmin,MVmax,PV,Ts,Kc,Ti,Td,alpha,E_init=0,method='EBD'):
    
    """
    The function "PID_RT" needs to be included in a "for or while loop".
    
    :SP: Setpoint
    :E: Control error vector
    :MV: Manipulated value vector
    :MVP: Proportional action vector
    :MVI: Integrating action vector
    :MVD: Derivative action vector
    :MVFF: Feed Forward action vector
    :MVMan: Manual value vector
    :MVmin: Minimum MV value
    :MVmax: Maximum MV value
    :PV: Process value vector
    :Ts: sampling period [s]
    :Kc: gain
    :Ti: Integral time constant
    :td: Derivative time constant
    :alpha: derivative filter smoothing factor
    :E_init: initial error value (optional: default value is 0)
    :man_mode: manual mode flag (optional: default value is False)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        TRAP: Trapezoïdal method
    
    The function "PID_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    

    ##### ERROR #####

    if len(PV) == 0:
        E.append(E_init)
    else:
        E.append(SP[-1] - PV[-1])
    
    ##### PID #####

    # PROPORTIONAL action
    MVP.append(Kc * E[-1])

    # INTEGRATING action
    if method == 'TRAP':
        if len(MV) == 0:
            MVI.append(((Kc * Ts) / (2 * Ti)) * E[-1])
        else:
            MVI.append(MVI[-1] + ((Kc * Ts) / (2 * Ti)) * (E[-1] + E[-2]))
    else:   # EBD
        if len(MV) == 0:
            MVI.append(((Kc * Ts) / Ti) * E[-1])
        else:
            MVI.append(MVI[-1] + ((Kc * Ts) / Ti) * E[-1])
    # no EFD because it introduces instability
    
    # DERIVATIVE action
    Tfd = alpha * Td
    if method == 'TRAP':
        if len(MV) == 0:
            MVD.append(((Kc * Td) / (Tfd + Ts / 2)) * E[-1])
        else:
            MVD.append(((Tfd - Ts / 2) / (Tfd + Ts / 2)) * MVD[-1] + ((Kc * Td) / (Tfd + Ts / 2)) * (E[-1] - E[-2]))
    else:
        if len(MV) == 0:
            MVD.append(((Kc * Td) / (Tfd + Ts)) * E[-1])
        else:
            MVD.append((Tfd / (Tfd + Ts)) * MVD[-1] + ((Kc * Td) / (Tfd + Ts)) * (E[-1] - E[-2]))
    

    # Manual mode integrating action reset
    if man_mode[-1]:
        MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] # - MVFF

    # Saturation
    MVtot = MVP[-1] + MVI[-1] + MVD[-1]
    if MVtot > MVmax:
        MVI[-1] = MVmax - MVP[-1] - MVD[-1]
    elif MVtot < MVmin:
        MVI[-1] = MVmin - MVP[-1] - MVD[-1]
    
    MVtot = MVP[-1] + MVI[-1] + MVD[-1]
    MV.append(MVtot)


def IMC_Tuning(K, T1, T2, gamma=0.5):
    """
    IMC_Tuning(K, T1, T2, theta, gamma=0.5)
    
    This function computes the optimized IMC PID tuning parameters for a SOPDT Process.
    
    :param K: process gain
    :param T1: first time constant [s]
    :param T2: second time constant [s]
    :param gamma: used to compute the closed-loop time constant TCLP [s] such as
                  TCLP = gamma * T1p with T1p = main time constant of the process.
                  (range for gamma is [0.2 ... 0.9], default value is 0.5)
    
    :returns: (Kc, Ti, Td) - the PID controller parameters
    """
    
    Tc = gamma * T1  

    KcK = (T1 + T2) / Tc
    
    Kc = KcK / K
    
    Ti = (T1 + T2)
    
    Td = (T1 * T2) / (T1 + T2)
    
    return (Kc, Ti, Td)
