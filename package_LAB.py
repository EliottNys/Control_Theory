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

def PID_RT(MV,MVP,MVI,MVD,MVFF,MVMan,Ts,Kc,Ti,Td,alpha,E,E_init=0,man_mode=False,method='EBD'):
    
    """
    The function "PID_RT" needs to be included in a "for or while loop".
    
    :MV: Manipulated value
    :MVP: Proportional action vector
    :MVI: Integrating action vector
    :MVD: Derivative action vector
    :MVFF: Feed Forward action vector
    :MVMan: Manual value vector
    :Ts: sampling period [s]
    :Kc: gain
    :Ti: Integral time constant
    :td: Derivative time constant
    :alpha: derivative filter smoothing factor
    :E: error vector
    :E_init: initial error value (optional: default value is 0)
    :man_mode: manual mode flag (optional: default value is False)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        TRAP: Trapezoïdal method
    
    The function "PID_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
    MVP.append(Kc * E[-1])

    if man_mode:
        MVI.append(MVMan - MVP - MVD - MVFF)
    else:
        if len(MV) == 0:
            MVI.append(((Kc * Ts) / Ti) * E[-1])
        elif method == 'TRAP':
            MVI.append(MVI[-1] + ((Kc * Ts) / (2 * Ti)) * (E[-1] + E[-2]))
        else:
            MVI.append(MVI[-1] + ((Kc * Ts) / Ti) * E[-1])
    
    Tfd = alpha * Td
    if method == 'TRAP':
        MVD.append(((Tfd - Ts / 2) / (Tfd + Ts / 2)) * MVD[-1] + ((Kc * Td) / (Tfd + Ts / 2)) * (E[-1] - E[-2]))
    else:
        MVD.append((Tfd / (Tfd + Ts)) * MVD[-1] + ((Kc * Td) / (Tfd + Ts)) * (E[-1] - E[-2]))
    
    # LL_RT(MVFF,Kp,Tlead,Tlag,Ts,PV,) ???

    MV.append(MVP[-1] + MVI[-1] + MVD[-1]) # + MVFF[-1] !!!