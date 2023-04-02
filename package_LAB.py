import numpy as np
import matplotlib.pyplot as plt

def LeadLag_RT(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0,method='EBD'):
    
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
                # PV.append(((Kp * K)/(1 + K))*(((Tlead / Ts) + 1) * MV[-1] + (1 - Tlead / Tlag) * MV[-1]))
            else:   #if typo -> EBD (default)
                PV.append((1 / (1 + K)) * PV[-1] + Kp * K * ((Tlead * MV[-1]) / Ts) + (1 - Tlead / Ts) * MV[-2])
    else:
        PV.append(Kp*MV[-1])

def PID_RT(SP,E,MV,MVP,MVI,MVD,MVFF,man_mode,MVMan,MVmin,MVmax,PV,Ts,Kc,Ti,Td,alpha,E_init=0,method='EBD'):
    
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
        MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]

    # Saturation
    MVtot = MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]
    if MVtot > MVmax:
        MVI[-1] = MVmax - MVP[-1] - MVD[-1] - MVFF[-1]
    elif MVtot < MVmin:
        MVI[-1] = MVmin - MVP[-1] - MVD[-1] - MVFF[-1]
    
    MVtot = MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1]
    MV.append(MVtot)


def IMC_Tuning(K, T1, T2=0, theta=0, gamma=0.5):
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

    KcK = (T1 + T2) / (Tc + theta)
    
    Kc = KcK / K
    
    Ti = (T1 + T2)
    
    Td = (T1 * T2) / (T1 + T2)
    
    return (Kc, Ti, Td)

class PID:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 1.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 0.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 0.0
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 1.0

def Margins(P,C,omega, Show = True):
    
    """
    :P: Process as defined by the class "Process".
        Use the following command to define the default process which is simply a unit gain process:
            P = Process({})
        
        A delay, two lead time constants and 2 lag constants can be added.
        
        Use the following commands for a SOPDT process:
            P.parameters['Kp'] = 1.1
            P.parameters['Tlag1'] = 10.0
            P.parameters['Tlag2'] = 2.0
            P.parameters['theta'] = 2.0
        
        Use the following commands for a unit gain Lead-lag process:
            P.parameters['Tlag1'] = 10.0        
            P.parameters['Tlead1'] = 15.0   
            
    :C: Controller as defined by the class "PID".
        Use the following command to define the default process which is simply a unit gain process:
            C = PID({})
        
        Use the following commands for a PID Controller:
            C.parameters['Kc'] = 1.1
            C.parameters['Ti'] = 10.0
            C.parameters['Td'] = 2.0
            C.parameters['alpha'] = 1
        
        Use the following commands for a unit gain PID:
            C.parameters['Ti'] = 10.0        
            C.parameters['Td'] = 15.0 
        
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
    :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise Ps (P(j omega)) (vector of complex numbers) is returned.
    
    The function "Margins" generates the Bode diagram of the process P*C and prints the margins of the system. 
    """     
    
    s = 1j*omega
    
    Ptheta = np.exp(-P.parameters['theta']*s)
    PGain = P.parameters['Kp']*np.ones_like(Ptheta)
    PLag1 = 1/(P.parameters['Tlag1']*s + 1)
    PLag2 = 1/(P.parameters['Tlag2']*s + 1)
    PLead1 = P.parameters['Tlead1']*s + 1
    PLead2 = P.parameters['Tlead2']*s + 1
    
    Ps = np.multiply(Ptheta,PGain)
    Ps = np.multiply(Ps,PLag1)
    Ps = np.multiply(Ps,PLag2)
    Ps = np.multiply(Ps,PLead1)
    Ps = np.multiply(Ps,PLead2)
    
    
    CGain = C.parameters['Kc']*np.ones_like(Ptheta)
    C1 = (1+1/(C.parameters['Ti']*s)+ (C.parameters['Td']*s/(C.parameters['alpha']*C.parameters['Td']*s+1)))
   
    
    Cs = np.multiply(CGain,C1)

    
    Ps = np.multiply(Ps,Cs)
    
    
    if Show == True:
    
        fig, (ax_gain, ax_phase) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        # Gain part
        ax_gain.semilogx(omega,20*np.log10(np.abs(Ps)),label='P(s)')
        ax_gain.semilogx(omega,20*np.log10(np.abs(np.ones_like(PGain))),label='0')
        
        gain = 20*np.log10(np.abs(Ps))
        
 
        gain_min = np.min(20*np.log10(np.abs(Ps)/5))
        gain_max = np.max(20*np.log10(np.abs(Ps)*5))
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |P| [db]')
        ax_gain.set_title('Bode plot of P')
        ax_gain.legend(loc='best')
    
        # Phase part
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ps)),label='P(s)')
        phase = (180/np.pi)*np.unwrap(np.angle(Ps))
        
        indexc = np.where(np.round(gain,2) == 0.00)[0][0]
        wc = omega[indexc]
        indexu = np.where(np.round(phase,1) == -180.0)[0][0]
        wu = omega[indexu]
        ax_phase.vlines(wc,-180,phase[indexc])
        print('Gain margin : ',np.round(gain[indexu],5),' at ', np.round(wu,5), ' rad/s')
        print('Phase margin : ',np.round(phase[indexc],5)+180,' at ', np.round(wc,5), ' rad/s')
        
        ax_gain.vlines(wu,0,gain[indexu])
        
        
        ax_phase.axhline(y=(-180))
        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ps))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ps))) + 10
        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle P$ [°]')
        ax_phase.legend(loc='best')
    else:
        return Ps