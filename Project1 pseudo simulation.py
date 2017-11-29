from brian2 import *
import numpy as np
from matplotlib import gridspec
import time
t0 = time.clock()

defaultclock.dt =0.01*ms
div=defaultclock.dt

El=-10*mV
gl=0.01*msiemens/cm**2 # start from 0.01 to 0.8
Cm=1.0*ufarad/cm**2
tau_ampa=5*ms # synaptic time constant = 5ms
area=20000*umetre**2

# Parameters of Izhikevich model
# LTS (0.02, 0.25, -65, 2)
a = 0.02; b = 0.25;
c = -65.0; d = 2.0;
g_synpk = 0.12*nsiemens # peak synaptic conductance
# maximum synaptic conductance
g_synmax = g_synpk/(tau_ampa/ms*exp(-1)) 

# NMDA neurotransmitter variables
Mg = 2 # unit mM. For ease of calculation, we do not give unit for this
       # found in Chapter 8 of "Foundations of Neuroscience" p. 163
betaN = 0.0066/ms
alphaN = 0.072/ms

v_th = 0*mV # reversal potential
transwdth = 1.0*ms

# GABA_b synaptic neuotransimtters
K3 = 0.18/ms
K4 = 0.034/ms
Kd = 100 # found in Chapter 8 of "Foundations of Neuroscience" in p.164
betab = 0.0012/ms
alphab = 0.09/ms
n = 4 # power of GABAb s
adjust = 100


eqs_Antenna = '''
# Izhikevich + alpha function
dv/dt = (0.04*(v)*(v) + 5*(v) + 140- u)/ms
                + (I/area/mV)/Cm : 1        
du/dt = a*(b*v - u)/ms : 1
I : amp
'''

eqs_motor = '''
# ANTENNA 1 (food agressor)
dg_ampaO/dt=-g_ampaO/tau_ampa + zO/ms: siemens
dzO/dt = -zO/tau_ampa: siemens

# ANTENNA 2 (positon coward)
dg_ampaT/dt=-g_ampaT/tau_ampa + zT/ms: siemens
dzT/dt = -zT/tau_ampa: siemens

# ANTENNA 3 (Female jitter)
gNMDA = g_synmax*B*sN*30000*bb : siemens
dsN/dt = alphaN*(1-sN)*Tr_pre - betaN*sN : 1
B = 1/(1+exp(-vN/(16.13*mV))*(Mg/3.57)) : 1
dvN/dt = gNMDA*(60*mV+vN)/(Cm*area) :  volt

Tr_pre=100*(tanh((t/ms-tspike/ms)/.005)-
            tanh((t/ms-(tspike/ms +transwdth/ms))/.005)):1

# ANTENNA 4 (Enemy GABAb excitatory)
gGABAb = (g_synmax*1000)*(sb**n)/(sb**n + Kd) : siemens
dr/dt = (alphab*Tr_pre_GABA*(1-r) - betab*r) : 1
dsb/dt = (K3*r - K4*sb) : 1

Tr_pre_GABA=aa*(tanh((t/ms-tspikeN/ms)/(.005*adjust))-
            tanh((t/ms-(tspikeN/ms +transwdth/
            (adjust*ms)))/(.005*adjust))):1
            
I: amp
tspike:second
tspikeN:second
aa:1
bb:1
'''

# Motor neurons
Motor1 = NeuronGroup(1,  eqs_motor,
                    clock=Clock(defaultclock.dt),method='euler')

Motor2 = NeuronGroup(1,  eqs_motor,
                    clock=Clock(defaultclock.dt),method='euler')


# Antenna neurons
Antenna1T = NeuronGroup(1,  eqs_Antenna,
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),method='euler')

Antenna2T = NeuronGroup(1,  eqs_Antenna,
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),method='euler')

Antenna3T = NeuronGroup(1,  eqs_Antenna,
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),method='euler')
Antenna4T = NeuronGroup(1,  eqs_Antenna,
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),method='euler')

Antenna1B = NeuronGroup(1,  eqs_Antenna,
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),method='euler')

Antenna2B = NeuronGroup(1,  eqs_Antenna,
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),method='euler')

Antenna3B = NeuronGroup(1,  eqs_Antenna,
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),method='euler')
Antenna4B = NeuronGroup(1,  eqs_Antenna,
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),method='euler')

# initialize voltage/conductance
# motor neurons voltage
# motor neurons conductance
Motor1.g_ampaO= 0*nS;Motor1.g_ampaT= 0*nS;
Motor1.sN= 0;Motor1.r= 0;Motor1.sb= 0;
Motor1.aa = 0; Motor1.bb = 0;Motor1.vN = 0

Motor2.g_ampaO= 0*nS;Motor2.g_ampaT= 0*nS;
Motor2.sN= 0;Motor2.r= 0;Motor2.sb= 0;
Motor2.aa = 0; Motor2.bb = 0;Motor2.vN = 0

# antenna neurons
Antenna1T.v= c;Antenna2T.v= c;Antenna3T.v= c;Antenna4T.v= c;
Antenna1B.v= c;Antenna2B.v= c;Antenna3B.v= c;Antenna4B.v= c;

# connect neurons

# Antenna1T, Antenna2T, Antenna3B, Antenna4B to Motor1
# Antenna1T -> Motor 1
SrM1A1 = Synapses(Antenna1T, Motor1, clock=Antenna1T.clock,
                model='''
        	w : siemens
        	''',
		pre='''
		g_ampaO += w	
		''')
SrM1A1.connect(i=[0],j=[0])
SrM1A1.w = g_synmax/(defaultclock.dt/ms)*3
SrM1A1.delay=5*ms
# Antenna2B -> Motor 1
SrM1A2 = Synapses(Antenna2B, Motor1, clock=Antenna2B.clock,
                model='''
        	w : siemens
        	''',
		pre='''
		g_ampaT += w		
		''')
SrM1A2.connect(i=[0],j=[0])
SrM1A2.w = g_synmax/(defaultclock.dt/ms)*3
SrM1A2.delay=5*ms
# Antenna3T -> Motor 1
SrM1A3 = Synapses(Antenna3T, Motor1, clock=Antenna3T.clock,
		pre='''
                tspike =t
                bb = 1
                vN = -65*mV 		
		''')
SrM1A3.connect(i=[0],j=[0])
SrM1A3.delay=5*ms
# Antenna4B -> Motor 1
SrM1A4 = Synapses(Antenna4B, Motor1, clock=Antenna4B.clock,
		pre='''
                tspikeN = t
                aa = 200
                bb = 0		
		''')
SrM1A4.connect(i=[0],j=[0])
SrM1A4.delay=5*ms

# Antenna1B, Antenna2B, Antenna3T, Antenna4T to Motor2
# Antenna1B -> Motor 2
SrM2A1 = Synapses(Antenna1B, Motor2, clock=Antenna1B.clock,
                model='''
        	w : siemens
        	''',
		pre='''
		g_ampaO += w		
		''')
SrM2A1.connect(i=[0],j=[0])
SrM2A1.w = g_synmax/(defaultclock.dt/ms)*3
SrM2A1.delay=5*ms
# Antenna2T -> Motor 2
SrM2A2 = Synapses(Antenna2T, Motor2, clock=Antenna2T.clock,
                model='''
        	w : siemens
        	''',
		pre='''
		g_ampaT += w		
		''')
SrM2A2.connect(i=[0],j=[0])
SrM2A2.w = g_synmax/(defaultclock.dt/ms)*3
SrM2A2.delay=5*ms
# Antenna3B -> Motor 2
SrM2A3 = Synapses(Antenna3B, Motor2, clock=Antenna3B.clock,
		pre='''
                tspike = t
                bb = 1
                vN = -65*mV  		
		''')
SrM2A3.connect(i=[0],j=[0])
SrM2A3.delay=5*ms
# Antenna4T -> Motor 2
SrM2A4 = Synapses(Antenna4T, Motor2, clock=Antenna4T.clock,
		pre='''
                tspikeN = t
                aa = 200
                bb = 0  		
		''')
SrM2A4.connect(i=[0],j=[0])
SrM2A4.delay=5*ms

# setup monitors
# motor neuron (voltage and conductance)
M1 = StateMonitor(Motor1,('g_ampaO','Tr_pre_GABA',
                 'g_ampaT', 'gNMDA','gGABAb'),record=True)
M2 = StateMonitor(Motor2,('g_ampaO', 'Tr_pre_GABA',
                 'g_ampaT', 'gNMDA','gGABAb'),record=True)
                 
# Antenna on top (only voltage)
A1T = StateMonitor(Antenna1T, ('I'), record = True)
A2T = StateMonitor(Antenna2T, ('I'), record = True)
A3T = StateMonitor(Antenna3T, ('I'), record = True)
A4T = StateMonitor(Antenna4T, ('I'), record = True)

# Antenna on Bottom (only voltage)
A1B = StateMonitor(Antenna1B, ('I'), record = True)
A2B = StateMonitor(Antenna2B, ('I'), record = True)
A3B = StateMonitor(Antenna3B, ('I'), record = True)
A4B = StateMonitor(Antenna4B, ('I'), record = True)

# give stimulus to antenna accordingly to the 
# scenairo
# run 200 msec for stablizing
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(200*ms,report='text')

# encounter 1st food
Antenna1T.I=10*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=10*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(20*ms,report='text')

# interval
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(60*ms,report='text')

# encounter 2nd food
Antenna1T.I=10*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=10*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(20*ms,report='text')

# interval
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(60*ms,report='text')

# encounter 1st poison
Antenna1T.I=0*nA; Antenna2T.I=10*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=3*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(20*ms,report='text')

# interval
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(60*ms,report='text')

# encounter 2nd poison
Antenna1T.I=0*nA; Antenna2T.I=3*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=10*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(20*ms,report='text')

# interval
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(60*ms,report='text')

# encounter female
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=10*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=10*nA; 
Antenna4B.I=0*nA;
run(20*ms,report='text')

# interval
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(60*ms,report='text')

# encounter enemy
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=2*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=10*nA;
run(20*ms,report='text')

# interval
Antenna1T.I=0*nA; Antenna2T.I=0*nA; Antenna3T.I=0*nA; 
Antenna4T.I=0*nA; 
Antenna1B.I=0*nA; Antenna2B.I=0*nA; Antenna3B.I=0*nA; 
Antenna4B.I=0*nA;
run(60*ms,report='text')

figure(1)
gs = gridspec.GridSpec(2,2, height_ratios = [2,1])
# conductance of motor 1
subplot(gs[0,0])
plot(M1.t/ms, M1.g_ampaO[0]/nS, '-g', lw=1, 
            label = 'g_Food Motor Neuron 0')
plot(M1.t/ms, M1.g_ampaT[0]/nS, '-c', lw=1, 
            label = 'g_Poison Motor Neuron 0')
plot(M1.t/ms, M1.gNMDA[0]/nS, '-r', lw=1, 
            label = 'g_Female Motor Neuron 0')
plot(M1.t/ms, M1.gGABAb[0]/nS, '-m', lw=1, 
            label = 'g_Enemy Motor Neuron 0')                      
xlabel('time (msec)', fontsize = 10)
ylabel('conductance (mS)', fontsize = 10)
xlim(150, 700); legend(loc = 0, prop={'size':10}); 
title('A', fontsize = 10)

# stimulus to antennas (top)
subplot(gs[1,0])
plot(A1T.t/ms, A1T.I[0]/nA, '-g', lw=1, 
            label = 'Stim to N1 (Top)')
plot(A2T.t/ms, A2T.I[0]/nA, '-c', lw=1, 
            label = 'Stim to N2 (Top)')
plot(A3T.t/ms, A3T.I[0]/nA, '-r', lw=1, 
            label = 'Stim to N3 (Top)')
plot(A4T.t/ms, A4T.I[0]/nA, '-m', lw=1, 
            label = 'Stim to N4 (Top)')
xlabel('time (msec)', fontsize = 10)
ylabel('Stimulus (nA)', fontsize = 10)
xlim(150, 700); ylim(0,11)
legend(loc = 0, prop={'size':10}) 
title('C', fontsize = 10)

# conductance of motor 2
subplot(gs[0,1])
plot(M2.t/ms, M2.g_ampaO[0]/nS, '-g', lw=1, 
            label = 'g_Food Motor Neuron 1')
plot(M2.t/ms, M2.g_ampaT[0]/nS, '-c', lw=1, 
            label = 'g_Poison Motor Neuron 1')
plot(M2.t/ms, M2.gNMDA[0]/nS, '-r', lw=1, 
            label = 'g_Female Motor Neuron 1')
plot(M2.t/ms, M2.gGABAb[0]/nS, '-m', lw=1, 
            label = 'g_Enemy Motor Neuron 1')                      
xlabel('time (msec)', fontsize = 10)
ylabel('conductance (mS)', fontsize = 10)
xlim(150, 700); legend(loc = 0, prop={'size':10}); 
title('B', fontsize = 10)

# stimulus to antennas (Bottom)
subplot(gs[1,1])
plot(A1B.t/ms, A1B.I[0]/nA, '-g', lw=1, 
            label = 'Stim to N1 (Bottom)')
plot(A2B.t/ms, A2B.I[0]/nA, '-c', lw=1, 
            label = 'Stim to N2 (Bottom)')
plot(A3B.t/ms, A3B.I[0]/nA, '-r', lw=1, 
            label = 'Stim to N3 (Bottom)')
plot(A4B.t/ms, A4B.I[0]/nA, '-m', lw=1, 
            label = 'Stim to N4 (Bottom)')
xlabel('time (msec)', fontsize = 10)
ylabel('Stimulus (nA)', fontsize = 10)
xlim(150, 700); ylim(0,11)
legend(loc = 0, prop={'size':10}) 
title('D', fontsize = 10)

print time.clock() - t0, "seconds process time"
show()