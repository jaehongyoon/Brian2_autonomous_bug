# import libraries/packages
from brian2 import *
import matplotlib.pyplot as plt

#----------------------------------------------------
# set map size (given in skeleton code)
map_size = 100
# set global variables
global poisonx, poisony, poison_plot, poison_count
global bug_plot, sr_plot, sl_plot
global foodx, foody, food_count, food_plot
global mate_plot, sr_mate_plot, sl_mate_plot, mate_count
global enemy_plot, sr_enemy_plot, sl_enemy_plot, enemy_count
global death_count

#----------------------------------------------------
#set coordinates for four objects
# 1. poison
poisonx = 50
poisony = 50

# 2.food
foodx = -50
foody = -50
food_count = 0

# 3.others
mate_count = 0
death_count = 0
poison_count = 0
enemy_count = 0

#----------------------------------------------------
# set bug sensor
# Sensor neurons = Izhichevich model
a = 0.02; b = 0.2;
c = -65; d = 0.5

# The virtual bug
taum = 4*ms
base_speed = 9.5
turn_rate = 5*Hz

# NMDA (mate) neurotransmitter variables
Mg = 2 # unit mM. For ease of calculation, we do not give unit for this
       # found in Chapter 8 of "Foundations of Neuroscience" p. 163
betaN = 0.0066/ms
alphaN = 0.072/ms
v_th = 0*mV # reversal potential
transwdth = 1*ms

I0 = 100
g_synpk = 0.12*nsiemens # peak synaptic conductance
# maximum synaptic conductance
g_synmax = g_synpk/(taum/ms*exp(-1))
Cm=1.0*ufarad/cm**2
area=20000*umetre**2

# GABA_b synaptic neuotransimtters
K3 = 0.18/ms
K4 = 0.034/ms
Kd = 100 # found in Chapter 8 of "Foundations of Neuroscience" in p.164
betab = 0.0012/ms
alphab = 0.09/ms
n = 4 # power of GABAb s
adjust = 50

#----------------------------------------------------
# bug motor neurons
# 1/2. poison/food = alpha synapse
# 3. mate = NMDA synapse
# 4. enemy = GABAb type excitatory synapse
bug_eqs = '''

# motor neurons
motorr=  (20*g_ampa2+10*g_ampa3+10*g_ampa5+20*g_ampa8)/nS: 1
motorl = (20*g_ampa1+10*g_ampa4+10*g_ampa6+20*g_ampa7)/nS: 1

# speed
speed = (motorl + motorr)/2 + base_speed : 1

# angle and displacement calculation
dangle/dt = (motorr - motorl)*turn_rate : 1
dx/dt = speed*cos(angle)*15*Hz : 1
dy/dt = speed*sin(angle)*15*Hz : 1

#----------------------------------------------------
# for I1 (right side sensors)
# poison right side (g_ampa1)
dg_ampa1/dt=-g_ampa1/taum + z1/ms: siemens
dz1/dt = -z1/taum: siemens

# food right side (g_ampa3)
dg_ampa3/dt=-g_ampa3/taum + z3/ms: siemens
dz3/dt = -z3/taum: siemens

# mate right side (g_ampa5)
g_ampa5 = 8000*g*B5*sN5*bb5 : siemens
dsN5/dt = alphaN*(1-sN5)*Tr_pre5 - betaN*sN5 : 1
B5 = 1/(1+exp(-vN5/(16.13*mV))*(Mg/3.57)) : 1
dvN5/dt = g_ampa5*(vN5+60*mV)/(Cm*area) :  volt

Tr_pre5=10*(tanh((t/ms-tspike5/ms)/.005)-
            tanh((t/ms-(tspike5/ms +transwdth/ms))/.005)):1
tspike5:second
bb5:1

# enemy right side (g_ampa7)
g_ampa7 = (g_synmax*2000)*(sb7**n)/(sb7**n + Kd) : siemens
dr7/dt = (alphab*Tr_pre_GABA7*(1-r7) - betab*r7) : 1
dsb7/dt = (K3*r7 - K4*sb7) : 1

Tr_pre_GABA7=aa7*(tanh((t/ms-tspikeN7/ms)/(.005*adjust))-
            tanh((t/ms-(tspikeN7/ms +transwdth/
            (adjust*ms)))/(.005*adjust))):1
tspikeN7:second
aa7:1

#----------------------------------------------------
# for I2 (left side sensors)
# poison left side (g_ampa2)
dg_ampa2/dt=-g_ampa2/taum + z2/ms: siemens
dz2/dt = -z2/taum: siemens

# food left side (g_ampa4)
dg_ampa4/dt=-g_ampa4/taum + z4/ms: siemens
dz4/dt = -z4/taum: siemens

# mate left side (g_ampa6)
g_ampa6 = 8000*g*B6*sN6*bb6 : siemens
dsN6/dt = alphaN*(1-sN6)*Tr_pre6 - betaN*sN6 : 1
B6 = 1/(1+exp(-vN6/(16.13*mV))*(Mg/3.57)) : 1
dvN6/dt = g_ampa6*(vN6+60*mV)/(Cm*area) :  volt

Tr_pre6=10*(tanh((t/ms-tspike6/ms)/.005)-
            tanh((t/ms-(tspike6/ms +transwdth/ms))/.005)):1
tspike6:second
bb6:1
g: siemens

# enemy left side (g_ampa8)
g_ampa8 = (g_synmax*2000)*(sb8**n)/(sb8**n + Kd) : siemens
dr8/dt = (alphab*Tr_pre_GABA8*(1-r8) - betab*r8) : 1
dsb8/dt = (K3*r8 - K4*sb8) : 1

Tr_pre_GABA8=aa7*(tanh((t/ms-tspikeN8/ms)/(.005*adjust))-
            tanh((t/ms-(tspikeN8/ms +transwdth/
            (adjust*ms)))/(.005*adjust))):1
tspikeN8:second
aa8:1

'''

bug = NeuronGroup(1, bug_eqs, clock=Clock(0.2*ms))
bug.angle = pi/2
bug.x = 0
bug.y = 0
bug.g = g_synmax

#----------------------------------------------------
# mate bug equation and neuron
mate_eqs = '''

# motor neurons
motorr=  (10*g_ampa5 + 20*g_ampa8)/nS: 1
motorl = (10*g_ampa6 + 20*g_ampa7)/nS: 1

# speed
speed = (motorl + motorr)/2 + base_speed : 1

# angle and displacement calculation
dangle/dt = (motorr - motorl)*turn_rate : 1
dx/dt = speed*cos(angle)*15*Hz : 1
dy/dt = speed*sin(angle)*15*Hz : 1

# mate right side (g_ampa5)
g_ampa5 = 8000*g*B5*sN5*bb5 : siemens
dsN5/dt = alphaN*(1-sN5)*Tr_pre5 - betaN*sN5 : 1
B5 = 1/(1+exp(-vN5/(16.13*mV))*(Mg/3.57)) : 1
dvN5/dt = g_ampa5*(vN5+60*mV)/(Cm*area) :  volt

Tr_pre5=10*(tanh((t/ms-tspike5/ms)/.005)-
            tanh((t/ms-(tspike5/ms +transwdth/ms))/.005)):1
tspike5:second
bb5:1

# enemy right side (g_ampa7)
g_ampa7 = (g_synmax*2000)*(sb7**n)/(sb7**n + Kd) : siemens
dr7/dt = (alphab*Tr_pre_GABA7*(1-r7) - betab*r7) : 1
dsb7/dt = (K3*r7 - K4*sb7) : 1

Tr_pre_GABA7=aa7*(tanh((t/ms-tspikeN7/ms)/(.005*adjust))-
            tanh((t/ms-(tspikeN7/ms +transwdth/
            (adjust*ms)))/(.005*adjust))):1
tspikeN7:second
aa7:1

# mate left side (g_ampa6)
g_ampa6 = 8000*g*B6*sN6*bb6 : siemens
dsN6/dt = alphaN*(1-sN6)*Tr_pre6 - betaN*sN6 : 1
B6 = 1/(1+exp(-vN6/(16.13*mV))*(Mg/3.57)) : 1
dvN6/dt = g_ampa6*(vN6+60*mV)/(Cm*area) :  volt

Tr_pre6=10*(tanh((t/ms-tspike6/ms)/.005)-
            tanh((t/ms-(tspike6/ms +transwdth/ms))/.005)):1
tspike6:second
bb6:1
g:siemens

g_ampa8 = (g_synmax*2000)*(sb8**n)/(sb8**n + Kd) : siemens
dr8/dt = (alphab*Tr_pre_GABA8*(1-r8) - betab*r8) : 1
dsb8/dt = (K3*r8 - K4*sb8) : 1

Tr_pre_GABA8=aa7*(tanh((t/ms-tspikeN8/ms)/(.005*adjust))-
            tanh((t/ms-(tspikeN8/ms +transwdth/
            (adjust*ms)))/(.005*adjust))):1
tspikeN8:second
aa8:1

'''

mate = NeuronGroup(1, mate_eqs, clock=Clock(0.2*ms))
mate.angle = pi/2
mate.x = -75
mate.y = 75
mate.g = g_synmax

#----------------------------------------------------
# enemy bug
# regard female and bug as both food
enemy_eqs = '''

# motor neurons
motorr=  10*(g_ampa1+g_ampa3)/nS: 1
motorl = 10*(g_ampa2+g_ampa4)/nS: 1

# speed
speed = (motorl + motorr)/2 + base_speed : 1

# angle and displacement calculation
dangle/dt = (motorr - motorl)*turn_rate : 1
dx/dt = speed*cos(angle)*15*Hz : 1
dy/dt = speed*sin(angle)*15*Hz : 1

#----------------------------------------------------
# for I1 (right side sensors)

# food1 = bug right side (g_ampa3)
dg_ampa3/dt=-g_ampa3/taum + z3/ms: siemens
dz3/dt = -z3/taum: siemens

# food2 = female right side (g_ampa1)
dg_ampa1/dt=-g_ampa1/taum + z1/ms: siemens
dz1/dt = -z1/taum: siemens

#----------------------------------------------------
# for I2 (left side sensors)

# food1 = bug left side (g_ampa4)
dg_ampa4/dt=-g_ampa4/taum + z4/ms: siemens
dz4/dt = -z4/taum: siemens

# food2 = female left side (g_ampa2)
dg_ampa2/dt=-g_ampa2/taum + z2/ms: siemens
dz2/dt = -z2/taum: siemens

'''
enemy = NeuronGroup(1, enemy_eqs, clock=Clock(0.2*ms))
enemy.angle = 3*pi/2
enemy.x = -100
enemy.y = -100

#----------------------------------------------------
# sensor (antenna) = Izhikevich
# 1. poison
sensor_eqs = '''
# Izhikevich + alpha function
dv/dt = (0.04*(v)*(v) + 5*(v) + 140- u)/ms
        + (I/mV)/(Cm*area):  1
du/dt = a*(b*v - u)/ms : 1

# stimulus from the objects to sensor
I = I0*nA/
    (sqrt((x-poisonx)**2 + (y-poisony)**2)): amp
    
# bug position
x : 1
y : 1
x_disp : 1
y_disp : 1

# poison coordinates
poisonx : 1
poisony : 1
'''

# 2. food
sensor_eqs2 = '''
# Izhikevich + alpha function
dv/dt = (0.04*(v)*(v) + 5*(v) + 140- u)/ms
        + (I/mV)/(Cm*area):  1
du/dt = a*(b*v - u)/ms : 1

# stimulus from the objects to sensor
I = I0*nA/
    (sqrt((x-foodx)**2 + (y-foody)**2)): amp
    
# bug position
x : 1
y : 1
x_disp : 1
y_disp : 1

# poison coordinates
foodx : 1
foody : 1
'''

# 3. mate
sensor_eqs3 = '''
# Izhikevich + alpha function
dv/dt = (0.04*(v)*(v) + 5*(v) + 140- u)/ms
        + (I/mV)/(Cm*area):  1
du/dt = a*(b*v - u)/ms : 1

# stimulus from the objects to sensor
I = 3*I0*nA/
    sqrt((x-matex)**2 + (y-matey)**2): amp
    
# bug position
x : 1
y : 1
x_disp : 1
y_disp : 1

# poison coordinates
matex : 1
matey : 1
'''

# 4. Enemy (predator)
sensor_eqs4 = '''
# Izhikevich + alpha function
dv/dt = (0.04*(v)*(v) + 5*(v) + 140- u)/ms
        + (I/mV)/(Cm*area):  1
du/dt = a*(b*v - u)/ms : 1

# stimulus from the objects to sensor
I = I0*nA/
    sqrt((x-enemyx)**2 + (y-enemyy)**2): amp
    
# bug position
x : 1
y : 1
x_disp : 1
y_disp : 1

# poison coordinates
enemyx : 1
enemyy : 1
'''

#----------------------------------------------------
# for each neuron group you will also need to 
# initialize parameters associated with membrane model, 
# and assign thresholds and resets as needed- sveral 
# neurons can be in each group

# set neuron for right side sensor
# 1. poison
sr_poison = NeuronGroup(1, sensor_eqs, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sr_poison.x_disp = 5; sr_poison.v = c;
sr_poison.y_disp = 5
sr_poison.x = sr_poison.x_disp
sr_poison.y = sr_poison.y_disp
sr_poison.poisonx = poisonx
sr_poison.poisony = poisony
# 2. food
sr_food = NeuronGroup(1, sensor_eqs2, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sr_food.x_disp = 5; sr_food.v = c;
sr_food.y_disp = 5
sr_food.x = sr_food.x_disp
sr_food.y = sr_food.y_disp
sr_food.foodx = foodx
sr_food.foody = foody
# 3. mate
sr_mate = NeuronGroup(1, sensor_eqs3, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sr_mate.x_disp = 10; sr_mate.v = c;
sr_mate.y_disp = 5
sr_mate.x = sr_mate.x_disp
sr_mate.y = sr_mate.y_disp
sr_mate.matex = mate.x
sr_mate.matey = mate.y

# 4. enemy
sr_enemy = NeuronGroup(1, sensor_eqs4, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sr_enemy.x_disp = 10; sr_enemy.v = c;
sr_enemy.y_disp = 5
sr_enemy.x = sr_enemy.x_disp
sr_enemy.y = sr_enemy.y_disp
sr_enemy.enemyx = enemy.x
sr_enemy.enemyy = enemy.y


# set neuron for left side sensor
# 1. poison
sl_poison = NeuronGroup(1, sensor_eqs, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sl_poison.x_disp = -5
sl_poison.y_disp = 5
sl_poison.x = sl_poison.x_disp; sl_poison.v = c;
sl_poison.y = sl_poison.y_disp
sl_poison.poisonx = poisonx
sl_poison.poisony = poisony
# 2. food
sl_food = NeuronGroup(1, sensor_eqs2, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sl_food.x_disp = -5
sl_food.y_disp = 5
sl_food.x = sl_food.x_disp; sl_food.v = c;
sl_food.y = sl_food.y_disp
sl_food.foodx = foodx
sl_food.foody = foody

# 3. mate
sl_mate = NeuronGroup(1, sensor_eqs3, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sl_mate.x_disp = -10
sl_mate.y_disp = 5
sl_mate.x = sl_mate.x_disp; sl_mate.v = c;
sl_mate.y = sl_mate.y_disp
sl_mate.matex = mate.x
sl_mate.matey = mate.y

# 4. enemy
sl_enemy = NeuronGroup(1, sensor_eqs4, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sl_enemy.x_disp = -10
sl_enemy.y_disp = 5
sl_enemy.x = sl_enemy.x_disp; sl_enemy.v = c;
sl_enemy.y = sl_enemy.y_disp
sl_enemy.enemyx = enemy.x
sl_enemy.enemyy = enemy.y


#----------------------------------------------------
# the mate bug
# only consist of enemy and mate

sr_m = NeuronGroup(1, sensor_eqs3, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sr_m.x_disp = 10; sr_m.v = c;
sr_m.y_disp = 5
sr_m.x = sr_m.x_disp
sr_m.y = sr_m.y_disp
sr_m.matex = bug.x
sr_m.matey = bug.y

sr_m_enemy = NeuronGroup(1, sensor_eqs4, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sr_m_enemy.x_disp = 10; sr_m_enemy.v = c;
sr_m_enemy.y_disp = 5
sr_m_enemy.x = sr_m_enemy.x_disp
sr_m_enemy.y = sr_m_enemy.y_disp
sr_m_enemy.enemyx = enemy.x
sr_m_enemy.enemyy = enemy.y

sl_m = NeuronGroup(1, sensor_eqs3, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sl_m.x_disp = -10
sl_m.y_disp = 5
sl_m.x = sl_m.x_disp; sl_m.v = c;
sl_m.y = sl_m.y_disp
sl_m.matex = bug.x
sl_m.matey = bug.y

sl_m_enemy = NeuronGroup(1, sensor_eqs4, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sl_m_enemy.x_disp = -10
sl_m_enemy.y_disp = 5
sl_m_enemy.x = sl_m_enemy.x_disp; sl_m_enemy.v = c;
sl_m_enemy.y = sl_m_enemy.y_disp
sl_m_enemy.enemyx = enemy.x
sl_m_enemy.enemyy = enemy.y

#----------------------------------------------------
# the enemy
# regard both bug (food1) and female (food2) as food

# right side
# food 1 = bug
sr_enemy_food1 = NeuronGroup(1, sensor_eqs2, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sr_enemy_food1.x_disp = 5; sr_enemy_food1.v = c;
sr_enemy_food1.y_disp = 5
sr_enemy_food1.x = sr_enemy_food1.x_disp
sr_enemy_food1.y = sr_enemy_food1.y_disp
sr_enemy_food1.foodx = bug.x
sr_enemy_food1.foody = bug.y
# food 2 = female
sr_enemy_food2 = NeuronGroup(1, sensor_eqs2, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sr_enemy_food2.x_disp = 5; sr_enemy_food2.v = c;
sr_enemy_food2.y_disp = 5
sr_enemy_food2.x = sr_enemy_food2.x_disp
sr_enemy_food2.y = sr_enemy_food2.y_disp
sr_enemy_food2.foodx = mate.x
sr_enemy_food2.foody = mate.y

# left side
# food 1 = bug
sl_enemy_food1 = NeuronGroup(1, sensor_eqs2, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sl_enemy_food1.x_disp = -5; sl_enemy_food1.v = c;
sl_enemy_food1.y_disp = 5
sl_enemy_food1.x = sl_enemy_food1.x_disp
sl_enemy_food1.y = sl_enemy_food1.y_disp
sl_enemy_food1.foodx = bug.x
sl_enemy_food1.foody = bug.y
# food 2 = female
sl_enemy_food2 = NeuronGroup(1, sensor_eqs2, 
                    threshold='v > 0',
                    reset = 'v = c; u = u + d',
                    refractory='v > -60',
                    clock=Clock(defaultclock.dt),
                    method='euler')
sl_enemy_food2.x_disp = -5; sl_enemy_food2.v = c;
sl_enemy_food2.y_disp = 5
sl_enemy_food2.x = sl_enemy_food2.x_disp
sl_enemy_food2.y = sl_enemy_food2.y_disp
sl_enemy_food2.foodx = mate.x
sl_enemy_food2.foody = mate.y

#----------------------------------------------------
# Synapses (sensors communicate with bug motor)
# right side
# 1. poison
# 2. food
# 3. mate
# 4. enemy
syn_r = Synapses(sr_poison, bug, 
                clock=sr_poison.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa2 += w
		''')
syn_r.connect(i=[0],j=[0])
syn_r.w = g_synmax/(defaultclock.dt/ms)*10

syn_rf = Synapses(sr_food, bug, 
                clock=sr_food.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa4 += w
		''')
syn_rf.connect(i=[0],j=[0])
syn_rf.w = g_synmax/(defaultclock.dt/ms)*10

syn_rm = Synapses(sr_mate, bug, 
                clock=sr_mate.clock,
		pre='''
                tspike6 =t
                bb6 = 1
                vN6 = -65*mV 	
		''')
syn_rm.connect(i=[0],j=[0])

syn_re = Synapses(sr_enemy, bug, clock=sr_enemy.clock,
		pre='''
                tspikeN = t
                aa8 = 200		
		''')
syn_re.connect(i=[0],j=[0])


# left side
# 1. poison
# 2. food
# 3. mate
syn_l = Synapses(sl_poison, bug, 
                clock=sl_poison.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa1 += w
		''')
syn_l.connect(i=[0],j=[0])
syn_l.w = g_synmax/(defaultclock.dt/ms)*10

syn_lf = Synapses(sl_food, bug, 
                clock=sr_food.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa3 += w
		''')
syn_lf.connect(i=[0],j=[0])
syn_lf.w = g_synmax/(defaultclock.dt/ms)*10

syn_lm = Synapses(sl_mate, bug, 
                clock=sr_mate.clock,
		pre='''
                tspike5 =t
                bb5 = 1
                vN5 = -65*mV
    		
		''')
syn_lm.connect(i=[0],j=[0])

syn_le = Synapses(sl_enemy, bug, clock=sl_enemy.clock,
		pre='''
                tspikeN = t
                aa7 = 200		
		''')
syn_le.connect(i=[0],j=[0])

#------------------------------------------------------------------
# Mate bug synapse
syn_mate_rm = Synapses(sr_m, mate, 
                clock=sr_mate.clock,
		pre='''
                tspike6 =t
                bb6 = 1
                vN6 = -65*mV 		
		''')
syn_mate_rm.connect(i=[0],j=[0])

syn_re = Synapses(sr_enemy, mate, clock=sr_enemy.clock,
		pre='''
                tspikeN = t
                aa8 = 200		
		''')
syn_re.connect(i=[0],j=[0])

syn_mate_lm = Synapses(sl_m, mate, 
                clock=sr_mate.clock,
		pre='''
                tspike5 =t
                bb5 = 1
                vN5 = -65*mV		
		''')
syn_mate_lm.connect(i=[0],j=[0])

syn_mate_le = Synapses(sl_enemy, mate, 
                clock=sl_enemy.clock,
		pre='''
                tspikeN = t
                aa7 = 200		
		''')
syn_mate_le.connect(i=[0],j=[0])

#------------------------------------------------------------------
# enemy bug synapse
syn_enemy_rf1 = Synapses(sr_enemy_food1, enemy, 
                clock=sr_enemy_food1.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa4 += w
		''')
syn_enemy_rf1.connect(i=[0],j=[0])
syn_enemy_rf1.w = g_synmax/(defaultclock.dt/ms)*8

syn_enemy_rf2 = Synapses(sr_enemy_food2, enemy, 
                clock=sr_enemy_food2.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa2 += w
		''')
syn_enemy_rf2.connect(i=[0],j=[0])
syn_enemy_rf2.w = g_synmax/(defaultclock.dt/ms)*5

syn_enemy_lf1 = Synapses(sl_enemy_food1, enemy, 
                clock=sl_enemy_food1.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa3 += w
		''')
syn_enemy_lf1.connect(i=[0],j=[0])
syn_enemy_lf1.w = g_synmax/(defaultclock.dt/ms)*8

syn_enemy_lf2 = Synapses(sl_enemy_food2, enemy, 
                clock=sl_enemy_food2.clock,model='''
        	w : siemens
        	''',
		pre='''
		g_ampa1 += w
		''')
syn_enemy_lf2.connect(i=[0],j=[0])
syn_enemy_lf2.w = g_synmax/(defaultclock.dt/ms)*5

#-------------------------------------------------------------------
# plot

f = figure(1)
bug_plot = plot(bug.x, bug.y, 'ko')
poison_plot = plot(poisonx, poisony, 'b^')
food_plot = plot(foodx, foody, 'r^')
mate_plot = plot(mate.x, mate.y, 'yo')
enemy_plot = plot(enemy.x, enemy.y, 'ro')

sr_plot = plot([0], [0], 'w')   # Just leaving it blank for now
sl_plot = plot([0], [0], 'w')
sr_mate_plot = plot([0], [0], 'w')
sl_mate_plot = plot([0], [0], 'w')
sr_enemy_plot = plot([0], [0], 'w')
sl_enemy_plot = plot([0], [0], 'w')

# array of x and y
x_plot = []; y_plot = []
x_mate_plot = []; y_mate_plot = []
x_enemy_plot = []; y_enemy_plot = []

# This block of code updates the position of the bug and 
# poison and makes the appropriate rotations
@network_operation()
def update_positions():
    global poisonx, poisony, poison_count
    global foodx, foody, food_count
    global matex, matey, mate_count
    global death_count, enemy_count


# <bug command>
# 1. poison position  
    sr_poison.x = bug.x + sr_poison.x_disp*cos(
                    bug.angle-pi/2
                    ) - sr_poison.y_disp*sin(
                    bug.angle-pi/2)
    sr_poison.y = bug.y + sr_poison.x_disp*sin(
                    bug.angle-pi/2
                    ) + sr_poison.y_disp*cos(
                    bug.angle-pi/2)

    sl_poison.x = bug.x + sl_poison.x_disp*cos(
                    bug.angle-pi/2
                    ) - sl_poison.y_disp*sin(
                    bug.angle-pi/2)
    sl_poison.y = bug.y + sl_poison.x_disp*sin(
                    bug.angle-pi/2
                    ) + sl_poison.y_disp*cos(
                    bug.angle-pi/2)

# 2. food position     
    sr_food.x = bug.x + sr_food.x_disp*cos(
                    bug.angle-pi/2
                    ) - sr_food.y_disp*sin(
                    bug.angle-pi/2)
    sr_food.y = bug.y + sr_food.x_disp*sin(
                    bug.angle-pi/2
                    ) + sr_food.y_disp*cos(
                    bug.angle-pi/2)

    sl_food.x = bug.x + sl_food.x_disp*cos(
                    bug.angle-pi/2
                    ) - sl_food.y_disp*sin(
                    bug.angle-pi/2)
    sl_food.y = bug.y + sl_food.x_disp*sin(
                    bug.angle-pi/2
                    ) + sl_food.y_disp*cos(
                    bug.angle-pi/2)

# 3. mate sensor position
    sr_mate.x = bug.x + sr_mate.x_disp*cos(
                    bug.angle-pi/2
                    ) - sr_mate.y_disp*sin(
                    bug.angle-pi/2)
    sr_mate.y = bug.y + sr_mate.x_disp*sin(
                    bug.angle-pi/2
                    ) + sr_mate.y_disp*cos(
                    bug.angle-pi/2)

    sl_mate.x = bug.x + sl_mate.x_disp*cos(
                    bug.angle-pi/2
                    ) - sl_mate.y_disp*sin(
                    bug.angle-pi/2)
    sl_mate.y = bug.y + sl_mate.x_disp*sin(
                    bug.angle-pi/2
                    ) + sl_mate.y_disp*cos(
                    bug.angle-pi/2)
                    
# 4. enemy sensor position
    sr_enemy.x = bug.x + sr_enemy.x_disp*cos(
                    bug.angle-pi/2
                    ) - sr_enemy.y_disp*sin(
                    bug.angle-pi/2)
    sr_enemy.y = bug.y + sr_enemy.x_disp*sin(
                    bug.angle-pi/2
                    ) + sr_enemy.y_disp*cos(
                    bug.angle-pi/2)

    sl_enemy.x = bug.x + sl_enemy.x_disp*cos(
                    bug.angle-pi/2
                    ) - sl_enemy.y_disp*sin(
                    bug.angle-pi/2)
    sl_enemy.y = bug.y + sl_enemy.x_disp*sin(
                    bug.angle-pi/2
                    ) + sl_enemy.y_disp*cos(
                    bug.angle-pi/2)                                
                    
                    
# < mate command>                                        
# 1. Mate position
    sr_m.x = mate.x + sr_m.x_disp*cos(
                    mate.angle-pi/2
                    ) - sr_m.y_disp*sin(
                    mate.angle-pi/2)
    sr_m.y = mate.y + sr_m.x_disp*sin(
                    mate.angle-pi/2
                    ) + sr_m.y_disp*cos(
                    mate.angle-pi/2)

    sl_m.x = mate.x + sl_m.x_disp*cos(
                    mate.angle-pi/2
                    ) - sl_m.y_disp*sin(
                    mate.angle-pi/2)
    sl_m.y = mate.y + sl_m.x_disp*sin(
                    mate.angle-pi/2
                    ) + sl_m.y_disp*cos(
                    mate.angle-pi/2)

# 2. enemy position
    sr_m_enemy.x = mate.x + sr_m_enemy.x_disp*cos(
                    mate.angle-pi/2
                    ) - sr_m_enemy.y_disp*sin(
                    mate.angle-pi/2)
    sr_m_enemy.y = mate.y + sr_m_enemy.x_disp*sin(
                    mate.angle-pi/2
                    ) + sr_m_enemy.y_disp*cos(
                    mate.angle-pi/2)

    sl_m_enemy.x = mate.x + sl_m_enemy.x_disp*cos(
                    mate.angle-pi/2
                    ) - sl_m_enemy.y_disp*sin(
                    mate.angle-pi/2)
    sl_m_enemy.y = mate.y + sl_m_enemy.x_disp*sin(
                    mate.angle-pi/2
                    ) + sl_m_enemy.y_disp*cos(
                    mate.angle-pi/2)                                

# <enemy command>                                        
    sr_enemy_food1.x = enemy.x + sr_enemy_food1.x_disp*cos(
                    enemy.angle-pi/2
                    ) - sr_enemy_food1.y_disp*sin(
                    enemy.angle-pi/2)
    sr_enemy_food1.y = enemy.y + sr_enemy_food1.x_disp*sin(
                    enemy.angle-pi/2
                    ) + sr_enemy_food1.y_disp*cos(
                    enemy.angle-pi/2)

    sl_enemy_food1.x = enemy.x + sl_enemy_food1.x_disp*cos(
                    enemy.angle-pi/2
                    ) - sl_enemy_food1.y_disp*sin(
                    enemy.angle-pi/2)
    sl_enemy_food1.y = enemy.y + sl_enemy_food1.x_disp*sin(
                    enemy.angle-pi/2
                    ) + sl_enemy_food1.y_disp*cos(
                    enemy.angle-pi/2)                    

    sr_enemy_food2.x = enemy.x + sr_enemy_food2.x_disp*cos(
                    enemy.angle-pi/2
                    ) - sr_enemy_food2.y_disp*sin(
                    enemy.angle-pi/2)
    sr_enemy_food2.y = enemy.y + sr_enemy_food2.x_disp*sin(
                    enemy.angle-pi/2
                    ) + sr_enemy_food2.y_disp*cos(
                    enemy.angle-pi/2)

    sl_enemy_food2.x = enemy.x + sl_enemy_food2.x_disp*cos(
                    enemy.angle-pi/2
                    ) - sl_enemy_food2.y_disp*sin(
                    enemy.angle-pi/2)
    sl_enemy_food2.y = enemy.y + sl_enemy_food2.x_disp*sin(
                    enemy.angle-pi/2
                    ) + sl_enemy_food2.y_disp*cos(
                    enemy.angle-pi/2) 
                                                                                                                            
# re-locate objects when encountered
                                                                                                                                                                                                                                                                                        
    if ((bug.x-poisonx)**2+(bug.y-poisony)**2) < 10:
        print 'The bug ate poison and died'
        # refresh NMDA learning behavior
        bug.g = g_synmax
        
        # count death and poison
        death_count += 1
        poison_count += 1
       
        # relocate poison
        relocatex = True
        relocatey = True
        while relocatex:             
        # relocate poison
            poisonx_candidate = randint(-map_size+10, map_size-10)
            if (abs(poisonx_candidate-foodx) > 30 and
    	       abs(poisonx_candidate-bug.x)>30):
    	       poisonx = poisonx_candidate
    	       relocatex = False
    	       print 'relocated poison in x-direction'
    	while relocatey:
    	    poisony_candidate = randint(-map_size+10, map_size-10)
    	    if (abs(poisony_candidate-foody) > 30 and
	       abs(poisony_candidate-bug.y)>30):
	       poisony = poisony_candidate
	       relocatey = False
	       print 'relocated poison in y-direction'

    if ((mate.x-enemy.x)**2+(mate.y-enemy.y)**2) < 10:
        print 'The predator ate the mate'
        # refresh NMDA learning behavior
        mate.g = g_synmax
        
        # count enemy catch number
        enemy_count += 1
        
        # relocate mate
        relocatex = True
        relocatey = True
        while relocatex:
            matex_candidate = randint(-map_size+10, map_size-10)
            if (abs(matex_candidate-poisonx) > 30 and 
                abs(matex_candidate-enemy.x)>30):
                mate.x = matex_candidate
                relocatex = False
                print 'relocated mate in x-direction'
                
       	
       	while relocatey:
       	    matey_candidate = randint(-map_size+10, map_size-10)  
   	    if (abs(matey_candidate-poisony) > 30 and
    	       abs(matey_candidate-enemy.y)>30):
   	        mate.y = matey_candidate
   	        relocatey = False
   	        print 'relocated mate in y-direction'
	    
    if ((bug.x-enemy.x)**2+(bug.y-enemy.y)**2) < 10:
        print 'The predator ate the bug'
        # refresh NMDA learning behavior
        bug.g = g_synmax
        
        # count death and enemy catch number
        death_count += 1
        enemy_count += 1
        
        # relocate mate
        relocatex = True
        relocatey = True
        while relocatex:
            bugx_candidate = randint(-map_size+10, map_size-10)
            if (abs(bugx_candidate-poisonx) > 30 and 
	       abs(bugx_candidate-enemy.x) > 30):
	       bug.x = bugx_candidate
	       relocatex = False
	       print 'relocated bug in x-direction'
	       
	while relocatey:
	    bugy_candidate = randint(-map_size+10, map_size-10)
	    if (abs(bugy_candidate-poisony) > 30 and 
	       abs(bugy_candidate-enemy.y) > 30):
	       bug.y = bugy_candidate
	       relocatey = False
	       print 'relocated bug in y-direction'      
	    
    if ((bug.x-foodx)**2+(bug.y-foody)**2) < 10:
        food_count += 1
	
	relocatex = True
	relocatey = True
	while relocatex:
	    foodx_candidate = randint(-map_size+10, map_size-10)
	    if (abs(foodx_candidate-poisonx) > 30):
	        foodx = foodx_candidate
	        relocatex = False
	        
	while relocatey:
	    foody_candidate = randint(-map_size+10, map_size-10)
	    if (abs(foody_candidate-poisony) > 30):
	        foody = foody_candidate
	        relocatey = False
        
    if ((bug.x-mate.x)**2+(bug.y-mate.y)**2) < 10:
        mate_count += 1
        relocatex = True
        relocatey = True
        
        bug.g += g_synmax*0.1 #LTP learning
        mate.g += g_synmax*0.1 #LTP learning
        
        while relocatex:
            matex_candidate = randint(-map_size+10, map_size-10)
            if (abs(matex_candidate-enemy.x) > 30):
                mate.x = matex_candidate
                relocatex = False
        
        while relocatey:
            matey_candidate = randint(-map_size+10, map_size-10)
            if (abs(matey_candidate-enemy.y) > 30):
                mate.y = matey_candidate 
                relocatey = False 	    

# rebound when hit margin
    if (bug.x < -map_size):
        bug.x = -map_size
        bug.angle = pi - bug.angle
    if (bug.x > map_size):
	bug.x = map_size
	bug.angle = pi - bug.angle
    if (bug.y < -map_size):
	bug.y = -map_size
	bug.angle = -bug.angle
    if (bug.y > map_size):
	bug.y = map_size
	bug.angle = -bug.angle
	
    if (mate.x < -map_size):
        mate.x = -map_size
        mate.angle = pi - mate.angle
    if (mate.x > map_size):
	mate.x = map_size
	mate.angle = pi - mate.angle
    if (mate.y < -map_size):
	mate.y = -map_size
	mate.angle = -mate.angle
    if (mate.y > map_size):
	mate.y = map_size
	mate.angle = -mate.angle

    if (enemy.x < -map_size):
        enemy.x = -map_size
        enemy.angle = pi - enemy.angle
    if (enemy.x > map_size):
	enemy.x = map_size
	enemy.angle = pi - enemy.angle
    if (enemy.y < -map_size):
	enemy.y = -map_size
	enemy.angle = -enemy.angle
    if (enemy.y > map_size):
	enemy.y = map_size
	enemy.angle = -enemy.angle

    sr_poison.poisonx = poisonx
    sr_poison.poisony = poisony
    sl_poison.poisonx = poisonx
    sl_poison.poisony = poisony
    sr_food.foodx = foodx
    sr_food.foody = foody
    sl_food.foodx = foodx
    sl_food.foody = foody
    sr_mate.matex = mate.x
    sr_mate.matey = mate.y
    sl_mate.matex = mate.x
    sl_mate.matey = mate.y
    sr_enemy.enemyx = enemy.x
    sr_enemy.enemyy = enemy.y
    sl_enemy.enemyx = enemy.x
    sl_enemy.enemyy = enemy.y
    
    sr_m.matex = bug.x
    sr_m.matey = bug.y
    sl_m.matex = bug.x
    sl_m.matey = bug.y    
    sr_m_enemy.enemyx = enemy.x
    sr_m_enemy.enemyy = enemy.y
    sl_m_enemy.enemyx = enemy.x
    sl_m_enemy.enemyy = enemy.y
    
    sr_enemy_food1.foodx = bug.x
    sr_enemy_food1.foody = bug.y
    sl_enemy_food1.foodx = bug.x
    sl_enemy_food1.foody = bug.y
    sr_enemy_food2.foodx = mate.x
    sr_enemy_food2.foody = mate.y
    sl_enemy_food2.foodx = mate.x
    sl_enemy_food2.foody = mate.y
    
# trace the path of enemy, bug, mate    
    x_plot.append(float(bug.x/1.0))
    y_plot.append(float(bug.y/1.0))
    
    x_mate_plot.append(float(mate.x/1.0))
    y_mate_plot.append(float(mate.y/1.0))
    
    x_enemy_plot.append(float(enemy.x/1.0))
    y_enemy_plot.append(float(enemy.y/1.0))

# this block of code updates the plots so 
# you can see the bug and poison move
@network_operation(dt=2*ms)
def update_plot():
    global poisonx, poisony, poison_plot
    global bug_plot, sr_plot, sl_plot
    global mate_plot, sr_mate_plot, sl_mate_plot
    global enemy_plot, sr_enemy_plot, sl_enemy_plot
    global food_plot, foodx, foody
        
    bug_plot[0].remove()
    poison_plot[0].remove()
    food_plot[0].remove()
    mate_plot[0].remove()
    enemy_plot[0].remove()
    sr_plot[0].remove()
    sl_plot[0].remove()
    sr_mate_plot[0].remove()
    sl_mate_plot[0].remove()
    sr_enemy_plot[0].remove()
    sl_enemy_plot[0].remove()
        
    bug_x_coords = [bug.x, bug.x-
                    2*cos(bug.angle), 
                    bug.x-4*cos(bug.angle)]
                    # ant-like body
    bug_y_coords = [bug.y, bug.y-
                    2*sin(bug.angle), 
                    bug.y-4*sin(bug.angle)]
                    
    mate_x_coords = [mate.x, mate.x-
                    2*cos(mate.angle), 
                    mate.x-4*cos(mate.angle)]
                    # ant-like body
    mate_y_coords = [mate.y, mate.y-
                    2*sin(mate.angle), 
                    mate.y-4*sin(mate.angle)]

    enemy_x_coords = [enemy.x, enemy.x-
                    2*cos(enemy.angle), 
                    enemy.x-4*cos(enemy.angle)]
                    # ant-like body
    enemy_y_coords = [enemy.y, enemy.y-
                    2*sin(enemy.angle), 
                    enemy.y-4*sin(enemy.angle)]                    
                                                            
    # Plot the bug's current position
    bug_plot = plot(bug_x_coords, bug_y_coords, 'ko')
    mate_plot = plot(mate_x_coords, mate_y_coords, 'yo')
    enemy_plot = plot(enemy_x_coords, enemy_y_coords, 'ro')
    
    # draws sensors. Since no need to draw all four sensors
    # this part remains unchanged  
    sr_plot = plot([bug.x, sr_poison.x], [bug.y, sr_poison.y], 'b')
    sl_plot = plot([bug.x, sl_poison.x], [bug.y, sl_poison.y], 'r')
    
    sr_mate_plot = plot([mate.x, sr_m.x], [mate.y, sr_m.y], 'b')
    sl_mate_plot = plot([mate.x, sl_m.x], [mate.y, sl_m.y], 'r')
    
    sr_enemy_plot = plot([enemy.x, sr_enemy_food1.x], 
                            [enemy.y, sr_enemy_food1.y], 'b')
    sl_enemy_plot = plot([enemy.x, sl_enemy_food1.x], 
                            [enemy.y, sl_enemy_food1.y], 'r')
    
    
    poison_plot = plot(poisonx, poisony, 'b^')
    food_plot = plot(foodx, foody, 'r^')

    axis([-100,100,-100,100])
    draw()
    #print "."
    pause(0.01)

run(3000*ms,report='text')

plt.clf()
plt.plot(x_plot, y_plot, 'k--')
plt.plot(x_mate_plot,y_mate_plot, 'b--')
plt.plot(poisonx, poisony, 'b^')
plt.plot(foodx, foody, 'b^')

axis([-100,100,-100,100])
title('Path')
show()
