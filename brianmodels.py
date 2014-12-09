# Equations for different Brian2 models. Can be combined by using string concatenation.

eqs_LIF = '''
dv/dt = (-g_l*(v-E_l) + I_syn_e + I_syn_i + I_rand + I_inject) / C : volt (unless refractory)
C : farad
g_l : siemens
E_l : volt
I_inject : amp
v_t : volt
v_r : volt
'''

eqs_ML = '''
dv/dt = (-g_l*(v-E_l) - g_Na*m_inf*(v-E_Na) - g_K*w*(v-E_K) + I_adapt + I_syn_e + I_syn_i + I_rand + I_inject) / C : volt
I_adapt = - g_adapt*z*(v-E_K) : amp
dw/dt = phi*(w_inf-w) / tau_w : 1
dz/dt = alpha * ((1/(1+exp((beta-v)/gamma))) - z) : 1
m_inf = 0.5 * (1 + tanh( (v-V_1)/V_2 )) : 1
w_inf = 0.5 * (1 + tanh( (v-V_3)/V_4 )) : 1
tau_w = ms / cosh( (v-V_3)/(2*V_4) ) : second
C : farad
g_l : siemens
E_l : volt
g_Na : siemens
E_Na : volt
g_K : siemens
E_K : volt
g_adapt : siemens
I_inject : amp
phi : 1
alpha : hertz
beta : volt
gamma : volt
V_1 : volt
V_2 : volt
V_3 : volt
V_4 : volt
'''

eqs_HHLS = '''
dv/dt = (-g_l*(v-E_l) - g_Na*(m*m*m)*h*(v-E_Na) - g_K*(n*n*n*n)*(v-E_K) + I_syn_e + I_syn_i + I_rand + I_inject) / C : volt
dm/dt = alpha_m*(1-m)-beta_m*m : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
alpha_m = ( (25.*mV-(v-E_l))/(10.*mV) ) / (exp((25.*mV-(v-E_l))/(10.*mV)) - 1.) / ms : Hz
beta_m = 4.*exp(-((v-E_l))/(18.*mV)) / ms : Hz
alpha_h = (0.07*exp( -(v-E_l)/(20*mV) )) / ms : Hz
beta_h = (1 / (exp( (30.*mV-(v-E_l))/(10.*mV) ) + 1.)) / ms : Hz
alpha_n = ( (10.*mV-(v-E_l))/(100.*mV) ) / (exp( (10.*mV-(v-E_l))/(10.*mV) ) - 1.) / ms : Hz
beta_n = 0.125*exp(-((v-E_l)) / (80*mV)) / ms : Hz
C : farad
g_l : siemens
E_l : volt
g_Na : siemens
E_Na : volt
g_K : siemens
E_K : volt
I_inject : amp
'''

eqs_I_exp = '''
dI_syn_e/dt = -I_syn_e/tau_syn_e : siemens
dI_syn_i/dt = -I_syn_i/tau_syn_i : siemens
tau_syn_e : second
tau_syn_i : second
'''

eqs_I_exp_refrac = '''
dI_syn_e/dt = -I_syn_e/tau_syn_e : siemens (unless refractory)
dI_syn_i/dt = -I_syn_i/tau_syn_i : siemens (unless refractory)
tau_syn_e : second
tau_syn_i : second
'''

eqs_g_exp = '''
dg_e/dt = -g_e/tau_syn_e : siemens
dg_i/dt = -g_i/tau_syn_i : siemens
I_syn_e = g_e*(E_syn_e - v) : amp
I_syn_i = g_i*(E_syn_i - v) : amp
tau_syn_e : second
tau_syn_i : second
E_syn_e : volt
E_syn_i : volt
'''

eqs_g_exp_refrac = '''
dg_e/dt = -g_e/tau_syn_e : siemens (unless refractory)
dg_i/dt = -g_i/tau_syn_i : siemens (unless refractory)
I_syn_e = g_e*(E_syn_e - v) : amp
I_syn_i = g_i*(E_syn_i - v) : amp
tau_syn_e : second
tau_syn_i : second
E_syn_e : volt
E_syn_i : volt
'''

eqs_I_rand = '''
I_rand = I_mean + I_std*zeta : amp
dzeta/dt = -(zeta/tau_rand) + ((2/tau_rand)**0.5) * xi : 1
I_mean : amp
I_std : amp
tau_rand : second
'''

eqs_I_rand_refrac = '''
I_rand = I_mean + I_std*zeta : amp
dzeta/dt = -(zeta/tau_rand) + ((2/tau_rand)**0.5) * xi : 1 (unless refractory)
I_mean : amp
I_std : amp
tau_rand : second
'''

eqs_g_rand = '''
I_rand = g_rand*(E_rand - v) : amp
g_rand = g_mean + g_std*zeta : siemens
dzeta/dt = -(zeta/tau_rand) + ((2/tau_rand)**0.5) * xi : 1
g_mean : siemens
g_std : siemens
tau_rand : second
E_rand : volt
'''

eqs_g_rand_refrac = '''
I_rand = g_rand*(E_rand - v) : amp
g_rand = g_mean + g_std*zeta : siemens
dzeta/dt = -(zeta/tau_rand) + ((2/tau_rand)**0.5) * xi : 1 (unless refractory)
g_mean : siemens
g_std : siemens
tau_rand : second
E_rand : volt
'''

eqs_g_rand_ei = '''
I_rand = g_rand_e*(E_rand_e - v) + g_rand_i*(E_rand_i - v) : amp
g_rand_e = g_mean_e + g_std_e*zeta_e : siemens
g_rand_i = g_mean_i + g_std_i*zeta_i : siemens
dzeta_e/dt = -(zeta_e/tau_rand_e) + ((2/tau_rand_e)**0.5) * xi : 1
dzeta_i/dt = -(zeta_i/tau_rand_i) + ((2/tau_rand_i)**0.5) * xi : 1
g_mean_e : siemens
g_mean_i : siemens
g_std_e : siemens
g_std_i : siemens
tau_rand_e : second
tau_rand_i : second
E_rand_e : volt
E_rand_i : volt
'''

eqs_g_rand_ei_refrac = '''
I_rand = g_rand_e*(E_rand_e - v) + g_rand_i*(E_rand_i - v) : amp
g_rand_e = g_mean_e + g_std_e*zeta_e : siemens
g_rand_i = g_mean_i + g_std_i*zeta_i : siemens
dzeta_e/dt = -(zeta_e/tau_rand_e) + ((2/tau_rand_e)**0.5) * xi : 1 (unless refractory)
dzeta_i/dt = -(zeta_i/tau_rand_i) + ((2/tau_rand_i)**0.5) * xi : 1 (unless refractory)
g_mean_e : siemens
g_mean_i : siemens
g_std_e : siemens
g_std_i : siemens
tau_rand_e : second
tau_rand_i : second
E_rand_e : volt
E_rand_i : volt
'''

eqs_threshold = 'not_refractory and (v > -5.*mV)'

eqs_refractory = 'v > -5.*mV'