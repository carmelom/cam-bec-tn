dynamic_expressions = (
    '(N_Na-N_K)/N_K',
    'N_Na+N_K',
    'Nbec_Na/N_Na',
    'Nbec_K/N_K',
    'm1x_K-m1x_Na',
    'mx_K-mx_Na',
    '(N0_K-N1_K)/N0_K',
    'N_Na *1e-3* (1 + (2*20/9.7946)**2)',
    "np.sum(np.array(Timestamp.split(':'), dtype=float)*np.array([1,60,3600]))",
    "0.5e-6*(sx_K**2 + sy_K**2)",     # avg s^2 [mm]
    "2.78e-9*(2*np.pi*8.87*sx_K)**2", #Temperature in situ from sx [uK]
    "2.78e-9*(2*np.pi*86.7*sy_K)**2",  #Temperature in situ from sy [uK]
	"np.sqrt(((8.87*sx_K)**2 - (85.45*sy_K)**2)/(8.87**2 - 85.45**2))", # sigma blur
    "float(re.search('t_exp = ([0-9]+\.[0-9]+)ms', params_K).groups()[0])",
    )
