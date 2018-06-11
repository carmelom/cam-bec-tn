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
    )
