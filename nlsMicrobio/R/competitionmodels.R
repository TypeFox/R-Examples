"jameson_buchanan" <- as.formula("
LOG10N~(flora==1)*((t<=lag_1)*LOG10N0_1+
((t>lag_1)&(t<tmax))*(LOG10N0_1 + mumax_1/log(10)*(t-lag_1))+
(t>=tmax)*(LOG10N0_1 +mumax_1/log(10)*(tmax-lag_1)))+
(flora==2)*((t<=lag_2)*LOG10N0_2+
((t>lag_2)&(t<tmax))*(LOG10N0_2 + mumax_2/log(10)*(t-lag_2))+
(t>=tmax)*(LOG10N0_2 +mumax_2/log(10)*(tmax-lag_2)))")

"jameson_baranyi" <- as.formula("
LOG10N ~ (flora==1)* ( (t<=tmax)*(LOG10N0_1 + mumax_1 * t/log(10) + log10(exp(-mumax_1 * t) * 
    (1 - exp(-mumax_1 * lag_1)) + exp(-mumax_1 * lag_1))) +
(t>tmax) * (LOG10N0_1 + mumax_1 * tmax/log(10) + log10(exp(-mumax_1 * tmax) * 
    (1 - exp(-mumax_1 * lag_1)) + exp(-mumax_1 * lag_1))) ) + 
(flora==2)* ( (t<=tmax)*(LOG10N0_2 + mumax_2 * t/log(10) + log10(exp(-mumax_2 * t) * 
    (1 - exp(-mumax_2 * lag_2)) + exp(-mumax_2 * lag_2))) +
(t>tmax) * (LOG10N0_2 + mumax_2 * tmax/log(10) + log10(exp(-mumax_2 * tmax) * 
    (1 - exp(-mumax_2 * lag_2)) + exp(-mumax_2 * lag_2))) ) ")

"jameson_without_lag" <- as.formula("
LOG10N~(flora==1)*( (t<tmax)*(LOG10N0_1 + mumax_1/log(10)*t)+
(t>=tmax)*(LOG10N0_1 +mumax_1/log(10)*tmax))+
(flora==2)*( (t<tmax)*(LOG10N0_2 + mumax_2/log(10)*t)+
(t>=tmax)*(LOG10N0_2 +mumax_2/log(10)*tmax))")
