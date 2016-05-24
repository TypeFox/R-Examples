"baranyi"	<- as.formula(LOG10N ~ LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0))))

"baranyi_without_Nmax"	<- as.formula(LOG10N ~ LOG10N0 + mumax * t / log(10) + log10(exp(-mumax * t) * (1 - exp(-mumax * lag)) + exp(-mumax * lag)))

"baranyi_without_lag"	<- as.formula(LOG10N ~ (LOG10Nmax - log10(1 + (10^(LOG10Nmax - LOG10N0) - 1) * exp(-mumax * t))))

"buchanan"	<- as.formula(LOG10N ~ LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) * log(10) / mumax)) * mumax * (t - lag) / log(10) + (t >= lag) * (t>(lag + (LOG10Nmax - LOG10N0) * log(10) / mumax)) * (LOG10Nmax - LOG10N0))

"buchanan_without_Nmax"	<- as.formula(LOG10N ~ (t <= lag) * LOG10N0 + (t > lag) * (LOG10N0 + mumax / log(10) * (t - lag)))

"buchanan_without_lag"	<- as.formula(LOG10N ~ LOG10N0 + (t <= ((LOG10Nmax - LOG10N0) * log(10) / mumax)) * mumax * t / log(10) + (t > ((LOG10Nmax - LOG10N0) * log(10) / mumax)) * (LOG10Nmax - LOG10N0))

"gompertzm"	<- as.formula(LOG10N ~ LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * (lag - t) / ((LOG10Nmax - LOG10N0) * log(10)) + 1)))
