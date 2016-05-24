"geeraerd" <- as.formula(LOG10N ~ LOG10Nres + log10((10^(LOG10N0 - LOG10Nres) - 1) * exp(kmax * Sl) / (exp(kmax * t)+(exp(kmax * Sl) - 1)) + 1))

"geeraerd_without_Nres" <- as.formula(LOG10N ~ LOG10N0 - kmax * t / log(10) + log10(exp(kmax * Sl) / (1 + (exp(kmax * Sl) - 1) * exp(-kmax * t))))

"geeraerd_without_Sl" <- as.formula(LOG10N ~ (LOG10Nres + log10(1 + (10^(LOG10N0 - LOG10Nres) - 1) * exp(-kmax * t))))

"mafart" <- as.formula(LOG10N ~ LOG10N0 - (t / delta)^p)

"albert" <- as.formula(LOG10N ~ LOG10Nres + log10((10^(LOG10N0 - LOG10Nres) - 1) * 10^(-(t / delta)^p) + 1))

"trilinear" <- as.formula(LOG10N ~ LOG10N0 - (t >= Sl) * (t <= (Sl + (LOG10N0 - LOG10Nres) * log(10) / kmax)) * kmax * (t - Sl) / log(10) + (t >= Sl) * (t > (Sl + (LOG10N0 - LOG10Nres) * log(10) / kmax)) * (LOG10Nres - LOG10N0))

"bilinear_without_Nres" <- as.formula(LOG10N ~ (t <= Sl) * LOG10N0 + (t > Sl) * (LOG10N0 - kmax / log(10) * (t - Sl)))

"bilinear_without_Sl" <- as.formula(LOG10N ~ LOG10N0 - (t <= ((LOG10N0 - LOG10Nres) * log(10) / kmax)) * kmax * t / log(10) + (t > ((LOG10N0 - LOG10Nres) * log(10) / kmax)) * (LOG10Nres - LOG10N0))
