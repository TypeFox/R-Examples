sersic = function(mag, re, n, e = 0, r = re){
    bn = qgamma(0.5,2*n)
    lumtot = 1*(re^2)*2*pi*n*((exp(bn))/(bn^(2*n)))*gamma(2*n)*(1-e)
    magtot = -2.5*log10(lumtot)
    Ie = 1/(10^(0.4*(mag-magtot)))
    x = bn*(r/re)^(1/n)
    lumr = Ie*lumtot*pgamma(x,2*n)
    a = r
    b = a*(1-e)
    intenr = Ie*exp(-bn*(((r/re)^(1/n))-1))
    lumtot = Ie*lumtot
    magtot = -2.5*log10(lumtot)
    magr = -2.5*log10(lumr)
    mur = -2.5*log10(intenr)
    muavgr = -2.5*log10(lumr/(pi*a*b))
    return(list(mag=magr, magdiff=magtot-magr, mu=mur, muavg=muavgr, inten=intenr, lum=lumr, lumtot=lumtot, lumfrac=lumr/lumtot))
}

