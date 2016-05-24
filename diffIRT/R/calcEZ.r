calcEZ = function(Pc, VRT, MRT, N=100)
{
    # adapted from Wagenmakers, Van der Maas, & Grasman (2007) JMP 
    
    s=s2=1
    # If Pc equals 0, .5, or 1, an edge-correction is conducted, see Wagenmakers et al (2007)
    if (Pc == 1) Pc = Pc - 1 / (2*N)
    if (Pc == 0.5) Pc = Pc + 1 / (2*N)
    if (Pc == 0) Pc = Pc + 1 / (2*N)

    L = qlogis(Pc)
    x = L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
    v = sign(Pc-0.5)*s*x^(1/4)
    a = s2*qlogis(Pc)/v
    y   = -v*a/s2
    MDT = (a/(2*v))*(1-exp(y))/(1+exp(y))
    Ter = MRT-MDT
    return(list(v, a, Ter))
}

