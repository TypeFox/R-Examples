planck = function(wave,temp) {
    bbflux = wave*0.
    w = wave / 1.e8                              # angstroms to cm
    c1 =  3.7417749e-5                # =2*!dpi*h*c*c
    c2 =  1.4387687                  # =h*c/k
    val =  c2/w/temp
    good = ( val< .Machine$double.max.exp)    #avoid floating underflow
    bbflux[ good ] =  c1 / ( w[good]^5 * ( exp( val[good])-1. ) )
    return (bbflux*1.e-8)              # convert to ergs/cm2/s/a
}
