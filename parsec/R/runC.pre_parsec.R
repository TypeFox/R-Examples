runC.pre_parsec <-
function(l) {
    polbool <- -as.double(l$inequality < 0)
    v <- .C("bd",
            linext = as.integer(l$linext - 1),
            n_scal = as.integer(l$n),
            nit_scal = as.integer(l$nit),
            zeta = as.integer(t(l$zeta)),
            rankfreq = rep(0L, l$n^2), # in questo modo il conteggio
                                       # parte sempre da 0, poi in A
            threshold = as.integer(l$threshold),
            thrfreq = rep(0L, l$n),
            loweqthr = rep(0L, l$n),
            weights = as.double(l$weights),
            distances = as.double(t(l$distances)),
            cumdist = as.double(rep(0, l$n)),
            gapAP = as.double(rep(0, l$n)),
            gapRP = as.double(rep(0, l$n)),
            gapAR = as.double(rep(0, l$n)),
            gapRR = as.double(rep(0, l$n)),
            polarization_scal = polbool
    )
    # aggiorno solo quello che serve
    l$linext <- v$linext + 1

    # A: il conteggio viene incrementato in R
    l$rankfreq <- l$rankfreq + matrix(v$rankfreq, l$n, l$n, byrow = TRUE)
    l$thrfreq <- l$thrfreq + v$thrfreq
    l$loweqthr <- l$loweqthr + v$loweqthr
    l$gapAP <- l$gapAP + v$gapAP
    l$gapRP <- l$gapRP + v$gapRP
    l$gapAR <- l$gapAR + v$gapAR
    l$gapRR <- l$gapRR + v$gapRR
    # la polarizzazione richiede molta memoria perche' coinvolge iterazioni
    # dell'ordine di N^2/2, pertanto do la possibilita' di non calcolarla
    # nel caso la si calcoli la media la faccio alla fine
    if (polbool >= 0)
        l$inequality <- l$inequality + v$polarization_scal * l$nit
    return(l)
}

runC.pre_parsec_simp <-
function(l) {
    #polbool <- -as.double(l$inequality < 0)
    v <- .C("bd_simp",
            linext = as.integer(l$linext - 1),
            n_scal = as.integer(l$n),
            nit_scal = as.integer(l$nit),
            zeta = as.integer(t(l$zeta)),
            rankfreq = rep(0L, l$n^2), # in questo modo il conteggio
            # parte sempre da 0, poi in A
            threshold = as.integer(l$threshold),
            thrfreq = rep(0L, l$n),
            loweqthr = rep(0L, l$n),
            weights = as.double(l$weights)#,
            #distances = as.double(t(l$distances)),
            #cumdist = as.double(rep(0, l$n)),
            #gapAP = as.double(rep(0, l$n)),
            #gapRP = as.double(rep(0, l$n)),
            #gapAR = as.double(rep(0, l$n)),
            #gapRR = as.double(rep(0, l$n)),
            #polarization_scal = polbool
    )
    # aggiorno solo quello che serve
    l$linext <- v$linext + 1
    
    # A: il conteggio viene incrementato in R
    l$rankfreq <- l$rankfreq + matrix(v$rankfreq, l$n, l$n, byrow = TRUE)
    l$thrfreq <- l$thrfreq + v$thrfreq
    l$loweqthr <- l$loweqthr + v$loweqthr
    #l$gapAP <- l$gapAP + v$gapAP
    #l$gapRP <- l$gapRP + v$gapRP
    #l$gapAR <- l$gapAR + v$gapAR
    #l$gapRR <- l$gapRR + v$gapRR
    # la polarizzazione richiede molta memoria perche' coinvolge iterazioni
    # dell'ordine di N^2/2, pertanto do la possibilita' di non calcolarla
    # nel caso la si calcoli la media la faccio alla fine
    #if (polbool >= 0)
    #    l$inequality <- l$inequality + v$polarization_scal * l$nit
    return(l)
}
