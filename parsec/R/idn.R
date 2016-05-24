idn <-
function(
    profiles = NULL,
    threshold, # non puo' assolutamente mancare perche' su di essa si basa tutta la metodologia
    error = 10^(-3),
    zeta = getzeta(profiles),
    weights = {
        if(!is.null(profiles))
            profiles$freq
        else
            rep(1, nrow(zeta))
    },
    #distances = {n <- nrow(zeta); matrix(1, n, n) - diag(1, n)},
    linext = lingen(zeta),
    nit = floor({n <- nrow(zeta); n^5*log(n)+n^4*log(error^(-1))}),
    maxint = 2^31-1#,
    #inequality = TRUE
)
{
    
    n <- nrow(zeta)
    
    
    if (is.numeric(threshold)) {
        stopifnot(max(threshold) <= n)
        stopifnot(min(threshold) >= 1)
        threshold <- 1:n %in% threshold
    }
    if (is.character(threshold)) {
        oldl <- length(threshold)
        threshold <- rownames(zeta) %in% threshold
        if (sum(threshold) != oldl)
            stop("not all thresold profiles can be found in the poset")
    }
    stopifnot(is.logical(threshold))
    
    lev <- levels.incidence(zeta)[threshold]
    if (any(lev==1))
        stop(paste("The element/s", paste(names(which(lev==1)), collapse=", "),
            "of the threshold can define all profiles are poor"))
    
    # frammenta le esecuzioni in modo tale da non passare a C numeri interi
    # pi? grandi di maxint
    nitot <- nit
    nit <- rep(maxint, nitot %/% maxint)
    resto <- nitot %% maxint
    if (resto > 0)
        nit <- c(nit, resto)
    
    pb <- txtProgressBar(style = 3, min = 0, max = nitot)
    cont <- 0
    
    l <- list(
        zeta = zeta,
        linext = linext,
        n = n,
        nit = 0,
        rankfreq = matrix(0, n, n, dimnames=list(rownames(zeta), 1:n)),
        threshold = threshold,
        thrfreq = rep(0, n),
        loweqthr = rep(0, n),
        weights = weights#,
        #distances = distances,
        #gapAP = rep(0, n),
        #gapRP = rep(0, n),
        #gapAR = rep(0, n),
        #gapRR = rep(0, n),
        #inequality = -(!inequality)
    )
    class(l) <- "pre_parsec_simp"
    
    for(j in nit) {
        
        l$nit <- j
        l <- runC(l)
        
        cont <- cont + j
        setTxtProgressBar(pb, cont)
        
    }
    
    close(pb)

    #l$gapAP <- l$gapAP/nitot
    #l$gapRP <- l$gapRP/nitot
    #l$gapAR <- l$gapAR/nitot
    #l$gapRR <- l$gapRR/nitot

    names(l$threshold) <- names(l$thrfreq) <-
    names(l$loweqthr) <- #names(l$gapAP) <- names(l$gapRP) <-
    #names(l$gapAR) <- names(l$gapRR) <-
    #rownames(l$distances) <- colnames(l$distances) <-
        rownames(zeta)
    
    names(l$linext) <- rownames(zeta)[l$linext][l$linext]
    
    l$nit <- nitot
    N <- sum(l$weights)
    #if (inequality) {
    #    maxpolar <- N^2/4*(n - 1)
    #    l$inequality <- l$inequality/nitot/maxpolar
    #} else {
    #    l$inequality <- NA
    #}
    
    #########################
    # CREAZIONE DELL'OUTPUT #
    #########################
    
    res <- list(
        profiles = profiles,
        number_of_profiles = l$n,
        number_of_variables = ncol(profiles$profiles),
        incidence = l$zeta,
        cover = incidence2cover(l$zeta),
        threshold = l$threshold,
        number_of_iterations = l$nit,
        rank_dist = l$rankfreq/l$nit,
        thr_dist = l$thrfreq/l$nit,
        prof_w = l$weights,
        edg_w = l$distances,
        idn_f = l$loweqthr/l$nit,
        svr_abs = NA, #l$gapAP,
        svr_rel = NA, #l$gapRP,
        wea_abs = NA, #l$gapAR,
        wea_rel = NA, #l$gapRR,
        poverty_gap = NA, #weighted.mean(l$gapRP, l$weights),
        wealth_gap = NA, #weighted.mean(l$gapRR, l$weights),
        inequality = NA #l$inequality
    )
    
    class(res) <- "parsec"
    
    return(res)
}
