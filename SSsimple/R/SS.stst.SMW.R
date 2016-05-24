SS.stst.SMW <-
function (F, H, Q, inv.R, P0, epsilon, verbosity = 0) {
	R <- inv.R ; rm("inv.R")  ######################### R IS INVERSE !!!!!!!!
    d <- ncol(H)
    n <- nrow(H)
    params <- .internal.chk.mod.params(F, H, Q, R, P0 = P0, d = d, n = n) ######################### R IS INVERSE !!!!!!!!
    F <- params$F
    H <- params$H
    Q <- params$Q
    R <- params$R ######################### R IS INVERSE !!!!!!!!
    P0 <- params$P0
    P.apos <- P0
    P.apri <- P0
    I <- diag(1, d)
    tHR <- crossprod(H, R)
    RH <- R %*% H
    tHRH <- tHR %*% H
    j <- 0
    stst <- FALSE
    while (!stst) {
        j <- j + 1
        if (verbosity > 0) {
            cat(paste("iteration: ", j), "\n")
        }
        P.apri.temp <- F %*% P.apos %*% t(F) + Q
		
        delta.apri <- sum( diag(P.apri - P.apri.temp)^2 )
		
		P.apri <- P.apri.temp
        
        ######## exucution order is critical !!!!!!!!!!!!!!!!!!!!!!
        SMW.HPHR <- R - RH %*% solve( solve( P.apri, tol=0 ) + tHRH,  tol=0 ) %*% tHR ###### R IS INVERSE !!!
		
        K <- P.apri %*% t(H) %*% SMW.HPHR
		
        P.apos.temp <- (I - K %*% H) %*% P.apri
		
        delta.apos <- sum( diag(P.apos - P.apos.temp)^2 )
		
        P.apos <- P.apos.temp
        if (delta.apri < epsilon & delta.apos < epsilon) {
            stst <- TRUE
        }
        if (verbosity > 1) {
            print(P.apri) ; print(P.apos)
            print(delta.apri) ; print(delta.apos)
        }
    }
    cat(paste("coverged at iteration:", j), "\n")
    return(list("P.apri" = P.apri, "P.apos" = P.apos, "tconverge"= j))
}

