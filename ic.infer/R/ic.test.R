ic.test <- function (obj, TP = 1, s2 = 1, df.error = Inf, ui0.11 = diag(rep(1, 
    length(obj$b.restr))), ci0.11 = NULL, meq.alt = 0, df = NULL, 
    wt = NULL, tol = sqrt(.Machine$double.eps), ...) 
{
    if (!"orest" %in% class(obj)) 
        stop("obj must be of class orest, e.g. result of ic.est")
    cov <- obj$Sigma
    ## make this case equal to fixed sigma case (later multiplied again)
    if ("orlm" %in% class(obj)) { 
        df.error <- obj$df.error 
        s2 <- obj$s2
        cov <- obj$Sigma/s2    ## calculation of test statistic requires
                               ## unscaled covariance matrix 
        }
    ui <- obj$ui
    ci <- obj$ci
    if (obj$meq == nrow(ui)) 
        stop("ic.test not applicable for obj with equality restrictions only.")
    b.unrestr <- obj$b.unrestr
    b.restr <- obj$b.restr
    b.alt <- NULL
    index <- obj$restr.index
    g <- length(b.restr)
    meq <- obj$meq
    ## prevent long runs with late aborts because of too low memory
    if (nrow(ui) - meq - 2 > 2){
      if (!is.numeric(try(matrix(0, floor((nrow(ui) - meq - 2)/2), 
        choose(nrow(ui) - meq, floor((nrow(ui) - meq - 2)/2))), silent = TRUE))) 
        stop(paste("ic.test does not work, too many inequality restrictions in obj, \n", 
            "interim matrix with ", floor((nrow(ui) - meq)/2) * 
                choose(nrow(ui) - meq, floor((nrow(ui) - meq)/2)), 
            "elements cannot be created", sep = ""))
    }
    ## initialize for all TPs
    ## will later be modified for TP11
    ui.extra <- NULL
    b0 <- NULL
    ## initialize for all TPs
    ## will later be modified for TP1
    b.eqrestr <- NULL

    ### check inputs
    if (!(is.matrix(cov) & nrow(cov) == ncol(cov))) 
        stop("cov must be a square matrix.")
    if (!(all(eigen(cov)$value > 0))) 
        stop("cov must be positive definite.")
    if (!length(b.unrestr) == length(b.restr)) 
        stop("b.unrestr and b.restr are of different length.")
    if (!length(b.unrestr) == nrow(cov)) 
        stop("b.unrestr and nrow(cov) must be identical.")
    if (!length(b.unrestr) == length(b.restr)) 
        stop("b.unrestr and b.restr are of different length.")
    if (s2 <= 0) 
        stop("test not possible with non-positive s2.")
    if (df.error == Inf) 
        cov <- cov * s2
    if (TP == 11) {
        if (!is.null(ci0.11)) 
            if (any(!ci0.11 == 0)) 
                stop("ci0.11 must be a vector of 0es.")
        if (is.vector(ui0.11)) 
            ui0.11 <- matrix(ui0.11, 1, length(ui0.11))
        if (is.null(ci0.11)) 
            ci0.11 <- rep(0, nrow(ui0.11))
        if (!ncol(ui0.11) == g) 
            stop(paste("ui0.11 must have", g, "columns."))
        if (!nrow(ui0.11) == length(ci0.11)) 
            stop(paste("mismatch between ui0.11 and ci0.11"))
    }

    ## check for full row rank
    ## show linear independent subset of rows of ui, if violated
    hilf <- RREF(t(ui), tol = tol)
    if (hilf$rank < nrow(ui)) 
        stop(paste("Matrix ui must have full row-rank (choose e.g. rows", 
            paste(hilf$pivot, collapse = " "), ")."))
    if (TP == 21 & meq.alt == 0) 
        TP <- 2
    if (TP == 21 & meq.alt > meq) 
        stop("meq.alt must not be larger than obj$meq.")
    
    ## initialize weights and df
    wt.bar <- NULL
    df.bar <- NULL

    ## expand ui with 0 columns to match dimension of b.unrestr
    uiw <- ui
    if (length(index) < nrow(cov)) {
        uiw <- matrix(0, nrow(ui), nrow(cov))
        uiw[, index] <- ui
    }
    if (TP == 11) {
        ## correct user errors / laziness in specifying matrices
        ui.extra <- ui0.11 %*% (diag(rep(1, g)) - t(uiw) %*% 
            solve(uiw %*% t(uiw), uiw))
        if (all(abs(ui.extra) < tol)) 
            TP <- 1
        else {
            ### eliminate 0 rows
            ui.extra <- ui.extra[!rowSums(abs(ui.extra) < tol) == 
                g, ]
            if (is.vector(ui.extra)) 
                ui.extra <- matrix(ui.extra, 1, length(ui.extra))
            else {
                hilf <- RREF(t(ui.extra), tol = tol)
                if (hilf$rank == 1) 
                  ui.extra <- matrix(ui.extra[1, ], 1, ncol(ui.extra))
                else {
                  ## extract linear independent subset of rows of ui.extra
                  if (hilf$rank < nrow(ui.extra)) {
                    ui.extra <- ui.extra[hilf$pivot, ]
                    if (is.vector(ui.extra)) 
                      ui.extra <- matrix(ui.extra, 1, length(ui.extra))
                  }
                }
            }
            ui0.11 <- rbind(ui.extra, uiw)
            ci0.11 <- c(rep(0, nrow(ui.extra)), ci)
        }
    }
    if (TP == 3 & obj$meq > 0) 
        stop("TP3 not applicable in case of equality restrictions.")
    if (!TP == 3) 
        if (!all(uiw %*% b.restr - ci >= -tol)) 
            stop("b.restr must fulfill restriction ui%*%b.restr>=ci component wise.")
        ## rounding error prevented

    ### determine weights
    ### different TPs are accomodated by different (orders of) degrees of freedom
    if (TP %in% c(1, 2, 11, 21)) {
        ## no equality restrictions in original set of restrictions: 
        ## --> weights are defined by (5.5) in Shapiro 
        ## --> or - in case of TP 11 - by (5.10) which leads to the same weights
        ##
        ## equality restrictions in original set of restrictions:
        ## --> weights are defined by (5.9) in Shapiro 
        if (is.null(wt)) {
            if (obj$meq == 0) 
                wt.bar <- ic.weights(uiw %*% cov %*% t(uiw), 
                  ...)
            else wt.bar <- ic.weights(solve(solve(uiw %*% cov %*% 
                t(uiw))[-(1:meq), -(1:meq)]), ...)
        }
        else wt.bar <- wt
    }
    ### restriction equal vs >= without =
    if (TP == 1) {
        b.eqrestr = solve.QP(Dmat = solve(cov), dvec = solve(cov, 
            b.unrestr), Amat = t(uiw), bvec = ci, meq = nrow(ui))$solution
        b.eqrestr[abs(b.eqrestr) < tol] <- 0
        names(b.eqrestr) <- names(b.restr)
        T <- t(b.restr - b.eqrestr) %*% solve(cov, b.restr - 
            b.eqrestr)
        if (is.null(df)) 
            df.bar <- (nrow(ui) - obj$meq):0
        else df.bar <- df
        if (df.error == Inf) 
            p.value <- 1 - pchibar(T, df = df.bar, wt = wt.bar)
        else {
            T <- T/(df.error * s2 + T)
            p.value <- 1 - pbetabar(T, df1 = df.bar/2, df2 = df.error/2, 
                wt = wt.bar)
        }
    }
    ### ui0.11%*%beta==0 vs >= without =
    ### only reasonable if ui0.11%*%beta==0 fulfills Restriction ui%*%beta==ci with "="
    ### this is guaranteed by sanity check
    if (TP == 11) {
        b0 <- ic.est(b.unrestr, ui = ui0.11, ci = ci0.11, Sigma = cov, 
            meq = length(ci0.11))$b.restr
          ## it is assumed that ui0.11 has full row rank and extends restrictions 
          ##     on ui; this has been ascertained by projecting uiw out, 
          ##     reducing rest to linearly independent rows and appending it afterwards!
        if (!all(abs(uiw %*% b0 - ci) < tol)) 
            stop("TP11 is not implemented for situations \n for which b0 violates restriction", 
                 "\n           ui%*%b0 == ci in some components.")
        T <- t(b.restr - b0) %*% solve(cov, b.restr - b0)
        ### degrees of freedom: maximum is dimension of space containing simply 
        ### restricted b.restr if all inequality restrictions are inactive, 
        ### minimum is dimension of this space if all inequality restrictions are 
        ### active
        
        ### corrected August 1 2014: changed length(b.restr) to length(ci0.11)
        if (is.null(df)) 
            df.bar <- (length(ci0.11) - obj$meq):(length(ci0.11) - 
                nrow(uiw))
        else df.bar <- df
        if (df.error == Inf) 
            p.value <- 1 - pchibar(T, df = df.bar, wt = wt.bar)
        else {
            T <- T/(df.error * s2 + T)
            ### T excluding s2, i.e. not vcov but cov.unscaled
            p.value <- 1 - pbetabar(T, df1 = df.bar/2, df2 = df.error/2, 
                wt = wt.bar)
        }
    }
    ### restricted vs. not
    if (TP == 2) {
        T <- t(b.unrestr - b.restr) %*% solve(cov, b.unrestr - 
            b.restr)
        if (is.null(df)) 
            df.bar <- meq:nrow(ui)
        else df.bar <- df
        if (df.error == Inf) 
            p.value <- 1 - pchibar(T, df = df.bar, wt = wt.bar)
        else {
            T <- T/(df.error * s2 + T)
            p.value <- 1 - pbetabar(T, df1 = df.bar/2, df2 = df.error/2, 
                wt = wt.bar)
        }
    }
    ### restricted vs. not, with some equality restrictions preserved
    if (TP == 21) {
        b.alt <- ic.est(b.unrestr, ui = uiw[1:meq.alt, ], ci = ci[1:meq.alt], 
            Sigma = cov, meq = meq.alt)$b.restr
        T <- t(b.restr - b.alt) %*% solve(cov, b.restr - b.alt)
        ### degrees of freedom: maximum is dimension of space containing simply 
        ### restricted b.restr if all inequality restrictions are inactive, 
        ### minimum is dimension of this space if all inequality restrictions are 
        ### active
        if (is.null(df)) 
            df.bar <- (meq - meq.alt):(nrow(uiw) - meq.alt)
        else df.bar <- df
        if (df.error == Inf) 
            p.value <- 1 - pchibar(T, df = df.bar, wt = wt.bar)
        else {
            T <- T/(df.error * s2 + T)
            ### T excluding s2, i.e. not vcov but cov.unscaled
            p.value <- 1 - pbetabar(T, df1 = df.bar/2, df2 = df.error/2, 
                wt = wt.bar)
        }
    }
    ### everything else vs. restricted
    if (TP == 3) {
        T <- min((uiw %*% b.unrestr - ci)/sqrt(diag(uiw %*% cov %*% 
            t(uiw))))
        if (df.error == Inf) 
            p.value <- 1 - pnorm(T)
        else {
            T <- T/sqrt(s2)
            ### T excluding s2, i.e. not vcov but cov.unscaled
            p.value <- 1 - pt(T, df.error)
        }
    }
    ### output results
    aus <- list(TP = TP, b.unrestr = b.unrestr, b.restr = b.restr, 
        ui = ui, ci = ci, restr.index = index, meq = obj$meq, 
        iact = obj$iact, ui.extra = ui.extra, b.eqrestr = b.eqrestr, 
        b.extra.restr = b0, b.alt = b.alt, T = T, p.value = p.value, 
        s2 = s2, cov = cov, df.error = df.error, df.bar = df.bar, 
        wt.bar = wt.bar)
    class(aus) <- "ict"
    aus
}
