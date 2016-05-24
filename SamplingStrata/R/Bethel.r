#
# ----------------------------------------------------------------------
# Function BETHEL definition Multivariate optimal
# allocation for different domains of interest in
# stratified sample design
# 
# Extension of Bethel methodology with Chromy Algorithm see
# Bethel(1989)'Sample Allocation in Multivarate Surveys' -
# Survey Methodology Author: Daniela Pagliuca
# <pagliuca@istat.it> with contributions of M. Teresa
# Buglielli <bugliell@istat.it> and Giulio Barcaroli
# <barcarol@istat.it>
# 
# ----------------------------------------------------------------------
bethel <- function(stratif, errors, minnumstrat = 2, maxiter = 200, 
    maxiter1 = 25, printa = FALSE, realAllocation = FALSE, epsilon = 1e-11) # Begin body
{
    # First input data frame
    colnames(stratif) <- toupper(colnames(stratif))
    
    # Second input data frame
    colnames(errors) <- toupper(colnames(errors))
    #
    checkData(strata = stratif, errors = errors)
    #
    ordina_variabili <- function(dati, prefisso, n_var) {
        if (!is.data.frame(dati)) 
            stop()
        as.matrix(dati[, paste(prefisso, 1:n_var, sep = ""), 
            drop = FALSE])
    }
    #
    # --------------------------------------------------------------
    # Initialization of parameters. Attribution of initial
    # values
    # ------------------------------------------------------------
    iter1 <- 0
    val <- NULL
    m <- NULL
    s <- NULL
    cv <- NULL
    # ------------------------------ Initialization of
    # parameters ------------------------------
    nstrat <- nrow(stratif)
    nvar <- length(grep("CV", names(errors)))
    ndom <- nrow(errors)
    
    varloop <- c(1:nvar)
    strloop <- c(1:nstrat)
    
    # ---------------------------------------------------------
    # Initial Data Structures - selection from input variables
    # ---------------------------------------------------------
    
    # means
    med <- ordina_variabili(stratif, "M", nvar)
    
    # variances estimates
    esse <- ordina_variabili(stratif, "S", nvar)
    
    if (ncol(med) != ncol(esse)) 
        stop(print("Error: Number of M variables don't match the number of S variables"))
    if (ncol(med) != nvar) 
        stop(print("Error: Number of variables don't match the number of planned CV"))
    # populations
    N <- as.vector(stratif$N)
    # vector cens
    cens <- as.vector(stratif$CENS)
    # strata with N < minnumstrata are set to take-all strata
    cens[N < minnumstrat] <- 1
    # vector cost
    cost <- as.vector(stratif$COST)
    # default for cost and cens
    if (is.null(cost)) 
        cost <- rep(1, nstrat)
    if (is.null(cens)) 
        cens <- rep(0, nstrat)
    nocens <- 1 - cens
    # check variable cens
    if (sum(cens) == length(cens)) {
        warning(print("Warning: Variable CENS always equal 1"))
    }
    # domains
    nom_dom <- sapply(1:ndom, function(i) paste("DOM", i, sep = ""))
    dom <- ordina_variabili(stratif, "DOM", ndom)
    
    # numbers of different domains of interest (numbers of
    # modalities/categories for each type of domain) if (ndom
    # ==1) (nvalues<-nlevels((as.factor(stratif$DOM1)))) if
    # (ndom>1) {nvalues<-sapply(nom_dom, function(vari)
    # {val<-c(val, nlevels(as.factor(dom[,vari])))})}
    nvalues <- sapply(nom_dom, function(vari) {
        val <- c(val, nlevels(as.factor(dom[, vari])))
    })
    
    # -------------------------------------------------
    # disjunctive matrix
    # -------------------------------------------------
    
    crea_disj <- function(data, vars) {
        out <- NULL
        sapply(vars, function(vari) {
            col <- as.factor(data[, vari])
            out <<- cbind(out, outer(col, levels(col), function(y, 
                x) ifelse(y == x, 1, 0)))
        })
        out
    }
    
    disj <- crea_disj(stratif, nom_dom)
    
    # -------------------------------------------------
    # Building means (m) and deviations matrices (s) for
    # different domains of interest
    # ------------------------------------------------- (m) and
    # (s)
    nc <- ncol(disj)
    for (i in 1:nc) {
        m <- cbind(m, disj[, i] * med)
        s <- cbind(s, disj[, i] * esse)
    }
    
    #
    # -------------------------------------------------------------
    # computation of the coefficients of variation CVs for
    # different domains of interest
    # -------------------------------------------------------------
    cvDom <- NULL  #Per stampa
    cvDom2 <- NULL  #per stampa
    for (k in 1:ndom) {
        cvx <- ordina_variabili(errors[k, ], "CV", nvar)
        ndomvalues <- c(1:nvalues[k])
        for (k1 in ndomvalues) {
            cv <- cbind(cv, cvx)
            cvDom <- c(cvDom, rep(nom_dom[k], length(cvx)))
            cvDom2 <- c(cvDom2, rep(levels(as.factor(dom[, k]))[k1], 
                length(cvx)))
        }
    }
    
    #
    # ------------------------------------------------------------
    # New definition of initial values
    # ------------------------------------------------------------
    
    nvar <- ncol(cv)  # new numbers of variables
    varloop <- c(1:nvar)
    
    varfin <- c(rep(0, nvar))
    totm <- c(rep(0, nvar))
    alfa2 <- c(rep(0, nvar))
    #
    # ------------------------------------------------------------
    # Calculation of aij - matrix of standardized precision
    # units
    # ------------------------------------------------------------
    crea_a <- function() {
        numA <- (N^2) * (s^2) * nocens
        
        denA1 <- colSums(t(t(N * m) * c(cv)))^2
        denA2 <- colSums(N * (s^2) * nocens)
        
        denA <- denA1 + denA2 + epsilon
        
        a <- t(t(numA)/denA)
        return(a)
    }
    
    ###
    ### -----------------------------------------------------------
    ### Computation of alfa's values - Chromy Algorithm
    ### Iteration
    ### -----------------------------------------------------------
    
    chromy <- function(alfatot, diff, iter, alfa, alfanext, x) {
        while (diff > epsilon && iter < maxiter) {
            iter <- iter + 1
            
            den1 <- sqrt(rowSums(t(t(a) * c(alfa))))
            den2 <- sum(sqrt(rowSums(t(t(a * cost) * c(alfa)))))
            
            x <- sqrt(cost)/(den1 * den2 + epsilon)
            
            alfatot <- sum(c(alfa) * (t(a) %*% x)^2)
            alfatot[alfatot == 0] <- epsilon  # substitution of zeroes with epsilon
            alfanext <- c(alfa) * (t(a) %*% x)^2/alfatot
            
            diff <- max(abs(alfanext - alfa))
            alfa <- alfanext
            alfa2 <<- alfanext  #For sensitivity
        }
        # Allocation vector
        if (realAllocation == FALSE) 
            n <- ceiling(1/x)
        if (realAllocation == TRUE) 
            n <- 1/x
        return(n)
    }
    #
    # -----------------------------------------------------------
    # Check results
    # -----------------------------------------------------------
    # check n<minnumstr
    check_n <- function() {
        n[n < minnumstrat] <- pmin(minnumstrat, N)[n < minnumstrat]
        n
    }
    
    #
    # -----------------------------------------------------------
    
    a <- crea_a()
    # n<-chromyF()
    n <- chromy(0, 999, 0, c(rep(1/nvar, nvar)), c(rep(0, nvar)), 
        array(0.1, dim = c(nstrat, 1)))
    
    #
    # -----------------------------------------------------------
    # check n>N
    contx <- sum(n > N)
    cens[n > N] <- 1
    nocens <- 1 - cens
    
    n <- check_n()
    
    #### ITERATIONS###
    while (contx > 0 && iter1 < maxiter1) {
        iter1 <- iter1 + 1
        #
        # -----------------------------------------------------------
        # Calculation of aij - matrix of standardized precision
        # units Computation of alfa's values - Chromy Algorithm
        # Iteration
        # -----------------------------------------------------------
        a <- crea_a()
        
        # n <- chromyF()
        n <- chromy(0, 999, 0, c(rep(1/nvar, nvar)), c(rep(0, 
            nvar)), array(0.1, dim = c(nstrat, 1)))
        # check n>N
        contx <- sum(n > N)
        cens[n > N] <- 1
        nocens <- 1 - cens
        
        n <- check_n()
    }
    
    # ------------------------- Definitive best allocation
    # -------------------------
    
    n <- (nocens * n) + (cens * N)
    # dyn.unload('chromy2DLL.dll')
    # -----------------------------------------------------------
    if (printa == TRUE) 
        {
            # --------------------- Printing allocation
            # ---------------------
            
            stampa_confronto <- function(n, N, strato) {
                nomi <- c("STRATUM", "POPULATION", "BETHEL", 
                  "PROPORTIONAL", "EQUAL")
                df <- NULL
                df <- cbind(df, as.character(strato), N, n, ceiling(sum(n) * 
                  N/sum(N)), ceiling(sum(n)/length(n)))
                tot <- apply(matrix(as.numeric(df[, 2:5]), ncol = 4), 
                  2, sum)
                df <- rbind(df, c("TOTAL", tot))
                colnames(df) <- nomi
                df
            }
            calcola_cv <- function() {
                
                NTOT <- c(rep(0, nvar))
                CVfin <- c(rep(0, nvar))
                #
                # -------------------------------------------------------------
                # Populations in strata for different domains of interest
                # NTOTj
                # -------------------------------------------------------------
                NTOT <- colSums((m > 0) * N)
                #
                # ------------------------------------------------------------
                # Computation of the CVs
                # ------------------------------------------------------------
                varfin <- rowSums(t((s * N)^2 * (1 - round(n)/N)/round(n))/NTOT^2)
                totm <- rowSums(t(m * N))
                
                CVfin <- round(sqrt(varfin/(totm/NTOT)^2), digits = 4)
                return(CVfin)
            }
            calcola_sensibilita <- function() {
                #
                # -------------------------------------------------------------
                # Sensitivity (10%)
                # ------------------------------------------------------------
                
                t <- g <- 0
                for (i in 1:nstrat) {
                  t <- sum(alfa2 * a[i, ])
                  g <- g + sqrt(cost[i] * t)
                }
                g <- g^2
                sens <- 2 * 0.1 * alfa2 * g
                return(sens)
            }
            
            # --------------------- Printing CV's ---------------------
            stampa_cv <- function() {
                CVfin <- calcola_cv()
                sens <- calcola_sensibilita()
                
                domcard <- c(rep(0, ndom))
                dom <- as.data.frame(dom)
                for (i in 1:ndom) domcard <- c(domcard, nlevels(dom[, 
                  i]))
                
                tit_cv <- c("TYPE", "DOMAIN/VAR.", "PLANNED CV ", 
                  "ACTUAL CV", "SENSITIVITY 10%")
                outcv <- cbind(as.vector(cvDom), paste(as.vector(cvDom2), 
                  "/V", c(1:ncol(med)), sep = ""), as.vector(cv), 
                  as.vector(CVfin), as.vector(ceiling(sens)))
                colnames(outcv) <- tit_cv
                return(outcv)
            }
            strato <- stratif$STRATO
            # default for cost and cens
            if (is.null(strato)) 
                strato <- paste("STR", 1:nstrat, sep = "")
            confr <- stampa_confronto(n, N, strato)
            outcv <- stampa_cv()
            attr(n, "confr") <- confr
            attr(n, "outcv") <- outcv
            
            # print (sum(n))
        }  #End Print
    
    return(n)
    
    # End body
}

checkData <- function(strata = NULL, errors = NULL) {
    # controls on strata dataframe
    if (!is.null(strata)) {
        if (sum(grepl("N", colnames(strata))) < 1) 
            stop("In strata dataframe the indication of population (N) is missing")
        
        if (sum(grepl("STRAT", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of stratum (STRATUM) is missing")
        if (sum(grepl("DOM1", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of at least Ine domain (DOM1) is missing")
        # if
        # (sum(grepl('CENS',toupper(colnames(strata)),fixed=TRUE))
        # < 1) stop('In strata dataframe the indication of strata
        # to be sampled or censused (CENS) is missing') if
        # (sum(grepl('COST',toupper(colnames(strata)),fixed=TRUE))
        # < 1) stop('In strata dataframe the indication of
        # interviewing cost in strata (COST) is missing')
        if (sum(grepl("M1", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of at least one mean (M1) is missing")
        if (sum(grepl("S1", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of at least one standard deviation (S1) is missing")
        # if
        # (sum(grepl('M+[0123456789]',toupper(colnames(strata)),perl=TRUE))
        # !=
        # sum(grepl('S+[0123456789]',toupper(colnames(strata)),perl=TRUE))+1)
        # stop('In strata dataframe the number of means (Mx)
        # differs from the number of standard deviations (Sx)')
    }
    # controls on errors dataframe
    if (!is.null(errors)) {
        if (sum(grepl("DOM", toupper(colnames(errors)), fixed = TRUE)) < 
            1) 
            stop("In errors dataframe the indication of domain (DOM) is missing")
        if (sum(grepl("CV", toupper(colnames(errors)), fixed = TRUE)) < 
            1) 
            stop("In errors dataframe the indication of at least one constraint (CV) is missing")
    }
    # crossed controls between errors and strata
    if (!is.null(errors) && !is.null(strata)) {
        # if
        # (sum(grepl('S+[0123456789]',toupper(colnames(strata)),perl=TRUE))
        # != sum(grepl('CV',toupper(colnames(errors)),fixed=TRUE)))
        # stop('In strata dataframe the number of means and std
        # deviations differs from the number of coefficient of
        # variations in errors dataframe')
        if (sum(grepl("DOM", toupper(colnames(strata)), fixed = TRUE)) != 
            nrow(errors)) 
            stop("The different domains (DOMx) in strata dataframe are not represented in errors dataframe")
    }
    
}
