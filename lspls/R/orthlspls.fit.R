###
### Orthogonal case, fitting of model
###

orthlspls.fit <- function(Y, X, Z, ncomp) {
    ## Parametres:
    nObs <- nrow(X)
    totNumComps <- sum(unlist(ncomp))
    totNumCoefs <- ncol(X) + totNumComps
    if(totNumCoefs > nObs) stop("Too many variables/components selected.")

    ## Containers:
    B <- matrix(nrow = totNumCoefs, ncol = ncol(Y)) # Regr. coefs.
    V <- matrix(nrow = nObs, ncol = totNumCoefs) # Regr. variables
    models <- list()                    # plsr models
    orthCoefs <- list()                 # Orthogonalisation matrices
    ## These two are not strictly neccessary:
    S <- list()                         # Scores
    L <- list()                         # Loadings


    ## Choose PLS fit algorithm.  FIXME: Support other algs?
    pls.fit <- oscorespls.fit

    ## Start with X:
    lsX <- lm.fit(X, Y)
    nVar <- ncol(X)
    ## Extract
    V[,1:nVar] <- X
    B[1:nVar,] <- lsX$coefficients
    res <- lsX$residuals

    ## For testing:
    ##Balt <- B

    ## Walk through Z:
    for (i in seq(along = Z)) {
        ##cat("i =", i, "\n")
        M <- Z[[i]]
        if (is.matrix(M)) {             # Single matrix
            Mo <- orth(M, V[,1:nVar])   # Orth. M against all used variables
            orthCoefs[[i]] <- Corth(M, V[,1:nVar]) # For pred
            models[[i]] <- pls.fit(Mo, res, ncomp[[i]])# Could use Y
            V[,nVar + (1:ncomp[[i]])] <- S[[i]] <- models[[i]]$scores
            L[[i]] <- models[[i]]$loadings
            ## FIXME: Testing:
            ##lmZ <- lm.fit(models[[i]]$scores, res)
            ##Balt[nVar + (1:ncomp[[i]]),] <- lmZ$coefficients
            ## FIXME: This depends on the pls algorithm (at least the scaling):
            B[nVar + (1:ncomp[[i]]),] <- t(models[[i]]$Yloadings)
            res <- models[[i]]$residuals[,,ncomp[[i]]]
            nVar <- nVar + ncomp[[i]]
        } else {                        # Parallell matrices
            Vadd <- matrix(nrow = nObs, ncol = sum(ncomp[[i]])) # The variables to be added in the present step
            added <- 0
            S[[i]] <- list()
            L[[i]] <- list()
            orthCoefs[[i]] <- list()
            models[[i]] <- list()
            for (j in seq(along = M)) {
                ##cat("j =", j, "\n")
                ## Walk through Z[[i]]
                Mo <- orth(M[[j]], V[,1:nVar])
                orthCoefs[[i]][[j]] <- Corth(M[[j]], V[,1:nVar]) # For pred
                models[[i]][[j]] <- pls.fit(Mo, res, ncomp[[i]][[j]])# Could use Y
                Vadd[,added + (1:ncomp[[i]][[j]])] <- S[[i]][[j]] <- models[[i]][[j]]$scores
                L[[i]][[j]] <- models[[i]][[j]]$loadings
                added <- added + ncomp[[i]][[j]]
            }
            V[,nVar + (1:sum(ncomp[[i]]))] <- Vadd
            ## Not strictly neccessary in orth. version:
            lmZ <- lm.fit(Vadd, res)
            B[nVar + (1:sum(ncomp[[i]])),] <- lmZ$coefficients
            ##Balt[nVar + (1:sum(ncomp[[i]])),] <- lmZ$coefficients
            res <- lmZ$residuals
            nVar <- nVar + sum(ncomp[[i]])
        } # if
    } # for
    list(coefficients = B, predictors = V, orthCoefs = orthCoefs,
         models = models, ncomp = ncomp, scores = S, loadings = L, residuals = res)
} # function

## Forenklingstanke: Gjoer om alle enkeltmatrisene i Z til lister med
## ett element.  Da blir algoritmene enklere (og burde ikke bli
## nevneverdig saktere).  Ved behov kan man teste paa length(M)
## (f.eks. ved beregning av B og ny res).

## FIXME:
## - Change name of 'models' component to 'plsmodels'?
## - Add fitted values?
