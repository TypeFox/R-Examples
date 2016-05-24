# aggiornata il 23-12-2009: modificata con il test del prof Pesarin
# funzione per la combinazione iterata, utilizza
# le tre funzioni di combinazioni di:
# > Tippet
# > Fisher
# > Liptak (Quantile Normale)
#------------------------------------------------------------
# Inputs:
# > P:   matrice, ogni colonna ha la distribuzione
#        di permutazione dei p-values di uno dei
#        test parziali da combinare
# > tol: tolleranza per il test di arresto, valore
#        maggiore di 1, rappresenta il rapporto tra
#        il massimo errore tollerato e la minima
#        precisione raggiungibile "1/B"
# > maxIter: massime iterazioni consentite
# > test: tipologia del test d'arresto,
#         > "SSQ" = numerico, Sum of SQuares, i.e. (n-1)*(s^2)
#         > "ABS" = logico, errore relativo assoluto delle
#                   differenze a coppie
#         > "NORM2" = numerico, norma euclidea del vettore
#                     delle differenze tra due iterazioni
#                     consecutive
#         > "EDF" = basato sulla ripartizione empirica,
#                   complemento a uno dei p-value, una
#                   sorta di t-test
# > plotIt: disegnare o no il grafico delle iterazioni?
#------------------------------------------------------------
# Outputs:
# > P.final: matrice contenente la distribuzione di permutaz.
#            dei p-values combinati secondo le tre combinazioni
# > P.iter: matrice contenente i p-values osservati di ogni
#           iterazione
#------------------------------------------------------------
# 
#' Iterated NonParametric Combination of the test statistic matrix, mostly for 
#' internal use.
#' 
#' It takes as input a matrix whose columns have to be combined with the 
#' iterated NPC procedure, default combining functions are 3: Fisher, 
#' Liptak (normal version), and Tippett (minP)
#' 
#' @title Iterated NonParametric Combination
#' @param P 
#'     input \code{matrix} containing the test statistic in the form of p.values 
#'     (permutation or asymptotic)
#' @param tol 
#'     \code{integer} representing the desired tolerance, the actual one being 
#'     \eqn{\frac{tol}{B}}{tol/B} where \eqn{B} is the number of permutations
#' @param maxIter 
#'     \code{integer} maximum number of iterations to be performed, default 10
#' @param plotIt 
#'     \code{logical}, if \code{TRUE} (default) plots the diagnostic grahp of 
#'     $p$-values obtained with each combining function \emph{vs.} iteration index
#' @param combFun1 
#'     first combining \code{function} needed for the algorithm, default is 
#'     \emph{Fisher}'s: \eqn{-2 \sum_i \log{p_i}}{-2 * sum(log(p_i))}
#' @param combFun2 
#'     second combining \code{function}, default is \emph{Liptak}:
#'     \eqn{\sum_i \Phi^{-1} \left( 1 - p_i \right)}{sum(\Phi^{-1} (1-p_i))}
#' @param combFun3 
#'     third combining \code{function}, default is \emph{Tippett}:
#'     \eqn{-\min_i p_i}{-min_i p_i}
#' @param test
#'     \code{character}, it is the stopping rule used to check for convergence, 
#'     each one of the 4 kinds currently implemented takes as input the vector with
#'     the result of the combination with the different combining functions for 
#'     one permutation. There are 4 choices for this argument:\cr
#'     \describe{
#'         \item{\code{"SSQ"}}{
#'             Sum of SQuares, the algorithm stops when 
#'             \eqn{\sqrt{(n-1)(s^2)}}{square root of (n-1)*(s^2)} is smaller 
#'             than the actual tolerance; here where \eqn{s} is the sample 
#'             variance of the vector.}
#'         \item{\code{"ABS"}}{
#'             The algorithm stops if not all pairwise absolute differences 
#'             between the elements are smaller than the actual tolerance}
#'         \item{\code{"NORM2"}}{
#'             The algorithm stops if the euclidean distance between two 
#'             consecutive iterations is smaller than the actual tolerance.}
#'         \item{\code{"EDF"}}{
#'             It is based on the Empirical Distribution Function of the p.values. 
#'             The algorithm stops if the standardized absolute difference  
#'             between the average of two consecutive iterations is smaller 
#'             than the actual tolerance. The standardization involves the variance
#'             of the numerator, it is a sort of t-test.}
#'     }
#' @param Pmat 
#'     \code{logical}, if \code{TRUE} returns the final matrix of combined $p$-values,
#'     default is \code{FALSE} 
#' @param onlyCombined 
#'     \code{logical}, if \code{TRUE} returns only the first column of the final matrix,
#'     in case only the distribution of combined $p$-values is needed
#' @return 
#'     The output is conditioned on some of the input argument.
#'     The default is a \code{list} containing only the element \code{"P.iter"}: 
#'     a \code{matrix} with 3 columns containing the observed p.values across 
#'     iterations and for all combining functions (to manually check convergence).
#'     If \code{Pmat} is \code{TRUE} than the list contains also the element 
#'     \code{"P.final"} that is the final permutation space of p.values obtained 
#'     with all combining function.
#'     If \code{onlyCombined} is \code{TRUE} than the resulting output is just the
#'     vector containing the first column of \code{"P.final"}. 
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @export
#' 
iterNPC <- function(P, tol = 1, maxIter = 10, plotIt = TRUE,
        combFun1 = function(x)        # Fisher
        {
            -2 * sum(log(x), na.rm = TRUE)
        }, 
        combFun2 = function(x)        # Liptak Normal
        {
#            x[x == 1] <- 1 - 2e-12
#            x[x == 0] <- 2e-12
            sum(qnorm(1 - x), na.rm = TRUE)
        },  
        combFun3 = function(x)        # Tippett
        {
            -min(x, na.rm = TRUE)     # Tippett is with <=
        }, 
        test = c("SSQ", "ABS", "NORM2", "EDF"), Pmat = FALSE, 
        onlyCombined = FALSE)
{
    
#stat2p <- function(x, B, cf = TRUE){
#    r <- 1 - rank(x[-1],ties.method="min")/B
#    if(cf)  r <- r + 1/B
#    z <- x[is.na(x) == FALSE]
#    return( c(mean(z[-1] >= z[1]), r) )
#}#end_stat2p
    
    ### names of the combining functions, if default then we know the names ;-)
    if (missing(combFun1) && missing(combFun2) && missing(combFun3))
    {
        nams <- c("Fisher", "Liptak", "Tippett")
    } else
    {
        nams <- paste0("combFun", 1:3L)
    }
    
    stat2p <- function(x, B, cf = FALSE)
    {
        r <- rank(-x, ties.method = "max", na.last = "keep")/B
#        z <- x[!is.na(x)]
#        return(c(mean(z[-1] >= z[1]), r))
        return(r)
    }# END:stat2p
    
    if(missing(test))
        test <- "SSQ"
    
    if(!test %in% c("SSQ", "ABS", "NORM2", "EDF"))
        stop("test dev'essere uno tra: 'SSQ','ABS','NORM2','EDF'")
    
    B <- NROW(P)          # number of permutations
    k <- NCOL(P)          # number of hypothesis
    
    eps <- tol * (1/(B - 1)) # actual tolerance
    
    #Fish  <- function(x)  -2*log(prod(x, na.rm = TRUE)) # Fisher
    #Tipp  <- function(x)  min(x, na.rm = TRUE) # Tippett
    #Lipt  <- function(x){
    #    x[x == 1] <- 1 - 2e-16 ; x[x == 0] <- 2e-16
    #    sum(qnorm(1 - x), na.rm = TRUE)
    #}# Liptak
        
    # -> matrici per le iterazioni:
    #    p-values combinati osservati, piu'
    #    loro distribuzione di permutazione
    g1 <- array(NA, dim = c(B, 3))
    if((test == "NORM2") || (test == "EDF"))
        g2 <- g1
    
    # -> matrice dei p-values ottenuti
    #    ad ogni iterazione:
    #    riga = num. delle iterazioni
    iterPVals <-  rep(NA, 3)
    names(iterPVals) <- nams 
    cont <- 0 # indicatore iterazione corrente
    stopTest <-  c() # vettore per test d'arresto
    
    
    #---> primo step: 
    #     combino e passo ai p-values
    g1[, 1] <- stat2p(apply(P, 1, FUN = combFun1), B = B)
    g1[, 2] <- stat2p(apply(P, 1, FUN = combFun2), B = B)
    g1[, 3] <- stat2p(apply(P, 1, FUN = combFun3), B = B) 
    #g1[B + 1,] <- g1[1,] # ultima permutazione = osservata
    
    #---> p-values oss. al primo step
    iterPVals <- g1[1, ]
    cont <- 1 # 1^ combinazione eseguita
    
    g2 <- g1 # per la prima iterazione
    
#---> primo test d'arresto
    if(test == "ABS")
    {
        aux <- tcrossprod(iterPVals)
        stopTest[1] <- max(aux[lower.tri(aux)]) > eps
#        stopTest[1] <- (
#                    (abs(iterPVals[1] - iterPVals[2]) > eps) ||
#                    (abs(iterPVals[1] - iterPVals[3]) > eps) ||
#                    (abs(iterPVals[2] - iterPVals[3]) > eps))
    }
    
    if(test == "SSQ")
        stopTest[1] <- sqrt(2 * var(iterPVals)) > eps
    
    if((test == "NORM2") || (test == "EDF"))
    {
        stopTest[1] <- TRUE
        g2[, 1] <- stat2p(apply(g1, 1, FUN = combFun1))
        g2[, 2] <- stat2p(apply(g1, 1, FUN = combFun2))
        g2[, 3] <- stat2p(apply(g1, 1, FUN = combFun3))
        # g2[B + 1,] <- g2[1,] # ultima permutazione = osservata
        
        iterPVals <- rbind(iterPVals, g2[1, ])
        cont <- cont+1
        stopTest[2] <- sqrt(crossprod(iterPVals[2, ] - iterPVals[1, ])) > eps
    }#end NORM2
    
    if(test == "EDF")
    {
        p.0 <- mean(iterPVals[1, ])
        p.1 <- mean(iterPVals[2, ])
        varP <- p.0 * (1 - p.0) * eps
        stopTest[2] <- abs(p.1 - p.0) > 2 * sqrt(varP)
    }#end EDF
    
    
    while((stopTest[cont]) && (cont <= maxIter))
    {
        #---> aggiorno la matrice di p-val
        g1 <- g2
        
        #---> combino e passo ai p-values
        g2[, 1] <- stat2p(apply(g1, 1, FUN = combFun1), B = B)
        g2[, 2] <- stat2p(apply(g1, 1, FUN = combFun2), B = B)
        g2[, 3] <- stat2p(apply(g1, 1, FUN = combFun3), B = B) 
        # g2[B + 1,] <- g2[1,] # ultima permutazione = osservata
        
        #---> p-values osservati
        iterPVals <- rbind(iterPVals, g2[1, ])
        cont <- cont + 1
        
        #---> test d'arresto
        if(test == "ABS")
        {
            stopTest[cont] <- (
                        (abs(iterPVals[cont, 1] - iterPVals[cont, 2]) > eps) ||
                        (abs(iterPVals[cont, 2] - iterPVals[cont, 3]) > eps) ||
                        (abs(iterPVals[cont, 2] - iterPVals[cont, 3]) > eps))
        }#end ABS
        
        if(test == "SSQ")
            stopTest[cont] <- sqrt(2 * var(iterPVals[cont, ])) > eps
        
        if(test == "NORM2")
            stopTest[cont] <- sqrt(crossprod(iterPVals[cont,] - iterPVals[cont-1,])) > eps
        
        if(test=="EDF")
        {
            p.0 <- p.1
            p.1 <- mean(iterPVals[cont, ])
            varP <- p.0 * (1 - p.0) * eps
            stopTest[cont] <- abs(p.1 - p.0) > 2 * sqrt(varP)
        }#end EDF
    }#end while
    
    if(!is.matrix(iterPVals))
        iterPVals <- t(iterPVals)
    
    
    if((plotIt == TRUE))
    {
        if (NROW(iterPVals) > 1 )
        {
            # yLim  <- range(pv.it)
            # yLim <- c(min(iterPVals, .001), max(iterPVals, 0.25))
            yLim <- c(0, max(iterPVals, 0.25))
            if(length(dev.list()) <= 1)
                dev.new()
            
            plot(iterPVals[, 1], type = 'l', ylim = yLim, lwd = 2, col = 1, lty = 2,
                xaxt = "n")
            
            title(main = "Iterations Behaviour")
            lines(iterPVals[, 2], lwd = 2, col = 2, lty = 1)
            lines(iterPVals[, 3], lwd = 2, col = 4, lty = 3)
        
            axis(side = 1, at = seq_len(nrow(iterPVals)), 
                    labels = seq_len(nrow(iterPVals)))
            
            
            abline(h = .01, lwd = 2, col = "gray30", lty = 4)
            text(x = (nrow(iterPVals) - 1), y = 0.01, pos = 3, col = "gray30",
                    offset = .33,
                    labels = expression(alpha == 0.01))
            abline(h = .05, lwd = 2, col = "gray50", lty = 4)
            text(x = (nrow(iterPVals) - 1), y = 0.05, pos = 3, col = "gray50",
                    offset = .33,
                    labels = expression(alpha == 0.05))
            abline(h = .1, lwd = 2, col = "gray70", lty = 4)
            text(x = (nrow(iterPVals) - 1), y = 0.1, pos = 3, col = "gray70",
                    offset = .33,
                    labels = expression(alpha == 0.1))
            
            legendPos <- ifelse(
                    abs(mean(iterPVals[nrow(iterPVals), ]) - yLim[2]) < .1 * diff(yLim), 
                    "right", "topright")
            legend(x = legendPos, legend = c("Fisher", "Liptak", "Tippett"),
                    col = c(1, 2, 4), lty = c(2, 1, 3), lwd = rep(2, 3))
        } else
        {
            message("Convergence at the first step, no plot produced")
#            title(main = "Only one iterations")
#            abline(h = iterPVals[, 1], lwd = 2, col = 1, lty = 2)
        }# END: ifelse - more than one iteration
    }# END:plotIt
    
    
    rownames(iterPVals) <- paste0("iter", seq_len(NROW(iterPVals)))
    colnames(iterPVals) <- colnames(g2) <- nams
    Res <- list("P.iter" = iterPVals)
    
    if(Pmat)
        Res[["P.final"]] <- g2
    
    if (onlyCombined)
        Res <- g2[, 1]
    
    return(invisible(Res))
}
