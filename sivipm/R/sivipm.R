##################################################################
# sivipm R package
# Copyright INRA 2016
# INRA, UR1404, Research Unit MaIAGE
# F78352 Jouy-en-Josas, France.
#
# URL: http://cran.r-project.org/web/packages/sivipm
#
# This file is part of sivipm R package.
# sivipm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################



##########################################################################
# Compute the sensitivity indices isivip and tsivip after regression PLS
##########################################################################
#  Total or individual sensitivity indices SIVIP
# @title compute several TSIVIP results 
# @sivipm
# @param Y outputs  data.frame 
# @param XIndic, an object of class 'polyX'
# @param nc number of components
# @param options outputs selection. String vector.
# Valid values: 'tsivip', 'isivip', 'signif', 'all'
# @param graph if TRUE, display graph (for tsivip when alea =false)
# @param alea alea is added (only for tsivip)
# @param fast (if no missing only): fast procedure
# @param output additional results are returned
# @return
# \itemize{
# \item{tsivip}{ - when \code{options} include 'tsivip',  total sensitivity indices for each input variable and their percentage.}
# \item{isivip}{ - when \code{options} include 'isivip',
# first order individual sensitivity indices for each monomial.
# When alea=TRUE, the non significant monomials are removed}
# \item{SIMCA}{ - when \code{options} include 'simca',
# the significative components (when no missing) calculated
# by SIMCA-P rule,
# \item{Lazraq}{ - when \code{options} include 'lazraq',
# the significative components
# (missing accepted if one response, only) calculated
# by the test de Lazraq).}
# \item{output}{ - required additional results.}
# }
# @export

sivipm <- function(Y, XIndic, nc = 2, options = c("fo.isivip", "tsivip", "simca", 
    "lazraq"), graph = FALSE, alea = FALSE, fast = FALSE,
                   output = NULL) {
    # Initialisation of the returned structures
    retour <- sivip() # creation of an object of class 'sivip'
    retisivip <- NULL
    rettsivip <- NULL
    retsimca <- NULL
    retlaz <- NULL
    regpls <- NULL

    
    # Access to the data
    data.exp <- XIndic@dataX.exp
    nmono <-ncol(data.exp)
    if (nc > nmono )
      stop(paste("The number of components, ", nc, ", should be less or equal to the number of monomials, ", nmono, "\n", sep=""))
    
    
    if (length(XIndic@Pindic@P) >0) {
        Indic <- XIndic@Pindic@indic
        nvar <- ncol(Indic)
    } else {
        Indic <- NULL
        nvar <- ncol(data.exp)
    }
    Y <- as.matrix(Y) # in case of a vector=> a matrix with one column
    
    # Access to the options
    # Check options
    if (is.null(options)) stop("At least one value in argument 'options' must be given")
    oo <- options
    oo[grepl("isivip", options, ignore.case = TRUE)]=NA
    oo[grepl("tsivip", options, ignore.case = TRUE)]=NA
    oo[grepl("simca", options, ignore.case = TRUE)]=NA
    oo[grepl("lazraq", options, ignore.case = TRUE)]=NA
    
    if (!all(is.na(oo))) {
      stop(paste("Bad  values in argument 'options'. Valid values are: 'fo.isivip', 'tsivip', 'simca', 'lazraq' ", sep=""))
    }
    
    # Set option indicators      
    opt.isivip <- any(grepl("isivip", options, ignore.case = TRUE)) || any(grepl("isivip", 
        output, ignore.case = TRUE))
    opt.tsivip <- any(grepl("tsivip", options, ignore.case = TRUE))
    opt.simca <- any(grepl("simca", options, ignore.case = TRUE))
    opt.laz <- any(grepl("lazraq", options, ignore.case = TRUE))
    
    
    if (!opt.tsivip && alea) {
        warning("sivipm : option alea ignored when tsivip not required")
    }
    if (!opt.tsivip && graph) {
        warning("sivipm : option graph ignored when tsivip not required")
    }
    
    # with/without missing values
    if (any(is.na(data.exp)) || any(is.na(Y))) 
        na.miss <- TRUE else na.miss <- FALSE
    
    # with/without additional results
    if (!is.null(output))  {
      if (!all(output %in% c("isivip", "VIP","RSS","PRESS","Q2","betaNat","PLS")))
        stop("Bad values in argument 'output'. Valid values are: 'isivip', 'VIP','RSS','PRESS','Q2','betaNat','PLS' ")
      routput <- TRUE
    }  else routput <- FALSE
    
    if (fast  && na.miss) {
        warning("sivipm : option fast ignored when there are missing values.")
        fast <- FALSE
      }

    if (opt.simca && na.miss) {
        warning("sivipm : option simca ignored when there are missing values.")
        opt.simca <- FALSE
    }
    if (opt.laz && na.miss && (ncol(as.matrix(Y)) > 1)) {
        warning("sivipm : option Lazraq ignored when there are missing values and several output variables.")
        opt.laz <- FALSE
    }
    
    if (opt.simca || opt.laz) {
        # Call regpls2 with all outputs because Q2 required
      
        regpls <- regpls2(Y, data.exp, nc,
                          na.miss= na.miss, fast=fast, output = TRUE)
        
        if (opt.simca) 
            retsimca <- fsimcarule(regpls$Q2)
        if (opt.laz) 
            retlaz <- rlaz(regpls$x.scores, Y, nc)
    }


    
    # significant components not required and
    # (isivip or (tsivip without alea))
    if (is.null(regpls) && (opt.isivip || (opt.tsivip && alea == FALSE))) {
        regpls <- regpls2(Y, data.exp, nc,
                          na.miss= na.miss, fast=fast, output = routput)
    }
    
    
    
    if (opt.isivip) {
      retisivip <- isivip(regpls$VIP[, nc], data.exp)
        # Consider only the VIP of the last component
    }

    if (opt.tsivip) 
        {
            if (is.null(Indic) || all(is.na(Indic))) {
                warning("sivipm : option tsivip ignored because there is no polynomial description")
            } else {
                if (alea) {
                  # case where there is alea
                  rett <- tsivipalea(Y, XIndic, nc,
                                     na.miss, fast=fast, output = routput)
                  if (!is.null(retisivip)) 
                    {
                      # add the isivip of the alea into retisivip
                      nvar <- ncol(Indic)
                      mono <- names(retisivip)
                      nmono <- length(mono)
                      bb <- c(retisivip[1:nvar], rett$isivipalea, retisivip[(nvar + 
                        1):nmono])
                      names(bb) <- c(mono[1:nvar], "X_alea", mono[(nvar + 1):nmono])
                      retisivip <- bb
                    }  # fin !is.null(retisivip)
    # Remove from the outputs what is about the alea
    rett$isivipalea <- NULL
    nrett <- length(rett$tsivip)
    rett$tsivip <- rett$tsivip[-nrett]
    rett$tsivip.percent <- rett$tsivip.percent[-nrett]
    nmono <- length(rett$monosignif)
    rett$monosignif <- rett$monosignif[-nmono]

                  
                  regpls <- rett$regpls
                  rett$regpls <- NULL
                  rettsivip <- rett
                  if (graph == TRUE) 
                    tsivipgraph(rettsivip$tsivip, retisivip, Indic, nc, alea)
                } else {
                  # case where there is no alea
                  rettsivip <- tsivipnoalea(regpls$VIP, Y, Indic, nc)
                  if (graph == TRUE) {
                    tsivipgraph(rettsivip$tsivip, retisivip, Indic, nc, alea, rettsivip$Evol)
                  }
                  rettsivip$Evol <- NULL
                }
            }  # end else
        }  # end options=='tsivip'
    
    # Prepare the output    
    if (any(grepl("isivip", options))) {
        # First order isivip
        retour@fo.isivip <- retisivip[1:nvar]
        if (!is.null(colnames(Indic))) {
          names(retour@fo.isivip) <-  colnames(Indic)
        }
        # Percentage fo.isivip
        som <- sum(retour@fo.isivip)
        retour@fo.isivip.percent <- sort( (retour@fo.isivip * 100)/som,
                                         decreasing=TRUE)
    }

     
    for (nom in names(rettsivip)) {
      slot(retour, nom) <- rettsivip[[nom]]
    }
    for (nom in names(retsimca)) {
      slot(retour, nom) <- retsimca[[nom]]
    }
    for (nom in names(retlaz)) {
      slot(retour, nom) <- retlaz[[nom]]
    }
   
    if (!is.null(regpls)) 
        {
          
            if (any(grepl("isivip", output, ignore.case = TRUE))) {
                retour@output$isivip <- retisivip
            }
            
            if (any(grepl("VIP", output))) {
                retour@output$VIP <- regpls$VIP
                retour@output$VIPind <- regpls$VIPind
                regpls$VIP <- NULL
                regpls$VIPind <- NULL
            }
          # Results when there are no missing, only
            # Ignored: PRESS and RSS are no more an option 
           if (any(grepl("RSS", output, ignore.case = TRUE))) {
             if (na.miss) {
               warning("RSS is not calculated when there are missing values")
             } else {
                 retour@output$RSS <- regpls$RSS
                 regpls$RSS <- NULL
               }
           }
             if (any(grepl("PRESS", output, ignore.case = TRUE))) {
               if (na.miss) {
               warning("PRESS is not calculated when there are missing values")
             } else {
                 retour@output$PRESS <- regpls$PRESS
                 regpls$PRESS <- NULL
             }
             }

            

            if (any(grepl("Q2", output, ignore.case = TRUE))) {
              if (na.miss) {
              warning("Q2 is not calculated when there are missing values")
            } else {
                retour@output$Q2 <- regpls$Q2
                retour@output$Q2cum <- regpls$Q2cum
                
                regpls$Q2 <- NULL
                regpls$Q2cum <- NULL
            }
            }
          
            if (any(grepl("betaNat", output, ignore.case = TRUE))) {
                retour@output$betaNat <- regpls$betaNat
                retour@output$betaNat0 <- regpls$betaNat0
                regpls$betaNat <- NULL
                regpls$betaNat0 <- NULL
             }
 

            if (any(grepl("PLS", output, ignore.case = TRUE))) {
# ne  mettre dans PLS que les sorties spécifiques au cas ou on
# ne serait pas passé par les tests ci-dessus 
# (PLS demandé seulement, par exemple)
                regpls$VIP <- NULL
                regpls$VIPind <- NULL
                regpls$RSS <- NULL
		regpls$PRESS <- NULL
		regpls$Q2 <- NULL
                regpls$Q2cum <- NULL
		 regpls$betaNat <- NULL
                regpls$betaNat0 <- NULL
                retour@output$PLS <- regpls
            }
        }  # fin !regpls
    
    if (any(grepl("isivip", output, ignore.case = TRUE))) 
        retour@output$isivip <- retisivip
    
    return(retour)
}  # end sivipm

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# @fsimcarule
# @title compute rule implemented in SIMCA-P to determine the number of significative components from Q2 (no missing values)
# Internal function
#@return  significative components

fsimcarule <- function(Q2) {
    
    #    a <- ncol(Q2)
    # signif <- as.integer(Q2[, a] >= 0.0975) s <- max(2, max(which(signif == 1)) +
    # 1) ret <- list(signifcomponents = s)
    signif <- (Q2 >= 0.0975)

    return(list(simca.signifcomponents=signif))
}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Compute the total (or 'global') indices tsivip
# @title compute TSIVIP indices
# @tsivipf
# @param vecin VIP vector
# @param Indic indic matrix
# @return tsivip sensitivity indices for each input variable
# Internal function

tsivipf <- function(vecin, Indic) {
    a <- isivip(vecin, t(Indic))
    tsivip <- t(Indic) %*% a
    
    tri <- sort(tsivip, index.return = TRUE, decreasing = TRUE)
    if (length(tri$x) == 0) {
        rank <- 1:length(tsivip)
    } else {
        rank <- tri$ix
        tsivip <- tri$x
    }
    
    # b<-barplot(tsivip,xlab=as.character(rank))
    ret <- list(tsivip, rank)
    return(ret)
}  # end tsivipf

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# @tsivipgraph
# internal
# TSIVIP graphs

tsivipgraph <- function(tsivip, isivip, Indic, nc, alea, Evol = NULL) {
    taillept <- 0.8  # taille des annotations
    d <- seq(1, nc)
    grcomp <- (!is.null(Evol) && (nc > 1))
    if (grcomp) {
        # we will do the graph of the components
        nombrecomposantes <- d
        for (i in 1:(ncol(Indic) - 1)) {
            nombrecomposantes <- cbind(nombrecomposantes, d)
        }
        par(mfrow = c(2, 1))
    } else {
        # we won't do the graph of the components
        par(mfrow = c(1, 1))
    }
    sortsivip <- sort(tsivip, decreasing = TRUE, index.return = TRUE)
    xsivip <- sortsivip$ix
    label <- names(tsivip)
    # space=0, pour que les barres soient jointes las=1, pour que les noms des
    # monomes soient horizontaux
    barplot(sortsivip$x, xlab = "TSIVIP (%)", names.arg = label[xsivip], horiz = TRUE, 
        space = 0, cex.axis = taillept, cex.names = taillept, cex.lab = taillept, 
        las = 1)
    if (!is.null(isivip)) 
        {
            # superpose the isivip of the unit monomials
            visivip <- isivip[xsivip]
            barplot(visivip, add = TRUE, names.arg = rep("", length(visivip)), col = 1, 
                angle = 45, density = 20, legend.text = "First order ISIVIP", space = 0, 
                horiz = TRUE, args.legend = list(bty = "n", cex = taillept), 
                yaxt = "n", xaxt = "n")
        }  # fin (!is.null(isivip))
    if (grcomp) 
        {
            # Evol est nul dans le cas alea pour mettre les ticks-marks sur l'axe des x en
            # entier, on les met 'a la main'
            
            matplot(nombrecomposantes, Evol, type = "l", ylab = "TSIVIP", xlab = "number of components", 
                xaxt = "n", cex.axis = taillept, cex.lab = taillept)
            axis(1, at = 1:nc, labels = 1:nc, cex.axis = taillept)
            # pie(Evol[nrow(Evol),])
        }  # fin (!is.null(Evol)
    return(invisible())
}  # fin tsivipgraph


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# @tsivipnoalea
# @title compute TSIVIP total sensitivity indices from VIP
# when no alea
# Internal function
# @param VIP component VIP returned by regpls2
# @param Y outputs  data.frame 
# @param nc number of components
# @param Indic matrix computed by indic function
# @param graph  display graph
# @return tsivip total indices for each input variable

tsivipnoalea <- function(VIP, Y, Indic, nc = 2, graph = FALSE) {
    
    if (!is.null(colnames(Indic))) 
        varnames <- colnames(Indic) else varnames <- paste("X", 1:ncol(Indic), sep = "")
    
    
    TSIVIP <- matrix(nrow = nc, ncol = ncol(Indic))
    TSIVIPr <- matrix(nrow = nc, ncol = ncol(Indic))
    for (j in 1:nc) {
        VIPj <- VIP[, j]
        
        aa <- tsivipf(VIPj, Indic)
        TSIVIP[j, ] <- aa[[1]]
        TSIVIPr[j, ] <- aa[[2]]

    }
    r <- matrix(nrow = ncol(Indic), ncol = nc)
    for (i in 1:ncol(Indic)) {
        for (j in 1:nc) {
            r[i, j] <- which(TSIVIPr[j, ] == i)
        }
    }
    Evol <- matrix(ncol = ncol(Indic), nrow = nc)
    for (i in 1:ncol(Indic)) {
        for (k in 1:nc) {
            Evol[k, i] <- TSIVIP[k, r[i, k]]
        }
    }
    dimnames(Evol) <- list(1:nc, varnames)
    # Prendre la derniere ligne de Evol, i.e ce qui correspond a la derniere
    # composante

    tsivip <- Evol[nrow(Evol), ]
    names(tsivip) <- varnames
    
    percent <- sort((tsivip/sum(tsivip)) * 100, decreasing = TRUE)
    ret <- list(tsivip = tsivip, tsivip.percent = percent, Evol = Evol)
    return(ret)
    
}  # end tsivipmnoalea
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# @tsivipalea
# @title compute TSIVIP total sensitivity indices from VIP
# with alea
# Internal function
# @param Y outputs  data.frame 
# @param XIndic, an object of class 'polyX'
# @param nc number of components
# @param Indic matrix computed by indic function
# @return tsivip total indices for each input variable
#  The non significant monomials are removed
# from the tsivip calculation.


tsivipalea <- function(Y, XIndic, nc,
                       na.miss, fast=FALSE, output = FALSE) {
    dataX.exp <- XIndic@dataX.exp
    Indic <- XIndic@Pindic@indic
    
    if (!is.null(colnames(Indic))) 
        varnames <- colnames(Indic) else varnames <- paste("X", 1:ncol(Indic), sep = "")
    
    X_alea <- runif(nrow(dataX.exp))
    # Add a random variable into dataX.exp
    Xalea <- cbind(dataX.exp, X_alea)
    # Add a column and a row into Indic
    Indicalea <- cbind(Indic, rep(0, nrow(Indic)))
    varnames <- c(varnames, "X_alea")
    colnames(Indicalea) <- varnames
    Indicalea <- rbind(Indicalea, rep(0, ncol(Indicalea)))
    Indicalea[nrow(Indicalea), ncol(Indicalea)] <- 1
    
    regpls <- regpls2(Y, Xalea, nc, na.miss, fast=fast, output = output)
    VIP <- regpls$VIP
    VIP <- VIP[, ncol(VIP)]
    isivipalea <- isivip(VIP, Xalea)
    leisivipalea <- isivipalea[length(isivipalea)]
    signif <- (isivipalea > leisivipalea)
    isivipalea <- isivipalea[signif]
    Indicalea <- Indicalea[signif, ]
    tsivipa <- t(Indicalea) %*% isivipalea
    tsivipa <- as.vector(tsivipa)
    names(tsivipa) <- varnames
    percentt <- sort((tsivipa/sum(tsivipa)) * 100, decreasing = TRUE)
    correlation <- cor(X_alea, Y)
    retour <- list(isivipalea = leisivipalea,
                   tsivip = tsivipa,
                   tsivip.percent = percentt, 
        monosignif = signif, correlalea = correlation)
    if (output) {
        retour$regpls <- regpls
    }
    
    return(retour)
    
    
}  # end tsivipalea
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# rlaz
# compute the Lazraq et Cléroux inferential test It is used to determine the
# significant components when there are missing values and a single response
# variable.  Input TT: x.scores. Matrix nc X nobs Y : response variable. Matrix
# nobs X 1 nc: number of components Internal function
rlaz <- function(TT, Y, nc) {
    YY <- scale(Y)  # centered-reduced response variable
    alpha <- 0.05  # alpha is fixed
    nobs <- nrow(YY)
    stu <- qt(1 - alpha, nobs)

   TTmat <- as.matrix(TT)
     mat <- cbind(YY, matrix(TTmat[1, ],  ncol = 1))
    res <- rep(FALSE, nc)
    for (irow in 1:nc) {
        mat[, ncol(mat)] <-TTmat[irow,]
        cr <- cor(mat, use = "na.or.complete")[1, 2]
         tlaz <- (sqrt(nobs) * cr)/sqrt(1 - cr * cr)
        
        if (tlaz > stu) {
            res[irow] <- TRUE
        }
    } # fin irow
    names(res) <- paste("c", 1:length(res), sep="")
    return(list(lazraq.signifcomponents = res))
}  # end rlaz 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calcul des indices individuels isivip
# @title compute isivip individual sensitivity indices for each monomial
# @isivip
# @param vecin VIP vector
# @param dataX.exp inputs data.frame
# @return isivip individual sensitivity indices for each monomial

isivip <- function(vecin, dataX.exp) {
  
    sivip <- rep(0, ncol(dataX.exp))
    if (!is.null(colnames(dataX.exp)))
      names(sivip) <- colnames(dataX.exp)
    else
      names(sivip) <- paste("X", 1:ncol(dataX.exp))
    for (i in 1:ncol(dataX.exp)) {
        sivip[i] <- (vecin[i] * vecin[i])/ncol(dataX.exp)
    }
    return(sivip)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# @title compute TSIVIP total sensitivity indices 
# @gosivip
# @param Y outputs  data.frame
# @param XIndic, an object of class 'polyX'
# @param nc number of components
# @param graph=FALSE display graph
# @return tsivip total indices for each input variable
# When alea=TRUE, the non significant monomials are removed
# from the tsivip calculation.
# It is a fast version of function of sivipm, because reduced to the
# case where only tsivip are required. Used in bootstrap.

gosivip <- function(Y, XIndic, na.miss, nc=2, graph = FALSE, alea = FALSE,
                    fast= FALSE) {

    dataX.exp <- XIndic@dataX.exp
    Indic <- XIndic@Pindic@indic
    if (!alea) 
        {
          regpls <- regpls2(Y, dataX.exp, nc, na.miss=na.miss,
                            fast=fast, output=FALSE)
          ret <- tsivipnoalea(regpls$VIP, Y,   Indic, nc, graph= graph)

        }  # end !alea
    
    if (alea) {
        ret <- tsivipalea(Y, XIndic, nc,  output=FALSE)
      } # end alea
    
    return(ret)  
    
} # end gosivip
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# @title realise des intervalles de confiance pour les tsivip par methode de bootstrap
# @sivipboot
# @param Y ouputs data.frame 
# @param XIndic, an object of class \code{\linkS4class{polyX}}
# @param nc number of components
# @param B number of bootstrap replicates
# @param alpha level of bootstrap confidence intervals
# @return tsivip percentile bootstrap confidence intervals

sivipboot <- function(Y, XIndic, B, nc=2, graph=FALSE, alea=FALSE,
                      fast=FALSE,
                      alpha = 0.05) {
    dataX.exp <- XIndic@dataX.exp
    Indic <- XIndic@Pindic@indic
    if (!alea)
      boot <- matrix(nrow = B, ncol = ncol(Indic))
    else
      boot <- matrix(nrow = B, ncol = ncol(Indic)+1)
    # in the case where there is alea, we add a variable
    
    # MAXTRY: maximal number of X generation at each loop
    MAXTRY = 5
    for (b in 1:B) {
        T <- sample(nrow(dataX.exp), nrow(dataX.exp), replace = TRUE)

        # AB: check scale is possible
        Xscale <- scale(as.matrix(dataX.exp[T, ]))
        ntry <- 1
        while (any(is.nan(Xscale)) && (ntry < MAXTRY) ) {
          T <- sample(nrow(dataX.exp), nrow(dataX.exp), replace = TRUE)
          Xscale <- scale(as.matrix(dataX.exp[T, ]))
          ntry <- ntry + 1
        } # end while

        
        if (any(is.nan(Xscale)))
            stop("sivipboot: unsuccessful generation of X")
        # End check scale

        XIndic@dataX.exp <- dataX.exp[T, ]
        ressivipm <- gosivip(Y[T, ], XIndic, na.miss=NA,
                             nc=nc, alea=alea,
                             fast=fast,
                            graph = FALSE)$tsivip
          
        boot[b, ] <- ressivipm
         
    }
    varnames <- names(ressivipm)
    IC.inf <- rep(0, ncol(boot))
    IC.sup <- rep(0, ncol(boot))
    
    for (i in 1:ncol(boot)) {
        IC.inf[i] <- sort(boot[, i])[ceiling(B * (alpha/2))]
        IC.sup[i] <- sort(boot[, i])[ceiling(B * (1 - alpha/2))]
    }
    IC <- cbind(IC.inf, IC.sup)
    if (graph) {
      taillept <- 0.8 # taille des annotations
      colnames(boot) <- colnames(XIndic@Pindic@indic)
      # search for the median
         bb <- boxplot(boot, plot=FALSE)
      # tri selon la mediane
      bootsort <- boot[, order(bb$stats[3,])]
      boxplot(bootsort,cex.axis=taillept, cex.names=taillept, cex.lab=taillept, las=3)
       }
    rownames(IC) <- varnames
    return(IC)
} # end sivipboot
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compute the PLS coefficents
# In case it is useful
# beta can be the natural beta or the centered-reducted ones
calccoefPLS <- function(beta, Indic, nvar, nY) {
  coefPLS <- matrix (NA, nrow=nvar, ncol=nY)
  for (j in 1:nY) {
    for(i in 1:nvar) {
      a <- Indic[,i] * beta[,j]
      coefPLS[i,j] <- sum(a*a)
    }
  }
  return(coefPLS)
} # end calccoefPLS 
