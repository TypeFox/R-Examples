"useFP" <- function
(
 obj,
 data
 )
### alternative idea to get plot fp terms etc ... missing a lot ...
### 2008-11
{
        browser()
        obj.full <- decomposeSurv(obj$formula, data, sort=FALSE)
        ## -> $fac,resp,mm1,NTDE,timedata,timeind,NFP,fpnames,fpind,PTcoefs,covnames,ind,offset.values

        ## take obj.full$mm1[, obj$ind]


        CODE <- paste("(I(PT(\\1)) + powM2(PT(\\1)) + powM1(PT(\\1)) + powM0.5(PT(\\1)) + log(PT(\\1)) + ",
                      "sqrt(PT(\\1)) + pow2(PT(\\1)) + pow3(PT(\\1)) + ",
                      "RI(PT(\\1)) + RpowM2(PT(\\1)) + RpowM1(PT(\\1)) + RpowM0.5(PT(\\1)) + ",
                      "Rlog(PT(\\1)) + Rsqrt(PT(\\1)) + Rpow2(PT(\\1)) + Rpow3(PT(\\1)) )")
        sub3 <- gsub("fp\\(([a-zA-Z0-9]*)\\)", CODE, obj$formula[3])
        formula <- as.formula(paste(as.character(obj$formula)[2], "~", sub3))

        ## define simple transformations
        powM2 <- function(z) z^(-2)
        powM1 <- function(z) z^(-1)
        powM0.5 <- function(z) z^(-0.5)
        pow2 <- function(z) z^2
        pow3 <- function(z) z^3

        ## define repeated powers
        RI <- function(z) z * log(z)
        RpowM2 <- function(z) z^(-2) * log(z)
        RpowM1 <- function(z) z^(-1) * log(z)
        RpowM0.5 <- function(z) z^(-0.5) * log(z)
        Rlog <- function(z) log(z) * log(z)
        Rsqrt <- function(z) sqrt(z) * log(z)
        Rpow2 <- function(z) z^2 * log(z)
        Rpow3 <- function(z) z^3 * log(z)

        ##
        browser()

}

fp.scale <- function
(
 x
 )
### taken from package <mfp>
### 2008-06
{
        scale <- 1
        shift <- 0
        if (min(x) <= 0) {
                z <- diff(sort(x))
                shift <- min(z[z > 0]) - min(x)
                shift <- ceiling(shift * 10)/10
        }
        range <- mean(x + shift)
        scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
        
        list(
             shift=shift,
             scale=scale
             )
}



"decomposeSurv" <- function
(
 formula,
 data,
 sort=FALSE,
 offset=NULL
 )
### decomposes complex survival formula
### trans: I(), powM2, powM1, powM0.5, sqrt, pow2, pow3, log,
### 2008-04
{
	orig.formula <- formula

	## expand formula if needed:
        repeat {
                terms <- terms(formula, "fp", data=data)
		fac <- attr(terms, "factors")

		needed <- rownames(fac)[!rownames(fac) %in% colnames(fac)][-1]

		if(length(needed) == 0) break

		formula <- as.formula(paste(as.character(formula)[2], "~",
                                            as.character(formula)[3], "+",
                                            paste(needed, sep="+")))
	}

        ## determine position of <fp> terms (of all orders) in
        ## not expanded formula
        spec <- untangle.specials(terms, "fp", order=1:3)
        fpnames <- gsub("fp\\(([a-zA-Z0-9]*)\\)", "\\1", spec$vars) # c("a", "b")
        fpiden <- paste("(", fpnames, ")", sep="") # c("(a)", "(b)")

        formulaOrig <- formula
        ## replace fp() by the 8x2 possible terms (in formula[3])
        ## attention: fp-variables must consist of a-z,A-Z,0-9, otherwise
        ## must be made modifications here (gsub) ...
        ## FURTHER: 1st term is the linear one
        ## FURTHER: 1st half (=8) are the usual terms, 2nd half (9-16) are the repeated powers
        ## FURTHER: problems possible if one variable name is part of the other, e.g. ABC & ABCDE
        ##CODE <- paste("(I(\\1) + powM2(\\1) + powM1(\\1) + powM0.5(\\1) + log(\\1) + sqrt(\\1) + pow2(\\1) + pow3(\\1) + ",
        ##              "  RI(\\1) + RpowM2(\\1) + RpowM1(\\1) + RpowM0.5(\\1) + Rlog(\\1) + Rsqrt(\\1) + Rpow2(\\1) + Rpow3(\\1) )")
        CODE <- paste("(I(PT(\\1)) + powM2(PT(\\1)) + powM1(PT(\\1)) + powM0.5(PT(\\1)) + log(PT(\\1)) + ",
                      "sqrt(PT(\\1)) + pow2(PT(\\1)) + pow3(PT(\\1)) + ",
                      "RI(PT(\\1)) + RpowM2(PT(\\1)) + RpowM1(PT(\\1)) + RpowM0.5(PT(\\1)) + ",
                      "Rlog(PT(\\1)) + Rsqrt(PT(\\1)) + Rpow2(PT(\\1)) + Rpow3(PT(\\1)) )")
        sub3 <- gsub("fp\\(([a-zA-Z0-9]*)\\)", CODE, formulaOrig[3])
        formula <- as.formula(paste(as.character(formula)[2], "~", sub3))

        ## define simple transformations
        powM2 <- function(z) z^(-2)
        powM1 <- function(z) z^(-1)
        powM0.5 <- function(z) z^(-0.5)
        pow2 <- function(z) z^2
        pow3 <- function(z) z^3

        ## define repeated powers
        RI <- function(z) z * log(z)
        RpowM2 <- function(z) z^(-2) * log(z)
        RpowM1 <- function(z) z^(-1) * log(z)
        RpowM0.5 <- function(z) z^(-0.5) * log(z)
        Rlog <- function(z) log(z) * log(z)
        Rsqrt <- function(z) sqrt(z) * log(z)
        Rpow2 <- function(z) z^2 * log(z)
        Rpow3 <- function(z) z^3 * log(z)

        ## pretransformation function
        PT <- function(z) {
                obj <- fp.scale(z)
                (z + obj$shift) / obj$scale }

	## construct 3-col response:
	resp <- model.extract(model.frame(formula, data = data), "response")
	if(is.null(dim(resp))) 	 resp<-cbind(resp-min(min(resp-0.01),0), rep(1,length(resp)))
	if(ncol(resp) == 2)      resp <- cbind(start=rep(0, nrow(resp)), resp)

	## sortieren nach STOPzeit und -Cens
	if(sort) {
                sort <- order(resp[, 2],  -resp[, 3])
                data <- data[sort, , drop=FALSE]
		resp <- resp[sort, ]
	}

	mm <- model.matrix(formula, data = data) ## Model-Matrix
        mm1 <- mm[, -1, drop=FALSE]	# w/o intercept

        ## *** RETREIVE PRETRANS COEFS
        sub3 <- gsub("fp\\(([a-zA-Z0-9]*)\\)", "PT(\\1)", formulaOrig[3])
        formulaPT <- as.formula(paste(as.character(formulaOrig)[2], "~", sub3))
        ## varied pretransformation function to get coefficients:
        PT <- function(z) {
                obj <- fp.scale(z)
                c(shift=obj$shift, scale=obj$scale, rep(-999, length(z) - 2)) }
        mmPT <- model.matrix(formulaPT, data = data)
        PTcoefs <- mmPT[1:2, !is.na(mmPT[1,]) & (mmPT[3,]==-999), drop=FALSE]
        rownames(PTcoefs) <- c("shift", "scale")
        ## *** END

        ## offset ..
        if(length(offset) != 0)
          offset.values <- offset
        else
          offset.values <- NA

        terms <- terms(formula, "fp", data=data)
	fac <- attr(terms, "factors")
	labels <- attr(terms, "term.labels")

	## splits by special chars
	f <- function(str)
          for(chars in c("(", ")", ":", " ", ",", "*", "^"))
            str <- unlist(strsplit(str, split=chars, fixed=TRUE))

	rowSplit <- sapply(rownames(fac), f, simplify=FALSE)	# splitted effects
	stopName <- tail(rowSplit[[1]], 2)[1]	# name of stoptime
	rowInter <- unlist(lapply(rowSplit[-1], function(z) any(z == stopName)))
        ##  rowOffset <- unlist(lapply(rowSplit[-1], function(z) any(z == "offset")))
        ##  rowOffset


	fac <- fac[-1, , drop=FALSE]	# omit Surv

	colSplit <- lapply(colnames(fac), f)
	colInter <- unlist(lapply(colSplit, function(z) any(z == stopName)))

	nTimes <- colSums(fac[rowInter, , drop=FALSE])
	nFac   <- colSums(fac[!rowInter, , drop=FALSE])

	inters <- (nFac>0) & (nTimes>0)
	NTDE <- sum(inters)


	timedata <- matrix(0, nrow(data), 0)
	timeind <- c()

	## loop for (time x effect)
	for(i in which(inters)) {
		## search pure time:
		ind <- (colSums(fac[rowInter, i] != fac[rowInter, , drop=FALSE]) == 0) & (nFac==0)
		timedata <- cbind(timedata, mm1[, ind, drop=FALSE])

		## search pure effect:
		ind <- (colSums(fac[!rowInter, i] != fac[!rowInter, , drop=FALSE]) == 0) & (nTimes == 0)
		timeind <- c(timeind, which(ind[!colInter]))
	}
	mm1 <- mm1[, !colInter, drop=FALSE]

	covnames <- c(colnames(mm1),
                      paste(colnames(timedata), colnames(mm1)[timeind], sep=":")
                      )
        alter <- c(colnames(mm1),
                   paste(colnames(mm1)[timeind], colnames(timedata), sep=":")
                   )

	## indicator to identify the original formula:
	ind <- covnames %in% colnames(attr(terms(orig.formula, "fp", data=data), "factors")) |
               alter %in% colnames(attr(terms(orig.formula, "fp", data=data), "factors"))


        ## FP indicator matrix (1:16)
        NFP <- length(fpnames)
        fpind <- matrix(0, NFP, length(covnames), dimnames=
                        list(fpnames, covnames))
        for(i in seq(length=NFP)) {
                inds <- grep(fpiden[i], covnames)
                ## when interactions occur, 1:16,1:16 instead of 1:32 is needed
                fpind[i, inds] <- ((seq(along=inds) - 1) %% 16) + 1
        }

        ## return object
	list(
             fac=fac,                    # factor matrix ..
             resp=resp,                  # N x 3 - response matrix
             mm1=mm1,                    # model matrix without time effects

             NTDE=NTDE,                  # number time dep. effects
             timedata=timedata,          # matrix with time functions as columns
             timeind=timeind, 		       # indicator of time-dependend effects

             NFP=NFP,                    # number of frac.polys
             fpnames=fpnames,            # names of the variables used as fractional poly's
             fpind=fpind,                # matrix with frac.polyn. in each row, numbered by 1 to 5
             PTcoefs=PTcoefs,            # coefficients of pretransformations

             covnames=covnames,          # names of covariates
             ind=ind,			               # indicator vector:
                                         # which terms are really part of the formula
             offset.values=offset.values # offset values
             )
}




`concreg` <- function                                                        #! ändern!!!
(
 formula=attr(data, "formula"),         # formula
 data=sys.parent(),                     #
 id=NULL,                               # identifier: numeric or character or factor
 normalize=TRUE,
 scale.weights=1,
 offset=NULL,
 alpha=0.05,                            # confidence limit
 maxit=50,                              # max. iterations
 maxhs=5,                               # half steps
 epsilon=1e-6,                          #
 maxstep=2.5,                           #
 x=TRUE,                                # für Ausgabe
 y=TRUE,                                # für Ausgabe
 print=TRUE,                            # print fitting information on screen
 c.risk=NULL,                           # Vektor mit Competing Risk 0, 1=event, 2=competing risk (status ist normal!!) NEU
 strata.var=NULL,                       # Namen der Stratifizierungsvariablen (Daten selbst sind bei data dabei) NEU
 trunc.weights=1,                       # quantile for weight truncation: all weights greater than that quantile will be truncated to that value
 npar=FALSE,                            # nonparametric estimation?
 # breslow=NA, # righthand formula, if breslow weighted terms, e.g. ~ A + A:C
 # prentice=NA, # righthand formula, if prentice weighted terms, e.g. ~ A + A:C
 # taroneware=NA, # righthand formula, if tarone-ware weighted terms, e.g. ~ A + A:C
 # censcorr=FALSE,
 # robust=FALSE,
 # jack=FALSE,
 ...                                    # N -> breslow; km -> prentice
 )
### by MP und GH, 2008
### Surv-Objekt entweder usual  (z.B. Surv(time,event)~A+B+C+D)
### event: 2=competing risk, 1=dead/event, 0=censored
{
        ## evt. use of abbrev. parameter names
#        l <- list(...)
#        if("N" %in% names(l)) breslow <- l$N
#        if("km" %in% names(l)) prentice <- l$km
 alpha.fp=c(0.20, 0.05)
 fp.iter=10                            # maximum number of iterations of large <fp> loop
       n <- nrow(data)
#        pl <- FALSE
#        robcov <- if(robust) 1 else 0

        ## generate or reorder id's such that values are within 1:n
        if (is.null(id)) { id <- 1:n } else id <- as.numeric(as.factor(id))
        maxid <- max(id)

        ## here only ONCE the full model matrix is spanned with all possible fp-effects
	      obj.full <- decomposeSurv(formula, data, sort=FALSE, offset)

  # change Daniela Competing Risk
  if (is.null(c.risk)) {
    crisk <- obj.full$resp[,3]
    obj.full$resp <- cbind(obj.full$resp, crisk)

  }

  if (is.null(c.risk)==FALSE) {
    crisk <- obj.full$resp[,3]
    c.risk[c.risk==1] <- 2
    crisk[crisk==0] <- c.risk[crisk==0]
    obj.full$resp <- cbind(obj.full$resp, crisk)
  }


  ## stratagruppen bilden
  if (is.null(strata.var)[1]) {
    obj.full$stratum <- rep(1, n)
  }

  if (is.null(strata.var)[1] == FALSE) {
    mult.strata <- 1
    if (length(strata.var) > 1)  {
      for (l in 1:length(strata.var)) { mult.strata <- c(mult.strata, mult.strata[length(mult.strata)]*10) }
    }

    daten.tmp <- data.frame(data[,strata.var])
    for (l in 1:length(strata.var)) {
      anz.strata <- length(table(data[,strata.var[l]]))

      for (l1 in 1:anz.strata)  {                                               # Strata neu durchnummerieren
        daten.tmp[data[,strata.var[l]] == as.numeric(names(table(data[,strata.var[l]])))[l1], l] <- l1
      }

     daten.tmp[, l] <- daten.tmp[, l]*mult.strata[l]
    }

    if (length(strata.var)==1)  { daten.tmp$tmp <- daten.tmp  }
    if (length(strata.var)> 1)  { daten.tmp$tmp <- apply(daten.tmp[, strata.var], 1, sum) }
    tab <- table(daten.tmp$tmp)

#    if (sum(tab==1)!=0) { stop("number of observations is in at least 1 stratum 1")}     #! ev. einfach diese Daten weglassen
    anz.strata <- length(tab)
    for (l2 in 1:anz.strata)  {
      daten.tmp$stratum[daten.tmp$tmp == as.numeric(names(tab)[l2])] <- l2
    }
    obj.full$stratum <- daten.tmp$stratum
  }

  max.strata <- max(obj.full$stratum)

  # change Daniela Competing Risk Weights (getrennt nach Strata)
  # geht nicht so einfach, Problem mit ties (gleiche Zeiten nur 1 mal in der Liste)
  if(!is.null(c.risk))   #### Georg 110517
  {
   vgl <- G <- data.frame(cbind(obj.full$resp[,2], obj.full$stratum, -99))
   G <- unique(G)
   for (i.strata in 1: max.strata) {
    g.tmp <- rep(NA, times=sum(obj.full$stratum==i.strata))
    g.tmp[obj.full$resp[obj.full$stratum==i.strata,3]==0] <- 1
    g.tmp[obj.full$resp[obj.full$stratum==i.strata,3]==1 | obj.full$resp[obj.full$stratum==i.strata,3]==2] <- 0

    G[G$X2==i.strata,3] <- survfit(Surv(obj.full$resp[obj.full$stratum==i.strata,2], g.tmp)~1)$surv
  }

  if (nrow(G) != nrow(vgl))  {                        # Bei Bindungen innerhalb der Strata
    G.n <- 1
    for (i.n in 1:n)  {
      if (sum(vgl[i.n, 1:2] == G[G.n, 1:2])==2) {
        vgl[i.n, 3] <- G[G.n,3]
        G.n <- G.n+1
      }
    }

    vgl[order(vgl[,3])[1:sum(vgl[,3] == -99)],3] <- vgl[order(vgl[,3])[1:sum(vgl[,3] == -99)]-1,3]
    G <- vgl
  }
  G <- G[,3]

  } else G<-rep(1,nrow(obj.full$resp))     #### Georg 110517
  # neue Time Variable für competing risk
  time.crisk <- obj.full$resp[,2]
  time.crisk[obj.full$resp[,"crisk"]==2] <- max(obj.full$resp[,2])+1
  obj.full$resp <- cbind(obj.full$resp, time.crisk)


  # calculate weights
  W <- concreg.wei(resp=obj.full$resp, max.strata=max.strata, stratum=obj.full$stratum, trunc.weights=trunc.weights)


        obj <- obj.full

        kk <- ncol(obj.full$mm1) # should be k2 - NTDE                          #! weg? und dafür k nach vorne holen
        obj$mm1 <- obj.full$mm1[, obj$ind[1:kk], drop=FALSE]
        obj$covnames <- obj.full$covnames[obj$ind]

        obj$timeind  <- obj.full$timeind[obj$ind[-(1:kk)]]                      #! weg?
        obj$timedata <- obj.full$timedata[, obj$ind[-(1:kk)], drop=FALSE]       #! weg?
        ## re-index $timeind
        obj$timeind <- match(obj$timeind, (1:kk)[obj$ind[1:kk]])                #! weg?
        NTDE <- obj$NTDE <- length(obj$timeind)                                 #! weg?


        ind.offset <- sum(length(offset) != 0)

        k <- ncol(obj$mm1)    # number covariates w/o time-dep effects          # Anzahl Variablen
        k2 <- k + NTDE                                                          #! weg?
        if (npar) nonpar<-1
        else nonpar<-0
#! weg?: W$NGV, NTDE
#! hinzufügen: obj.full$stratum, W, G
        PARMS <- c(n, k, 0, maxit, maxhs, maxstep, epsilon,  0, 0, 0, 0, 0, 0, max.strata, maxid, ind.offset, nonpar)


  # alles sortieren (obj, W, G, id ist sortiert)
  ord <- order(obj.full$stratum, obj.full$resp[,"time.crisk"],(1-obj.full$resp[,3]))
  obj$mm1     <- obj$mm1[ord,]                                                  # ist jetzt keine Matrix mehr. stört das?
  obj$resp    <- obj$resp[ord,]
  obj$stratum <- obj$stratum[ord]
  G           <- G[ord]
  W           <- W[ord]
  id          <- id[ord]

        ##   if (offset) {
        ##    IOARRAY[1,1]<-0    # first variable is offset
        ##    IOARRAY[2,1]<-Z.sd[1]    # first variable is offset
        ##   }
        ## **************** fit model *********************
        value0 <- concreg.fit(obj=obj, id=id, W=W, G=G, PARMS=PARMS, npar=npar)            #! Aufruf der Funktion mit dem Fortran Aufruf

#        cov.ls <- value0$cov.ls
#        cov.lw <- cov.j <- dfbeta.resid <- NULL

#        if(robust) {
#                cov.lw <- matrix(value0$outtab[(k2+4):(3+2*k2), ], ncol=k2) / value0$ZxZ
#                dfbeta.resid <- value0$dfbetaresid / matrix(value0$Z.sd, maxid, k2, byrow=TRUE)
#                covs <- cov.lw
#                cov.method <- "Lin-Wei"
#
#        } else if(jack) {
#                dfbeta.resid <- matrix(0, maxid, k2)
#                nc <- ncol(value0$cards)
#                for(iid in 1:maxid) {
#                        cardsJ <- value0$cards[id != iid,]
#                        cardsJ[id[id!=iid]>iid,nc] <- cardsJ[id[id!=iid]>iid,nc]-1
#                        n.jack <- nrow(cardsJ)
#                        PARMS <- c(n.jack, k, 0, maxit, maxhs,maxstep,epsilon, 0, 0, 0, 0, 0, W$NGV, NTDE, maxid-1, ind.offset)
#                        valueJ <- coxphw.fit(obj, id, W$weights, PARMS, cardsJ) # fit model jack(i) *******
#                        dfbeta.resid[iid,] <- value0$coef.orig - valueJ$coef.orig
#                }
#                cov.j <- (maxid-1) / maxid * crossprod(dfbeta.resid) / value0$ZxZ
#                dfbeta.resid <- dfbeta.resid / matrix(value0$Z.sd, maxid, k2, byrow=TRUE)
#                covs <- cov.j
#                cov.method <- "Jackknife"
#
#        } else {
#                cov.method <- "Lin-Sasieni"
#                covs <- cov.ls
#        }

        vars <- as.matrix(diag(value0$cov.rob))
#        dimnames(vars) <- list(obj$covnames, obj$covnames)

        ## DECIDE THE MODEL: USE PVALS ##############
        probs <- 1 - pchisq((value0$coefs^2/vars), 1)

        ## ########## NOW ONLY FINAL MODEL IS CONSIDERED ############
        names(value0$coefs) <- obj$covnames
        if(value0$outpar[10]>=maxit)
          cat("No convergence attained in ", value0$outpar[10], " iterations.\n", sep="")

        Means <- colMeans(value0$mmm)

        ## return object
        fit <- list(coefficients = value0$coefs,     # coefficients of the fit
                    cards    = value0$cards,         #
                    parms    = value0$outpar,
                    ioarray  = value0$outtab,
#                    dfbeta.resid = dfbeta.resid,
                    alpha    = alpha,                # significance level
                    var      = value0$cov.rob,                 # covariance matrix
                    df       = k2,                   # degrees of freedom (k + NTDE)          #! ändern
                    iter     = value0$outpar[10],
                    method.ties = "no",         #
                    n = n,                           # number observations
                    ##terms = terms(formula),
                    y = obj$resp,                    # responses
                    formula = formula,               # original formula
                    exit.code=value0$outpar[8],

  #                  fpind=obj$fpind,
  #                  PTcoefs=obj.full$PTcoefs,                                   #! weg?
  #                  ind=obj$ind,                                                #! was ist das?

                    call    = match.call(),
                    cov.mb  = value0$cov.mb,
#                    cov.j   = cov.j,              # jackknife covariance matrix
#                    cov.lw  = cov.lw,             #
#                    cov.ls  = cov.ls,             #
#                    cov.method=cov.method,
#                    w.matrix= W$w.matrix,         # weight matrix
#                    Wald    = value0$outpar[9],
                    Wald =  (t(value0$coefs) %*% solve(value0$cov.rob)) %*% value0$coefs,
                    means   = Means,               # means of <X> (model matrix)
                    linear.predictors= as.vector(scale(value0$mmm, Means, scale=FALSE) %*% value0$coefs),
                    method  = "Weighted Estimation",
                    method.ci= "Wald",             #
                    ci.lower= exp(value0$coefs + qnorm(alpha/2) * vars^0.5), #
                    ci.upper= exp(value0$coefs + qnorm(1 - alpha/2) * vars^0.5), #
                    prob    = probs,               # p-values
                    G       = G,
#                    W       = W,
                    W       = cbind(obj$stratum,obj$resp[,2],W)[obj$resp[,3]==1,],
                    offset.values= obj$offset.values, #
                    x       = if(x) obj$mm1 else NA,   # return original model matrix if requested
                    npar = npar
                    )

#        ## if all weights are the same ...
#        if(W$const) {
#                fit$loglik <- value0$outpar[12:11]
#                fit$score <- value0$outpar[7]
#        }

        ##        if(offset) {
        ##          obj$covnames[1]<-as.character(paste("offset(",obj$covnames[1],")"))
        ##          fit$Wald <- t(coefs[2:k2]) %*% solve(covs[2:(k2),2:(k2)]) %*% coefs[2:(k2)]
        ##          if(x) {
        ##           fit$x <- mm1.orig[,2:k2]
        ##           fit$offset <- mm1.orig[,1]
        ##          }
        ##        }
        names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- obj$covnames
        attr(fit, "class") <- c("concreg")                 #! ändern
        fit
}

coef.concreg<-function(object, ...) object$coefficients

vcov.concreg<-function(object, ...) object$var

cindex<-function(object, confint=FALSE, level=0.95) {
  OC<-exp(coef(object))
  cindex<-OC/(1+OC)
  if(confint) {
    ci<-confint(object, level=level, what="cindex")
    cindex<-rbind(cindex,ci)
    }
  return(cindex)
  }
  
confint.concreg<-function(object, parm, level=0.95, what="coefficients", ...){   
  if(what != "coefficients" & what !="OC" & what !="cindex") {
    what<-"coefficients"
    warning("Will produce confidence intervals for regression coefficients.\n")
    }
  if(missing(parm)) parm<-1:length(coef(object))
  se<-sqrt(diag(vcov(object)))
  ci.coef<-(coef(object)+se*cbind(qnorm((1-level)/2), qnorm(1-(1-level)/2))  )[parm,]
  if(what=="coefficients") return(ci.coef)
  else if (what=="OC") return (exp(ci.coef))
  else if (what=="cindex") return (exp(ci.coef)/(1+exp(ci.coef)))
  }


"plotw" <- function                                                             #! noch anpassen
(
 x,    # object of class coxphw
 rank=FALSE,
 log=FALSE,
 xlim=c(0,max(time)),
 
  ...            # dummy
 )
 {
   if(rank) {
    time<-order(x$W[,2])
    label<-"Ranked time"
    }
   else {
    time<-x$W[,2]
    label<-"Time"
    }
   if(log) {
    weights<-log(x$W[,3])
    wlabel<-"Log of weight"
   }
   else {
    weights<-x$W[,3]
    wlabel<-"Weight"
    }
#   if (is.na(xlim)) xlim=c(0,max(time))
   ltyi<-1
   for (i in unique(x$W[,1])) {
    if (ltyi==1) plot(time[x$W[,1]==i], weights[x$W[,1]==i],type="l",lty=ltyi,ylim=c(min(weights[time<=xlim[2]]),max(weights[time<=xlim[2]])), xlab=label, ylab=wlabel, xlim=xlim)
    else lines(time[x$W[,1]==i], weights[x$W[,1]==i],lty=ltyi,ylim=c(min(weights[time<=xlim[2]]),max(weights[time<=xlim[2]])), xlab=label, ylab=wlabel, xlim=xlim)
    ltyi<-ltyi+1
    }
#   lines(time, weights,lty=2)
#   lines(time, weights,lty=3)
   legend(min(time),0.95*max(weights[time<=xlim[2]]),unique(x$W[,1]),lty=1:(ltyi-1), lwd=1)
 }


"print.concreg" <- function                                                   #! noch anpassen
(
  x,     # object of class coxphw
  ...            # dummy
 )
### MP and GH
### overall test is score test
### 2007-07
{
        print(x$call)
        cat("Model fitted by", x$method, "\n\n")
        if(x$npar) cat("Nonparametric estimation \n")
        se<- diag(x$var)^0.5
        out <- cbind(x$coefficients, se, exp(x$coefficients),
                     x$ci.lower, x$ci.upper, x$coefficients/se, x$prob)
        dimnames(out) <- list(names(x$coefficients),
                              c("coef", "se(coef)", "exp(coef)",
                                paste(c("lower", "upper"), 1 - x$alpha), "z", "p"))

        if (x$method.ci != "Wald")
          dimnames(out)[[2]][6] <- "Chisq"
        print(out)

#        if("loglik" %in% names(x)) cat("\nScore test=", x$score, " on ", x$df,
#            " df, p=", 1 - pchisq(x$score, x$df), ", n=", x$n, "", sep = "")
        cat("\nWald test=", x$Wald, " on ", x$df, "df, p=", 1 - pchisq(x$Wald, x$df), ", n=", x$n, "\n\n", sep = "")

        invisible(x)
}


"summary.concreg" <- function                                                 #! noch anpassen
(
 object,              # object of class coxphf
 ...                  # dummy
 )
### MP and GH
### 2007-07
{
        print(object$call)
        cat("\nModel fitted by", object$method, "\n\n")
        if(object$npar) cat("Nonparametric estimation \n")
        ##cat("Confidence intervals and p-values by", object$method.ci, "\n\n")
        se<-diag(object$var)^0.5
        c.index<-exp(object$coefficients)/(1+exp(object$coefficients))
        c.index.lower<-object$ci.lower/(1+object$ci.lower)
        c.index.upper<-object$ci.upper/(1+object$ci.upper)
        out <- cbind(object$coefficients, se,
                     exp(object$coefficients),
                     object$ci.lower, object$ci.upper,
                      c.index, c.index.lower, c.index.upper, object$coefficients/se, object$prob)
        dimnames(out) <- list(names(object$coefficients),
                              c("coef", "se(coef)", "exp(coef)",
                                paste(c("lower", "upper"), 1 - object$alpha), "c", "c (lower)", "c (upper)", "z", "p"))
        if (object$method.ci != "Wald")
          dimnames(out)[[2]][6] <- "Chisq"

        print(out)

        if("loglik" %in% names(object)) {
                LL <- 2 * diff(object$loglik)
                cat("Likelihood ratio test=", LL, " on ", object$df,
                    " df, p=", 1 - pchisq(LL, object$df), ", n=", object$n, "\n", sep = "")
                cat("\nScore test=", object$score, " on ", object$df,
            " df, p=", 1 - pchisq(object$score, object$df), ", n=", object$n, "\n", sep = "")
        }
#        wald.z <- t(coef(object)) %*% solve(object$var) %*% coef(object)
        cat("Wald test =", object$Wald, "on", object$df, " df, p =",
            1 - pchisq(object$Wald, object$df))
        cat("\n\nCovariance-Matrix:\n")
        print(object$var)

        invisible(object)
}


concreg.fit <- function(obj, id, W, G, PARMS, npar)                                   #! für was braucht man CARDS=NULL?
{
  # fitter function

#        k <- ncol(obj$mm1)    # number covariates w/o time-dep effects         #! in PARMS
        k <- PARMS[2]
        k2 <- k + obj$NTDE                                                      #! weg?
#        maxid <- max(id)                                                       #! in PARMS
        maxid <- PARMS[15]

        ## standardize model matrix, but only in semiparametric mode
        if(!npar) {
         sd1 <- apply(as.matrix(obj$mm1),2,sd)
         sd2 <- apply(as.matrix(obj$timedata),2,sd)                                                 #! weg?
         Z.sd <- c(sd1, sd2 * sd1[obj$timeind])                                  #! weg?
         ZxZ <- as.matrix(Z.sd) %*% t(as.matrix(Z.sd)) # used often to restandardize ...    #! weg?
         obj$mm1 <- scale(obj$mm1, FALSE, sd1)
         } else {
          sd1<-1
          sd2 <- 1
          Z.sd<-1
          ZxZ <-1
          }
         

        ##if(ind.offset)
        obj$mm1o <- if(PARMS[16] != 0) cbind(obj$offset.values, obj$mm1) else obj$mm1

#        if(is.null(CARDS))
           CARDS <- cbind(obj$mm1o, obj$resp[,c(2, 4)], W, G, id, obj$stratum)
#          CARDS <- cbind(obj$mm1o, obj$resp, weights, obj$timedata, id)        #! was ist weights?

        if(!npar) obj$timedata <- scale(obj$timedata, FALSE, sd2)                         # weg?
        mmm <- cbind(obj$mm1, obj$timedata) # model matrix inc. time data       #! ohne timedata?
        ##   if (offset) {
        ##    IOARRAY[1,1]<-0    # first variable is offset
        ##    IOARRAY[2,1]<-Z.sd[1]    # first variable is offset
        ##   }
        DFBETA <- matrix(0, maxid, k2)                                          #! k2 durch k ersetzen?
        IOARRAY <- rbind(rep(1, k2), matrix(0, 2+2*k2, k2))                     #! k2 durch k ersetzen?
        if(obj$NTDE >0)                                                         #! weg?
          IOARRAY[4, (k+1):k2] <- obj$timeind                                   #! weg?

        ## --------------- Aufruf Fortran-Routine ????????? ------------
        storage.mode(CARDS) <- storage.mode(PARMS) <- storage.mode(IOARRAY) <- storage.mode(DFBETA) <- "double"
#       dyn.load("concreg_dll.dll")
#        dyn.load("D:\\WORK\\concreg_dll.dll")
        value <- .Fortran("CONCREG",                                            #! anpassen
                          cards=CARDS,
                          outpar = PARMS,
                          outtab = IOARRAY)
#                          PACKAGE=concreg.fp)
        if(value$outpar[8])
          warning("Error in Fortran routine concreg; parms8 <> 0")

        coefs <- value$coefs / Z.sd                                             #! value$coefs gibt es nicht
        cov.mb <- matrix(value$outtab[4:(k2+3), ], ncol=k2) / ZxZ            # "model-based" varianz
        cov.rob <- matrix(value$outtab[(3+k2+1):(3+2*k2), ], ncol=k2) / ZxZ  # "robuste" varianz (standard)

        res <- list(                                                            #! anpassen
                    cards=value$cards,
                    outpar=value$outpar,
                    outtab=matrix(value$outtab, nrow=3+2*k2),
 #                   dfbetaresid=value$dfbetaresid,
                    coef.orig=value$outtab[3,  ],
                    coefs=value$outtab[3,  ] / Z.sd, # coefficients
                    cov.rob=cov.rob,                 # covariances
                    cov.mb=cov.mb,
                    Z.sd=Z.sd,
                    ZxZ=ZxZ,                                                    #! weg?
                    mmm=mmm                          # model matrix
                    )
        res
}


#! Berechne das Gewicht für die logistische Regression: W
                                     
concreg.wei <- function(resp, max.strata, stratum, trunc.weights)
{
  #! function to calculate the weigth for concreg
  #! rechnet mit time.crisk

  data.tmp <- data.frame(resp[,c(5,3,2,4)])                                         #! achtung, falls start bei resp wegkommt, hier ändern
  dimnames(data.tmp)[[2]][2] <- "status"
  dimnames(data.tmp)[[2]][3] <- "time.orig"
  dimnames(data.tmp)[[2]][4] <- "status012"
  data.tmp$statusfu<-(data.tmp$status==0)+0
  
  W <- rep(0, nrow(data.tmp))
  if(all(data.tmp[,2]==1)) return(rep(1, nrow(data.tmp)))

  for (i.max.strata in 1:max.strata)  {
    data.strata <- data.tmp[stratum==i.max.strata,]

    n.strata <- nrow(data.strata)
    data.strata$pat <- 1:n.strata

#    attach(data.strata, warn.conflicts=FALSE)         # survfit geht sonst nicht aus irgendeinem Grund
    fit <- survfit(Surv(time.crisk, status)~1, data=data.strata)           #pseudo survival to reconstruct number of pairs
    fit.fu<-survfit(Surv(time.orig, statusfu)~1, data=data.strata)         #follow-up KM
#    detach(data.strata)
    fsurv<-fit$surv

#   if (length(fsurv) == n.strata)  {                        # keine Bindungen
#      fsurv <- c(1, fit$surv)[-(length(fit$surv)+1)]             
#     data.strata$surv <- fsurv
#      data.strata$n.risk <- fit$n.risk
#    }

 #   if (length(fsurv) != n.strata)  {                        # Bei Bindungen innerhalb der Strata
      data.strata$surv <- -99
      data.strata$n.risk <- -99
      data.strata$gsurv <- -99

      fit.n <- 1
      fsurv<-c(1,fit$surv)
      fn.risk<-c(n.strata, fit$n.risk)
      gsurv<-c(1,fit.fu$surv)
      gtime<-c(0,fit.fu$time)
      ftime<-c(0,fit$time)
      for (i.n in 1:n.strata)  {
        if (data.strata$status[i.n]==1)  indices<-sum(ftime<data.strata[i.n, "time.crisk"])
        if (data.strata$status[i.n]==0)  indices<-sum(ftime<=data.strata[i.n, "time.crisk"])
        index.fu<-sum(gtime<data.strata[i.n, "time.crisk"])
        data.strata[i.n, c("surv", "n.risk", "gsurv")] <- c((fsurv[indices]), (fn.risk[indices]), gsurv[index.fu])
#        if (data.strata[i.n, "time.crisk"] == fit$time[fit.n]) {
#          fit.n <- fit.n+1
        }
#      }

#    data.strata[data.strata[order(data.strata[,"surv"])[1:sum(data.strata[,"surv"] == -99)],3], c("surv", "n.risk")] <-
#    data.strata[data.strata[order(data.strata[,"surv"])[1:sum(data.strata[,"surv"] == -99)]-1,3], c("surv", "n.risk")]
    

#    W <- c(W, (n.strata*data.strata[,"surv"] / (data.strata[,"n.risk"]-1))* (data.strata[,"surv"] / (data.strata[,"n.risk"])))
#    W <- c(W, ((n.strata*data.strata[,"surv"]-1) / (data.strata[,"n.risk"]-1))* (n.strata*data.strata[,"surv"] / (data.strata[,"n.risk"])))
    W.strata <- ((n.strata*data.strata[,"surv"]-1) / (data.strata[,"n.risk"]-1))/ data.strata[,"gsurv"]
    W[stratum==i.max.strata] <- W.strata

  }
#  W <- W[-1]
  W[W==Inf] <- 1
  W[is.na(W)] <- 1
  truncat<-quantile(W, probs=trunc.weights)
  W[W>truncat]<-truncat
  W[abs(W)==Inf]<-0
  W
}




