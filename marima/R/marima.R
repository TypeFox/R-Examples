##' @title marima
##'
##' @description Estimate multivariate arima and arima-x models.
##' Setting up the proper model for (especially) arima-x estimation
##' can be accomplished using the routine 'define.model' that can
##' assist in setting up the necessary autoregressive and moving average
##' patterns used as input to 'marima'.
##'
##' A more elaborate description of 'marima' and how it is used
##' can be downloaded from:
##' 
##' http://www.imm.dtu.dk/~hspl/marima.use.pdf
##' 
##' @param DATA time series matrix, dim(DATA) = c(kvar,n),
##' where 'kvar' is the dimension of the time series and 'n' is the
##' length of the series. If DATA is organized (n,kvar) (as a data.frame
##' e.g.) it is automatically transposed in marima, and the user need
##' not care about it. Also, and consequently, the output residuals and
##' fitted values matrices are both organised c(kvar,n) at return from marima.
##' The DATA is checked for completeness. Cases which include 'NA's or 'NaN's
##' are initially left out. A message is given (on the console) and the active
##' cases are given in the output object (...$used.cases).
##' @param ar.pattern autoregressive pattern for model (see define.model).
##' If ar.pattern is not specified a pure ma-model is estimated.
##' @param ma.pattern moving average pattern for model (see define.model).
##' If ma.pattern is not specified a pure ar-model is estimated. In this case
##' the estimation is carried by regression analysis in a few steps.
##' @param max.iter max. number of iterations in estimation (max.iter=50
##' is default which, generally, is more than enough).
##' @param means 0/1 indicator vector of length kvar, indicating
##' which variables in the analysis should be means adjusted or not.
##' Default: means=1 and all variables are means adjusted.
##' If means=0 is used, no variables are means adjusted.
##' @param weight weighting factor for smoothing the repeated
##' estimation procedure. Default is weight=0.33 which often works
##' well. If weight>0.33 (e.g. weight=0.66) is specified more damping
##' will result. If a large damping factor is used, the successive
##' estimations are more cautious, and a slower (but safer)
##' convergence (if possible) may result (max.iter may have to be
##' increased to, say, max.iter=75.
##' @param Plot 'none' or 'trace' or 'log.det' indicates a plot that shows
##' how the residual covariance matrix (resid.cov) develops with the
##' iterations.
##' If Plot= 'none' no plot is generated.
##' If Plot= 'trace' a plot of the trace of the residual
##' covariance matrix versus iterations is generated.
##' If Plot='log.det' the log(determinant) of the residual
##' covariance matrix (resid.cov) is generated. Default is Plot= 'none'.
##' @param Check (TRUE/FALSE) results (if TRUE) in a printout of some
##' controls of the call to arima. Useful in the first attemp(s) to use
##' marima. Default=FALSE.
##' @param penalty parameter used in the R function 'step' for
##' stepwise model reduction. If penalty=2, the conventional
##' AIC criterion is used. If penalty=0, no stepwise reduction of model is
##' performed. Generally 0<=penalty<=2 works well (especially penalty=1).
##' The level of
##' significance of the individual parameter estimates in the final
##' model can be checked by considering the (approximate) 'ar.pvalues'
##' and the 'ma.pvalues' calculated by marima.
##' 
##'
##' @return Object of class marima containing:
##'
##'  N            = N length of analysed series
##'
##'  kvar         = dimension of time series (all random and
##'    non-random variables.
##'
##'  ar.estimates = ar-estimates
##'
##'  ma.estimates = ma-estimates
##'
##'  ar.fvalues   = ar-fvalues (approximate)
##'
##'  ma.fvalues   = ma-fvalues (approximate)
##'
##'  ar.pvalues   = ar-p-values (approximate)
##'
##'  ma.pvalues   = ma-p-values (approximate)
##'
##'  residuals = estimated residuals (for used.cases)
##'
##'  fitted    = estimated/fitted values for all data
##' (including non random variables) (for used.cases)
##'
##'  resid.cov = covariance matrix of (all) residuals
##' (including non random variables) (for used.cases)
##'
##'  data.cov  = covariance matrix of (all)
##'  input data (for used.cases)
##'
##'  averages  = averages of input variables
##'
##'  Constant  = estimated model constant = (sum_i(ar[,,i])) x averages
##'
##'  call.ar.pattern = calling ar.pattern
##'
##'  call.ma.pattern = calling ma.pattern
##'
##'  out.ar.pattern = resulting ar.pattern
##'             (after possible model reduction)
##'
##'  out.ma.pattern = resulting ar.pattern
##'             (after possible model reduction)
##'
##'  max.iter = max no. of iterations in call
##'
##'  penalty = factor used in AIC model reduction
##'
##'  weight  = weighting of successive residuals
##'                     updating (default=0.33)
##'
##'  used.cases = cases in input which are analysed
##'
##'  trace   = trace(random part of resid.cov)
##'
##'  log.det = log(det(random part of resid.cov))
##'
##'  randoms = which are random variables in problem?
##'
##' @source The code is an R code which is based on the
##' article (below) by Spliid (1983). A repeated (socalled) pseudo
##' regression procedure is used in order to estimate the multivariate
##' arma model.
##'
##' @examples
##' # Example 1:
##' library(marima)
##' # Generate a 4-variate time series (in this example):
##' #
##' kvar<-4 ; set.seed(4711)
##' y4<-matrix(round(100*rnorm(4*1000,mean=2.0)),nrow=kvar)
##' # If wanted define differencing of variable 4 (lag=1)
##' # and variable 3 (lag=6), for example:
##' y4.dif<-define.dif(y4,difference=c(4,1,3,6))
##' # The differenced series will be in y4.dif$y.dif, the observations
##' # lost by differencing being excluded.
##' #
##' y4.dif.analysis<-y4.dif$y.dif
##' # Give lags the be included in ar- and ma-parts of model:
##' #
##' ar<-c(1,2,4)
##' ma<-c(1)
##' # Define the multivariate arma model using 'define.model' procedure.
##' # Output from 'define.model' will be the patterns of the ar- and ma-
##' # parts of the model specified.
##' #
##' Mod <- define.model(kvar=4,ar=ar,ma=ma,reg.var=3)
##' arp<-Mod$ar.pattern
##' map<-Mod$ma.pattern
##' # Print out model in 'short form':
##' #
##' short.form(arp)
##' short.form(map)
##' # Now call marima:
##' Model <- marima(y4.dif.analysis,ar.pattern=arp,ma.pattern=map,
##'                 penalty=0.0)
##' # The estimated model is in the object 'Model':
##' #
##' ar.model<-Model$ar.estimates
##' ma.model<-Model$ma.estimates
##' dif.poly<-y4.dif$dif.poly  # = difference polynomial in ar-form.
##' # Multiply the estimated ar-polynomial with difference polynomial
##' # to compute the aggregated ar-part of the arma model:
##' #
##' ar.aggregated <- pol.mul(ar.model,dif.poly,L=12)
##' # and print everything out in 'short form':
##' #
##' short.form(ar.aggregated,leading=FALSE)
##' short.form(ma.model,leading=FALSE)
##'
##' @references
##'
##' Jenkins,G.M. & Alavi,A.: Some aspects of modelling and forecasting
##'     multivariate time series, Journal of Time Series Analysis,
##'     Vol. 2, Issue 1, Jan. 1981, pp. 1-47.
##'
##' Madsen, H. (2008) Time Series Analysis, Chapmann \& Hall (in particular
##' chapter 9: Multivariate time series).
##'
##' Reinsel G.C. (2003) Elements of Multivariate Time Series Analysis,
##' Springer Verlag, 2$^{nd}$ ed. pp. 106-114.
##'
##' Spliid, H.: A Fast Estimation Method for the Vector
##'    Autoregressive Moving Average Model With Exogenous Variables, Journal
##'    of the American Statistical Association, Vol. 78, No. 384, Dec. 1983,
##'    pp. 843-849.
##'
##' www.itl.nist.gov/div898/handbook/pmc/section4/pmc45.htm
##'
##' @importFrom stats arima complete.cases cov lm step
##' 
##' @export

marima <- function(DATA = NULL, ar.pattern = NULL,
                   ma.pattern = NULL, means = 1,
                   max.iter = 50, penalty = 0, weight = 0.33,
                   Plot = 'none', Check = FALSE) {
    
    Test <- FALSE
    kvar <- min(dim(DATA))
    till = 6
    
    ## Check if ar- and/or ma-pattern specified when calling marima. If not
    ## set pattern to the correct unity-matrix (kvar by kvar): 
    if (is.null(ar.pattern) & is.null(ma.pattern)) {
        cat("Neither ar- nor ma-part of model specified. \n")
    }
    if (is.null(ar.pattern)) {
        ar.pattern = array(c(diag(kvar), 0 * diag(kvar)),
            dim = c(kvar, kvar, 2))
    }
    if (is.null(ma.pattern)) {
        ma.pattern = array(c(diag(kvar), 0 * diag(kvar)),
            dim = c(kvar, kvar, 2))
        max.iter <- 3
        till <- 1
        weight <- 0
    }
    ##
    
    AR <- ar.pattern
    MA <- ma.pattern
    max.iter <- round(max.iter)
    D <- dim(DATA)
    if (D[1] < D[2]) {
        DATA <- t(DATA)
    }
    
    case = which(complete.cases(DATA))
    All <- 1:max(dim(DATA))
    DATA <- DATA[case, ]
    Leftout <- setdiff(All, case)
    if (length(Leftout) <= 0) {
        cat("All cases in data, ", min(case), " to ", max(case),
            " accepted for completeness.\n")
    }
    
    if (length(Leftout) > 0) {
        cat("Cases ", Leftout, " left out for marima analysis. \n")
    }
    
    D <- dim(DATA)
    N <- D[1]
    kvar <- D[2]
    rownames(DATA) <- as.character(1:N)
    
    ## Determine which variables are random and which are not:
    randoms <- which(rowSums(matrix(abs(c(AR, MA)), nrow = kvar)) > 2)
    non.randoms <- which(rowSums(matrix(abs(c(AR, MA)), nrow = kvar)) == 2)
    ## 

    ## Check how averages for variables is to be handled:
    xmean <- means
    if (length(xmean) == 1) {
        if (xmean == 1) {
            means <- rep(1, kvar)
        }
    }
    if (length(xmean) == 1) {
        if (xmean != 1) {
            means <- rep(0, kvar)
        }
    }
    ##

    ## If specified (by setting Check=TRUE) if control output is wanted:  
    if (Check) {
        cat("Control printout (Check=TRUE (default)). \n")
        cat("----------------------\n")
        cat("Calling marima. Data dimensions (kvar,N) = ", kvar, N, "\n")
        cat("Coefficient polynomial dimensions = ", dim(AR), " and ",
            dim(MA), "\n")
        cat("including leading unity matrix. \n")
        
        if (is.null(ar.pattern)) {
            cat("no ar.pattern specified \n")
        }
        if (is.null(ma.pattern)) {
            cat("no ma.pattern specified \n")
        }
        
        if (!is.null(AR)) {
            cat("ar.pattern= \n")
            print(AR)
        }
        if (!is.null(MA)) {
            cat("ma.pattern= \n")
            print(MA)
        }
        cat("Start of data (5 first observations): \n")
        print(DATA[1:5, ])
        cat(" ... \n")
        cat("End of data (5 last observations): \n")
        print(DATA[(N - 4):N, ])
        cat( " \n")
        cat(" Calling parameters in use: \n")
        cat(" max.iter= ", max.iter, ", means =", means, ", penalty =",
            penalty, ", \n")
        cat(" weight =", weight, ", Plot = ", Plot,
            ", Check = ", Check, " \n")
        cat(" The above printout can be suppressed ",
         "by calling with Check=FALSE. \n ")
    }
    # End-of-control output

    ## Construct ar- and ma-polynomial arrays:
    ARmodel <- matrix(0, nrow = kvar, ncol = (kvar * dim(AR)[3]))
    MAmodel <- matrix(0, nrow = kvar, ncol = (kvar * dim(MA)[3]))
    for (i in 1:kvar) {
        ARmodel[i, which(AR[i, , ] > 0)] <- which(AR[i, , ] > 0)
        MAmodel[i, which(MA[i, , ] > 0)] <- which(MA[i, , ] > 0)
    }
    if(Test) {
             cat("ARmodel & MAmodel as constructed: \n")
             print(ARmodel)
             print(MAmodel)
             readline("0.4: Press <return to continue")
         }
    
    ## Construct proper names:
    arnam <- matrix(paste("ar", ARmodel, sep = ""), nrow = kvar)
    manam <- matrix(paste("ma", MAmodel, sep = ""), nrow = kvar)
    ##
    if(Test) {
             cat("arnam & manam as constructed: \n")
             cat("dim(arnam),dim(manam) \n")
             print(arnam)
             print(manam)
             readline("0.5: Press >return to continue ") 
         }

    On <- rep(0, kvar)
    for (i in 1:kvar) {
        if (sum(AR[i, , ]) > 1) {
            On[i] <- i
        }
        if (sum(AR[, i, ]) > 1) {
            On[i] <- i
        }
    }
    
    On <- On[On > 0]
    cases1 <- which(complete.cases(DATA[, On]))
    rownames(DATA) <- as.character(1:dim(DATA)[1])
    DATA <- DATA[cases1, ]
    D <- dim(DATA)
    X <- lagged.data(DATA = DATA, AR = AR, MA = MA,
                     init = TRUE, means = xmean)
    averages <- X$averages

    if(Test) { cat("About X \n")
               cat(names(X)," \n")
               cat("Averages =",X$averages,"\n")
               cat("X$ardata og X$madata = \n")
               na<-dim(X$ardata)[1]
               nm<-dim(X$madata)[1]
               print(X$ardata[c((1:8),((na-7):na)),])
               print(X$madata[c((1:8),((nm-7):nm)),])
               readline("1.5: Press <return to continue")
           }
    
    N <- dim(X$ardata)[1]
    
    SAR <- short.form(AR, "", tail = TRUE)
    ARlags <- as.numeric(dimnames(SAR)[[3]][])
    L <- dim(SAR)[3]
    SAM <- short.form(MA, "", tail = TRUE)
    MAlags <- as.numeric(dimnames(SAM)[[3]][])
    M <- dim(SAM)[3]

    if(Test) { 
        cat("SAR , ARlags , SAM , MAlags: \n")
        print(SAR)
        print(ARlags)
        print(SAM)
        print(MAlags)
        readline("2.0: Press >return to continue")
    }
        
    stt <- ARlags * kvar + 1
    stp <- ARlags * kvar + kvar
    for (i in 1:length(stt)) {
        if (i == 1) {
            arcols <- c(stt[i]:stp[i])
        }
        if (i > 1) {
            arcols <- c(arcols, c(stt[i]:stp[i]))
        }
    }
    
    stt <- MAlags * kvar + 1
    stp <- MAlags * kvar + kvar
    for (i in 1:length(stt)) {
        if (i == 1) {
            macols <- c(stt[i]:stp[i])
        }
        if (i > 1) {
            macols <- c(macols, c(stt[i]:stp[i]))
        }
    }

    if (Test) {
        cat("arcols = ",arcols,"\n")
        cat("macols = ",macols,"\n")
        readline("2.5: Press <-return to continue")
    }

    ARDATA <- fill.out(DATAarray = array(X$ardata,
                           dim = c(N, kvar, L)), SAR = SAR)
    MADATA <- fill.out(DATAarray = array(X$madata,
                           dim = c(N, kvar, M)), SAR = SAM)
    
    ARDATA <- matrix(ARDATA, nrow = N)
    MADATA <- matrix(MADATA, nrow = N)
    
    rownames(ARDATA) <- rownames(DATA)
    rownames(MADATA) <- rownames(DATA)
    
    colnames(ARDATA) <- paste("ar", arcols, sep = "")
    colnames(MADATA) <- paste("ma", macols, sep = "")

    if( Test ) {
             cat("SAR and SAM \n")
             cat(dim(SAR),"..",dim(SAM)," \n")
             print(SAR)
             print(SAM)
             readline( "3.14: Press <return to continue")           
             cat("ARDATA & MADATA \n")
             na <- dim(ARDATA)[1]
             nm <- dim(MADATA)[1]
             print(ARDATA[c((1:8),((na-7):na)),])
             cat("..... \n")
             print(MADATA[c((1:8),((nm-7):nm)),])
             readline( "3.15: Press <return to continue")
             }

    On2 <- matrix(rep(On, L), ncol = L)
    if (L > 1) {
        for (i in 2:L) {
            On2[, i] <- On2[, 1] + (i - 1) * kvar
        }
    }
    cases <- which(complete.cases(ARDATA[, c(On2)]))
    
    ARDATA <- ARDATA[cases, ]
    MADATA <- MADATA[cases, ]
    N <- length(cases)
    
    if( Test) {
        cat('colnames ARDATA', colnames(ARDATA), '\n')
        cat('colnames MADATA', colnames(MADATA), '\n')
        readline( "4.0: Press <return to continue")
    }
    
    DA <- get.models(ARDATA = ARDATA, AR = AR)
    DM <- get.models(ARDATA = MADATA, AR = MA)

    if( Test) {
        cat( "DA  and  DM \n")
        print(DA)
        print(DM)
        readline( "4.5: Press <return to continue")
    }
    
    randoms <- randoms
    trace <- 0
    log.det <- 0
    NAM = ""
    Iter1 <- round(max.iter/3)
    Iter3 <- Iter1 + till
    for (iterate in 1:max.iter) {
        for (i in 1:kvar) {
            # iterate=1 i=2
            colar <- colnames(ARDATA)
            colma <- colnames(MADATA)
            varx <- DA[i, ]
            varx <- varx[varx > 1]
            vare <- DM[i, ]
            vare <- vare[vare > 1]
            MO <- lm.form(i, colnames(ARDATA), zero = T)
            if ( Test & iterate==3 ) {
            cat(as.character(MO), "\n")
            readline("4.7: Press <return to continue")
            }
            MO <- lm.adds(MO, varx, colar)
            if ( Test & iterate==3 ) {
            cat(as.character(MO), "\n")
            readline("4.80: Press <return to continue")
            }
            MO <- lm.adds(MO, vare, colma)
            if ( Test & iterate==3 ) {
               cat(as.character(MO), "\n")
            readline("4.90: Press <return to continue")
            }
            # cat(i,M=,'\n')
            if ( iterate > Iter3 & penalty > 0 ) {
                MO <- NAM[i]
            }
                if( Test & iterate== Iter3 + 3 ){
                cat("NAM[i]= ",as.character(NAM[i]), "\n")
                cat("MO=     ",as.character(MO)    , "\n")
                readline("4.92: Press <return to continue")
            }
            
            MOD <- lm(MO, data = data.frame(ARDATA, MADATA))

            if ( Test & iterate==3 ) {
                    cat("ar-colnames=", colar, "\n")
                    cat("ma-colnames=", colma, "\n")
                    cat("iterate,i", iterate, i, "\n")
                    cat(as.character(MO), "\n")
                    readline("5: Press <return to continue")
                }
            summary(MOD)
            if (penalty > 1e-04) {
                if (iterate > Iter1) {
                  # cat('Stepwise Regression')
                  if (iterate >= Iter1 & iterate <= Iter3 & penalty > 0) {
                    MOD <- step(MOD, direction = "backward",
                                k = penalty, trace = FALSE)
                  }
                  if (iterate == Iter3) {
                    # cat('iterate and i: ',iterate,' , ',i,'\n')
                    NAM[i] <- as.character(MOD$call[2])
                    # cat('Model form i = ',NAM[i],' \n')
                  }
                }
            }
            
            if (iterate == max.iter) {
                if (i == 1) {
                  ar.estimates <- matrix(AR * 0, nrow = kvar)
                  ar.estimates[1:kvar, 1:kvar] <- diag(kvar)
                  ma.estimates <- matrix(MA * 0, nrow = kvar)
                  ma.estimates[1:kvar, 1:kvar] <- diag(kvar)
                  ar.fvalues <- ar.estimates
                  ma.fvalues <- ma.estimates
                  ar.pvalues <- ar.estimates
                  ma.pvalues <- ma.estimates
                }
                
                varia <- names(MOD$coefficients)
                if(Test){
                    cat('varia=names(MOD$coefficients=',varia,'\n')
                    readline("6: Press <return to continue")
                      }
                if (length(varia) > 0) {
                  types <- substr(varia, 1, 2)
                  numbers <- as.numeric(substr(varia, 3, 10))
                  ars <- which(types == "ar")
                  mas <- which(types == "ma")
                  ar.estimates[i, numbers[ars]] <- -MOD$coefficients[ars]
                  ma.estimates[i, numbers[mas]] <- +MOD$coefficients[mas]
                  f.values <- summary(MOD)[[4]][, 3]
                  f.values <- f.values^2
                  ar.fvalues[i, numbers[ars]] <- f.values[ars]
                  ma.fvalues[i, numbers[mas]] <- f.values[mas]
                  p.values <- summary(MOD)[[4]][, 4]
                  ar.pvalues[i, numbers[ars]] <- p.values[ars]
                  ma.pvalues[i, numbers[mas]] <- p.values[mas]
                }
                
                if (iterate == max.iter) {
                  # if(i==1){MOD1<-MOD} if(i==2){MOD2<-MOD}
                  # if(i==3){MOD3<-MOD} if(i==4){MOD4<-MOD}
                  if (i == kvar) {
                    ar.estimates <- array(ar.estimates, dim = dim(AR))
                    ma.estimates <- array(ma.estimates, dim = dim(MA))
                    ar.fvalues <- array(ar.fvalues, dim = dim(AR))
                    ma.fvalues <- array(ma.fvalues, dim = dim(MA))
                    ar.pvalues <- array(ar.pvalues, dim = dim(AR))
                    ma.pvalues <- array(ma.pvalues, dim = dim(MA))
                  }
                }
            }
            
            res <- MOD$residual
            # cat('i=',i,' ,length(MADATA[,i])=',length(MADATA[,i]),
            # 'length(residual)=',length(res),'\n')
            
            MADATA[, i] <- weight * MADATA[, i] + (1 - weight) * res
            # cat('Compute Trace', iterate,'\n')
            log.det[iterate] <- log(det(cov(MADATA[, randoms])))
            trace[iterate] <- sum(diag(cov(MADATA[, 1:kvar])))
            MADATA <- fill.out(DATAarray = array(MADATA,
                                   dim = c(N, kvar, M)), SAR = SAM)
            MADATA <- matrix(MADATA, nrow = N)
            colnames(MADATA) <- colma
        }
    }
    
    rownames(MADATA) <- rownames(ARDATA)
    # cat('Data analysed: .... \n') #
    # print(ARDATA[1:10,1:kvar])
    # print(MADATA[1:10,1:kvar])
    # cat('Dimensions are:',dim(ARDATA[,1:kvar]),' og ',
    # dim(MADATA[,1:kvar]),'\n')
    
    fitted <- ARDATA[, 1:kvar] - MADATA[, 1:kvar]
    
    # cat('Fitted = (1) \n') print(fitted[1:10,])
    # cat('averages= ',averages,'\n')
    # cat('means = ',means,'\n')
    
    for (i in 1:kvar) {
        fitted[, i] <- fitted[, i] + averages[i] * means[i]
    }
    colnames(fitted) <- paste("y", 1:kvar, sep = "")
    
    # cat('Fitted = (2) ... \n') print(fitted[1:10,])
    
    # cat('Iteration shift values',Iter1,' and ',Iter3,'.\n')
    
    fitted <- t(fitted)
    residuals <- t(MADATA[, 1:kvar])
    rownames(residuals) <- paste("u", 1:kvar, sep = "")
    
    covu <- cov(MADATA[, 1:kvar])
    covy <- cov(ARDATA[, 1:kvar])
    colnames(covu) <- rownames(residuals)
    rownames(covu) <- rownames(residuals)
    colnames(covy) <- paste("y", c(1:kvar), sep = "")
    rownames(covy) <- colnames(covy)
    
    mains = "Residual covariance matrix trace"
    if (Plot == "trace" & iterate >= 3) {
        plot(trace[1:max.iter], main = mains,
             xlab = "No. of iterations",
             ylab = "Computed trace", type = "l")
        grid(col = "blue")
    }
    
    # cat('Iteration shift at :',Iter1,' and ',Iter3, '\n')
    
    mains = "log(det(residual covariance matrix))"
    if (Plot == "log.det" & iterate >= 3) {
        plot(log.det[1:max.iter],
             main = mains,
             xlab = "No. of iterations",
             ylab = "log(det(residual covariance matrix))",
             type = "l")
        grid(col = "blue")
    }
    
    used.cases <- as.numeric(colnames(residuals))
    out.ar.pattern <- abs(ar.estimates)
    out.ar.pattern[which(out.ar.pattern > 0)] <- 1
    out.ma.pattern <- abs(ma.estimates)
    out.ma.pattern[which(out.ma.pattern > 0)] <- 1

    am<-averages*means
    Constant <- am
    lar <- length(c(ar.estimates))/(kvar*kvar)
    if (lar > 1) {
        for (i in 2:lar) {
        Constant <- Constant + ar.estimates[,,i]%*%am
#       cat("Constant =",Constant,"\n")
        }
    }   
    
    obj <- list(N = N,
                kvar = kvar,
                ar.estimates = ar.estimates,
                ma.estimates = ma.estimates,
                Constant     = Constant,
                ar.fvalues = ar.fvalues,
                ma.fvalues = ma.fvalues,
                ar.pvalues = ar.pvalues,
                ma.pvalues = ma.pvalues,
                residuals = residuals,
                fitted = fitted, 
                resid.cov = covu,
                data.cov = covy,
                averages = averages,
                mean.pattern = means,
                call.ar.pattern = AR,
                call.ma.pattern = MA,
                out.ar.pattern = out.ar.pattern,
                out.ma.pattern = out.ma.pattern,
                max.iter = max.iter,
                penalty = penalty, 
                weight = weight,
                used.cases = used.cases,
                trace = trace,
                log.det = log.det,
                randoms = randoms)
    class(obj) <- "marima"
    obj
}

########### END OF MARIMA ################

fill.out <- function(DATAarray = NULL, SAR = NULL) {
    Lags <- as.numeric(dimnames(SAR)[[3]])
    # cat("Input=\n")
    # print(DATAarray[1:10,,])
    # print(SAR)
    # cat(dim(SAR),dim(DATAarray),"\n")
    # cat('Lags=',Lags,'\n')
    N <- dim(DATAarray)[1]
    # cat('fill.out: dim=',dim(DATAarray),'\n')
    if (length(Lags) > 1) {
        for (i in 2:length(Lags)) {
            use <- (N + 1 - (1:(N - Lags[i])))
            DATAarray[use, , i] <- DATAarray[use - Lags[i], , 1]
        }
    }
    # cat("Output=\n")
    # cat(dim(SAR),dim(DATAarray),"\n")
    # print(DATAarray[1:10,,])
    return(DATAarray)
}

get.models <- function(ARDATA = NULL, AR = NULL) {
    SAR <- short.form(AR, "", tail = TRUE)
    # cat('ARDATA= \n') print(ARDATA)
    kvar <- dim(AR)[1]
    DA <- matrix(0, nrow = kvar, ncol = dim(ARDATA)[2])
    # print(DA) cat('kvar=',kvar,'\n') print(SAR)
    for (i in 1:kvar) {
        P <- SAR[i, , ]
        # cat('i,P',i,P,'\n')
        LP <- dim(P)[2]
        SP <- sum(P)
        if (SP > 1) {
            order = LP - 1
        }
        if (SP == 1) {
            order = 0
        }
        # cat('order=',order,'\n')
        P <- c(P)
        # print(P)
        DA[i, ] <- P
        # cat('length(P)=',length(P),'\n')
        for (j in (kvar + 1):length(P)) {
            DA[i, j] = j * P[j]
        }
    }
    return(DA)
}

lm.adds <- function(formula, varx, arnames) {
    formul <- formula
    if (length(varx) > 0) {
        for (i in 1:length(varx)) {
            formul <- paste(formul, "+", arnames[varx[i]], sep = "")
        }
    }
    if (length(varx) > 0) {
        formula <- formul
    }
    return(formula)
}

lm.form <- function(y, arnames, zero) {
    if (zero) {
        formula <- paste(arnames[y], " ~ -1 ", sep = "")
    }
    if (!zero) {
        formula <- paste(arnames[y], " ~    ", sep = "")
    }
    return(formula)
}

lagged.data <- function(DATA = NULL, AR = NULL, MA = NULL,
                        init = FALSE, means = 1)
{
    xmean <- means
    N <- dim(DATA)[1]
    kvar <- dim(DATA)[2]
    
    # cat('N,kvar=',N,kvar,'\n')
    
    if (length(xmean) == 1) {
        if (xmean == 1) {
            means <- rep(1, kvar)
        }
    }
    if (length(xmean) == 1) {
        if (xmean != 1) {
            means <- rep(0, kvar)
        }
    }
    
    averages <- rep(0, kvar)
    
    for (i in 1:kvar) {
        t <- which(!is.na(DATA[, i]))
        averages[i] <- mean(DATA[t, i])
        if (means[i] == 1) {
            DATA[, i] <- DATA[, i] - averages[i]
        }
        # cat('i',i,'mean',averages[i],'\n')
    }
    
    AR <- check.one(AR)
    MA <- check.one(MA)
    L <- dim(AR)
    M <- dim(MA)
    # cat('L=',L,'M=',M,'\n')
    ar.pos <- rep(0, L[3])
    j <- 0
    for (i in 1:L[3]) {
        coef <- sum(abs(AR[, , i]))
        # cat('i=',i,'ar.coef.sum=',coef,'\n')
        if (coef > 0) {
            j <- j + 1
            ar.pos[i] <- j
        }
    }
    
    ma.pos <- rep(0, M[3])
    j <- 0
    for (i in 1:M[3]) {
        coef <- sum(abs(MA[, , i]))
        # # cat('i=',i,'ma.coef.sum=',coef,'\n')
        if (coef > 0) {
            j <- j + 1
            ma.pos[i] <- j
        }
    }
    
    pr <- dim(short.form(AR, tail = TRUE))[3]
    qr <- dim(short.form(MA, tail = TRUE))[3]
    ardata <- array(rep(DATA, pr), dim = c(N, kvar, pr))
    madata <- array(rep(DATA, qr), dim = c(N, kvar, qr))
    rownames(ardata) <- rownames(DATA)
    rownames(madata) <- rownames(DATA)
    
    # cat('dim(ardata):',dim(ardata),'\n')
    if (pr > 1) 
        ardata[, , 2:pr] <- NA
    madata[, , 1:qr] <- 0
    
    ardata <- matrix(ardata, nrow = N)
    madata <- matrix(madata, nrow = N)
    rownames(ardata) <- rownames(DATA)
    rownames(madata) <- rownames(DATA)
    
    # cat('rownames(DATA)',rownames(DATA)[1:10])
    # cat('pr=',pr,'qr=',qr,'\n')
    # cat('ar-data & ma-data \n')
    # print(ardata) print(madata)
    
    # Start by initialising residuals for 'active' variables:
    # print(AR)
    
    if (init) {
        vars <- 1:kvar
        # cat('variables:',vars,'\n')
        for (i in 1:kvar) {
            if (sum(AR[i, , ]) == 1) {
                vars[i] = 0
            }
        }
        for (i in 1:kvar) {
            S <- sum(AR[i, , ])
            E <- rep(0, N)
            if (S > 1) {
                Y <- vars
                # cat('Y',Y,'\n')
                X <- Y[-i]
                X <- X[which(X != 0)]
                # cat('X',X,'\n') cat('i=',i,' Extra=',X,'\n')
                E <- arima(ardata[, i],
                      order = c(3, 0, 0),
                            xreg = ardata[, X])$residuals
                madata[, i] <- E
            }
        }
    }
    
    # cat('pr=',pr,'qr=',qr,'\n')
    # cat('ar-data \n')
    # print(ardata)
    # cat('ma-data \n')
    # print(madata)
    
    return(list(ardata = ardata,
                madata = madata,
                averages = averages))
} 
