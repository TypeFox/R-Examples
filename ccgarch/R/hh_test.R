#************************************************************************************************
hh.test <- function(dvar){
   dvar <- as.matrix(dvar)              # transforming a variable into a matrix object. effective, eg, for zoo object.
   nobs <- dim(dvar)[1]                 # number of observations
   ndim <- dim(dvar)[2]                 # number of dimensions
   ind <- 1:ndim
   sqy <- dvar^2
   sqy.mean <- colMeans(sqy)

     if(ndim==1){
        stop("the argument must be a matrix")
     }

    uni.garch <- function(dvar){           # a function for fitting a GARCH(1,1) model to each series
       nobs <- dim(dvar)[1]; ndim <- dim(dvar)[2]; 

        # estimating a univariate GARCH(1,1)
        uni.garch.est <- function(dvar, ini.par=c(0.01, 0.1, 0.85)){
            # a loglikelihood function of data with GARCH(1,1) cond. varaince.
                garch.ll <- function(para, dvar){
                    para <- para^2 # positivity condition
                    h <- uni.vola(para, dvar)
                    lf <- dnorm(dvar, mean=0, sd=sqrt(h), log=T)
                    sum(-lf)
                }
            est <- optim(sqrt(ini.par), garch.ll, dvar=dvar, method="BFGS",control=list(reltol=1e-15,maxit=10000))
            cond.var <- uni.vola(est$par^2, dvar)
            std.resid <- dvar/sqrt(cond.var)
            list(par.est=est$par^2, cond.var=cond.var, std.resid=std.resid)
        }

	   std.resid <- matrix(0, nobs, ndim)
	   cond.var <- std.resid
	   param <- matrix(0, 3, ndim)
	   for(i in 1:ndim){
	      tmp.garch <- uni.garch.est(dvar[,i])
	      std.resid[,i] <- tmp.garch$std.resid      # Standardised Residuals
	      cond.var[,i] <- tmp.garch$cond.var        # Conditional Variance
	      param[,i] <- tmp.garch$par.est            # Estimates of GARCH(1,1) Parameters
	   }
	   list(std.resid=std.resid, cond.var=cond.var, param=param)
	}

	tmp.garch <- uni.garch(dvar)	# fitting GARCH(1,1) to each series
   	h <- tmp.garch$cond.var         # estimated conditional variances
	para <- tmp.garch$param         # parameter estimates

    zeta <- sqy/h-1                 # denoted by zeta^2 - 1 in HH
    LM.HH <- numeric(ndim)
      for(j in 1:ndim){
            v.t <- cbind(1, sqy[,j], h[,j])
            v.t <- rbind(c(1, sqy.mean[j], sqy.mean[j]), v.t)
            v.t <- v.t[1:nobs, ] 
            dhdw <- matrix(0, nobs, 3)                              # partial derivatives of GARCH equation with respect to parameters
            dhdw[1,] <- v.t[1,]
            for(i in 2:nobs)  dhdw[i,] <- v.t[i,] + para[3,j]*dhdw[(i-1),]
            x.it <- dhdw/h[,j]                                      # denoted by x_{it} in HH 
            z.k <- cbind(sqy[,ind[ind!=j]], h[,ind[ind!=j]])
            z.k <- rbind(c(sqy.mean[ind[ind!=j]], sqy.mean[ind[ind!=j]]), z.k);  z.k <- z.k[1:nobs,] # denoted by z_{j} in HH
            reg <- summary(lm(zeta[,j]~x.it+z.k-1))                 # T*R^2 in this regression gives HH_{jk} test statistic 
         LM.HH[j] <- nobs*reg$r.squared                             # this follows the chi-squared dist with df<-2.
      }
   LM.HH <- sum(LM.HH)
   p.HH <- pchisq(LM.HH, 2*ndim*(ndim-1),lower.tail=F)              # p-values for HH test
   out <- c(LM.HH, p.HH)
   names(out) <- c("Test Stat", "p-Value")
   out
   }
#********************************************************************************************




