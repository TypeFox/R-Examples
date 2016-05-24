#********************************************************************************************
nt.test <- function(dvar){
   dvar <- as.matrix(dvar)              # transforming a variable into a matrix object. effective, eg, for zoo object.
   nobs <- dim(dvar)[1]
   ndim <- dim(dvar)[2]
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

   v <- cbind(1, sqy, h)
   v <- rbind(c(1, sqy.mean, sqy.mean), v)
   v <- v[1:nobs,]
   LM <- numeric(ndim)
   LMrob <- numeric(ndim)
   one <- rep(1, nobs)
   zeta <- sqy/h-1 #std.resid^2-1
   for(j in 1:ndim){
      # standard version of the test statistic
      dhdw <- matrix(0, nobs, 2*ndim+1)
      dhdw[1,] <- v[1,]
      for(i in 2:nobs){
         dhdw[i,] <- v[i,] + para[3,j]*dhdw[(i-1),]
      }
      dhdw <- dhdw/h[,j]
      reg <- summary(lm(zeta[,j]~dhdw-1))
      LM[j] <- nobs*reg$r.squared
      # Robust version of the test statistic
      ind <- numeric(ndim); ind[j] <- 1
      ind <- c(1, ind, ind)
      x1 <- dhdw[,ind==1]
      x2 <- dhdw[,ind!=1]
      reg1 <- lm(x2~x1-1)
      nd <- zeta[,j]*reg1$residuals
      reg2 <- lm(one~nd-1)
      LMrob[j] <- nobs - sum(reg2$residuals^2)
   }
   out <- c(sum(LM), sum(LMrob))
   out <- cbind(out, pchisq(out, 2*ndim*(ndim-1), lower.tail=F))
   colnames(out) <- c("Test stat", "p-value")
   rownames(out) <- c("NT","robNT")
   out
}
#********************************************************************************************
