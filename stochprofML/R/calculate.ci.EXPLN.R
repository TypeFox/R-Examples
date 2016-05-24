calculate.ci.EXPLN <-
function(alpha,parameter,prev.result,dataset,n,TY,fix.mu=F,fixed.mu) {
# Calculates approximate marginal maximum likelihood confidence intervals with
# significance level alpha for all parameters in the model. The intervals are 
# approximate because the function uses the formula
#    [ theta_i +/- q_{1-alpha/2} * sqrt(H_ii) ],
# where
# * theta_i is the i.th parameter,
# * q_{1-alpha/2} is the 1-alpha/2 quantile of the standard normal distribution,
# * H is the inverse Hessian of the negative log likelihood function evaluated 
#   at the maximum likelihood estimate; H_ii is the i.th diagonal element of H.
# This approximation implicitly assumes that the log likelihood function is unimodal.
# The confidence interval is first calculated on the transfomed, unrestricted parameter 
# space and then back-transformed to the original one.
# 
# Parameters:
# - alpha is the significance level
# - parameter is the maximum likelihood estimate around which the confidence interval 
#   is centered; if this value is missing, it is determined from prev.result. This 
#   parameter has to be on the original scale, not on the logit-/log-transformed scale
#   as used during the estimation procedure. It has to be (TY*(m+2)-1)-dimensional, even for
#   fix.mu==T.
# - prev.result is a list of previous results as returned by stochprof.results. This
#   variable is only used when "parameter" is missing.
# - dataset is a matrix which contains the cumulated expression data over all cells in a tissue sample.
#   Columns represent different genes, rows represent different tissue samples.
# - n is the number of cells taken from each tissue sample.
# - TY is the number of types of cells that is assumed in the stochastic model.
# - fix.mu (default=F) indicates whether the log-means are kept fixed in the estimation procedure or whether
#   they are to be estimated.
# - fixed.mu (no default, needs to be specified only when fix.mu==T) is a vector containing the values to which
#   the log-means should be fixed if fix.mu==T. The order of components is as follows:
#   (mu_type_1_gene_1, mu_type_1_gene_2, ..., mu_type_2_gene_1, mu_type__gene_2, ...)

   
   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)
   d.sum.of.mixtures <- NULL
   backtransform.par <- NULL
   stochprof.results <- NULL
   transform.par <- NULL   
   rm(d.sum.of.mixtures)
   rm(backtransform.par)
   rm(stochprof.results)
   rm(transform.par)   
   
   # number of genes   
   m <- ncol(dataset)

   ############################
   ## define target function ##
   ############################
   
   # log likelihood for one gene
   loglikeli <- function(y,p,mu,sigma,lambda) {
      # p is of length TY
      # mu and sigma are of length TY-1 (so possibly NULL)
      # lambda is of length 1
      max(-10^7,sum(d.sum.of.mixtures(y,n,p,mu,sigma,lambda,logdens=T)))
   }
      
   # this function will be minimized (the function "to.minimize" below just 
   # changes the parameterisation)
   target.function <- function(p,mu,sigma,lambda) {  
      # p is of length TY
      # mu is of length (TY-1)*m
      # sigma is of length 1 (if TY>1) or NULL
      # lambda is of length m
      
      # build mu.matrix such that the g.th row contains the values for gene g
      if (!is.null(mu)) {
         mu.matrix <- matrix(mu,byrow=F,ncol=TY-1)
      }
      
      # consider negative log likelihood because the target function will be minimized
      this.sum <- 0
      for (g in 1:m) {
         if (TY>1) {
            this.sum <- this.sum - loglikeli(dataset[,g],p,mu.matrix[g,,drop=T],rep(sigma,TY-1),lambda[m])
         }
         else {
            this.sum <- this.sum - loglikeli(dataset[,g],p,NULL,NULL,lambda[m])
         }
      }
      return(this.sum)
   }
      
   to.minimize <- function(theta) {
      # theta is the transformed parameter.
      # No penalization included here.

      # When a parameter is already at the boundary, e.g. log(sigma)=-Inf,
      # then the hessian function will pass an NaN to this function. In that case,
      # an error would occur. Hence, replace the NaNs by 
      na.positions <- which(is.na(theta))
      if (sum(na.positions)>0) {
         theta[na.positions] <- parameter.trans[na.positions]
      }

      back.theta <- backtransform.par(this.par=theta,m=m,fix.mu=fix.mu,fixed.mu=fixed.mu)
      
      if (TY>1) {
         p <- back.theta[1:(TY-1)]
         p <- c(p,1-sum(p))

         mu <- back.theta[TY:((m+1)*(TY-1))]
         sigma <- back.theta[(m+1)*(TY-1)+1]

         lambda <- back.theta[(m+1)*(TY-1)+1+(1:m)]         
      }
      else {
         p <- 1
         mu <- NULL
         sigma <- NULL
         lambda <- back.theta
      }            
          
      a <- target.function(p,mu,sigma,lambda)
      a <- min(10^7,a)
      
      return(a)
   }
     
   # if "parameter" is missing, it has to be determined as the optimum from the previous results
   if (missing(parameter)) {
      res <- stochprof.results(prev.result=prev.result,TY=TY,show.plots=F)
      parameter <- res[1,-ncol(res)]
   }
   
   # when fix.mu==T: check whether fixed.mu (to which mu should be fixed) agrees with the mu part
   # of "parameter". If not, something went wrong. If fixed.mu is missing, determine it from "parameter"
   if ((fix.mu) && (TY>1)) {
      mu.in.parameter <- parameter[TY:((m+1)*(TY-1))]
      if ((missing(fixed.mu)) || (is.null(fixed.mu))) {
         fixed.mu <- mu.in.parameter
      }
      else {
         differ <- fixed.mu - mu.in.parameter
         if (sum(round(differ,4))>0) {
            stop("calculate.ci: wrong mu fixed.")
         }         
      }
   }  
   
   # transform parameter to unrestricted parameter space;
   # if fix.mu==T, parameter.trans becomes lower-dimensional than parameter
   parameter.trans <- transform.par(parameter,m,fix.mu)
      
   # Calculate Hessian matrix of to.minimize, i.e. of *negative* log-likelihood.
   # This is the Hessian matrix in the unrestricted parameterization, not in the original one.
   hesse <- hessian(func=to.minimize,x=parameter.trans,method="Richardson")  
   hesse[is.na(hesse)] <- 0  

   # inverse Hessian
   Sigma <- ginv(hesse)

   # 1-alpha/2 quantile of standard normal distribution
   coeff <- qnorm(1-alpha/2)

   # check whether there are negative diagonal elements in the Hessian
   if (sum(diag(Sigma)<0)>0) {
      print("Hessian not positive semi-definite; return NULL.")
      return(NULL)   
   }
   
   # confidence interval for transformed parameter     
   CI.trans.lower <- parameter.trans - coeff * sqrt(diag(Sigma))
   CI.trans.upper <- parameter.trans + coeff * sqrt(diag(Sigma))   
   CI.trans <- cbind(CI.trans.lower,CI.trans.upper)

   # confidence interval for original parameter
   if (TY>1) {
      nrows <- TY*(m+1) 
   }
   else {
      nrows <- m
   }
   CI <- matrix(NA,ncol=2,nrow=nrows)   
     
   CI[,1] <- backtransform.par(this.par=CI.trans[,1],m=m,fix.mu=fix.mu,fixed.mu=fixed.mu)
   CI[,2] <- backtransform.par(this.par=CI.trans[,2],m=m,fix.mu=fix.mu,fixed.mu=fixed.mu)   
   colnames(CI) <- c("lower","upper")
   
   # bounds for p might not be in correct order (i.e. lower <= upper) due to transformation;
   # in that case, swap
   if (TY>1) {
      for (i in 1:(TY-1)) {
         if (CI[i,1]>CI[i,2]) {
            a <- CI[i,1]
            CI[i,1] <- CI[i,2]
            CI[i,2] <- a
         }
      }
   }
          
   return(CI)
}
