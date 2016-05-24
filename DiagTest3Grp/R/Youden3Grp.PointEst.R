Youden3Grp.PointEst <-
function(x,y,z, method="Normal",randomStart.N=1,optim.method=NULL,t.minus.start=NULL,t.plus.start=NULL,lam.minus=1/3,lam0=1/3,lam.plus=1/3,typeIerror=0.05,margin=0.05,FisherZ=FALSE,...)
  {
    
    
    ### This function finds the optimal cutoff point for Youden index under Normal/TN/EMP/KS approach by theoretical results or optimization and calculates the resulting Youden Index            
    ###Inputs:
    ###x,y,z  diagnostic test marker measurements in the D- (healthy), D0,D+(disease) groups,
    ###NOTE:!!!! x, y, z must have increasing mean, i.e., mean(x)<mean(y)<mean(z)
    ###method: "Normal"/"TN"/"EMP"/"KS"/"KS-SJ"
      ###For EMP/KS, we have to do numerical optimization for estimates of Youden index and cutoff points. Thus it will depend on initial values when data have a wide range
      ####If using TN, need all values to be positive,if original input data are not all positive, a small positive constant will be added, which will not
      ####               affect the youden index estimation but will on optimal cutoff, we adjust this back at the end for cut off estimation
    
    ###optim.method: Optimization methods have these options,method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),by default l-BFGS-B"
    ###                and BFGS needs gradient function provided, In KS, the gradient is easy to be provided, thus may be better use BFGS by providing gr=KS.gradient()
    
    ####???Should we do random start (multiple starting t.minus and t.plus for EMP and KS method in case of bumpy goal function
    
    ####randomStart.N=1, only 1 starting value for t.minus and t.plus , set randomStart.N=10 if want to get the t.minus and t.plus at the maximum value
    ####start.t.minus,start.t.plus:starting values for t.minus and t.plus for use in the optimization steps for emp and KS methods,
    ###                                              it missing or NULL, will use normal solutions or multiple start (if randomStart.N>1) points
    ####...: can take the other arguments as in optim() function

    ####Note:In very skewed distributions, starting points will be super important, if use random starting points, then results may be far off
    ###      Thus, we will do normal estimates anyway first, use that as an estimation
    ## lam.minus, lam0,lam.plus: for sample size calculation, the expected proportion of samples in the D-, D0 and D+ group, which can be equal or not

     ######2. Ouput: a list consisting of 3 components:
     #### (1)dat: a list of x, y, z, recording the raw data for D-,D0, D+
     #### (2)dat.summary: a data frame of 3 rows and 3 columns # of obs, mu and sd for each group (at row)
     ####(3) est: a data frame of 1 row with the point estimates at column
    
    
    ##the following three output will have values if only TN method is used
    lambda.est <- NA
    t.minus.TN <- NA
    t.plus.TN <- NA
    
    EMP.YoudenFun <- function(beta,x0=x,y0=y,z0=z)
      {
        ###maximization goal function: Youden Index under empirical cdf
        ###beta: a vector of length 2, 1st is t.minus, and 2nd is t.plus
        ###return the youden index using empirical cdf
        t.minus <- beta[1]
        t.plus <- beta[2]

        0.5*(empirical.cdf(x0,t.minus)-empirical.cdf(y0,t.minus)+empirical.cdf(y0,t.plus)-empirical.cdf(z0,t.plus))
        
      }
    
    KS.YoudenFun <- function(beta,x0=x,y0=y,z0=z,KS.method="KS-SJ")
      {
        ###maximization goal function: Youden Index under Kernel smoothing
        ###beta: a vector of length 2, 1st is t.minus, and 2nd is t.plus

        ###if KS.method="KS-SJ", bw use Sheather and Jones plug-in estimator,else use quick and dirty "normal reference rule"
        
        t.minus <- beta[1]
        t.plus <- beta[2]
        
        ###use either SJ (Sheather and Jones) estimator or normal reference rule as bw
        
        #if(KS.method=="KS-SJ") bw.method <- "SJ" else bw.method <- "Reference"
        
        bw.minus <- BW.ref(x0,KS.method)
        bw0 <- BW.ref(y0,KS.method)
        bw.plus <- BW.ref(z0,KS.method)

        0.5*(KernelSmoothing.cdf(x0,t.minus,bw.minus)-KernelSmoothing.cdf(y0,t.minus,bw0)+KernelSmoothing.cdf(y0,t.plus,bw0)-KernelSmoothing.cdf(z0,t.plus,bw.plus))
           
      }
    
        
    ####Before starting Youden Index and cutoff: make sure that x,y,z are arranged in increasing mean order
    #temp.xyz <- list(x=x,y=y,z=z)
    #temp.xyz <- OrderXYZ(temp.xyz)
    #x <- temp.xyz$x
    #y <- temp.xyz$y
    #z <- temp.xyz$z


    if(!(method%in%c("Normal","TN","EMP","KS","KS-SJ"))) stop("Method options are: Normal/TN/EMP/KS/KS-SJ!!")
    ###In bootstrap CI computation,sometimes if a test marker has only few unique values in one group (e.g, Alzheimer marker zbentd has only 5 unique values in D- group with most subjects taking one of the values,some bootstrap sample may only sample the mode value and the TN method will have -Inf (as ss=0) log-likelihood in lambda calculation and lambda is estimated as the upper bound.

    x <- na.exclude(x)
    y <- na.exclude(y)
    z <- na.exclude(z)
    
    ####Sample size
    n.minus <- length(x)
    n0 <- length(y)
    n.plus <- length(z)

    raw.dat <- list(x=x,y=y,z=z)
    
    ###data summary
    dat.summary <- data.frame(n=c(n.minus,n0,n.plus),mu=c(mean(x),mean(y),mean(z)),sd=c(sd(x),sd(y),sd(z)),row.names=c("D-","D0","D+"))
    
    ####If using TN, need all values to be positive, prepare data to be positive first
    ###step 1: estimate lambda first
    ###step 2:transform the data by the lambda
    ###step 3:apply normal
    
    if(method=="TN")
      {
        #about 1.4 unit user time per bootstrap
        
        #####Changed 08/06/2009, before estimate lambda.est first, now tranform data to be positive first then estimate lambda
        
        ###constant to be added to each value before transformation
        min0 <- min(c(x,y,z))

        if(min0<=0) add.const <- abs(min0)+1 else add.const <- 0

        x <- x+add.const
        y <- y+add.const
        z <- z+add.const

        lambda.est <- Youden3Grp.BoxCox.lambda(xx=x,yy=y,zz=z,nboot=50,lambda0=0.00001)
        ####!!!DEBUG:03/25/2010
        #if(is.na(lambda.est)) browser()
        
        #cat("lambda.est=",lambda.est,"\n")

        x <- bcPower(x,lambda=lambda.est)
        y <- bcPower(y,lambda=lambda.est)
        z <- bcPower(z,lambda=lambda.est)

        
       ###Note: when we simulate data, mean(x)<mean(y)<mean(z) but after Box-cox, the relationship might change!!
      }            


    #######Estimate under normal assumption first anyway, for TN, after transforming data to be positive, the same as Normal 
    mu.minus <- mean(x,na.rm=TRUE)##these are the mean and SD after Box-Cox transformation
    mu0 <- mean(y,na.rm=TRUE)
    mu.plus <- mean(z,na.rm=TRUE)
        
    
    ###check whether mean are monotonically increasing
    #if(!(mu.minus<=mu0 && mu0<=mu.plus)) {print("group means are not monotonically increasing!!");browser()}
    
    s.minus <- sd(x,na.rm=TRUE)
    s0 <- sd(y,na.rm=TRUE)
    s.plus <- sd(z,na.rm=TRUE)
    

    
    var.minus <- s.minus^2
    var0 <- s0^2
    var.plus <- s.plus^2

    if(method=="TN")
      {
        transform.dat <- list(x=x,y=y,z=z,lambda=lambda.est)
         ###data summary
        transform.dat.summary <- data.frame(n=c(n.minus,n0,n.plus),mu=c(mu.minus,mu0,mu.plus),sd=c(s.minus,s0,s.plus),row.names=c("D-","D0","D+"))
      }
    else
      {
        transform.dat <- NA
        transform.dat.summary <- NA
      }

    if(s.minus==s0)
      {
        t.minus <- 0.5*(mu.minus+mu0)
      }else
    {
      t.minus <- (mu0*var.minus-mu.minus*var0)-s.minus*s0*sqrt((mu.minus-mu0)^2+(var.minus-var0)*log(var.minus/var0))
      t.minus <- t.minus/(var.minus-var0)
    }
        
    if(s0==s.plus)
      {
        t.plus <- 0.5*(mu0+mu.plus)
      }else
    {        
      
      t.plus <- (mu.plus*var0-mu0*var.plus)-s0*s.plus*sqrt((mu0-mu.plus)^2+(var0-var.plus)*log(var0/var.plus))
      t.plus <- t.plus/(var0-var.plus)
    }
    
    youden <- Youden3Grp.Normal(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus)
    Se <- youden$Se
    Sp <- youden$Sp
    Sm <- youden$Sm
    #youden <- youden$youden                              

    
   if (method%in%c("EMP","KS","KS-SJ"))
     {

       #if(is.null(optim.method) |missing(optim.method)) optim.method <- "Nelder-Mead"##heuristic 
       #if(is.null(optim.method) |missing(optim.method)) optim.method <- "BFGS"##by approximation of hessian matrix
       if(is.null(optim.method) |missing(optim.method)) optim.method <- "L-BFGS-B"##by approximation of hessian matrix
       
       #browser()
       min0 <- min(c(x,y,z),na.rm=TRUE)
       max0 <- max(c(x,y,z),na.rm=TRUE)

       ####evaluate with Normal solution as starting point
       if(method=="EMP") value0 <- EMP.YoudenFun(beta=c(t.minus,t.plus),x0=x,y0=y,z0=z)##for L-BFGS-B, need lower and upper
       if (method=="KS" |method=="KS-SJ") value0 <- KS.YoudenFun(beta=c(t.minus,t.plus),x0=x,y0=y,z0=z,KS.method=method)

       start.t.minus <- t.minus
       start.t.plus <- t.plus
       value <- value0

       ####evaluate with specified starting point
       if(!missing(t.minus.start) & ! missing(t.plus.start)&!is.null(t.minus.start)&!is.null(t.plus.start))
         {
             ####use specified starting values eg. from TN solutions
             #start.t.minus <- t.minus.start
             #start.t.plus <- t.plus.start
           if(method=="EMP") value1 <- EMP.YoudenFun(beta=c(t.minus.start,t.plus.start),x0=x,y0=y,z0=z)##for L-BFGS-B, need lower and upper
           if (method=="KS" |method=="KS-SJ") value1 <- KS.YoudenFun(beta=c(t.minus.start,t.plus.start),x0=x,y0=y,z0=z,KS.method=method)
           if(value1-value0>=0.001)
             {
               start.t.minus <- t.minus.start
               start.t.plus <- t.plus.start
               value <- value1
               t.minus <- start.t.minus
               t.plus <- start.t.plus
             }
         }

       
       for(random.try in 1:randomStart.N)
          {     
            t.minus.try <- runif(1,min0,min0+(max0-min0)/3)
            t.plus.try <- runif(1,t.minus.try,t.minus.try+(max0-min0)/2)
                        
            if(method=="EMP")
              {                
                opt.res <- optim(par=c(t.minus.try,t.plus.try),fn=EMP.YoudenFun,x0=x,y0=y,z0=z,control=list(trace=F,fnscale=-1),method=optim.method,lower=min0,upper=max0,...)##fnscale=-1 means get maximum
              }
            if(method=="KS" |method=="KS-SJ")
              {                
                opt.res <- optim(par=c(start.t.minus,start.t.plus),fn=KS.YoudenFun,x0=x,y0=y,z0=z,KS.method=method,control=list(trace=F,fnscale=-1),method=optim.method,lower=min0,upper=max0,...)##!!need to work on this
              }
            
            value.try <- opt.res$value
            if(value.try-value>=0.001)
              {
                ##use the current output as solution
                start.t.minus <- opt.res$par[1]
                start.t.plus <- opt.res$par[2]

                ##update the starting point for next iteration
                #start.t.minus <- t.minus
                #start.t.plus <- t.plus
                value <- value.try
                t.minus <- start.t.minus
                t.plus <- start.t.plus

              }
            
          }
        
       if(method=="EMP")
         {
           Se <-1-empirical.cdf(z,t.plus)
           Sp <- empirical.cdf(x,t.minus)
           Sm <- empirical.cdf(y,t.plus)-empirical.cdf(y,t.minus)
         }
       
       if (method=="KS" | method=="KS-SJ")
         {
           bw.minus <- BW.ref(x,method)
           bw0 <- BW.ref(y,method)
           bw.plus <- BW.ref(z,method)
            
           Se <-1-KernelSmoothing.cdf(z,t.plus,bw.plus)
           Sp <- KernelSmoothing.cdf(x,t.minus,bw.minus)
           Sm <- KernelSmoothing.cdf(y,t.plus,bw0)-KernelSmoothing.cdf(y,t.minus,bw0)

         }
     }
    
    
    youden <- 0.5*(Se+Sp+Sm-1)
    
    #######If using "TN" method, data go through box-cox transformation first, t.minus and t.plus directly estimated are in transformed scale,
    
    if(method=="TN")
      {
        t.minus.TN <- t.minus ###record the optimal cutoff at the scale after box-cox transformation in TN
        t.plus.TN <- t.plus
        
        #browser()
        ###Now obtain the optimal cutoff point back to the original scale before box-cox transformation 
        if(lambda.est==0)
          {
            t.minus <- exp(t.minus.TN)
            t.plus <- exp(t.plus.TN)
          }
        else
          {
                
            t.minus <- (t.minus.TN*lambda.est+1)^(1/lambda.est)###if lambda.est>0, 
            t.plus <- (t.plus.TN*lambda.est+1)^(1/lambda.est)##will have problem if lambda.est=-0.9 and the first bracket produced a negative number!!
          }

        ###get back to original scale before adding the constant to make x, y, z all >0
        t.minus <- t.minus-add.const
        t.plus <- t.plus-add.const
      }
        

    ####add 07/22/2009, calculate Fisher's Z transformation on youden and the variance of z
    youden.z <- FisherZ(youden)
    ######

    #####Sample Size Calculation for all methods (Normal/TN/EMP/KS/KS-SJ) though the sample size equation is under normality assumption
    if(method=="TN") J.partial <- PartialDeriv.Youden(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus.TN,t.plus.TN) else J.partial <- PartialDeriv.Youden(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus)
    
    Mj <- s.minus^2*lam0/lam.minus*(J.partial$Y.mu.minus^2+0.5*J.partial$Y.s.minus^2)+s0^2*(J.partial$Y.mu0^2+0.5*J.partial$Y.s0^2)+s.plus^2*lam0/lam.plus*(J.partial$Y.mu.plus^2+0.5*J.partial$Y.s.plus^2)
    z0 <- qnorm(typeIerror/2,lower.tail=F)
    
    sampleSize <- z0^2*Mj/margin^2
    
    if(FisherZ)
      {
        sampleSize <- sampleSize/((1-youden^2)^2)
      }
      
    sampleSize <- ceiling(sampleSize)
    
    ###########

    
    out0 <- data.frame(mu.minus=mu.minus,mu0=mu0,mu.plus=mu.plus,s.minus=s.minus,s0=s0,s.plus=s.plus,t.minus=t.minus,t.plus=t.plus,Se=Se,Sp=Sp,Sm=Sm,youden=youden,youden.z=youden.z,lambda.est=lambda.est,t.minus.TN=t.minus.TN,t.plus.TN=t.plus.TN,sampleSize=sampleSize)##outdated output before the use of DiagTest3Grp class
    ###Note:In out0, for TN, the mu and sd are for the data after Box-Cox transformation while t.minus.TN,t.plus.TN are optimal cut-points on the transformed scale with t.minus, t.plus being the optimal cut-points that are already transformed back to original data scale; for others, they are the same as in dat.summary

    return(list(dat=raw.dat,dat.summary=dat.summary,est=out0))

  }

