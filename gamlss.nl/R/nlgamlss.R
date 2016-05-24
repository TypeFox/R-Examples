nlgamlss <- function (
                y = NULL, 
       mu.formula = ~1, 
    sigma.formula = ~1, 
       nu.formula = ~1, 
      tau.formula = ~1,
           mu.fix = FALSE, 
        sigma.fix = FALSE, 
           nu.fix = FALSE, 
          tau.fix = FALSE, 
          all.fix = FALSE, 
         mu.start = NULL, 
      sigma.start = NULL, 
         nu.start = NULL, 
        tau.start = NULL, 
           family = NO(),
          weights = 1, 
            exact = FALSE, 
            delta = 1, 
             data = parent.frame(),
  #       typsize = abs(p0),
  #        stepmax = sqrt(p0 %*% p0),
          control = NL.control(),         
      llik.output = FALSE )
{
.gamlss.multin.list<-c("MULTIN", "MN3", "MN4", "MN5")
#----------------------------------------------------------------------------------------
rqres <- function (pfun = "pNO", 
                   type = c("Continuous", "Discrete", "Mixed"),
               censored = NULL,  
                   ymin = NULL, 
                 mass.p = NULL, 
                prob.mp = NULL,
                      y = y,
                         ... )
{ }
body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
# this function is needed to get the infromation for the parameter formula
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
getParameters <- function(what = "mu", 
                       formula = mu.formula, 
                           fix = mu.fix, 
                         start = mu.start,
                            np = npmu,
                     start.at  = "1",
                     start.v   = 1)
 {
  if (inherits(formula, "formula")) # if formula
   {
        fo2 <- finterp(formula, .envir = envir, .start = start.v, .name = envname)
       npt1 <- length(attr(fo2, "parameters"))
      if (is.character(attr(fo2, "model"))) 
        {
            if (length(attr(fo2, "model")) == 1) 
            {
             fo1 <- function(p) eval(parse(text= paste(paste(paste(paste(what,
                 "h",sep="."),"(p[",sep=""), start.at, sep=""),"]*rep(1,N))",sep="")))                                                                          
             attributes(fo1) <- attributes(fo2)
             fo2 <- NULL
            }
        }
      else 
        {
            if (np  != npt1) 
            {
                cat("\nParameters are ")
                cat(attr(fo2, "parameters"), "\n")
                stop(paste(start, "should have", npt1, "estimates"))
            }
            if (is.list(start)) 
            {
                if (!is.null(names(start))) 
                {
                  o <- match(attr(fo2, "parameters"), names(start))
              start <- unlist(start)[o]
                  if (sum(!is.na(o)) != length(start)) 
                    stop("invalid estimates for", what, " - probably wrong names")
                }
                else start <- unlist(start)
            }
        }
      if (!is.null(fo2)) 
        {
               fo1 <- function(p) #mu.h(fo2(p))
                   { eval(parse(text=paste(paste(what,"h",sep="."),"(fo2(p))",sep="")))} 
   attributes(fo1) <- attributes(fo2)
        }
   } # end if formula
 else if (is.function(formula))  
   {# if function 
        fo1 <- formula
   }
   # end if
  if (!is.null(fo1) && is.null(attr(fo1, "parameters"))) 
    {
 attributes(fo1) <- if (is.function(formula)) 
                     {
                       if (!inherits(formula, "formulafn")) 
                        {
                        attributes(fnenvir(formula, .envir = envir))
                        }
                       else attributes(formula)
                     }
                    else 
                     {
            attributes(fnenvir(fo1, .envir = envir))
                     }
    }
    
    nlp <- if (is.function(fo1)) length(attr(fo1, "parameters"))
           else if (is.null(fo1))    NULL
           else npt1 #
    # end if else
   if (!is.null(nlp) && nlp != np) 
        stop(paste(what , "should have", nlp, "initial estimates"))

   if (inherits(formula, "formula") || is.function(formula)) 
    {
        if (!fix) 
        {
            if (is.numeric(start) && length(start) !=np) 
                stop(paste(what, "start must be of size ",np))
            if (!is.numeric(start)) 
                start <- rep(0,np)
                 fnfo <- fo1
        }
        else 
        {
            if (!is.numeric(start)) 
                stop("Missing initial conditions for mu")
            else if (length(start) !=np) 
                stop(paste(what, "start must be of size ",np))
            else fnfo <- function(p) fo1(start)
        }
    }
   else if (!fix) 
    {
        fnfo <- function(p) {function(p) eval(parse(text= paste(paste(paste(paste(what,"h",sep="."),"(p[",sep=""), start.at, sep=""),"]*rep(1,N))",sep="")))  }
        #mu.h(rep(p[1], N))
        if (!is.numeric(start)) 
        {
         stop(paste(what, "start must be of numeric"))
        }
    }
   else 
    {
        if (length(start) == N) 
            fnfo <- function(p) start
        else fnfo <- function(p) rep(start[1], N)
       np <- 1
    }
 list(func = fnfo, np = np , start = start)
 }
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
#  getParameters function ends here
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
## there are here to be consistent with gamlss but they have to change 
## gamlss.rc.list<-c("EX.rc","Exponential.rc") # the right censoring distribution list 
#gamlss.bi.list<-c("BI", "Binomial", "BB", "Beta Binomial") # binomial denominators
#----------------------------------------------------------------------------------------
      nlcall <- sys.call()
##-----------------------------------------------------
## this again meyby has to change if nlgamlss() is callled from gamlss()
      envir <- parent.frame()
    respenv <- FALSE  
    envname <-  deparse(substitute(envir))
##-----------------------------------------------------
      if(!missing(data)) 
        if(any(is.na(data)))stop("The data contains NA's, use data = na.omit(mydata)") 
# get the data (I am not sure if this is the best way)
    if (is.data.frame(data)) { attach(data); on.exit(detach(data))}
##-----------------------------------------------------
##     get the family
##-----------------------------------------------------
      family <- as.gamlss.family(family)         
       fname <- family$family[1]
        dfun <- paste("d",fname,sep="")
        pfun <- paste("p",fname,sep="")
        lpar <- length(family$parameters)   
        mu.h <- family$mu.linkinv
     sigma.h <- family$sigma.linkinv 
        nu.h <- family$nu.linkinv
       tau.h <- family$tau.linkinv
#------------------------------------------------------
    if (all.fix) 
     mu.fix <- sigma.fix <- nu.fix <- tau.fix <- FALSE
       npmu <- length(mu.start)
    npsigma <- length(sigma.start)
       npnu <- length(nu.start)
      nptau <- length(tau.start)
        mu1 <- sigma1 <- nu1 <- tau1 <- NULL
       type <- "unknown"
#------ y variable
##-----------------------------------------------------------------------------------------
## This part deals with the response variable 
     if (any(is.na(y)))     stop("NAs in y - use na.omit()")
         Y <- y       
         if(is.null(dim(Y)))                       # if y not matrix
         N <- length(Y) else N <- dim(Y)[1]   # calculate the dimension for y  
## extracting now the y and the binomial denominator in case we use BI or BB
    if(any(family$family%in%gamlss:::.gamlss.bi.list)) 
    { 
       if (NCOL(Y) == 1) 
            {
            y <- if (is.factor(Y))  Y != levels(Y)[1] else Y
            bd <- rep(1, N)
            if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
            } 
       else if (NCOL(Y) == 2) 
            {
            if (any(abs(Y - round(Y)) > 0.001)) {
            warning("non-integer counts in a binomial GAMLSS!")
                                                }
            bd <- Y[,1] + Y[,2]
            y <-  Y[,1]
            } 
       else stop(paste("For the binomial family, Y must be", 
            "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes", 
            "and col 2 is no. failures"))
     }
     # multinomial checking
     else if(any(family$family%in%.gamlss.multin.list))
            {
               y <- if(is.factor(Y))   unclass(Y)
                    else Y
            } 
     else if(is.Surv(Y))
          { 
           ## checking that the family is censored
           if (length(grep("censored",family$family[[2]]))==0) 
            stop(paste("the family in not a censored distribution, use cens()"))
           ## checking compatability of Surv object and censored distribution
           if (length(grep(attr(Y,"type"),family$family[[2]]))==0) 
            stop(paste("the Surv object and the censored distribution are not of the same type"))
           y <- Y    
           }     
     else {y <- Y }
     if (!family$y.valid(y))  stop( "response variable out of range")
##---------------------------------------------------------------------------------------
##checking the permissible y values      
   if (!family$y.valid(y)) # MS Thursday, June 20, 2002 at 16:30 
       stop( "response variable out of range")
        censor <- FALSE
    if (length(weights) == 1)    weights <- rep(weights, N)
##---------------------mu----------------------------------------------------------------
if ("mu"%in%names(family$parameters))
{ 
  all <- getParameters(what = "mu", formula = mu.formula, fix = mu.fix, 
                       start = mu.start,  np = npmu, start.at = "1", start.v=1)
  fnmu <- all$func
  npmu <- all$np
  mu.start <- all$start

    npl2 <- if (!sigma.fix)   npmu *(!mu.fix) + 1
            else 1
            
}
#-----------------------sigma------------------------------------------------------------
if ("sigma"%in%names(family$parameters)) 
 {
     all <- getParameters(what = "sigma", formula = sigma.formula, fix = sigma.fix, 
                         start = sigma.start,  np = npsigma, 
                         start.at = "npl2", start.v = npl2)
     fnsigma <- all$func
     npsigma <- all$np
 start.sigma <- all$start
        npl3 <- if (!nu.fix)  npmu *(!mu.fix) +npsigma * (!sigma.fix) + 1
               else 1
 }

#----------------------------nu----------------------------------------------------------
 if ("nu"%in%names(family$parameters))
 { 
       all <- getParameters(what = "nu", formula = nu.formula, fix = nu.fix, 
                         start = nu.start,  np = npnu, 
                         start.at = "npl3", start.v = npl3)
     fnnu <- all$func
     npnu <- all$np
 start.nu <- all$start
    npl4 <- if (!tau.fix) 
       npmu *(!mu.fix) +npsigma * (!sigma.fix) +npnu * (!nu.fix) + 1
    else 1
  }
    #------------------------------tau---------------------------------------------------
  if ("tau"%in%names(family$parameters))
  {  
      all <- getParameters(what = "tau", formula = tau.formula, fix = tau.fix, 
                         start = tau.start,  np = nptau, 
                         start.at = "npl4", start.v = npl4)
     fntau <- all$func
     nptau <- all$np
 start.tau <- all$start
   }
    #------------------------------------------------------------------------------------
    ## definition of the log-likelihood function
    #------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------
   logLikelihood <- function(p) 
       {
       
       if(lpar==1) 
         {   
                fmu <- fnmu(p)      
            if (exact) 
            {
                    if(any(family$family%in%gamlss:::.gamlss.bi.list))
                     {
                tamp1 <- call(pfun, q = y + delta/2, bd=bd, mu = fmu )
                tamp2 <- call(pfun, q = y - delta/2, bd=bd, mu = fmu )
                     }   
                     else    
                     {
                tamp1 <- call(pfun, q = y + delta/2, mu = fmu )
                tamp2 <- call(pfun, q = y - delta/2, mu = fmu )
                     }
               
                 tamp <- eval(tamp1)-eval(tamp2)
             llikcomp <- -(2*log(tamp)) * weights
            }
            else 
            {
             ctamp <-if(any(family$family%in%gamlss:::.gamlss.bi.list))   call(dfun, x = y, bd=bd, mu = fmu )  
                     else    call(dfun, x = y, mu = fmu )   
                tamp <-eval(ctamp)
               # llikcomp <- -(log(tamp) + log(delta)) * weights
                llikcomp <- -2*log(tamp)* weights
            }
         }
       if(lpar==2)
         {
              fmu <- fnmu(p)
            fsigma <- fnsigma(p)      
             if (exact) 
            {         
             if(any(family$family%in%gamlss:::.gamlss.bi.list))
                     {
                tamp1 <- call(pfun, q = y + delta/2, bd=bd, mu = fmu, sigma = fsigma )
                tamp2 <- call(pfun, q = y - delta/2, bd=bd, mu = fmu, sigma = fsigma )
                     }   
                     else    
                     {
                 tamp1 <- call(pfun, q = y + delta/2, mu = fmu )
                 tamp2 <- call(pfun, q = y - delta/2, mu = fmu , sigma = fsigma)
                     }      
                 tamp <- eval(tamp1)-eval(tamp2)
             llikcomp <- -(2*log(tamp)) * weights
            }
            else 
            {
             ctamp <-if(any(family$family%in%gamlss:::.gamlss.bi.list))    call(dfun, x = y, bd=bd, mu = fmu , sigma = fsigma)  
                     else   call(dfun, x = y, mu = fmu , sigma = fsigma)  
                 tamp <-eval(ctamp)
               # llikcomp <- -(log(tamp) + log(delta)) * weights
                llikcomp <- -2*log(tamp)* weights
            }
         }     
       if(lpar==3)
         {
          fmu <- fnmu(p)
        fsigma <- fnsigma(p)
           fnu <- fnnu(p)   
             if (exact) 
            {
            if(any(family$family%in%gamlss:::.gamlss.bi.list))
                     {
                tamp1 <- call(pfun, q = y + delta/2, bd=bd, mu = fmu, sigma = fsigma, nu= fnu )
                tamp2 <- call(pfun, q = y - delta/2, bd=bd, mu = fmu, sigma = fsigma, nu= fnu )
                     }   
                     else    
                     {
               tamp1 <- call(pfun, q = y + delta/2, mu = fmu , sigma = fsigma, nu = fnu)
               tamp2 <- call(pfun, q = y - delta/2, mu = fmu , sigma = fsigma, nu = fnu)
                     }      
              tamp <- eval(tamp1)-eval(tamp2)
          llikcomp <- -(2*log(tamp)) * weights
            }
            else 
            {
             ctamp <-if(any(family$family%in%gamlss:::.gamlss.bi.list))    
                              {call(dfun, x = y, bd=bd, mu = fmu , sigma = fsigma, nu = fnu)}  
                        else  {call(dfun, x = y,        mu = fmu , sigma = fsigma, nu = fnu)}  
                tamp <-eval(ctamp)
               # llikcomp <- -(log(tamp) + log(delta)) * weights
                llikcomp <- -2*log(tamp)* weights
            }
           }
       if(lpar==4)
         {
            fmu <- fnmu(p)
        fsigma <- fnsigma(p)
           fnu <- fnnu(p)
          ftau <- fntau(p)
             if (exact) 
            {
             if(any(family$family%in%gamlss:::.gamlss.bi.list))
                     {
                tamp1 <- call(pfun, q = y + delta/2, bd=bd, mu = fmu, sigma = fsigma, nu= fnu, tau = ftau )
                tamp2 <- call(pfun, q = y - delta/2, bd=bd, mu = fmu, sigma = fsigma, nu= fnu, tau = ftau )
                     }   
                     else    
                     {
               tamp1 <- call(pfun, q = y + delta/2, mu = fmu , sigma = fsigma, nu = fnu, tau = ftau)
               tamp2 <- call(pfun, q = y - delta/2, mu = fmu , sigma = fsigma, nu = fnu, tau = ftau)
                     }   
         tamp <- eval(tamp1)-eval(tamp2)
     llikcomp <- -(2*log(tamp)) * weights
            }
            else 
            {
             ctamp <-if(any(family$family%in%gamlss:::.gamlss.bi.list))    
                              {call(dfun, x = y, bd=bd, mu = fmu , sigma = fsigma, nu = fnu, tau = ftau)}  
                        else  {call(dfun, x = y,        mu = fmu , sigma = fsigma, nu = fnu, tau = ftau)}  
                tamp <-eval(ctamp)
               # llikcomp <- -(log(tamp) + log(delta)) * weights
                llikcomp <- -2*log(tamp)* weights
            }
         }     
        llik <- sum(llikcomp)
        if (llik.output) 
        {
            if (length(p) == 0) 
                cat("-LogLik: ", sum(llikcomp), "\n")
            else cat("-LogLik: ", sum(llikcomp), " ", p, "\n")
        }
               
    z <-  if (lpar==1) list(llik = llik, llikcomp = llikcomp, mu = fmu)
          else if (lpar==2) list(llik = llik, llikcomp = llikcomp, mu = fmu, sigma = fsigma)
          else if (lpar==3) list(llik = llik, llikcomp = llikcomp, mu = fmu, sigma = fsigma, nu = fnu)#
          else  list(llik = llik, llikcomp = llikcomp, mu = fmu, sigma = fsigma, nu = fnu, tau = ftau)
    z
       }
    #------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------
    # the function to be mimimized by nlm()
    #------------------------------------------------------------------------------------
    optFunction <- function(p) 
    {
        tamp <-logLikelihood(p)$llik
        if (llik.output) 
            cat("-LogLik: ", tamp, " (", p, ")", "\n")
        if (is.na(tamp)) 
            1e+20
        else tamp
    }
    #----------------------------------------------------------------------------------
    p0 <- c()
    if (!mu.fix) {
        p0 <- c(mu.start)
        names(p0) <- c(rep("mu.start", length(mu.start)))
    }
    if (!sigma.fix) {
        tamp <- names(p0)
        p0 <- c(p0, sigma.start)
        names(p0) <- c(tamp, rep("sigma.start", length(sigma.start)))
    }
    if (!nu.fix) {
        tamp <- names(p0)
        p0 <- c(p0, nu.start)
        names(p0) <- c(tamp, rep("nu.start", length(nu.start)))
    }
    if (!tau.fix) {
        tamp <- names(p0)
        p0 <- c(p0, tau.start)
        names(p0) <- c(tamp, rep("tau.start", length(tau.start)))
    }
    if (llik.output) {
        cat("No. of parameters: ",npmu, "",npsigma, "",npnu, 
            "",nptau, "\n")
        if (!mu.fix || !sigma.fix || !nu.fix || !tau.fix) {
            cat("Vector of initial conditions on IR^p:", "\n")
            print(p0)
        }
    }
    #browser()
   # llik0 <-logLikelihood(p = p0)
    np0 <- length(p0)
    
    #----basic nlm run
    # getting the control parameters
    typsize <- if(is.null(control$typsize)) abs(p0) else control$typsize
    if(any(typsize <= 0)) 
      {
      warning("the value of typsize supplied is zero or negative the default value of abs(p0) was used instead")
      typsize <- abs(p0)+0.001 #MSThursday, June 1, 2006 at 10:13
      }
    stepmax <- if(is.null(control$stepmax)) sqrt(p0 %*% p0) else control$stepmax
    if(stepmax < 0) 
      {
    warning("the value of stepmax supplied is zero or negative the default value of  sqrt(p0 %*% p0) was used instead")
         stepmax <-  sqrt(p0 %*% p0)
      }
    hessian <- control$hessian
    fscale  <- control$fscale
print.level <- control$print.level
     ndigit <- control$ndigit
    gradtol <- control$gradtol
    steptol <- control$steptol
    iterlim <- control$iterlim
    #------------------------------------------------------------------------------------
    # the actual non linear fitting here                                                |
    #------------------------------------------------------------------------------------
    if (np0 > 0) {
        p.opt <- nlm(optFunction, p = p0, hessian = hessian, fscale = fscale, 
            typsize = rep(1, length(p0)), print.level = print.level, 
            ndigit = ndigit, gradtol = gradtol, steptol = steptol, 
            iterlim = iterlim, stepmax = stepmax)
        z <-logLikelihood(p.opt$estimate)
    }
    else z <-logLikelihood(p0)
    #------------------------------------------------------------------------------------
 coefficients <- p.opt$estimate
# browser()
   #  gradient <- p.opt$gradient
   #     error <- p.opt$error
   #np <- if (lpar==1) npmu   
   #     else if (lpar==2) npmu +npsigma  
   #     else if (lpar==3)npmu +npsigma +npnu
   #     else npmu +npsigma +npnu +nptau
         nobs <- sum(as.numeric(weights))
           mu <- as.vector(z$mu)
        sigma <- as.vector(z$sigma)
           nu <- as.vector(z$nu)
          tau <- as.vector(z$tau)
#----------------------------------------------------------------------------------------
  if (np0 > 0) 
   {
 #   browser()
        cov <- diag(np0)
        if (hessian) {
            if (np0 == 1) 
                cov <- 1/(2*p.opt$hessian)
            else {
                a <- if (any(is.na(p.opt$hessian)) || any(abs(p.opt$hessian) == Inf)) 0
                     else qr(p.opt$hessian)$rank
                if (a == np0) cov <- solve(p.opt$hessian/2)
                else cov <- matrix(NA, ncol = np0, nrow = np0)
            }
        }
        se <- if (hessian) sqrt(diag(cov))
              else NA
        corr <- if (hessian) cov2cor(cov) # cov/(se %o% se) 
                else NA
        
   }
  else coefficients <- se <- cov <- corr  <-  NULL
#------------------ output --------------------------------------------------------------
    out <- list(family = family$family , 
            parameters = names(family$parameters), 
                  call = nlcall, 
                     y = y, 
               control = control, 
               weights = weights, 
            G.deviance = z$llik, 
     #      P.deviance = z$llik, 
                     N = nobs, 
                 rqres = family$rqres,  
                  type = family$type, 
                   aic = z$llik + 2*np0, 
                   sbc = z$llik + log(nobs)*np0,
           df.residual = nobs - np0,
                df.fit = np0,
             converged = p.opt$code, 
                  iter = p.opt$iter, 
      #            pen = 0,
             residuals = eval(family$rqres), 
                method = "JL()",
          coefficients = coefficients,
                    se = se,
                   cov = cov,
                  corr = corr
                    )
   # if(any(family$family%in%gamlss.bi.list))  out$bd <- bd                           
#========================================================================================
##---------------------------------------------------------------------------------------
## this function is used in for outputing the parameters
##---------------------------------------------------------------------------------------
##=======================================================================================
parameterOut <- function(what="mu")
 {
  out <- list() 
  if(family$parameter[[what]]==TRUE && eval(parse(text=paste(what,".fix",sep="")))==FALSE)
     {              
  out$fv <- eval(parse(text=paste("z", what ,sep="$")))
      out$lp <-  eval(parse(text=paste(paste(paste("family", what ,sep="$"),"linkfun(",sep="."),"out$fv)", sep="")))
    out$link <- eval(parse(text=paste(paste("family", what ,sep="$"),"link",sep=".")))
 out$formula <- eval(parse(text=paste(paste("attr(",paste("fn",what,sep=""),sep=""),"\"formula\")",sep=",")))
       names <- if(is.expression(eval(parse(text=paste(paste("attr(",paste("fn",what,sep=""),sep=""),"\"model\")",sep=",")))))
                 {
               eval(parse(text=paste(paste("attr(",paste("fn",what,sep=""),sep=""),"\"parameters\")",sep=",")))
                 }
                else  
                 { 
               eval(parse(text=paste(paste("attr(",paste("fn",what,sep=""),sep=""),"\"model\")",sep=",")))
                 }
out$coefficients <- coefficients[paste(what,"start", sep=".")==names(p0)]
names(out$coefficients)<-names
      out$se <- se[paste(what,"start", sep=".")==names(p0)]
names(out$se)<- names   
      out$df <- eval(parse(text=(paste("np", what, sep="")))) 
   #out$nl.df <- 0
   #out$pen  <- 0 
     }
   else
     { out$fix <- eval(parse(text=paste(what,".fix",sep="")))
        out$df <- 0
        out$fv <- eval(parse(text=paste("z", what ,sep="$")))
     }
#if (is.data.frame(data)) detach(data)
  out
  }
##----------------------------------------------------------------------------------------
##========================================================================================
##  Output for mu model: -----------------------------------------------------------------
 if ("mu"%in%names(family$parameters))   out <- c(out, mu = parameterOut(what="mu") )
 else            out$mu.df <- 0
##  Output for sigma model: --------------------------------------------------------------
 if ("sigma"%in%names(family$parameters))out <- c(out, sigma = parameterOut(what="sigma"))
 else            out$sigma.df <- 0
## Output for nu model: ------------------------------------------------------------------
 if ("nu"%in%names(family$parameters))   out <- c(out, nu = parameterOut(what="nu") )
 else              out$nu.df <- 0
## output for tau model ------------------------------------------------------------------
 if ("tau"%in%names(family$parameters))  out <- c(out, tau = parameterOut(what="tau") )
else               out$tau.df <- 0
    class(out) <- list("nlgamlss", "gamlss")
#if (is.data.frame(data)) detach(data)
    out
}
