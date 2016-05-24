### This version of lpred() function now performs (for a given gamlss parameter)
### everything that the equivalent predict.lm() or predict.glm() performs 
### except for the newdata argument 
### for newdata the predict.gamlss() must be used
### MS Tuesday, June 15, 2004 at 11:23
lpred <-function (obj, 
                  what = c("mu","sigma","nu","tau"), 
                  parameter= NULL,
                  type = c("link", "response", "terms"), 
                 terms = NULL, 
                se.fit = FALSE, 
                  ... ) 
{
 ## internal function 1
 ## calculates the s.e.'s for a given parameters mu, sigma, nu or tau 
 ## based on predict.lm() function 
 # --------cal.se  start here --------------------------------------------------
 cal.se <- function(obj, what ) 
  {   
  # for some additive terms the X contains zero columns 
      X <- obj[[paste(what,"x",sep=".")]]
     # n <- obj$N
      p <- obj[[paste(what,"qr",sep=".")]]$rank # 
     p1 <- seq(len = p)
    piv <- obj[[paste(what,"qr",sep=".")]]$pivot[p1]
  XRinv <- X[, piv] %*% qr.solve(qr.R(obj[[paste(what,"qr",sep=".")]])[p1, p1]) 
     vl <- drop(XRinv^2 %*% rep(1, p))
smo.var <- obj[[paste(what,"var",sep=".")]] 
     # if smoothing
   if (!is.null(smo.var)) 
    {     vl <- vl +  rowSums(smo.var) }  
  sqrt(vl)             
  }
 #---------cal.se ends here ----------------------------------------------------
 #------------------------------------------------------------------------------
 ## internal function 2
 ## calculates the terms for a given parameters mu, sigma, nu or tau 
 ## based on predict.lm() function 
 #--------cal.terms starts here 
 cal.terms <- function(obj, what, se.fit, terms)
  {
     tt <- obj[[paste(what,"terms",sep=".")]]
      n <- obj$N
      p <- obj[[paste(what,"qr",sep=".")]]$rank 
     p1 <- seq(len = p)
    piv <- obj[[paste(what,"qr",sep=".")]]$pivot[p1] 
     mm <- obj[[paste(what,"x",sep=".")]]
      X <- obj[[paste(what,"x",sep=".")]]
      w <- obj[[paste(what,"wt",sep=".")]]
   beta <- obj[[paste(what,"coefficients",sep=".")]]
smo.mat <- obj[[paste(what,"s",sep=".")]]
smo.var <- obj[[paste(what,"var",sep=".")]]      
     aa <- attr(mm, "assign")
     ll <- attr(tt, "term.labels")
   if (attr(tt, "intercept") > 0) 
     ll <- c("(Intercept)", ll)
    aaa <- factor(aa, labels = ll)
   asgn <- split(order(aa), aaa)
  hasintercept <- attr(tt, "intercept") > 0
   if (hasintercept) 
        {
       asgn$"(Intercept)" <- NULL
                      avx <- colMeans(mm)
               termsconst <- sum(avx[piv] * beta[piv])
        }
 nterms <- length(asgn)
        if (nterms > 0) 
        {
           predictor <- matrix(ncol = nterms, nrow = NROW(X))
           dimnames(predictor) <- list(rownames(mm), names(asgn))
            if (se.fit ) 
            {
                ip <- matrix(ncol = nterms, nrow = NROW(X))
                dimnames(ip) <- list(rownames(X), names(asgn))
                Rinv <- qr.solve(qr.R(obj[[paste(what,"qr",sep=".")]])[p1, p1])
            }
            if (hasintercept) 
                X <- sweep(X, 2, avx)
            unpiv <- rep.int(0, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq(1, nterms, length = nterms)) 
            {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0] <- 0
                predictor[, i] <- if (any(iipiv > 0)) # ms Thursday, May 1, 2008 at 10:13
                  X[, iipiv, drop = FALSE] %*% beta[iipiv]
                else 0
                if (se.fit ) 
                  ip[, i] <- if (any(iipiv > 0))  # ms Thursday, May 1, 2008 at 10:13
                    as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii, 
                      , drop = FALSE])^2 %*% rep.int(1,p)
                  else 0
            }
            # if smoothing 
             if (!is.null(smo.mat)) 
                  {
                       colnameS <- colnames(smo.mat)
                   nofsmoothers <- length(colnameS)
                   for (iN in seq(1,nofsmoothers, length = nofsmoothers))
                     {
                     predictor[,colnameS[iN]] <- predictor[,colnameS[iN]]+smo.mat[,colnameS[iN]]
                      if (se.fit ) 
                       {
                        ip[,colnameS[iN]] <- ip[,colnameS[iN]]+smo.var[,colnameS[iN]]
                       }
                     }  
                  }
            if (!is.null(terms)) 
            {   
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit) 
                  ip <- ip[, terms, drop = FALSE]
            }
        }
        else 
        {
            predictor <- ip <- matrix(0, n, 0)
        }
        attr(predictor, "constant") <- if (hasintercept)  termsconst
                                       else 0
   if (se.fit) 
        list(fit = predictor, se.fit = sqrt(ip) )
   else predictor
  }
#----------- cal.terms stops here  ---------------------------------------------
#-------------------------------------------------------------------------------
## the proper function lp starts here 
#-------------------------------------------------------------------------------
    if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", "")) 
#what <-  as.character(what) ##problems when lpred() is calling other functions
     what <- if (!is.null(parameter))  {
    match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
    type <- match.arg(type)
if (!se.fit)
  {
  pred <- switch(type, 
                 link = obj[[paste(what,"lp", sep = ".")]],
             response = obj[[paste(what,"fv", sep = ".")]], 
                terms = cal.terms(obj = obj, what = what, 
                                  se.fit = se.fit, terms = terms) 
                 )
  }
  else 
  {
  pred <- switch(type, 
                 link = list(fit = obj[[paste(what,"lp",sep=".")]],
                             se.fit = cal.se(obj = obj, what = what) ),              # I need se(mu)=(dmu/deta)*se(eta)
             response = {
             dmudeta <- try( abs(gamlss.family(eval(parse(text=paste(family(obj)[1],"(",what,".link=",# ms 
        eval(parse(text=(paste("obj$",what,".link", sep="")))),")", sep=""))
                  ))[[paste(what,"dr",sep=".")]](obj[[paste(what,"lp",sep=".")]])) , silent = TRUE) 
if (any(class(dmudeta)%in%"try-error")) 
         {
 dmudeta <- abs(eval(parse(text=paste(obj$family[1],"$",what, ".dr","(", "obj$",what,".lp",")", sep=""))))
         }                       
             list(fit = obj[[paste(what,"fv",sep=".")]],                  # eta:  obj[[paste(what, "lp", sep=".")]]
                             se.fit = cal.se(obj = obj, what = what)*dmudeta)          
                         }, 
                terms = cal.terms(obj = obj, what = what , se.fit = se.fit, terms = terms) 
                 )                
  }  
  if (!is.null(obj$na.action))
   {   
    pred <- napredict(obj$na.action, lp)
    if (se.fit) vl <- napredict(obj$na.action, vl)
   }
pred
}
