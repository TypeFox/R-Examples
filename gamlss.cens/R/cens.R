#----------------------------------------------------------------------------------------
cens <- function (
              family = "NO", 
                type = c("right", "left",  "interval"),
                name = "cens", 
               local = TRUE,
               delta = NULL, 
                ...)
{
#------------------------------------------
## dummy name
     TEST <- "TEST" 
     type <- match.arg(type)   
##   fname <- family
##   if (mode(family) != "character" && mode(family) != "name")
##    fname <- as.character(substitute(family))  
##     fam1 <- eval(parse(text=fname)) 
     fam  <- as.gamlss.family(family) # family 
    fname <- fam$family[[1]] # family name
   family <- c("None", "None") 
   dorfun <- paste("d",fname,sep="") # say dNO
   porfun <- paste("p",fname,sep="") # say pNO
     dfun <- paste(paste("d",fname,sep=""), substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="") # say dNOrc
     pfun <- paste(paste("p",fname,sep=""), substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="") # say pNOrc
if (local)
 {
#--trying to get gamlss sys.frame--  
     rexpr<-regexpr("gamlss",sys.calls())
for (i in 1:length(rexpr)){ 
    position <- i 
    if (rexpr[i]==1) break}
gamlss.environment <- sys.frame(position)      
#--end here------------------------
 }
 else gamlss.environment <- sys.frame(0)
#   generate d within gamlss
    eval(dummy <- cens.d(family = fname, type = type, ...))
    eval(call("<-",as.name(dfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
# generate p within gamlss
    eval(dummy <- cens.p(family = fname, type = type, ...))
    eval(call("<-",as.name(pfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
# rename the family 
   family[[1]] <- paste(fname, substr(type,start=1,stop=1),substr(name,start=1,stop=1), sep="")
   family[[2]] <- paste(type, "censored",fam$family[[2]])
    fam$family <- family
# Global deviance increment  
           sGD <- gsub(dorfun, dfun, deparse(body(fam$G.dev.incr)))
  body(fam$G.dev.incr) <- parse(text=sGD)
# get the no of parameters  
        nopar <- fam$nopar
# check for the delta
 if (length(delta)==0) delta <- rep(NA,nopar) 
 if (length(delta)==1) delta <- rep(delta,nopar)
 if (length(delta)!=nopar)  stop("delta should be the same length the parameters in the family ") 
# now change the first derivatives
  switch(nopar,  
          {
          # 1 parameter 
      fam$dldm <- function(y,mu) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, log=TRUE), "mu", delta=NULL), "gradient")) 
           # mu
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU) 
body(fam$dldm) <- parse(text=sMU[length(sMU)])
          # residuals
      cen <- paste("censored=", "\"", type ,sep="") # i.e.censored= "right"
          sres <- gsub(porfun, paste(pfun,"\" ,", cen),  deparse(fam$rqres)) 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)
         # initial mu
           inimu <- gsub("y", "y[,1]",  deparse(fam$mu.initial))
           inimu <- gsub("expression", "",  inimu)
           fam$mu.initial <- parse(text=inimu)
          #y.valid
           yval <- gsub("y", "y[,1]",  deparse(body(fam$y.valid)))
           yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$y.valid) <- parse(text=yval)[[1]] 
          }, 
          { 
          # 2 parameters  
      fam$dldm <- function(y,mu,sigma) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function(y,mu,sigma) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, log=TRUE), "sigma", delta=NULL), "gradient"))
          # mu 
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU) 
body(fam$dldm) <- parse(text=sMU[length(sMU)]) 
         # sigma  
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA)
body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])
          #d2ldmdd 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdd)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldmdd) <- parse(text=yval)[[1]] 
           #d2ldm2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldm2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldm2) <- parse(text=yval)[[1]] 
          #d2ldd2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldd2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldd2) <- parse(text=yval)[[1]] 
          #residuals
         #cen <- paste("censored=", "\"", type, "censored" , "\"")
         # sres <- gsub(porfun, paste(pfun,cen),  deparse(fam$rqres)) 
           cen <- paste("censored=", "\"", type ,sep="") # i.e.censored= "right"
          sres <- gsub(porfun, paste(pfun,"\" ,", cen),  deparse(fam$rqres)) 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)
          # initial mu
           inimu <- gsub("y", "y[,1]",  deparse(fam$mu.initial))
           inimu <- gsub("expression", "",  inimu)
           fam$mu.initial <- parse(text=inimu)
         # initial sigma
           inisigma <- gsub("y", "y[,1]",  deparse(fam$sigma.initial))
           inisigma <- gsub("expression", "",  inisigma)
           fam$sigma.initial <- parse(text=inisigma) 
         #y.valid
           yval <- gsub("y", "y[,1]",  deparse(body(fam$y.valid)))
           body(fam$y.valid) <- parse(text=yval)[[1]]
           }, 
           {
           # 3 parameters   
      fam$dldm <- function(y,mu,sigma,nu) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function(y,mu,sigma,nu) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "sigma", delta=NULL), "gradient"))
      fam$dldv <- function(y,mu,sigma,nu) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "nu", delta=NULL), "gradient"))
           # mu
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1]))sMU <- sub("NULL",  as.character(delta[1]), sMU)          
body(fam$dldm) <- parse(text=sMU[length(sMU)]) 
           # sigma 
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])
           # nu   
           sNU <- sub("TEST", dfun, body(fam$dldv))
if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)
body(fam$dldv) <- parse(text=sNU[length(sNU)])
           # d2ldmdd 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdd)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldmdd) <- parse(text=yval)[[1]]
           #d2ldmdv 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdv)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldmdv) <- parse(text=yval)[[1]]
           #d2ldddv 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldddv)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldddv) <- parse(text=yval)[[1]] 
            #d2ldm2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldm2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldm2) <- parse(text=yval)[[1]] 
          #d2ldd2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldd2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldd2) <- parse(text=yval)[[1]] 
           #d2ldv2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldv2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldv2) <- parse(text=yval)[[1]]   
           #residuals
      cen <- paste("censored=", "\"", type ,sep="") # i.e.censored= "right"
          sres <- gsub(porfun, paste(pfun,"\" ,", cen),  deparse(fam$rqres)) 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres) 
           # initial mu
           inimu <- gsub("y", "y[,1]",  deparse(fam$mu.initial))
           inimu <- gsub("expression", "",  inimu)
           fam$mu.initial <- parse(text=inimu)
          # initial sigma
           inisigma <- gsub("y", "y[,1]",  deparse(fam$sigma.initial))
           inisigma <- gsub("expression", "",  inisigma)
           fam$sigma.initial <- parse(text=inisigma) 
           # initial nu
           ininu <- gsub("y", "y[,1]",  deparse(fam$nu.initial))
           ininu <- gsub("expression", "",  ininu)
           fam$nu.initial <- parse(text=ininu)
          #y.valid
           yval <- gsub("y", "y[,1]",  deparse(body(fam$y.valid)))
           body(fam$y.valid) <- parse(text=yval)[[1]]
           },
           {
           # 4 parameters
      fam$dldm <- function(y,mu,sigma,nu,tau) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function(y,mu,sigma,nu,tau) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=NULL), "gradient"))
      fam$dldv <- function(y,mu,sigma,nu,tau) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=NULL), "gradient"))
      fam$dldt <- function(y,mu,sigma,nu,tau) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "tau", delta=NULL), "gradient"))
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU)      
body(fam$dldm) <- parse(text=sMU[length(sMU)])   
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])  
           sNU <- sub("TEST", dfun, body(fam$dldv))
if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)           
body(fam$dldv) <- parse(text=sNU[length(sNU)]) 
          sTAU <- sub("TEST", dfun, body(fam$dldt))
if (!is.na(delta[4])) sTAU <- sub("NULL",  as.character(delta[4]), sTAU)
body(fam$dldt) <- parse(text=sTAU[length(sTAU)]) 
           # d2ldmdd 1
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdd)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldmdd) <- parse(text=yval)[[1]]
           # d2ldmdv 2
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdv)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldmdv) <- parse(text=yval)[[1]]
           # d2ldmdt 3 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldmdt)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldmdt) <- parse(text=yval)[[1]]
           # d2ldddv 4
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldddv)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldddv) <- parse(text=yval)[[1]]   
           # d2ldddt 5 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldddt))) 
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldmdt) <- parse(text=yval)[[1]]
           #d2ldvdt 6 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldvdt)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldvdt) <- parse(text=yval)[[1]]
           #d2ldm2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldm2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldm2) <- parse(text=yval)[[1]] 
           #d2ldd2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldd2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T) 
           body(fam$d2ldd2) <- parse(text=yval)[[1]] 
           #d2ldv2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldv2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldv2) <- parse(text=yval)[[1]] 
           #d2ldt2 
          yval <- gsub("y", "y[,1]",  deparse(body(fam$d2ldt2)))
          yval <-  gsub('any[,1]', 'any', yval, fixed=T)
           body(fam$d2ldt2) <- parse(text=yval)[[1]]
          # residuals 
      cen <- paste("censored=", "\"", type ,sep="") # i.e.censored= "right"
          sres <- gsub(porfun, paste(pfun,"\" ,", cen),  deparse(fam$rqres)) 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)
           # initial mu
           inimu <- gsub("y", "y[,1]",  deparse(fam$mu.initial))
           inimu <- gsub("expression", "",  inimu)
           fam$mu.initial <- parse(text=inimu)
          # initial sigma
           inisigma <- gsub("y", "y[,1]",  deparse(fam$sigma.initial))
           inisigma <- gsub("expression", "",  inisigma)
           fam$sigma.initial <- parse(text=inisigma) 
           # initial nu
           ininu <- gsub("y", "y[,1]",  deparse(fam$nu.initial))
           ininu <- gsub("expression", "",  ininu)
           fam$nu.initial <- parse(text=ininu)
          # initial tau
           initau <- gsub("y", "y[,1]",  deparse(fam$tau.initial))
           initau <- gsub("expression", "",  initau)
           fam$tau.initial <- parse(text=initau)
          #y.valid
           yval <- gsub("y", "y[,1]",  deparse(body(fam$y.valid)))
           body(fam$y.valid) <- parse(text=yval)[[1]] 
           })
##          nfam <- function() fam
## formals(nfam) <- formals(fam1) 
 fam
}
