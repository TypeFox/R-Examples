## This function creates an gamlss.family object which can be used for fitting a 
## a model
## the problem arise when we need different link function 
## the solution at the moment is to as for the different link in the beginning
#----------------------------------------------------------------------------------------
trun <-function (par = c(0), 
              family = "NO", 
                type = c("left", "right", "both"),
                name = "tr", 
               local = TRUE,
               delta = NULL,
             varying = FALSE,     
                ...)
{
#------------------------------------------
     TEST <- "TEST" # dummy name
     type <- match.arg(type)
    #fname <- family
   #if (mode(family) != "character" && mode(family) != "name")
    #fname <- as.character(substitute(family))
    # fam1 <- eval(parse(text=fname))
     fam  <- as.gamlss.family(family)
    fname <- fam$family[[1]] 
   family <- c("None", "None")  
   dorfun <- paste("d",fname,sep="")
   porfun <- paste("p",fname,sep="")
     dfun <- paste(paste("d",fname,sep=""), name, sep="")
     pfun <- paste(paste("p",fname,sep=""), name, sep="")
   #  qfun <- paste(paste("q",fname,sep=""), name, sep="")
   #  rfun <- paste(paste("r",fname,sep=""), name, sep="")
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

if (varying) assign("PAR_", par, envir=gamlss.environment)
#   generate d within gamlss
    eval(dummy <- trun.d(par, family = fname, type = type, varying = varying, ...))
    eval(call("<-",as.name(dfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
# generate p within gamlss
    eval(dummy <- trun.p(par, family = fname, type = type, varying = varying, ...))
    eval(call("<-",as.name(pfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
# rename the family 
   family[[1]] <- paste(paste(fname, name, sep=""))
   family[[2]] <- paste(type, "truncated",fam$family[[2]])
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
          # dldm
      fam$dldm <- function(y,mu) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, log=TRUE), "mu", delta=NULL), "gradient")) 
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU) 
body(fam$dldm) <- parse(text=sMU[length(sMU)])
          # residuals
          sres <- gsub(porfun, pfun,  deparse(fam$rqres))
          if  (fam$type == "Discrete")
             {
              if (varying==FALSE)
              {
                sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",  par[1]+1),  sres) 
                else sres  
              }
              else
              {
                if (type=="left")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
                if (type=="both")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
              }  
             }
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)  
          },
          {
          # 2 parameters   
      # dldm and dldd    
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
          # residuals
          sres <- gsub(porfun, pfun,  deparse(fam$rqres))
         if  (fam$type == "Discrete")
          {
           if (varying==FALSE)
            {
           sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
           else sres  
            }
           else
            {
              if (type=="left")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
              if (type=="both")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
            }  
          }
#         if  (fam$type == "Discrete")
#             {
#               sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
#                        else sres
#             } 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)    
           }, # 3 parameters
          # dldm dldd dldv 
           {   
      fam$dldm <- function(y,mu,sigma,nu) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function(y,mu,sigma,nu) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "sigma", delta=NULL), "gradient"))
      fam$dldv <- function(y,mu,sigma,nu) as.vector(attr(gamlss:::numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "nu", delta=NULL), "gradient"))
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1]))sMU <- sub("NULL",  as.character(delta[1]), sMU)          
body(fam$dldm) <- parse(text=sMU[length(sMU)])  
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])  
           sNU <- sub("TEST", dfun, body(fam$dldv))
if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)
body(fam$dldv) <- parse(text=sNU[length(sNU)])
          sres <- gsub(porfun, pfun,  deparse(fam$rqres))
      if  (fam$type == "Discrete")
      {
        if (varying==FALSE)
        {
          sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
          else sres  
        }
        else
        {
          if (type=="left")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
          if (type=="both")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
        }  
      }
#          if  (fam$type == "Discrete")
#             {
#               sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
#                        else sres
#             }
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)    
           }, # 4 paramers
           # dldm dldd dldv dldt
           {
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
          sres <- gsub(porfun, pfun,  deparse(fam$rqres))
      if  (fam$type == "Discrete")
      {
        if (varying==FALSE)
        {
          sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
          else sres  
        }
        else
        {
          if (type=="left")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
          if (type=="both")  sres <- gsub("ymin = 0",  paste("ymin = PAR_+1"),  sres)
        }  
      }
#         if  (fam$type == "Discrete")
#             {
#               sres <-  if (type=="left"|type=="both")  gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
#                        else sres
#             }
     fam$rqres <- parse(text=sres)    
           })
fam
}
 
