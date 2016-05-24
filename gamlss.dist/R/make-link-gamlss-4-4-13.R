#----------------------------------------------------------------------------------------
# introducing own link MS Wednesday, September 21, 2005
# the user is allowed to used its own link now by defining 
# the functions 
# own.linkfun, own.linkin, own.mu.eta and own.valideta
# the distribution also has to change to include "own" in the available
# link for the parameter 
# intoducing the shifted log and logit links  MS Sunday, February 20, 2005 
# this will replace the make.link used for GLM's and GAM's
# the add link functions are 
# (1) lofshifted needing one parameters the left shift
# (2)logitshift.5 needing two parameters the left and rirght shift
# logit, 
# probit, 
# cauchit, 
# cloglog, 
# identity, 
# log, 
# sqrt, 
# "1/mu^2, 
# logshiftto1,
# logshiftto2
# logshiftto2 - Slog
# own
# inverse


 make.link.gamlss<-function (link) 
{
if (is.character(link) && length(grep("^power", link) > 0)) 
    {
        warning("calling make.link(\"power(z)\") is deprecated", 
            domain = NA)
        return(eval(parse(text = link)))
    }
    else if (!is.character(link) && !is.na(lambda <- as.numeric(link))) 
    {
        warning("calling make.link(number) is deprecated", domain = NA)
        return(power(lambda))
    }
    else switch(link, logit = 
    {
        linkfun <- function(mu)  qLO(mu) #.Call("logit_link", mu, PACKAGE = "stats")
        linkinv <- function(eta)         #.Call("logit_linkinv", eta, PACKAGE = "stats")
       {
            thresh <- -qLO(.Machine$double.eps)
               eta <- pmin(thresh, pmax(eta, -thresh))
            pLO(eta)
           }  
         mu.eta <- function(eta) pmax(dLO(eta), .Machine$double.eps)#.Call("logit_mu_eta", eta, PACKAGE = "stats")
       valideta <- function(eta) TRUE
    }, probit = 
    {
        linkfun <- function(mu) qnorm(mu)
        linkinv <- function(eta) 
           {
            thresh <- -qnorm(.Machine$double.eps)
               eta <- pmin(thresh, pmax(eta, -thresh))
            pnorm(eta)
           }
        mu.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
      valideta <- function(eta) TRUE
    }, cauchit = 
    {
        linkfun <- function(mu) qcauchy(mu)
        linkinv <- function(eta) 
        {
            thresh <- -qcauchy(.Machine$double.eps)
            eta <- pmin(pmax(eta, -thresh), thresh)
            pcauchy(eta)
        }
        mu.eta <- function(eta) pmax(dcauchy(eta), .Machine$double.eps)
        valideta <- function(eta) TRUE      
     }, cloglog = 
     {
        linkfun <- function(mu) log(-log(1 - mu))
        linkinv <- function(eta) pmax(pmin(-expm1(-exp(eta)), 
            1 - .Machine$double.eps), .Machine$double.eps)
        mu.eta <- function(eta) 
        {
            eta <- pmin(eta, 700)
            pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
        }
        valideta <- function(eta) TRUE
    }, identity = 
    {
       linkfun <- function(mu) mu
       linkinv <- function(eta) eta
        mu.eta <- function(eta) rep(1, length(eta))
      valideta <- function(eta) TRUE
    }, log = 
    {
        linkfun <- function(mu) log(mu)
        linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
        mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
        valideta <- function(eta) TRUE
    }, sqrt = 
    {
       linkfun <- function(mu) mu^0.5
       linkinv <- function(eta) eta^2
        mu.eta <- function(eta) 2 * eta
       valideta <- function(eta) all(eta > 0)
    }, "1/mu^2" = 
    {
       linkfun <- function(mu) 1/mu^2
       linkinv <- function(eta) 1/eta^0.5
        mu.eta <- function(eta) -1/(2 * eta^1.5)
      valideta <- function(eta) all(eta > 0)
     }, "mu^2" = 
    {
       linkfun <- function(mu)  mu^2
       linkinv <- function(eta) eta^0.5
        mu.eta <- function(eta)  0.5 * (eta^-0.5)
      valideta <- function(eta) all(eta > 0)
     },logshiftto1 = 
     {    # renamed Thursday, March 27, 2008 MS Saturday, February 19, 2005   
       linkfun <- function(mu)  log(mu-1+0.00001)
       linkinv <- function(eta) 1+pmax(.Machine$double.eps, exp(eta)) 
        mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
      valideta <- function(eta) TRUE   
     },  logshiftto2 = 
     {    # 6=7-2012  
       linkfun <- function(mu)  log(mu-2+0.00001) # changed 6-7-12
       linkinv <- function(eta) 2+pmax(.Machine$double.eps, exp(eta)) 
        mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
      valideta <- function(eta) TRUE   
     },
        logshiftto0  =
     {
       linkfun <- function(mu) { log(mu - 1e-05)}
       linkinv <- function(eta){ 1e-05 + pmax(.Machine$double.eps, exp(eta))}
        mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
      valideta <- function(eta) TRUE
     }, 
    Slog  = # identical to logshiftto0 
     {
       linkfun <- function(mu) { log(mu - 1e-05)}
       linkinv <- function(eta){ 1e-05 + pmax(.Machine$double.eps, exp(eta))}
        mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
      valideta <- function(eta) TRUE
}, 
       #logitshift.5 = { # MS Saturday, February 19, 2005 depreciated 
       #linkfun <- function(mu, shift = par )           
       #                  log((mu-shift[1])/(shift[2]-mu))
       #linkinv <- function(eta,  shift = par) 
       #     {
       #     thresh <- -log(.Machine$double.eps)
       #        eta <- pmin(thresh, pmax(eta, -thresh))
       #               shift[2]-(shift[2]-shift[1])/(1 + exp(eta))
       #     } 
       # mu.eta <- function(eta, shift = par ) 
       #     {
       #     thresh <- -log(.Machine$double.eps)
       #        res <- rep(.Machine$double.eps, length(eta))
       #     res[abs(eta) < thresh] <- ((shift[2]-shift[1])*exp(eta)/(1 + exp(eta))^2)[abs(eta) < thresh]
       #     res
       #     }
       #valideta <- function(eta) TRUE       
       #}, 
         own     = 
     {
       linkfun <- function(mu)   eval(body(get("own.linkfun", envir=globalenv())))  
       linkinv <- function(eta)  eval(body(get("own.linkinv", envir=globalenv())))  
        mu.eta <- function(eta)  eval(body(get("own.mu.eta", envir=globalenv())))   
      valideta <- function(eta)  eval(body(get("own.valideta", envir=globalenv())))      
      # linkfun <- function(mu) 
      #                    { 
      #             if (!exists("own.linkfun")) stop("own.linkfun is not defined") else own.linkfun(mu) 
      #                    }
      # linkinv <- function(eta) 
      #                    { 
      #             if (!exists("own.linkinv")) stop("own.linkinv is not defined") else own.linkinv(eta) 
      #                    }
      #  mu.eta <- function(eta) 
      #                    { 
      #             if (!exists("own.mu.eta")) stop("own.mu.eta is not defined") else own.mu.eta(eta) 
      #                    }
      #valideta <- function(eta) 
      #                    { 
      #             if (!exists("own.valideta")) stop("own.valideta is not defined") else own.valideta(eta) 
      #                    }
     }, inverse = {
       linkfun <- function(mu) 1/mu
       linkinv <- function(eta) 1/eta
        mu.eta <- function(eta) -1/(eta^2)
      valideta <- function(eta) all(eta != 0)
    }, stop(sQuote(link), " link not recognised"))
   structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
        valideta = valideta, name = link), class = "link-gamlss")
}
#----------------------------------------------------------------------------------------
#own.linkfun <- function(mu) { stop("own.linkfun is not defined") }
#own.linkinv <- function(eta){ stop("own.linkinv is not defined")}
#own.mu.eta <- function(eta){ stop("own.mu.eta is not defined")} 
#own.valideta <- function(eta){ stop("own.valideta is not defined")} 

# own.linkfun <- get(" own.linkfun", envir=globalenv()) 
#own.linkinv <- get("own.linkinv", envir=globalenv()) 
#own.mu.eta <- get("own.mu.eta", envir=globalenv()) 
#own.valideta <- get("own.valideta", envir=globalenv()) 
#---------------------------------------------------------------------------------------
show.link <- function(family = "NO")
 {
     what <-  c("mu", "sigma", "nu", "tau")
  family1 <- as.gamlss.family(family)  
link.list <- list()  
     npar <- family1$nopar
for (i in 1:npar)
   {
              name <- what[i]
 link.list[[name]] <- body(family)[[1+i]][[3]][[5]]
   }
link.list
 }
#---------------------------------------------------------------------------------------- 
