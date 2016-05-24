# this is an attempt to create 2 way intraction plot fo factors withing GAMLSS
# the problem arrised  analysing the BAT (illicit ciggaret packs in Greece)
# to DD
# #-----------------------------------------------------------------------------
plot2way <- function(obj,
                              terms = list(), # I think should be a list
                               what = c("mu","sigma","nu","tau"), 
                               parameter = NULL,
                              #type = c("link", "response"),
                        show.legend = TRUE,
                               ... )
{
        what <- if (!is.null(parameter))  {
       match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
        if (!is.gamlss(obj))  stop(paste("This is not a gamlss object", "\n", ""))
        if (!is.character(terms))  stop("the terms argument should be a character vector with two elements" )
#       type <- match.arg(type)
    theTerms <- lpred(obj, what=what, type="terms")
## check whether the names are in the list
      nfact1 <- as.character(terms[1])
      if (nfact1=="NULL") stop("the name of the first factor is not given")
     nfact2 <- as.character(terms[2])
      if (nfact2=="NULL") stop("the name of the secondt factor is not given")
wheretofind <- union(grep(c(nfact1),colnames(theTerms)),grep(c(nfact2),colnames(theTerms)))
wheretofind <-  wheretofind[order(wheretofind)]  
##      exists(as.character(obj$call$data)
     if   (!is.null(obj$call$data)) 
      {
       DaTa <- eval(obj$call$data)
      } else stop("the fitted model needs the data argument")
      fact1 <-  factor(with(DaTa, get(nfact1)))         
      fact2 <-  factor(with(DaTa, get(nfact2)))         
#       with(as.character(obj$call$data), nfact1) 
        SYM <- tapply(rowSums(theTerms[,wheretofind]), list(fact1, fact2) , mean)
  namelevf1 <- levels(fact1)
  namelevf2 <- levels(fact2)
     nlevf1 <- nlevels(fact1)
     nlevf2 <- nlevels(fact2)
    matplot( SYM,, type="n", xlab=nfact1,  ylab="partial contribution", xlim=c(.5, nlevf1+0.5))
    matpoints(SYM,pch=19)
    matlines(SYM)
    if (show.legend)
    {
      legend("topleft", legend =namelevf2, col = seq(1,nlevf2), lty = 1:nlevf2, ncol = 1) 
    }
}
