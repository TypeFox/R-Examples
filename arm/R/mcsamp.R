# mcsamp function (wrapper for mcmcsamp in lmer())
# Quick function to run mcmcsamp() [the function for MCMC sampling for
# lmer objects) and convert to Bugs objects for easy display

mcsamp.default <- function (object, n.chains=3, n.iter=1000, n.burnin=floor(n.iter/2), 
    n.thin=max(1, floor(n.chains * (n.iter - n.burnin)/1000)), 
    saveb=TRUE, deviance=TRUE, make.bugs.object=TRUE)
{
  cat("mcsamp() used to be a wrapper for mcmcsamp() in lme4.\nCurrently, mcmcsamp() is no longer available in lme4.\nSo in the meantime, we suggest that users use sim() to get\nsimulated estimates.\n")
}



#mcsamp.default <- function (object, n.chains=3, n.iter=1000, n.burnin=floor(n.iter/2), 
#    n.thin=max(1, floor(n.chains * (n.iter - n.burnin)/1000)), 
#    saveb=TRUE, deviance=TRUE, make.bugs.object=TRUE)
#{
#  
#  if (n.chains<2) stop ("n.chains must be at least 2")
#  n.keep <- n.iter - n.burnin
#  first.chain <- mcmcsamp (object, n.iter, saveb=saveb, trans=TRUE, deviance=deviance)[(n.burnin+1):n.iter,]
#  n.parameters <- ncol(first.chain)
#    
#  if (deviance) {
#    sims <- array (NA, c(n.keep, n.chains, n.parameters+1))
#  }
#  if (!deviance){
#    sims <- array (NA, c(n.keep, n.chains, n.parameters))
#  }
#
#  pred.names <- attr(terms(object), "term.labels")
#  par.names <- dimnames(first.chain)[[2]]
#  par.names <- gsub("b.", "b@", par.names, ignore.case = FALSE, # Su: rename "b.*" to ""
#                    extended = TRUE, perl = FALSE,
#                    fixed = TRUE, useBytes = FALSE)   
#  par.names <- gsub("b@.*", "", par.names, ignore.case = FALSE, 
#                    extended = TRUE, perl = FALSE,
#                    fixed = FALSE)    
#  par.names <- par.names[is.na(match(par.names,""))] 
#  name.chk.idx <- as.logical(match(par.names, pred.names, nomatch=0))
#  par.names[name.chk.idx] <- paste("beta", par.names[name.chk.idx], sep=".")
#
#  if (saveb){
#    b.hat <- se.coef (object)                   # Su: use se.coef() 
#    n.groupings <- length(b.hat) - 1
#    J <- NA
#    K <- NA
#    for (m in 1:n.groupings){
#      J[m] <- dim(b.hat[[m+1]])[1]
#      K[m] <- dim(b.hat[[m+1]])[2]
#      var.names <- paste (abbreviate(names(b.hat)[m+1],4), ".",
#                          unlist (dimnames(b.hat[[m+1]])[2]), sep="") ##sep="."
#      par.names <- c (par.names,
#        paste (rep(var.names,J[m]), "[", rep(1:J[m],each=K[m]), "]", sep=""))
#    }
#  }
#  sims[,1,1:n.parameters] <- first.chain
#
#  for (k in 2:n.chains){
#    sims[,k,1:n.parameters] <- mcmcsamp (object, n.iter, saveb=saveb, trans=TRUE, deviance=deviance)[(n.burnin+1):n.iter,]
#  }
#  
#  select <- c(rep(FALSE, n.thin-1),TRUE)
#  sims <- sims[select,,]
#  
#  for (j in 1:n.parameters){
#    if (pmatch("log(sigma^2)", par.names[j], nomatch=0)){#=="log(sigma^2)"){
#      par.names[j] <- "sigma.y"
#      sims[,,j] <- exp (sims[,,j]/2)
#    }
#    else if (pmatch("log(", par.names[j], nomatch=0)){#(substr(par.names[j],1,4)=="log("){
#      par.names[j] <- paste ("sigma.", substr(par.names[j], 5, nchar(par.names[j])-1), sep="")
#      sims[,,j] <- exp (sims[,,j]/2)
#    }
#    else if (pmatch("atanh(", par.names[j], nomatch=0)){#(substr(par.names[j],1,6)=="atanh("){
#      par.names[j] <- paste ("rho.", substr(par.names[j], 7, nchar(par.names[j])-1), sep="")
#      sims[,,j] <- tanh (sims[,,j])
#    }
#    #else if (substr(par.names[j],1,4)=="eta."){#(pmatch("eta.", par.names[j], nomatch=0)){#(substr(par.names[j],1,4)=="eta."){
#    #  par.names[j] <- paste ("", substr(par.names[j], 5, nchar(par.names[j])), sep="")
#    #  par.names[j] <- par.names[j]
#    #}
#    else if (pmatch("deviance", par.names[j], nomatch=0)){#(par.names[j]=="deviance"){          # Su: keep par.names for "deviance"
#        sims[,,n.parameters+1] <- sims[,,j] 
#        sims <- sims[,,-j]                          # Su: delete deviance value from sims
#    } 
##    else {  
##    }
#  } 
#  par.names <- gsub("(", "", par.names, ignore.case = FALSE,
#                    extended = TRUE, perl = FALSE,
#                    fixed = TRUE, useBytes = FALSE)   
#  par.names <- gsub(")", "", par.names, ignore.case = FALSE,
#                    extended = TRUE, perl = FALSE,
#                    fixed = TRUE, useBytes = FALSE) 
# # par.names <- gsub(".Intercept", ".Int", par.names, ignore.case = FALSE,
##                    extended = TRUE, perl = FALSE,
##                    fixed = TRUE, useBytes = FALSE) 
#  par.names <- gsub("rescale", "z.", par.names, ignore.case = FALSE,
#                    extended = TRUE, perl = FALSE,
#                    fixed = TRUE, useBytes = FALSE) 
#  
#  par.names <- par.names[is.na(match(par.names,"deviance"))] # Su: delete par.names for "deviance"
#  
#  if (deviance){
#      dimnames(sims) <- list (NULL, NULL, c(par.names,"deviance"))
#  }
#  if (!deviance){
#    dimnames(sims) <- list (NULL, NULL, par.names)
#  }
#  if (make.bugs.object){
#    return (as.bugs.array (sims, program="lmer", n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=deviance))
#  }
#  else {
#    return (sims)
#  }
#}
#
#
#
setMethod("mcsamp", signature(object = "merMod"),
    function (object, ...)
{
    mcsamp.default(object, deviance=TRUE, ...)
}
)
#
#setMethod("mcsamp", signature(object = "glmer"),
#    function (object, ...)
#{
#    mcsamp.default(object, deviance=FALSE, ...)
#}
#)
