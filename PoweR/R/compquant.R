compquant <- function(n,law.index,stat.index,probs=NULL,M=10^5,law.pars=NULL,stat.pars=NULL,model=NULL,Rlaw=NULL,Rstat=NULL,center=FALSE,scale=FALSE) {

  if ((stat.index == 0) & is.null(Rstat)) stop("'Rstat' should be a function when 'stat.index' is equal to 0.")
  if ((stat.index != 0) & is.function(Rstat)) stop("'stat.index' should be set to 0 when 'Rstat' is a function.")
  
  if (is.function(Rlaw) & (law.index != 0)) stop("You should set 'law.index' to 0 when 'Rlaw' is a (random generating) function.")
    
  nbparlaw <- length(law.pars)
  if (nbparlaw > 4) stop("The maximum number of law parameters is 4. Contact the package author to increase this value.")
  # This is a technical requirement for the C++ lawxxxx function call below,
  # because parlaw (arg. params of the C function) should always be a 4-length vector of double.
    law.pars.save <- law.pars
  law.pars <- c(law.pars,rep(0,4-nbparlaw))

  if (!is.function(Rstat) & (is.null(stat.pars) || is.na(stat.pars))) {
    stat.pars <- rep(0,getnbparstats(stat.index)) # C++ technical requirement as above.
    nbparstat <- 0 # The default values will be used by the C++ function.
  } else {
    nbparstat <- length(stat.pars)
  }
  
# TO BE CHECKED WHEN I WILL REALLY ADD THE MODEL FEATURE!!
  if (is.double(model) || is.integer(model)) {
    modelnum <- model
    funclist <- list(function(){})
    thetavec <- 0
    xvec <- 0
    p <- length(thetavec)
    np <- length(xvec)
  } else {
    if (is.null(model)) {
      modelnum <- 1
      funclist <- list(function(){})
      thetavec <- 0
      xvec <- 0
      p <- length(thetavec)
      np <- length(xvec)
    } else { # model should be a list (function(x,thetavec,xvec),theta,xvec)
      modelnum <- 0
      funclist <- list(model[[1]])
      thetavec <- model[[2]]
      xvec <- model[[3]]
      p <- length(thetavec)
      np <- length(xvec)     
    }
  }

  if (is.null(Rstat)) Rstat <- function(){}

  if ((law.index == 0) | (stat.index == 0)) {

      out <- .Call("compquantRcpp",n=as.integer(n),law=as.integer(law.index),stat=as.integer(stat.index),M=as.integer(M),statvec=as.double(rep(0,M)),nbparlaw=as.integer(nbparlaw),law.pars=as.double(law.pars),nbparstat=as.integer(nbparstat),stat.pars=as.double(stat.pars),as.integer(modelnum), as.list(funclist), as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),as.function(Rlaw), as.function(Rstat),as.integer(center), as.integer(scale),PACKAGE="PoweR")
      tmp <- paste(text=match.call()$Rlaw)
      tmp2 <- unlist(formals(eval(parse(text=tmp)))[-1])
      tmp3 <- c(paste(names(tmp2),"=",if (is.null(law.pars.save)) tmp2 else law.pars.save,sep=""))
        lawname <- paste(tmp,"(",paste(tmp3,collapse=","),")",sep="")
      if (is.null(law.pars)) {out$law.pars <- tmp2;out$nbparlaw <- length(out$law.pars)}

      statname <- "stat0"

  } else {
      out <- .C("compquant",n=as.integer(n),law=as.integer(law.index),stat=as.integer(stat.index),M=as.integer(M),statvec=as.double(rep(0,M)),nbparlaw=as.integer(nbparlaw),law.pars=as.double(law.pars),nbparstat=as.integer(nbparstat),stat.pars=as.double(stat.pars),as.integer(modelnum), funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),as.integer(center), as.integer(scale),PACKAGE="PoweR")
      lawname <- law.cstr(law.index,out$law.pars[1:getnbparlaws(law.index)])$name

      statname <- stat.cstr(stat.index)$name
  }



  return(list(stat=out$statvec,quant=quantile(out$statvec,if (is.null(probs)) c(0.025,0.05,0.1,0.9,0.95,0.975) else probs),lawname=lawname,law.pars=out$law.pars[1:out$nbparlaw],statname=statname,stat.pars=out$stat.pars[1:out$nbparstat]))

}
  
