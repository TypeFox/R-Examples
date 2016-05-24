powcomp.easy <- function(params,M=10^5,model=NULL,Rlaws=NULL,Rstats=NULL,center=FALSE,scale=FALSE) {

  if (is.vector(params)) params <- t(as.matrix(params)) else params <- as.matrix(params) # params n'a qu'une ligne
  nsim <- nrow(params)

  Rcpp <- any(params[,2] == 0)

  if (any(params[,"stat"] == 0) & is.null(Rstats)) stop("'Rstats' should be a list whose components are R functions.")
  if (any(params[,"stat"] == 0)) {
    if (!is.list(Rstats)) stop("'Rstats' should be a list whose components are R functions.")
    for (i in 1:nsim) if ((params[i,"stat"] == 0) & !is.function(Rstats[[i]])) stop(paste("The ",i,"th component of 'Rstats' should be an R function",sep=""))
  }
  if(is.null(Rstats)) Rstats <- list(NULL)
  
  if (Rcpp) {
      if (length(Rlaws) != nsim) stop("When some law indices in second column of 'params' are equal to 0, this means that you will be using some R random generators. In that case, you should provide the names of the random generation functions in the corresponding components of 'Rlaws' list, the other components should be set to NULL.")
      tmp <- gsub(" ","",paste(text=match.call()$Rlaws))
      tmp <- strsplit(substr(tmp,6,nchar(tmp)-1),",")[[1]]
      for (i in 1:nsim) {
          if (!is.null(Rlaws[[i]])) {
              if (class(Rlaws[[i]]) != "function") stop("Each non-null compoment of the list 'Rlaws' should be a (random generation) function.")
              if (params[i,2] != 0) stop(paste("params[",i,",2] should be set to 0.",sep=""))
              if (!all(is.na(params[i,8:11]))) {
                  npartmp <- length(unlist(formals(eval(parse(text=tmp[i])))[-1]))
                  if(sum(!is.na(params[i,8:11])) !=  npartmp) stop(paste("The number of law parameters set for row ",i," of 'params' should be ",npartmp,sep=""))
              }
          } else {
              Rlaws[[i]] <- function(){}
          }
      }
  }

  # We remove the parstatx columns that contain only NA values.
  if (ncol(params) > 11) params <-  params[,c(rep(TRUE,11),apply(params[,-(1:11),drop=FALSE],FUN=function(x) !all(is.na(x)),MARGIN=2)),drop=FALSE]

  params.save <- params
  paramsnew <- matrix(NA,nrow=nsim,ncol=ncol(params)+2)
  nbparstatmax <- ncol(params) - 11
  
# We extend the 'params' matrix so that it will contain the correct number of parameters
# for the laws and for the test statistics considered:
  for (i in 1:nsim) {
    beginvec <- params[i,(1:7)]
    endvec <- na.omit(params[i,(8:11)])
    nbparlaw <- length(endvec)
    
    # add parstat and nbparstat 
    # if parstat = NA or NULL, nbparstat = 0
    parstat <- params[i,-(1:11)]
    if (!all(is.na(parstat)) && !is.null(parstat)) {
      nbparstat <- length(na.omit(parstat))
    } else {
      nbparstat <- 0
    }

    if (nbparlaw > 4) stop("The maximum number of law parameters is 4. Contact the package author to increase this value!")
    if (nbparlaw == 0) paramsnew[i,] <- c(beginvec,nbparlaw,rep(0,4),nbparstat,parstat)
    else if (nbparlaw == 1) paramsnew[i,] <- c(beginvec,nbparlaw,endvec,rep(0,3),nbparstat,parstat)
    else if (nbparlaw == 2) paramsnew[i,] <- c(beginvec,nbparlaw,endvec,rep(0,2),nbparstat,parstat)
    else if (nbparlaw == 3) paramsnew[i,] <- c(beginvec,nbparlaw,endvec,rep(0,1),nbparstat,parstat)
    else paramsnew[i,] <- c(beginvec,nbparlaw,endvec,nbparstat,parstat)
  }

  params <- paramsnew

  colnames(params) <- rep("",ncol(params))
  colnames(params)[3] <- "stat"
  
  usecrit <- rep(1,nsim)
  # if cL and cR are NA, then usecrit <- 0
  for (i in 1:nsim) if (is.na(params[i,5]) && is.na(params[i,6])) usecrit[i] <- 0
  
  # if nsim = 1, we test for only one stat and one law
  if (nsim == 1) {
    params <-  t(as.matrix(c(params[,1:7],usecrit,params[,-(1:7)])))
  } else {
    params <- cbind(params[,1:7],usecrit,params[,-(1:7)])
  }

  params[,5][is.na(params[,5])] <- 0
  params[,6][is.na(params[,6])] <- 0

  decision.len <- nsim
  decision <- rep(0,decision.len)

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

  if (Rcpp | any(params[,"stat"] == 0)) {
      out <- .Call("powcompeasyRcpp",as.integer(M),params=as.double(as.vector(t(params))),as.integer(ncol(params)),decision=as.integer(decision),as.integer(decision.len),
                   as.integer(modelnum),as.list(funclist), as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),as.list(Rlaws),Rstats,as.integer(center), as.integer(scale),NAOK=TRUE,PACKAGE="PoweR")
  } else {
  # call function powcompeasy in C  
      out <- .C("powcompeasy",as.integer(M),params=as.double(as.vector(t(params))),as.integer(ncol(params)),decision=as.integer(decision),as.integer(decision.len),
                as.integer(modelnum), funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np), as.integer(center), as.integer(scale), NAOK=TRUE,PACKAGE="PoweR")
  }
  
  decision <- out$decision/M

  params <- matrix(out$params,nrow=nsim,byrow=TRUE)
  chlaws <- chstats <- rep("",nsim)
  for (i in 1:nsim) {
    indlaw <- params[i,2]
    nbparlaw <- params[i,9]
    if (params[i,2]!= 0) {
        chlaws[i] <- law.cstr(indlaw,params[i,10 + 0:(nbparlaw-1)])$name
    } else {
        tmp2 <- unlist(formals(eval(parse(text=tmp[i])))[-1])
        tmp3 <- c()
        for (j in 1:length(tmp2)) tmp3 <- c(tmp3,paste(names(tmp2[j]),"=",if (is.na(params.save[i,8+j-1])) tmp2[j] else params.save[i,8+j-1],sep=""))
        chlaws[i] <- paste(tmp[i],"(",paste(tmp3,collapse=","),")",sep="")
    }

    indstat <- params[i,3]
    chstats[i] <- stat.cstr(indstat)$name
  }
    
  res <- as.data.frame(matrix(NA,nrow=nsim,ncol=6+nbparstatmax))
  res[,1] <- params[,1]
  res[,2] <- chlaws
  res[,3] <- chstats
  if (nbparstatmax > 0) res[,4+0:(nbparstatmax-1)] <- params[,15+0:(nbparstatmax-1)]
  res[,3+nbparstatmax+1:2] <- params[,c(4,7)]
  res[,6+nbparstatmax] <- round(decision*100,2)
  colnames(res) <- c("n","law","stat",if (nbparstatmax > 0) paste("parstat",1:nbparstatmax,sep=""),"level","alter","power")

  
  return(res)
  
}
