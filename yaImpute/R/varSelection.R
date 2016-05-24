varSelection <- 
function (x,y,method="addVars",yaiMethod="msn",wts=NULL,
         nboot=20,trace=FALSE,
         useParallel=if (.Platform$OS.type == "windows") FALSE else TRUE,...)
{
  if (missing(x)) stop ("x must be specified.")
  if (missing(y)) stop ("y must be specified.")

  okMethods <- c("addVars","delVars")
  if (!(method %in% okMethods)) 
    stop("method=\"",method,"\" must be one of: \"",
      paste0(okMethods,collapse="\", \""),"\"")
  if (is.null(wts)) wts <- rep(1,ncol(y))
  if (useParallel && .Platform$OS.type != "Windows" && 
      requireNamespace ("parallel"))
  {
    myapply <- parallel::mclapply
  } else { 
    if (useParallel) 
      warning("package parallel was not loaded and is not being used")
    myapply <- lapply
  }
   
  cl <- match.call()
  bootstrap <- nboot > 0
  
  # load required packages...this is done here so that forked 
  # processes will have the required packages...different logic is needed
  # to support parallel on windows.
  if (yaiMethod == "gnn") # (GNN), make sure we have package vegan loaded
  {
    if (!requireNamespace ("vegan")) stop("install vegan and try again")
  }
  if (yaiMethod == "ica") # (ica), make sure we have package fastICA loaded
  {
    if (!requireNamespace ("fastICA")) stop("install fastICA and try again")
  }
  if (yaiMethod == "randomForest") # make sure we have package randomForest loaded
  {
    if (!requireNamespace ("randomForest")) stop("install randomForest and try again")
  }     
  if (yaiMethod == "msnPP") # make sure we have package ccaPP loaded
  {
    if (!requireNamespace ("ccaPP")) stop("install ccaPP and try again")
  }
  
  # single variable elimination logic:
  if (method=="delVars")    
  {
    allErr <- unlist( myapply(1:max(1,nboot), function (i,x,y,wts,yaiMethod,...)
               suppressWarnings(grmsd(one=suppressWarnings(yai(x=x,y=y,
               method=yaiMethod,bootstrap=bootstrap,bootstrap,...)),
               ancillaryData=y,wts=wts)) ,x,y,wts,yaiMethod,bootstrap,...) )
    if (trace) cat ("With all vars, mean grmsd (over boostraps) = ",mean(allErr),
            "; stddev=",sd(allErr),"; Num cols = ",ncol(x),"\n",sep="")
     
    xa <- x
    selvars <- list(None=allErr)
    while (ncol(xa) > 1)
    {
      err <- list()
      for (var in 1:ncol(xa)) err[[var]] <- unlist(myapply(1:max(1,nboot), 
          function (i,xa,y,wts,var,yaiMethod,bootstrap,...)
             suppressWarnings(grmsd(one=suppressWarnings(yai(x=xa[,-var, 
                drop=FALSE],y=y, method=yaiMethod,bootstrap=bootstrap,...)),
                ancillaryData=y,wts=wts)), xa,y,wts,var,yaiMethod,bootstrap,...) )
      names(err) <- names(xa)

      # drop the variable that creates the least error by not including it.
      del <- which.min(unlist(lapply(err,mean)))   
      selvars[[names(del)]] <- as.vector(unlist(err[del]))
      xa <- xa[,-del,drop=FALSE]
      remVars <- colnames(xa)
      if (trace) cat ("Delete var= ",names(del),
         "; mean grmsd (over boostraps) = ",mean(err[[del]]),
         "; stddev=",sd(err[[del]]),"; Num cols remaining= ",ncol(xa),"\n",sep="")
    }  
  } else if (method=="addVars") {
    remVars <- NULL  
    selvars <- list()
    keep <- NULL
    while (length(keep) < ncol(x))
    {    
      err <- list()
      for (var in setdiff(names(x),keep)) 
      {
        xa <- x[,c(keep,var),drop=FALSE]
        err[[var]] <- unlist(myapply(1:max(1,nboot), 
          function (i,xa,y,wts,yaiMethod,bootstrap,...) 
            suppressWarnings(grmsd(one=suppressWarnings(yai(x=xa,y=y,
              method=yaiMethod,bootstrap=bootstrap,...)),
          ancillaryData=y,wts=wts)), xa,y,wts,yaiMethod,bootstrap,...) )
      }

      # keep the variable that reduces the error the most when it is included
      add <- names(which.min(unlist(lapply(err,mean))))
      selvars[[add]] <- as.vector(unlist(err[add]))
      keep <- c(keep,add)  
      if (trace) cat ("Added var= ",add,
        "; mean grmsd (over boostraps) = ",mean(err[[add]]),
        "; stddev=",sd(err[[add]]),"; Num cols being used= ",
        ncol(xa),"\n",sep="")
    }  
  }
  err <- lapply(selvars,function (x) c (mean(x),sd(x)))
  rn <- names(err)
  err <- matrix(unlist(err),ncol=2,byrow=TRUE)
  rownames(err) <- rn
  colnames(err) <- c("mean","sd")
  rtn <- list(call=cl,grmsd=err,allgrmsd=selvars,method=method)
  if (!is.null(remVars)) rtn$remVars <- remVars
  class(rtn) <- c("varSel","list")
  rtn
}

plot.varSel <- function (x,main=NULL,nbest=NULL,arrows=TRUE,...)
{
  if (missing(x)) stop ("x must be present")
  x = x
  if (!inherits(x,"varSel")) stop ("class of x must be varSel")
  if (is.null(main)) main <- 
    switch(x$method, addVars="Mean distance as variables are added",
                     delVars="Mean distance as variables are removed",
                     stop("method '",x$method,"' not found in x"))
  if (is.null(nbest)) nbest <- length(bestVars(x))
  bcc <- rep(gray(.35),length(x$allgrmsd))
  if (is.null(x$remVars)) bcc[1:nbest] <- "black" else 
    bcc[(length(bcc)-nbest+2):length(bcc)] <- "black"
  orgmar <- par()$mar
  par(mar=c(6,3,2,1))
  boxplot(x$allgrmsd,border=bcc,horizontal= FALSE, las=2,main=main,...)
  lines(x$grmsd[,1]            ,x=1:nrow(x$grmsd),col=2)
  lines(x$grmsd[,1]+x$grmsd[,2],x=1:nrow(x$grmsd),col=3,lty=2)
  lines(x$grmsd[,1]-x$grmsd[,2],x=1:nrow(x$grmsd),col=3,lty=2)

  if (is.null(x$remVars))
  {
    if (arrows & nbest>2)
    {
      ytop <- par()$usr[4]-par()$cxy[2]
      arrows(x0=nbest,y0=ytop,x1=1+par()$cxy[1]*.5,y1=ytop,length=par()$cin[2]*.5) 
      text (paste0("Best ",nbest),x=nbest,y=ytop-par()$cxy[2],pos=2)
    }
  } else {
    ytop <- par()$usr[4]-par()$cxy[2]
    txt <- paste0("Remaining variable",
          if (length(x$remVars) > 1) "s" else "",": ",
          paste0(x$remVars,collapse=", "))
    text (txt,x=par()$usr[2],y=par()$usr[3]+(par()$cxy[2]*.5),pos=2) 
    x0 <- length(x$allgrmsd)-nbest+2   
    if (arrows & nbest>2)
    {
      arrows(x0=x0,y0=ytop,x1=length(x$allgrmsd)-par()$cxy[1]*.5,y1=ytop,
        length=par()$cin[2]*.5) 
      text (paste0("Best ",nbest),x=x0,y=ytop-par()$cxy[2],pos=4)
    }
  }
  par(mar=orgmar)
  invisible(NULL)
}

bestVars <- function (obj,nbest=NULL)
{
  if (missing(obj)) stop ("obj must be present")
  if (!inherits(obj,"varSel")) stop ("class of obj must be varSel")
  grmsd <- switch(obj$method,
    addVars=obj$grmsd[,1],
    delVars=rev(obj$grmsd[2:nrow(obj$grmsd),1]),
    stop("method '",obj$method,"' not found in obj"))
  if (is.null(nbest)) 
  {
    le <- length(grmsd)
    nbest <- if (le > 2) 
    {
      s <- (grmsd[le]-grmsd[1])/le
      ss <- unlist(lapply(2:(le-1), function (i,ss) (ss[i-1]+ss[i])/2,diff(grmsd)))
      sb <- abs(s) > abs(ss)
      if (any(sb)) min(le,which.max(sb)+1) else le
    } else le  
  }
  vars <- if (!is.null(obj$remVars)) c(obj$remVars,names(grmsd)) else names(grmsd)
  vars[1:min(nbest,length(vars))]
}


