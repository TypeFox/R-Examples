notablyDifferent <- function (object,vars=NULL,threshold=NULL,p=.05,...)
{
   if (missing(object)) stop ("object required.")
   if (class(object)[1] != "yai")  stop ("object must be class yai")
   refIds <- rownames(object$neiDstRefs)
   if (length(refIds) == 0) stop ("references are required")
   cl <- match.call()
   trgIds <- rownames(object$neiDstTrgs)
   
   if (is.null(vars)) vars=xvars(object)
   impObj <- impute.yai(object,vars=vars,observed=TRUE,...)
   if (is.null(impObj)) stop ("no imputations found using this object")
   nuke <- unlist(lapply(impObj,function (x) all(is.na(x))))
   nuke=nuke[nuke]
   if (length(nuke) > 0) impObj <- 
       impObj[,-match(names(nuke),names(impObj)),drop=FALSE]
   nuke <- unlist(lapply(impObj,function (x) is.factor(x)))
   nuke <- nuke[nuke]
   if (length(nuke) > 0) impObj <- 
        impObj[,-match(names(nuke),names(impObj)),drop=FALSE]
   impObj <- na.omit(impObj)
   if (is.null(vars)) vars <- names(impObj)
   vi <- paste(unique(strsplit(vars,".o",fixed=TRUE)))
   vi <- intersect(vi,names(impObj))
   notFound <- setdiff(vars,names(impObj))
   if (length(notFound)>0) warning ("variables not found or had missing values: ",
                                     paste(notFound,collapse=", "))
   if (length(vi) == 0) stop("nothing to compute")
   vo <- paste(vi,"o",sep=".")
   notFound <- setdiff(vo,names(impObj))
   if (length(notFound)>0) warning ("variables not found or had missing values: ",
                                     paste(notFound,collapse=", "))
   vo <- intersect(vo,names(impObj))
   both <- intersect(paste(unique(strsplit(vo,".o",fixed=TRUE))),vi)  
   if (length(both) == 0) stop("no variables with observed and imputed values")
   vo <- paste(both,"o",sep=".")

   var <- unlist(lapply(impObj[,vo,drop=FALSE],var)) 
   names(var) = sub(".o","",names(var),fixed = TRUE)  
   diff <- impObj[,both,drop=FALSE]-impObj[,vo,drop=FALSE]
   for (nam in both) diff[,nam] = (diff[,nam]*diff[,nam])/var[nam]

   rmsd <- sqrt(apply(diff,1,mean))
   if (is.null(threshold)) threshold <- quantile(rmsd[refIds],1-p)
   ans <- list(call=cl,vars=both,threshold=threshold,
               notablyDifferent.refs=sort(rmsd[refIds][rmsd[refIds]>threshold]),
               notablyDifferent.trgs=sort(rmsd[trgIds][rmsd[trgIds]>threshold]),
               rmsdS.refs=sort(rmsd[refIds]),
               rmsdS.trgs=sort(rmsd[trgIds]))
   class(ans) <- "notablyDifferent"
   ans
}


plot.notablyDifferent <- function (x,add=FALSE,...)
{ 
  if (missing(x)) stop ("x required")
  if (class(x) == "list")
  {
    if (!all(unlist(lapply(x,function (xx) class(xx)=="notablyDifferent"))))
       stop ("all members in the x list must be class notablyDifferent") 
    ans <- matrix(unlist(lapply(x,function (xx) 
      {
        all <- c(xx$rmsdS.refs,xx$rmsdS.trgs)
        c(max (all), length(all))
      })),length(x),2,byrow=TRUE)
    xlim <- c(1,max(ans[,2]))
    ylim <- c(0,max(ans[,1]))    
    for (i in 1:length(x))
    {
      plot.notablyDifferent(x[[i]],xlim=xlim,ylim=ylim,add=i>1,col=i,...)
      myusr = par()$usr
      par(usr=c(0,1,0,1))
      cy <- par()$cxy[2]
      cy <- cy*1.1
      if (i == 1) text(x=.05,y=.95,pos=4,"Legend")
      text(x=.05,y=.95-(i+2)*cy,pos=4,
          if (is.null(names(x)[i])) paste("Case",i) else
                      names(x)[i], col=i)
      par(usr=myusr)
    }
  }
  else
  { 
    if (class(x) != "notablyDifferent")  
        stop ("x must be class notablyDifferent")
    all <- c(x$rmsdS.refs,x$rmsdS.trgs)
    pch <- c(rep(1,length(x$rmsdS.refs)),rep(2,length(x$rmsdS.trgs)))
    names(pch) <- names(all)
    all <- sort(all)
    pch[names(x$notablyDifferent.refs)] <- 19
    pch[names(x$notablyDifferent.trgs)] <- 17
    xx <- 1:length(all)
    if (add) points(x=xx,y=all,pch=pch[names(all)],...) else 
               plot(x=xx,y=all,pch=pch[names(all)],
               main="Imputation Error Profile",
               ylab="Scaled RMSD",xlab="Observation",...)
    abline(h=x$threshold,...)
    if (!add)
    {
      myusr = par()$usr
      par(usr=c(0,1,0,1))
      cxy = par()$cxy*1.1
      bxy = c(.05,.95)
      points(x=cxy[1]  +bxy[1], y=bxy[2]-  cxy[2], pch=1)
      points(x=cxy[1]*2+bxy[1], y=bxy[2]-  cxy[2], pch=19)
      text(  x=cxy[1]*3+bxy[1], y=bxy[2]-  cxy[2], pos=4,"References")
      points(x=cxy[1]  +bxy[1], y=bxy[2]-2*cxy[2], pch=2)
      points(x=cxy[1]*2+bxy[1], y=bxy[2]-2*cxy[2], pch=17)
      text(  x=cxy[1]*3+bxy[1], y=bxy[2]-2*cxy[2], pos=4,"Targets")
      par(usr=myusr)
    }
  }
}
  
  
