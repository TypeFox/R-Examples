
buildConsensus <- function (reps, noTrgs=FALSE, noRefs=FALSE, k=NULL)
{
  if (!all(unlist(lapply(reps,function (x) class(x) == "yai")))) 
    stop("class must be yai for all reps")

  cl=match.call()
  
  if (length(reps) == 1) 
  {
    warning ("only one rep, nothing to do.")
    return (reps[[1]])
  }
  
  mink = min(unlist(lapply(reps,function (x) x$k)))
  if (!is.null(k)) 
  {
    if (k > mink) 
    {
      warning ("k=",k," ignored, replaced by k=", mink)
      k = mink  
    }
  } else { k = mink }

  rowNT = NULL
  rowNR = NULL
  
  if (!noTrgs)
  {
    rowNT = lapply(reps,function (x) if (is.null(x$neiIdsTrgs)) NULL else rownames(x$neiIdsTrgs))
    if (all(is.null(unlist(rowNT)))) noTrgs = TRUE else 
    {
      if (all(unlist(lapply(rowNT,function(x,y) identical(x,y), rowNT[[1]])))) rowNT = rowNT[[1]] else
      {
        rids = NULL
        for (r in rowNT) rids = union(rids,r)
        rowNT = rids
      }
    }
  }

  if (!noRefs)
  {
    rowNR = lapply(reps,function (x) if (is.null(x$neiIdsRefs)) NULL else rownames(x$neiIdsRefs))
    if (all(is.null(unlist(rowNR)))) noRefs = TRUE else 
    {
      if (all(unlist(lapply(rowNR,function(x,y) identical(x,y), rowNR[[1]])))) rowNR = rowNR[[1]] else
      {
        rids = NULL
        for (r in rowNR) rids = union(rids,r)
        rowNR = rids
      }
    }
  }
  if (!noTrgs) 
  {
    rowNT = setdiff(rowNT,rowNR)
    if (length(rowNT) == 0) 
    {
      rowNT = NULL
      noTrgs = TRUE
    }
  } else rowNT = NULL
  
  if (noTrgs & noRefs) stop("Can't find neighbors in any objects")  
 
  #build bootstrap sample weights
  
  cnts = table(unlist(lapply(reps,function (x) unique(x$bootstrap))))
  if (length(cnts) == 1) wts = NULL else 
  {
    wts  = length(reps)/cnts
    names(wts) = names(cnts)
  }

  # define an internal function to do the mergers
  mkMerge <- function (kIds,kDst,rown,nreps,wts)
  {
    if (is.null(rown)) # assume all rows line up which is much faster
    {
      kIds = matrix(unlist(kIds),ncol=nreps)
      kDst = matrix(unlist(kDst),ncol=nreps)
    } else {            # need to merge the rows to create the matrix.
      kid = matrix("",ncol=length(reps),nrow=length(rown))
      kds = matrix(NA,ncol=length(reps),nrow=length(rown))
      rownames(kid) = rown
      rownames(kds) = rown
      for (i in 1:nreps) 
      {
        idx = match(rown,names(kIds[[i]]))       
        kid[,i] = kIds[[i]][idx]
        kds[,i] = kDst[[i]][idx]
      }
      kid[kid == "" ] = NA
      kIds = kid
      kDst = kds
    }
    newIds = apply(kIds,1,function (x,wts) 
      {
        cnts = table(x)
        if (is.null(wts)) names(which.max(cnts)) else names(which.max(cnts*wts[names(cnts)]))       
      }, wts)
    newDst = vector("numeric",length(newIds))

    for (i in 1:length(newIds))
    {
      inc = kIds[i,] ==  newIds[i]
      inc[is.na(inc)] = FALSE        
      newDst[i] = mean(kDst[i,inc])
    }
    list(ids=newIds,dst=newDst)
  } ############## end of function definition

  idsT = if (noTrgs) NULL else list()
  dstT = if (noTrgs) NULL else list()
  idsR = if (noRefs) NULL else list()
  dstR = if (noRefs) NULL else list()
  
  mxk = k
  for (k in 1:mxk)
  {
    if (!noTrgs) 
    {
      kIds = lapply(reps,function (x,k) x$neiIdsTrgs[,k], k)
      kDst = lapply(reps,function (x,k) x$neiDstTrgs[,k], k)
      tmp = mkMerge(kIds,kDst,rowNT,length(reps),wts)
      idsT[[k]] = tmp$ids
      dstT[[k]] = tmp$dst
    } 
    if (!noRefs)
    {      
      kIds = lapply(reps,function (x,k) x$neiIdsRefs[,k], k)
      kDst = lapply(reps,function (x,k) x$neiDstRefs[,k], k)
      tmp = mkMerge(kIds,kDst,rowNR,length(reps),wts)
      idsR[[k]] = tmp$ids
      dstR[[k]] = tmp$dst
    }
  }
  
  if (!noTrgs) 
  {
    idsT = do.call(cbind,idsT)  
    rownames(idsT) = rowNT
    colnames(idsT) = paste0("Id.k",1:mxk)
    dstT = do.call(cbind,dstT)
    rownames(dstT) = rowNT
    colnames(dstT) = paste0("Dst.k",1:mxk)
  }    

  if (!noRefs) 
  {
    idsR = do.call(cbind,idsR)  
    rownames(idsR) = rowNR
    colnames(idsR) = paste0("Id.k",1:mxk)
    dstR = do.call(cbind,dstR)
    rownames(dstR) = rowNR
    colnames(dstR) = paste0("Dst.k",1:mxk)
  } 
  
  # build a merged list of yRefs, and xRefs
  # find out if all the column names are the same.
  
  clsI = TRUE
  for (i in 1:(length(reps)-1)) 
  {
    if (!identical(colnames(reps[[i]]$yRefs),colnames(reps[[i+1]]$yRefs)) ||
        !identical(colnames(reps[[i]]$xRefs),colnames(reps[[i+1]]$xRefs)))
    {
      clsI = FALSE
      break
    }
  }  
  if (clsI) 
  {
    idx = if (is.null(rowNR)) NULL else na.omit(match(rowNR,rownames(reps[[1]]$xRefs))) 
    yRefs = if (is.null(idx)) reps[[1]]$yRefs else reps[[1]]$yRefs[idx,]
    xRefs = if (is.null(idx)) reps[[1]]$xRefs else reps[[1]]$xRefs[idx,]
    for (i in 2:length(reps))
    {
      rowNR = setdiff(rowNR,rownames(xRefs))
      if (length(rowNR) == 0) break
      idx = na.omit(match(rowNR,rownames(reps[[i]]$xRefs))) 
      if (length(idx) > 0) 
      {
        yRefs = rbind(yRefs,reps[[i]]$yRefs[idx,,drop=FALSE])
        xRefs = rbind(xRefs,reps[[i]]$xRefs[idx,,drop=FALSE])
      }
    }
  } else {
    yRefs = reps[[1]]$yRefs
    xRefs = reps[[1]]$xRefs
    for (i in 2:length(reps))
    {
      if (!identical(yRefs,reps[[i]]$yRefs)) yRefs = unionDataJoin(yRefs,reps[[i]]$yRefs,warn=FALSE)
      if (!identical(xRefs,reps[[i]]$xRefs)) xRefs = unionDataJoin(xRefs,reps[[i]]$xRefs,warn=FALSE)
    }
  }
    
  
  # build a merged list of xall
  xall = reps[[1]]$xall
  for (i in 2:length(reps))
  {
    if (!identical(xall,reps[[i]]$xall)) xall = unionDataJoin(xall,reps[[i]]$xall,warn=FALSE)
  }    

  mymean = function(x)
  {
     if (is.null(ncol(x)))
     {
        ans = if (is.factor(x)) NA else mean(x)
     }
     else
     {
        ans=as.numeric(rep(NA,ncol(x)))
        names(ans)=colnames(x)
        for (i in 1:ncol(x)) if (!is.factor(x[,i])) ans[i]=mean(x[,i])
     }
     ans
  }
  mysd = function(x)
  {
     if (is.null(ncol(x)))
     {
        ans = if (is.factor(x)) NA else sd(x)
     }
     else
     {
        ans=as.numeric(rep(NA,ncol(x)))
        names(ans)=colnames(x)
        for (i in 1:ncol(x)) if (!is.factor(x[,i])) ans[i]=sd(x[,i])
     }
     ans
  }

  xScale=list(center=mymean(xRefs),scale=mysd(xRefs))
  yScale=list(center=mymean(yRefs),scale=mysd(yRefs))
      
  out=list(call=cl,yRefs=yRefs,xRefs=xRefs,obsDropped=NULL,yDrop=NULL,bootstrap=FALSE,
           xDrop=NULL,trgRows=rowNT,xall=xall,cancor=NULL,theFormula=NULL,
           ftest=NULL,yScale=yScale,xScale=xScale,ccaVegan=NULL,ranForest=NULL,
           ICA=NULL,k=mxk,projector=NULL,nVec=NULL,pVal=NULL,method="consensus",ann=FALSE,
           neiDstTrgs=dstT,neiIdsTrgs=idsT,
           neiDstRefs=dstR,neiIdsRefs=idsR)

  class(out)="yai"
  out
}

