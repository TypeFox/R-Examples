errorStats <- function(mahal,...,scale=FALSE,pzero=0.1,plg=0.5,seeMethod="lm")
{
   obsMinusImp = function (object,...,vars=NULL)
   {
      if (missing(object)) stop ("object required.")
      if (class(object)[1] == "yai") object = impute(object,vars=vars,observed=TRUE,...)
      if (is.null(object)) stop ("no imputations found using this object")
      object=na.omit(object)
      if (is.null(vars)) vars=names(object)
      vi=paste(unique(strsplit(vars,".o",fixed=TRUE)))
      vi=intersect(vi,names(object))
      notFound=setdiff(vars,names(object))
      if (length(notFound)>0) warning ("variables not found: ",paste(notFound,collapse=", "))
      if (length(vi) == 0) stop("nothing to compute")
      vo=paste(vi,"o",sep=".")
      notFound=setdiff(vo,names(object))
      if (length(notFound)>0) warning ("variables not found: ",paste(notFound,collapse=", "))
      vo=intersect(vo,names(object))
      both=intersect(paste(unique(strsplit(vo,".o",fixed=TRUE))),vi)
      if (length(both) == 0) stop("nothing to compute")
      vo=paste(both,"o",sep=".")
      obsMinusImp=matrix(NA,nrow(object),length(both))
      colnames(obsMinusImp)=both
      rownames(obsMinusImp)=rownames(object)
      for (i in 1:length(both))
      {
         if (!is.factor(object[,both[i]])) obsMinusImp[,i]=object[,both[i]]-object[,vo[i]]
      }
      obsMinusImp
   }

   frmmsd0 = function(x,label,pzero)
   {
      if (missing(x)) stop ("x required.")
      if (class(x) != "yai") stop ("class must be yai")
      if (x$method != "mahalanobis") stop ("method must be mahalanobis")
      if (is.null(x$neiDstRef)) stop ("reference neighbors must be present")
      xr=obsMinusImp(x)^2
      if (ncol(xr) == 0) stop ("nothing to compute")
      xdstt=notablyDistant(x,p=1-pzero)$threshold
      obs=x$neiDstRef[,1]<xdstt
      if (sum(obs)<20) warning ("number of observations used is ",sum(obs)," -- consider increasing pzero")
      rmmsd0=matrix(NA,ncol(xr),1)
      rownames(rmmsd0)=colnames(xr)
      colnames(rmmsd0)=paste(label,"rmmsd0",sep=".")
      for (var in rownames(rmmsd0)) rmmsd0[var,1]=mean(xr[obs,var])
      rmmsd0[rmmsd0<=0] = NA
      sqrt(rmmsd0)
   }
   frmsdlg = function(x,label,plg)
   {
      if (missing(x)) stop ("x required.")
      if (class(x) != "yai") stop ("class must be yai")
      if (is.null(x$neiDstRef)) stop ("reference neighbors must be present")
      xr=obsMinusImp(x)^2
      if (ncol(xr) == 0) stop ("nothing to compute")
      xdstt=notablyDistant(x,p=plg)$threshold
      obs=x$neiDstRef[,1]>xdstt
      if (sum(obs)<20) warning ("number of observations used is ",sum(obs)," -- consider increasing plg")
      rmsdlg=matrix(NA,ncol(xr),1)
      rownames(rmsdlg)=colnames(xr)
      colnames(rmsdlg)=paste(label,"rmsdlg",sep=".")
      for (var in rownames(rmsdlg)) rmsdlg[var,1]=mean(xr[obs,var])
      rmsdlg[rmsdlg<0] = NA
      sqrt(rmsdlg)
   }
   frmsd2 = function(x,label)
   {
      rmsd=rmsd.yai(x,scale=FALSE)
      colnames(rmsd)=paste(label,"rmsd",sep=".")
      rmsd
   }
   fsee = function(x,label,method)
   {
      if (missing(x)) stop ("x required.")
      if (class(x) != "yai") stop ("class must be yai")
      mfg = paste(xvars(x),collapse=",5)+s(")
      mfl = paste(xvars(x),collapse="+")
      data=cbind(x$yRefs,x$xRefs)
      see=matrix(NA,ncol(x$yRefs),1)
      rownames(see)=yvars(x)
      colnames(see)=paste(label,"see",sep=".")
      for (var in rownames(see))
      {
         if (method=="gam")
         {
            if (!requireNamespace(gam))
            {
              stop ("install package gam and try again") 
              # the purpose of this line of code is to suppress CRAN check notes
              gam <- function (...) NULL
            }
            g = try(gam(as.formula(paste(var,"~s(",mfg,")",sep="")),data=data))
         }
         else g = try( lm(as.formula(paste(var,"~",mfl,sep="")),data=data))
         if (is.element("lm",class(g))) see[var,1]=mean(resid(g)^2)
     }
     see[see<0] = NA
     sqrt(see)
   }
   fsei = function(rmsd,rmmsd0)
   {
      sei = rmsd^2 - (rmmsd0^2/2)
      sei[sei<0] = NA
      colnames(sei)=paste(strsplit(colnames(rmsd),".",fixed=TRUE)[[1]][1],"sei",sep=".")
      sqrt(sei)
   }
   fdstc = function(rmsd,rmmsd0)
   {
      dstc =rmsd^2 - rmmsd0^2
      dstc[dstc<0] = NA
      colnames(dstc)=paste(strsplit(colnames(rmsd),".",fixed=TRUE)[[1]][1],"dstc",sep=".")
      sqrt(dstc)
   }
   fmlf = function(see,rmmsd0)
   {
      mlf = see^2 - (rmmsd0^2/2)
      mlf[mlf<0] = NA
      colnames(mlf)=paste(strsplit(colnames(see),".",fixed=TRUE)[[1]][1],"mlf",sep=".")
      sqrt(mlf)
   }

##############

#  compute basic error statistics for the mahalanobis method

   if (mahal$method != "mahalanobis") stop ("method for first argument must be mahalanobis")

   label = deparse(substitute(mahal))
   rmmsd0= frmmsd0(mahal,label,pzero)
   rmsdlg= frmsdlg(mahal,label,plg)
   rmsd  = frmsd2(mahal,label)
   see   = fsee(mahal,label,seeMethod)
   sei   = fsei(rmsd,rmmsd0)
   dstc  = fdstc(rmsd,rmmsd0)
   mlf   = fmlf(see,rmmsd0)
   out   = cbind(see,rmmsd0,mlf,rmsd,rmsdlg,sei,dstc)
   if (scale)
   {
      scl = 1/mahal$yScale$scale[rownames(rmsd)]
      out = out*scl
   }

#  if there are other objects, compute error statistics for them...

   args = list(...)
   if (length(args)==0)
   {
      ans = vector("list",2)
      names(ans) = c("common",label)
      ans[[1]] = out[,1:3]
      ans[[2]] = out[,4:7]
   }
   else
   {
      names(args) = as.list(substitute(list(...)))[-1]
      ans = vector("list",length(args)+2)
      names(ans) = c("common",label,names(args))
      ans[[1]] = out[,1:3]
      ans[[2]] = out[,4:7]
      i = 2
      for (object in args)
      {
         i = i+1
         label = names(ans)[i]
         rmsd  = frmsd2(object,label)[rownames(rmmsd0),,drop=FALSE]
         rmsdlg= frmsdlg(object,label,plg=plg)[rownames(rmmsd0),,drop=FALSE]
         sei   = fsei(rmsd,rmmsd0)
         dstc  = fdstc(rmsd,rmmsd0)
         out   = cbind(rmsd,rmsdlg,sei,dstc)
         if (scale)
         {
            scl = 1/object$yScale$scale[rownames(rmsd)]
            out = out*scl
         }
         ans[[i]] <- out
      }
   }
   ans
}
