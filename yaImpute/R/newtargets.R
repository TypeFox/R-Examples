# Arguments:
#   object is of class "yai" created by function yai.
#   newdata is a matrix or dataframe that contains data for x variables
#   k is the number of neighbors desired (see function yai).
#
# Value:  A list of class "yai", that is a copy of the input object with
#         the following members replaced:
#
#   call: the call
#
#   obsDropped: a list of the rownames for observations dropped for various
#               reasons (missing data).
#
#   trgRows: a list of the rownames for target observations as a subset
#            of observations in xall (normally, all rows).
#
#   xall: the x variables for all observations.
#
#   neiDstTrgs: A data frame of distances between a target (identified by it's
#      rowname) and the k references. There are k columns.
#
#   neiIdsTrgs: A data frame of reference identifications that correspond to
#       neiDstTrgs.
#
#   k: new value of k, if null, the value is taken from object.
#
#   ann: use ann or not...if null, the value is taken from object.
#

newtargets=function(object,newdata,k=NULL,ann=NULL)
{
   if (class(object) != "yai") stop ("object must be class yai")
   if (object$method == "ensemble") 
     stop ("newtargets can not be found for objects with method 'ensemble'.")
   if (is.null(newdata) | nrow(newdata)==0) stop ("newdata is required")
   if (object$method == "gnn") 
      if (!requireNamespace ("vegan")) stop("install vegan and try again")
   if (object$method == "randomForest") 
      if (!requireNamespace ("randomForest")) 
        stop("install randomForest and try again")

   sumSqDiff=function(x,y) { d=x-y; sum(d*d) }
   factorMatch = get("factorMatch",asNamespace("yaImpute"))

   if (is.null(ann)) ann=object$ann
   if (!is.null(k)) object$k=k

   object$call=match.call()

   obsDropped=NULL

   # don't redo the factor matching for objects that come already done.
   if (is.null(attr(newdata,"illegalLevelCounts")) &&
       length(intersect(xvars(object),names(object$xlevels))) > 0)
   {
      newdata = factorMatch(newdata,object$xlevels)
      if (is.list(attr(newdata,"illegalLevelCounts")))
      {
         warning ("NA's generated due to illegal level(s).")
         cat ("Illegal levels\n")
         print(attr(newdata,"illegalLevelCounts"))
      }
   }
   if (is.null(object$theFormula))
   {
      have=intersect(colnames(object$xRefs),colnames(newdata))
      if (length(have) != length(colnames(object$xRefs)))
      {
          missing = setdiff(colnames(object$xRefs),colnames(newdata))     
          stop(paste("required column(s) missing:", paste (missing, collapse=", ")))
      }
      xall=na.omit(newdata[,have])
      obsDropped=names(attributes(na.omit(xall))$na.action)
      if (length(obsDropped)>0) warning (nrow(newdata)-nrow(xall)," observation(s) removed")
   }
   else
   {
      xall=model.frame(object$theFormula$x,newdata)
      if (!is.null(object$xDrop)) xall=xall[,!object$xDrop,drop=FALSE]
      obsDropped=setdiff(rownames(newdata),rownames(xall))
      if (length(obsDropped)) warning (length(obsDropped)," observation(s) removed")
   }
   if (nrow(xall) == 0) stop ("no observations")

   trgs=setdiff(rownames(xall),rownames(object$xRefs))
   if (nrow(xall) != length(trgs))
   {
      obsDropped=union(obsDropped,intersect(rownames(object$xRefs),rownames(xall)))
      warning (nrow(xall)-length(trgs)," row(s) in newdata are original references and ignored")
   }

   theCols = colnames(object$xRefs)  # may be changed for reduced rank, depending on method.

   if (object$method %in% c("msn","msn2","msnPP","mahalanobis","ica"))
   {
      theCols = rownames(object$projector)
      xcvRefs=scale(object$xRefs,center=object$xScale$center,scale=object$xScale$scale)
      if (length(theCols)<ncol(xcvRefs)) xcvRefs=xcvRefs[,theCols,drop=FALSE]
   }

   xTrgs=as.data.frame(xall[trgs,theCols,drop=FALSE]) # this is needed by randomForest and gnn unscalled.
   if (nrow(xTrgs)==0) stop("no observations")
   if (object$method == "gnn") # gnn
   {
      # create a projected space for the reference observations
      xcvRefs=predict(object$ccaVegan,type="lc",rank="full")
      xcvRefs=xcvRefs %*% diag(sqrt(object$ccaVegan$CCA$eig/sum(object$ccaVegan$CCA$eig)))

      # create a projected space for the unknowns (target observations)
      xcvTrgs=scale(xTrgs,center=object$xScale$center,scale=object$xScale$scale)
      xcvTrgs=predict(object$ccaVegan,newdata=as.data.frame(xcvTrgs),type="lc",rank="full")
      xcvTrgs=xcvTrgs %*% diag(sqrt(object$ccaVegan$CCA$eig/sum(object$ccaVegan$CCA$eig)))
      nVec = ncol(xcvRefs)
   }
   else if (object$method == "randomForest") # randomForest
   {
      nodes=NULL
      predObs = if (is.null(attr(object$ranForest,"rfRefNodeSort"))) rbind(object$xRefs,xTrgs) else xTrgs
      for (i in 1:length(object$ranForest))
      {
         nodeset=attr(predict(object$ranForest[[i]],predObs,
                      proximity=FALSE,nodes=TRUE),"nodes")
         if (is.null(nodeset)) stop("randomForest did not return nodes")
         colnames(nodeset)=paste(colnames(nodeset),i,sep=".")
         nodes=if (is.null(nodes)) nodeset else cbind(nodes,nodeset)
      }
      if (is.null(attr(object$ranForest,"rfRefNodeSort")))
      {
        INTrefNodes=as.integer(nodes[rownames(object$xRefs),])
        INTnrow=as.integer(nrow(object$xRefs))
        INTncol=as.integer(ncol(nodes))
        INTsort = INTrefNodes
        dim(INTsort) = c(INTnrow,INTncol)
        INTsort=apply(INTsort,2,function (x) sort(x,index.return = TRUE, decreasing = FALSE)$ix-1)
        attributes(INTsort)=NULL
        INTsort = as.integer(INTsort)
        nodes = nodes[rownames(xTrgs),]
      }
      else
      { 
        INTrefNodes = attr(object$ranForest,"rfRefNodeSort")[["INTrefNodes"]]
        INTnrow     = attr(object$ranForest,"rfRefNodeSort")[["INTnrow"]]
        INTncol     = attr(object$ranForest,"rfRefNodeSort")[["INTncol"]]
        INTsort     = attr(object$ranForest,"rfRefNodeSort")[["INTsort"]]
      }
   }  
   else if (object$method == "random")
   { 
      xcvRefs=data.frame(random=runif(nrow(object$xRefs)),row.names=rownames(object$xRefs))
      xcvTrgs=data.frame(random=runif(length(trgs)),row.names=trgs)
   }
   else if (object$method %in% c("msn","msn2","msnPP","mahalanobis","ica"))
   {  
      xcvRefs=as.matrix(xcvRefs[,theCols,drop=FALSE]) %*% object$projector
      xcvTrgs=scale(xTrgs,center=object$xScale$center,scale=object$xScale$scale)
      xcvTrgs=as.matrix(xcvTrgs[,theCols,drop=FALSE]) %*% object$projector
   }
   else if (object$method == "euclidean") 
   {      
      xcvRefs=scale(object$xRefs,center=object$xScale$center,scale=object$xScale$scale)
      xcvRefs=as.matrix(xcvRefs[,theCols,drop=FALSE])
      xcvTrgs=scale(xTrgs,center=object$xScale$center,scale=object$xScale$scale)         
      xcvTrgs=as.matrix(xcvTrgs[,theCols,drop=FALSE]) 
   }
   else # method is raw
   {
      xcvRefs=as.matrix(object$xRefs[,theCols,drop=FALSE])
      xcvTrgs=as.matrix(xTrgs[,theCols,drop=FALSE])   
   }

   neiDstTrgs=matrix(data=NA,nrow=length(trgs),ncol=object$k)
   rownames(neiDstTrgs)=trgs
   colnames(neiDstTrgs)=paste("Dst.k",1:object$k,sep="")
   neiIdsTrgs=neiDstTrgs
   colnames(neiIdsTrgs)=paste("Id.k",1:object$k,sep="")

   if (object$method %in%  c("msn","msn2","msnPP","mahalanobis","ica","euclidean","gnn","raw"))
   {
      if (ann & nrow(xcvTrgs)>0)
      {
          k=object$k
          ann.out=ann(xcvRefs, xcvTrgs, k, verbose=FALSE)$knnIndexDist
          neiDstTrgs[TRUE]=sqrt(ann.out[,(k+1):ncol(ann.out)])
          for (i in 1:k)
             neiIdsTrgs[,i]=rownames(xcvRefs)[ann.out[,i]]
          rownames(neiDstTrgs)=rownames(neiIdsTrgs)
      }
      else
      {
         for (row in rownames(xcvTrgs))
         {
            d=sqrt(sort(apply(xcvRefs,MARGIN=1,sumSqDiff,xcvTrgs[row,])))[1:object$k]
            neiDstTrgs[row,]=d
            neiIdsTrgs[row,]=names(d)
         }
      }
   }
   else if (object$method == "randomForest")
   {
     prox=lapply(apply(nodes,1,as.list),function (x) {
           prx=.Call("rfoneprox", INTrefNodes, INTsort, INTnrow, INTncol,
                     as.integer(x), vector("integer",INTnrow),dup=FALSE) 
           if (object$k > 1) px=sort(prx,index.return = TRUE, decreasing = TRUE)$ix[1:object$k]
           else              px=which.max(prx)
           c(prx[px],px)  # counts followed by pointers to references
           })
     for (i in 1:object$k)
     {
       neiDstTrgs[,i]=unlist(lapply(prox,function (x,i) (INTncol-x[i])/INTncol,i))
       neiIdsTrgs[,i]=unlist(lapply(prox,function (x,i,k,Rnames) 
              Rnames[x[k+i]],i,object$k,rownames(object$xRefs)))
     } 
   }
   else if (object$method == "random")
   {
      l=k+1
      d = matrix(unlist(lapply(xcvTrgs[[1]],function (x, xcv, l) 
            {
              sort((xcv-x)^2,index.return=TRUE)$ix[2:l]
            },xcvRefs[[1]],l)),nrow=nrow(xcvTrgs),ncol=k,byrow=TRUE)
      for (ic in 1:ncol(d))
      {
        neiDstTrgs[,ic]=abs(xcvTrgs[,1]-xcvRefs[d[,ic],1])
        neiIdsTrgs[,ic]=rownames(xcvRefs)[d[,ic]]
      }
   }
   else  # default
   {
      stop("no code for specified method")   
   }

   # if bootstrap, then modify the reference ID's in the result ID tables. 
   if (length(object$bootstrap) > 1) 
      neiIdsTrgs[] = sub("\\.[0-9]$","",neiIdsTrgs[])
   
   object$obsDropped=obsDropped
   object$trgRows=trgs
   addX = setdiff (rownames(object$xRefs),rownames(xall))
   if (length(addX) > 0) xall = rbind(xall,object$xRefs[addX,])
   object$xall=xall
   object$neiDstTrgs=neiDstTrgs
   object$neiIdsTrgs=neiIdsTrgs
   noRefs=TRUE
   object$neiDstRefs=NULL
   object$neiIdsRefs=NULL
   object$ann=ann
   object
}
