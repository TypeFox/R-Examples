yai <- function(x=NULL,y=NULL,data=NULL,k=1,noTrgs=FALSE,noRefs=FALSE,
                nVec=NULL,pVal=.05,method="msn",ann=TRUE,mtry=NULL,ntree=500,
                rfMode="buildClasses",bootstrap=FALSE,ppControl=NULL,
                sampleVars=NULL,rfXsubsets=NULL)
{
   # define functions used internally.
   sumSqDiff=function(x,y) { d=x-y; sum(d*d) }

   findFactors =  get("findFactors",asNamespace("yaImpute"))

   ftest.cor = function (p,q,N,cor)
   {  
      s=min(p,q)
      if (s==0) stop ("p and q must be > 0")
      if (length(cor) < s) stop ("cor too short")
      lamda=array(dim=s)
      k=1:s
      for (i in k) lamda[i]=prod(1-cor[i:s]^2)
      r=(N-s-1)-((abs(p-q)+1)/2)
      Ndf=(p-k+1)*(q-k+1)
      u=(Ndf-2)/4
      xx=((p-k+1)^2+(q-k+1)^2)-5
      t=vector(mode="numeric",length=s)
      for (i in k) if (xx[i]>0) t[i]=sqrt(((p-k[i]+1)^2*(q-k[i]+1)^2-4)/xx[i])
      lamda.invt=lamda^(1/t)
      Ddf=(r*t)-(2*u)
      setNA = Ddf < 1 | Ndf < 1
      firstNA = which(setNA == TRUE)
      if (length(firstNA) == 0) firstNA=0
      if (length(firstNA) > 1)  firstNA=firstNA[1]
      if (firstNA > 0) setNA[firstNA:length(setNA)] = TRUE
      F=((1-lamda.invt)/lamda.invt)*(Ddf/Ndf)      
      F[setNA] = NA
      pgF = F
      firstNA = if (firstNA > 1) firstNA else length(F)
      {
        pgF[1:firstNA]=pf(F[1:firstNA],Ndf[1:firstNA],Ddf[1:firstNA],
                          lower.tail=FALSE)
        pgF[setNA] = NA
      } 
      list(F=F,pgF=pgF)
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

   #===============================================
   # ARGUMENT and DATA screening

   methodSet=c("msn","msn2","msnPP","mahalanobis","ica","euclidean","gnn",
               "randomForest","raw","random")
               
   if (!(method %in% methodSet))
      stop (paste("method not one of:",paste(methodSet,collapse =", ")))
      
   if (method == "gnn") # (GNN), make sure we have package vegan loaded
   {
      if (!requireNamespace ("vegan")) 
      {
        stop("install vegan and try again")
        # the purpose of this line of code is to suppress CRAN check notes
        cca <- rda <- function (...) NULL
      } else {
        cca <- vegan::cca
        rda <- vegan::rda
      }
   }
   if (method == "ica") # (ica), make sure we have package fastICA loaded
   {
      if (!requireNamespace ("fastICA")) 
      {
        stop("install fastica and try again")
        # the purpose of this line of code is to suppress CRAN check notes
        fastICA <- function (...) NULL        
      } else {
        fastICA <- fastICA::fastICA
      }

   }
   if (method == "randomForest") # make sure we have package randomForest loaded
   {
      if (!requireNamespace ("randomForest")) 
      {
        stop("install randomForest and try again")
        # the purpose of this line of code is to suppress CRAN check notes
        randomForest <- function (...) NULL
      } else {
        randomForest <- randomForest::randomForest
      }

   }     
   if (method == "msnPP") # make sure we have package ccaPP loaded
   {
      if (!requireNamespace ("ccaPP")) 
      {
        stop("install ccaPP and try again")
        # the purpose of this line of code is to suppress CRAN check notes
        fastMAD <- ccaGrid <- ccaProj <- function (...) NULL
      } else {
        fastMAD <- ccaPP::fastMAD
        ccaGrid <- ccaPP::ccaGrid
        ccaProj <- ccaPP::ccaProj
      }

   }     

   cl=match.call()
   obsDropped=NULL
   theFormula=NULL
   yall=NULL
   if (is.data.frame(x) | is.matrix(x))
   {
      if (mode(rownames(x)) != "character") rownames(x)=as.character(rownames(x))
      xall=na.omit (as.data.frame(x))
      if (nrow(xall) != nrow(x))
      {
         warning (nrow(x)-nrow(xall)," x observation(s) removed")
         obsDropped=names(attributes(na.omit(x))$na.action)
      }
      if (!is.null(y))
      {
         if (is.null(dim(y)))
         {
           if (length(y) == nrow (x)) y=data.frame(y,row.names=rownames(x), 
                                                   stringsAsFactors = TRUE)
           else stop(paste0("when formulas are not used,",
                     " y must be a matrix or dataframe,",
                     " or a vector the same length of rows in x"))
         } 
         if (is.matrix(y) | is.data.frame(y))
         {
            if (mode(rownames(y)) != "character") rownames(y)=as.character(rownames(y))
            yall=na.omit(as.data.frame(y))
            if (nrow(yall) != nrow(as.data.frame(y)))
            {
               warning (nrow(y)-nrow(yall)," y observation(s) removed")
               obsDropped=union(obsDropped,names(attributes(na.omit(y))$na.action))
            }
         }
         theFormula=NULL
      }
   }
   else if (class(x) == "formula")
   {
      if (class(y) == "formula") yall=model.frame(y,data=data)
      xall=model.frame(x,data=data)
      obsDropped=setdiff(rownames(data),rownames(xall))
      if (length(obsDropped)) warning (length(obsDropped)," observation(s) removed")
      theFormula=list(x=x,y=y)
   }
   else stop ("x is missing or not a matrix nor dataframe")
   if (is.null(yall) & (method %in% c("mahalanobis","ica",
                        "euclidean","randomForest","raw")))
   {
      ydum=TRUE
      yall=data.frame(ydummy=rep(1,nrow(xall)),row.names=rownames(xall))
   }
   else ydum=FALSE   
  
   if (is.null(yall)) stop("y is missing")
   if (nrow(xall) == 0) stop ("no observations in x")
   if (! (method %in% c("random","randomForest")))
   {
      fy=0
      if (!(method %in% c("mahalanobis","ica","euclidean","raw"))) 
        fy=sum(findFactors(yall))
      if (fy+sum(findFactors(xall)>0)>0) 
        stop("factors allowed only for methods randomForest or random")
   }
   refs=intersect(rownames(yall),rownames(xall))
   if (length(refs) == 0) stop ("no reference observations.")

   # if X variable subsets are used, make sure we are using method="randomForest" 
   if (!is.null(rfXsubsets) && method != "randomForest") 
   {
     warning ("specification of rfXsubsets is ignored when method is not randomForest.")
     rfXsubsets = NULL
   }
   if (!is.null(rfXsubsets))
   {
     vtok = match(unique(unlist(rfXsubsets)),names(xall))
     if (any(is.na(vtok))) stop("one or more variables in rfXsubsets are not present in x.")
     xall = xall[,vtok]
   }

   # if sampling variables, set up xRefs and yRefs accordingly

   if (!is.null(sampleVars))
   {
     if (length(sampleVars) == 1 && is.null(names(sampleVars))) 
       sampleVars=rep(sampleVars,2)
     names(sampleVars) = if (is.null(names(sampleVars))) c("X","Y") else 
                             toupper(names(sampleVars))
     nx = match("X",names(sampleVars))
     ny = match("Y",names(sampleVars))
     nx = if (!is.na(nx)) sampleVars[nx] else 0
     ny = if (!is.na(ny)) sampleVars[ny] else 0
     if (nx > 0)
     {
       nx = if (nx < 1.) max(1, ncol(xall)*nx) else min(nx, ncol(xall))
       nxn = sample(1:ncol(xall),nx)
       xall = xall[,nxn,drop=FALSE]
     }
     if (ny > 0)
     {
       ny = if (ny < 1.) max(1, ncol(yall)*ny) else min(ny, ncol(yall))
       nyn = sample(1:ncol(yall),ny)
       yall = yall[,nyn,drop=FALSE]
     }
   }
   
   
   # if this is a bootstrap run, draw the sample.
   if (bootstrap)
   { 
     if (length (grep ("\\.[0-9]$",rownames(xall))) > 0) 
         stop ("rownames must not end with .[0-9] when bootstrap is true.")

     bootsamp <- sort(sample(x=refs, size=length(refs), replace=TRUE))
     yRefs=yall[bootsamp,,drop=FALSE]
     xRefs=xall[bootsamp,,drop=FALSE]
     refs = bootsamp
   } else {
     yRefs=yall[refs,,drop=FALSE]
     xRefs=xall[refs,,drop=FALSE]
   }
   trgs=setdiff(rownames(xall),refs)
   
  
   if (method == "gnn") # remove rows with zero sums or vegan will error off...
   {
      zero = apply(yRefs,1,sum) <= 0
      ndrop=sum(zero)
      if (ndrop>0)
      {
         warning (ndrop,paste0(" rows have y-variable row sums <= 0 were ",
           "converted to target observations for method gnn"))
         if (ndrop==nrow(yRefs)) stop ("all references were deleted")
         obsDropped=union(obsDropped,refs[zero])
         refs=refs[!zero]
         yRefs=yall[refs,,drop=FALSE]
         xRefs=xall[refs,,drop=FALSE]
         trgs=setdiff(rownames(xall),refs)
      }
   }

   if (method == "gnn") # remove columns with zero sums.
   {
      yDrop=apply(yRefs,2,sum) <= 0
      if (sum(yDrop) > 0) warning ("y variables with zero sums: ",
                                    paste(colnames(yRefs)[yDrop],collapse=","))
      if (sum(yDrop) == ncol(yRefs)) stop("no y variables")
      if (sum(yDrop) > 0) yRefs=yRefs[,!yDrop,drop=FALSE]
   }

   # initial scale values (maybe reduced by some methods).
   if (method != "msnPP")
   {
     xScale=list(center=mymean(xRefs),scale=mysd(xRefs))
     yScale=list(center=mymean(yRefs),scale=mysd(yRefs))
   } else {
     msn3cntr = function (x) 
     {
       if (is.factor(x)) return (list(NULL,NULL))
       cM = fastMAD(x)  # uses fastMAD for scaling.
       if (cM$MAD == 0)
       {
         cM$MAD = sd(x)
         cM$center = mean(x)
       }
       cM
     }
     ce=matrix(unlist(apply(xRefs,2,msn3cntr)),ncol=2,byrow=TRUE)
     rownames(ce) = colnames(xRefs)
     xScale=list(center=ce[,1],scale=ce[,2])
     ce=matrix(unlist(apply(yRefs,2,msn3cntr)),ncol=2,byrow=TRUE)
     rownames(ce) = colnames(yRefs)
     yScale=list(center=ce[,1],scale=ce[,2])
   }
   # for all methods except randomForest, random, and raw, 
   # variables with zero variance are dropped.
   if (!(method %in% c("randomForest","random","raw")))
   {
      xDrop=xScale$scale < 1e-10
      if (sum(xDrop) > 0) warning ("x variables with zero variance: ",
                                    paste(colnames(xRefs)[xDrop],collapse=","))
      if (sum(xDrop) == ncol(xRefs)) stop("no x variables")
      if (sum(xDrop) > 0)
      {
         xRefs=xRefs[,!xDrop,drop=FALSE]
         xScale$scale=xScale$scale[!xDrop]
         xScale$center=xScale$center[!xDrop]
      }
   }
   else xDrop=NULL

   # for this method, xRefs must be a matrix.
   if (method != "randomForest" && !is.matrix(xRefs)) xRefs=as.matrix(xRefs)

   # define these elements as NULL, some will be redefined below.

   cancor=NULL
   ftest=NULL
   yDrop=NULL
   projector=NULL
   ccaVegan=NULL
   ranForest=NULL
   xTrgs=NULL
   xcvTrgs=NULL
   ICA=NULL

   #======= Define projector (if used), scale the variables, and project the
   # reference space. Also project the target space if it is being used.

   if (method %in% c("msn","msn2","msnPP")) # msn (all kinds)
   {
      yDrop=yScale$scale < 1e-10
      if (sum(yDrop) > 0) warning ("y variables with zero variance: ",
                                    paste(colnames(yRefs)[yDrop],collapse=","))
      if (sum(yDrop) == ncol(yRefs)) stop("no y variables")
      if (sum(yDrop) > 0)
      {
         yRefs=yRefs[,!yDrop,drop=FALSE]
         yScale$scale=yScale$scale[!yDrop]
         yScale$center=yScale$center[!yDrop]
      }
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      ycvRefs=scale(yRefs,center=yScale$center,scale=yScale$scale)

      if (method %in% c("msn","msn2")) # msn and msn2 
      {
        cancor=cancor(xcvRefs,ycvRefs,xcenter = FALSE, ycenter = FALSE)                    
        theCols = rownames(cancor$xcoef)

        # scale the coefficients so that the cononical vectors will have unit variance.
        cscal = 1/apply(xcvRefs[,theCols,drop=FALSE] %*% 
                        cancor$xcoef[,1,drop=FALSE],2,sd)
        cancor$ycoef = cancor$ycoef * cscal
        cancor$xcoef = cancor$xcoef * cscal
      } else {                         # msnPP
        meth="spearman"
        ppfunc=ccaGrid
        if (!is.null(ppControl))
        {
          if (is.null(names(ppControl))) stop ("ppControl must have named strings.")
          for (ppn in names(ppControl))
          {
            if (ppn == "method") meth = ppControl[[ppn]]
            else if (ppn == "search") ppfunc = 
              if (ppControl[[ppn]] == "data" || 
                  ppControl[[ppn]] == "proj") ccaProj else ccaGrid
            else stop ("ppControl named element ",ppn," is invalid")
          }
        }
        # solve the canoncial correlation analysis via projection pursuit 
        cancor=ppfunc(xcvRefs,ycvRefs,method=meth,fallback=TRUE,
                       k=min(ncol(xcvRefs),ncol(ycvRefs)),nVec)
        # save the results using names and attributes that correspond 
        # to the cancor results
        cancor$ycoef = cancor$B
        rownames(cancor$ycoef) = colnames(ycvRefs)
        cancor$xcoef = cancor$A
        rownames(cancor$xcoef) = colnames(xcvRefs)
        theCols = rownames(cancor$xcoef)
        cancor$A = NULL
        cancor$B = NULL
        class(cancor) = "list"
      }        
      ftest=ftest.cor(p=nrow(cancor$ycoef),q=nrow(cancor$xcoef),
                      N=nrow(yRefs),cancor$cor)
      if (is.null(nVec))
      {
        fcheck = ftest$pgF[!is.na(ftest$pgF)]
        if (length(fcheck)> 0) nVec=max(1,length(fcheck)-sum(fcheck>pVal))
      }
      if (is.null(nVec)) nVec=1
      nVec=min(nVec,length(cancor$cor))
      nVec=max(nVec,1)
      if (method %in% c("msn","msnPP")) projector = 
        cancor$xcoef[,1:nVec,drop=FALSE] %*% 
          diag(cancor$cor[1:nVec,drop=FALSE],nVec,nVec)
                                        
      if (method == "msn2") 
      {
        if (any (1/(1-cancor$cor[1:nVec,drop=FALSE]^2) < 
            .Machine$double.eps*10000)) nVec=1
        if (any (1/(1-cancor$cor[1:nVec,drop=FALSE]^2) < 
            .Machine$double.eps*10000)) 
          stop("msn2 can not be run likely because there are too few obesrvations.")
        projector = cancor$xcoef[,1:nVec,drop=FALSE] %*%
                         diag(cancor$cor[1:nVec,drop=FALSE],nVec,nVec) %*%
                         diag(sqrt(1/(1-cancor$cor[1:nVec,drop=FALSE]^2)),nVec,nVec)
        if (any (projector == -Inf | projector == Inf | is.na(projector) | 
          is.nan(projector))) stop ("msn2 can not be run.")   
      }                                        
      if (length(theCols)<ncol(xRefs))
      {
         #just get the names and create a logical
         if (is.null(xDrop)) xDrop=xScale$center==0
         remove=setdiff(colnames(xRefs),theCols)
         xDrop[remove]=TRUE
         warning ("x variables with colinearity: ",paste(remove,collapse=","))
         xRefs=xRefs[,theCols,drop=FALSE]
         xScale$center=xScale$center[theCols]
         xScale$scale=xScale$scale[theCols]
         xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      }
      xcvRefs=xcvRefs %*% projector
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,theCols,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
         xcvTrgs=xcvTrgs %*% projector
      }
   }
   else if (method == "mahalanobis")
   {
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      qr = qr(xcvRefs)  # maybe we are not at full rank
      xcvRefs=xcvRefs[,qr$pivot[1:qr$rank],drop=FALSE]
      projector = solve(chol(cov(xcvRefs)))
      theCols = colnames(projector)
      if (length(theCols)<ncol(xRefs))
      {
         #just get the names and create a logical
         if (is.null(xDrop)) xDrop=xScale$center==0 
         remove=setdiff(colnames(xRefs),theCols)
         xDrop[remove]=TRUE
         warning ("x variables with colinearity: ",paste(remove,collapse=","))
         xRefs=xRefs[,theCols,drop=FALSE]
         xScale$center=xScale$center[theCols]
         xScale$scale=xScale$scale[theCols]
      }
      nVec = ncol(projector)  # same as qr$rank
      xcvRefs=xcvRefs %*% projector
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,theCols,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
         xcvTrgs=xcvTrgs %*% projector
      }
   }
   else if (method == "ica")
   {
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      qr = qr(xcvRefs)  # maybe we are not at full rank
      xcvRefs=xcvRefs[,qr$pivot[1:qr$rank],drop=FALSE]
      a=fastICA(xcvRefs,ncol(xcvRefs),method="C",)
      ICA=list(S=a$S,K=a$K,A=a$A,W=a$W)
      projector = a$K %*% a$W

      colnames(projector)=colnames(xcvRefs)
      rownames(projector)=colnames(xcvRefs)
      theCols = colnames(xcvRefs)
      if (length(theCols)<ncol(xRefs))
      {
         if (is.null(xDrop)) xDrop=xScale$center==0 
         remove=setdiff(colnames(xRefs),theCols)
         xDrop[remove]=TRUE
         warning ("x variables with colinearity: ",paste(remove,collapse=","))
         xRefs=xRefs[,theCols,drop=FALSE]
         xScale$center=xScale$center[theCols]
         xScale$scale=xScale$scale[theCols]
      }
      nVec = ncol(projector)  # same as qr$rank
      xcvRefs=xcvRefs %*% projector
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,theCols,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
         xcvTrgs=xcvTrgs %*% projector
      }
   }
   else if (method == "euclidean")
   {
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      nVec = ncol(xRefs)
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,!xDrop,drop=FALSE] 
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
      }
   }
   else if (method == "raw")
   {
      xcvRefs=xRefs
      nVec = ncol(xRefs)
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,,drop=FALSE]
         xcvTrgs=as.matrix(xTrgs)
      }
   }
   else if (method == "gnn") # GNN
   {
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      ccaVegan = cca(X=yRefs, Y=xcvRefs)
      if (is.null(ccaVegan$CCA) | 
          ccaVegan$CCA$rank == 0) 
      {
        warning (paste("cca() in package vegan failed, likely cause is",
               "too few X or Y variables.\nAttemping rda(),",
               "which is not well tested in the yaImpute package."))
        ccaVegan = rda(X=yRefs, Y=xcvRefs)
      }

      # create a projected space for the reference observations
      xcvRefs=predict(ccaVegan,type="lc",rank="full")
      xcvRefs=xcvRefs %*% diag(sqrt(ccaVegan$CCA$eig/sum(ccaVegan$CCA$eig)))

      # create a projected space for the unknowns (target observations)
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
         xcvTrgs=predict(ccaVegan,
                 newdata=as.data.frame(xcvTrgs),type="lc",rank="full")
         xcvTrgs=xcvTrgs %*% diag(sqrt(ccaVegan$CCA$eig/sum(ccaVegan$CCA$eig)))
      }
      nVec = ncol(xcvRefs)
   }
   else if (method == "randomForest")
   {  
      rfBuildClasses=NULL
      xTrgs=xall[trgs,1,drop=FALSE]
      rfVersion=packageDescription("randomForest")[["Version"]]  
      if (compareVersion(rfVersion,"4.5-22") < 0) 
          stop("Update your version of randomForest.")
      if (is.null(ntree)) ntree=500
      if (ydum)
      {
         if (!is.null(rfXsubsets)) warning(paste0("Specification of rfXsubsets",
           " ignored when unsupervised randomForest is run."))
         yone=NULL
         mt = if (is.null(mtry)) max(floor(sqrt(ncol(xRefs))),1) else 
                                 min(mtry, ncol(xRefs))
         ranForest=randomForest(x=xRefs,y=yone,proximity=FALSE,importance=TRUE,
                                keep.forest=TRUE,mtry=mt,ntree=ntree)
         ranForest$type="yaImputeUnsupervised"
         ranForest=list(unsupervised=ranForest)
      }
      else
      { 
         ranForest=vector("list",ncol(yRefs))
         if (length(ntree) < ncol(yRefs)) ntree=rep(max(50,
                                          floor(ntree/ncol(yRefs))),ncol(yRefs))
         for (i in 1:ncol(yRefs))
         {
            xN = names(xRefs)
            if (!is.null(rfXsubsets))
            {
              yV = names(yRefs)[i]
              if (!is.null(rfXsubsets[[yV]])) xN = intersect(rfXsubsets[[yV]],xN)
              if (length(xN)==0) stop ("rfXsubsets is empty for Y-variable ",yV)
            }
            yone=yRefs[,i]
            if (!is.factor(yone))
            { 
              if (is.null(rfBuildClasses) && rfMode=="buildClasses") 
                rfBuildClasses=TRUE
              if (is.null(rfBuildClasses))
              {
                  # if the version is prior to 19
                 if (compareVersion(rfVersion,"4.5-19") < 0) 
                 {
                   warning(paste0("yaImpute directly supports regression for ",
                     "continuous y's for randomForest version 4.5-19 and later."))
                   rfBuildClasses=TRUE
                 }
                 else rfBuildClasses=FALSE
              }
              if (rfBuildClasses)
              {
                yone=as.numeric(yone)
                breaks <- pretty(yone, n = min(20,nclass.Sturges(yone)),min.n = 1)
                div <- diff(breaks)[1]
                yone=as.factor(floor(yone/div))
              }
            }
            mt = if (is.null(mtry)) max(floor(sqrt(length(xN))), 1) else 
                                    min(mtry, length(xN))
            ranForest[[i]]=randomForest(x=xRefs[,xN,FALSE],
              y=yone,proximity=FALSE,importance=TRUE,keep.forest=TRUE,
              mtry=mt,ntree=ntree[i])
         }
         names(ranForest)=colnames(yRefs)
      }
      nodes=NULL
      for (i in 1:length(ranForest))
      {
         nodeset=attr(predict(ranForest[[i]],xall,
           proximity=FALSE,nodes=TRUE),"nodes")
         if (is.null(nodeset)) stop("randomForest did not return nodes")
         colnames(nodeset)=paste(colnames(nodeset),i,sep=".")
         nodes=if (is.null(nodes)) nodeset else cbind(nodes,nodeset)
      }
  
      if (bootstrap) 
      {
        rn = sub("\\.[0-9]$","",rownames(xRefs))
        refNodes = nodes[rn,]
        rownames(refNodes) = rownames(xRefs)
      } else refNodes = nodes[rownames(xRefs),]

      INTrefNodes=as.integer(refNodes)
      INTnrow=as.integer(nrow(xRefs))
      INTncol=as.integer(ncol(nodes))
      INTsort = INTrefNodes
      dim(INTsort) = c(INTnrow,INTncol)
      INTsort=apply(INTsort,2,function (x) sort(x,index.return = TRUE, 
        decreasing = FALSE)$ix-1)
      attributes(INTsort)=NULL
      INTsort = as.integer(INTsort)
      attr(ranForest,"rfRefNodeSort") = list(INTrefNodes=INTrefNodes, 
                         INTnrow=INTnrow, INTncol=INTncol, INTsort=INTsort)
   }
   else if (method == "random")
   { 
      nVec = 1
      ann=FALSE
      xcvRefs=data.frame(random=runif(nrow(xRefs)),row.names=rownames(xRefs))
      if (!noTrgs && length(trgs) > 0) xcvTrgs=
        data.frame(random=runif(length(trgs)),row.names=trgs)
   }
   else # default
   {
      stop("no code for specified method")
   }

   # if bootstrap, then modify the reference list essentually removing the
   # duplicate samples.  
   if (bootstrap) 
   {
     unq = unique(bootsamp)
     xcvRefs = xcvRefs[unq,,drop=FALSE]
     xRefs   = xRefs  [unq,,drop=FALSE]
   }

   k=min(k,nrow(xRefs))

   # ======= find neighbors for TARGETS
   if (noTrgs || length(trgs) == 0)
   {
      neiDstTrgs=NULL
      neiIdsTrgs=NULL
   }
   else
   {
      neiDstTrgs=matrix(data=NA,nrow=length(trgs),ncol=k)
      rownames(neiDstTrgs)=trgs
      colnames(neiDstTrgs)=paste("Dst.k",1:k,sep="")
      neiIdsTrgs=matrix(data="",nrow=length(trgs),ncol=k)
      rownames(neiIdsTrgs)=trgs
      colnames(neiIdsTrgs)=paste("Id.k",1:k,sep="")
      if (method %in%  c("msn","msn2","msnPP","mahalanobis",
                         "ica","euclidean","gnn","raw"))
      {
         if (ann)                                  
         { 
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
               d=sqrt(sort(apply(xcvRefs,MARGIN=1,sumSqDiff,xcvTrgs[row,]))[1:k])
               neiDstTrgs[row,]=d
               neiIdsTrgs[row,]=names(d)
            }
         }
      }
      else if (method == "random")
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
      else if (method == "randomForest")
      {
        prox=lapply(apply(nodes[rownames(xTrgs),,drop=FALSE],1,as.list),function (x) 
          {
             prx=.Call("rfoneprox", INTrefNodes, INTsort, INTnrow, INTncol,
                       as.integer(x), vector("integer",INTnrow),dup=FALSE) 
             if (k > 1)  px=sort(prx,index.return = TRUE, decreasing = TRUE)$ix[1:k]
             else        px=which.max(prx)
             c(prx[px],px)  # counts followed by pointers to references
          })
        for (i in 1:k)
        {
          neiDstTrgs[,i]=unlist(lapply(prox,function (x,i) (INTncol-x[i])/INTncol,i))
          neiIdsTrgs[,i]=unlist(lapply(prox,function (x,i,k,Rnames) 
                 Rnames[x[k+i]],i,k,rownames(xRefs)))
        } 
      }
      
      else # default
      {
         stop("no code for specified method")
      }
   }

   # ======= find neighbors for REFERENCES
   if (noRefs)
   {
      neiDstRefs=NULL
      neiIdsRefs=NULL
   }
   else
   {
      neiDstRefs=matrix(data=NA,nrow=nrow(xRefs),ncol=k)
      rownames(neiDstRefs)=rownames(xRefs)
      colnames(neiDstRefs)=paste("Dst.k",1:k,sep="")
      neiIdsRefs=matrix(data="",nrow=nrow(xRefs),ncol=k)
      rownames(neiIdsRefs)=rownames(xRefs)
      colnames(neiIdsRefs)=paste("Id.k",1:k,sep="")
      l=k+1
      if (method %in%  c("msn","msn2","msnPP","mahalanobis","ica","euclidean","gnn","raw"))
      {
         if (ann & nrow(xcvRefs)> 0)
         {
             ann.out=ann(xcvRefs, xcvRefs, l, verbose=FALSE)$knnIndexDist
             neiDstRefs[TRUE]=sqrt(ann.out[,(l+2):ncol(ann.out)])
             # check for a second neighbor being the reference itself (can happen 
             # if the first and second neighbors both have distance of zero).
             fix = ann.out[,1] != 1:nrow(ann.out)
             if (any(fix)) ann.out[fix,2] = ann.out[fix,1]             
             for (i in 2:l)
             {
                neiIdsRefs[,(i-1)]=rownames(xcvRefs)[ann.out[,i]]
             }
             rownames(neiDstRefs)=rownames(neiIdsRefs)
         }
         else
         {
            for (row in 1:nrow(xcvRefs))
            {
               d=sort(apply(xcvRefs,MARGIN=1,sumSqDiff,xcvRefs[row,])[-row])[1:k]
               neiDstRefs[row,]=d
               neiIdsRefs[row,]=names(d)
            }
         }
      }
      else if (method == "randomForest")
      {
        prox=lapply(apply(refNodes,1,as.list),function (x) 
          {
             prx=.Call("rfoneprox", INTrefNodes, INTsort, INTnrow, INTncol,
                       as.integer(x), vector("integer",INTnrow),dup=FALSE) 
             if (k > 1) px=sort(prx,index.return = TRUE, decreasing = TRUE)$ix[2:l]
             else
             { 
               px=which.max(prx)
               prx[px]=-1
               px=which.max(prx)
             }
             c(prx[px],px)  # counts followed by pointers to references
           })
        for (i in 1:k)
        {
           neiDstRefs[,i]=unlist(lapply(prox,function (x,i) (INTncol-x[i])/INTncol,i))
           neiIdsRefs[,i]=unlist(lapply(prox,function (x,i,k,Rnames) 
                  Rnames[x[k+i]],i,k,rownames(xRefs)))
        } 
      }
      else if (method == "random")
      {
         l=k+1
         d = matrix(unlist(lapply(xcvRefs[[1]],function (x, xcv, l) 
               {
                 sort((xcv-x)^2,index.return=TRUE)$ix[2:l]
               },xcvRefs[[1]],l)),nrow=nrow(xcvRefs),ncol=k,byrow=TRUE)
               
         for (ic in 1:ncol(d))
         {
           neiDstRefs[,ic]=abs(xcvRefs[,1]-xcvRefs[d[,ic],1])
           neiIdsRefs[,ic]=rownames(xcvRefs)[d[,ic]]
         }
      }
      else # default
      {
         stop("no code for specified method")
      }
   }

   xlevels=NULL
   fa=findFactors(xRefs)
   if (sum(fa)>0)
   {
      xlevels=vector(mode="list",length=sum(fa))
      kk=0
      for (i in 1:length(fa))
      {
         if (fa[i])
         {
            kk=kk+1
            xlevels[[kk]]=levels(xRefs[,i])
            names(xlevels)[[kk]]=names(xRefs)[i]
         }
      }
   }
   

   out=list(call=cl,yRefs=yRefs,xRefs=xRefs,obsDropped=obsDropped,yDrop=yDrop,
       bootstrap= if (bootstrap) bootsamp else bootstrap,
       xDrop=xDrop,trgRows=trgs,xall=xall,cancor=cancor,theFormula=theFormula,
       ftest=ftest,yScale=yScale,xScale=xScale,ccaVegan=ccaVegan,ranForest=ranForest,
       ICA=ICA,k=k,projector=projector,nVec=nVec,pVal=pVal,method=method,ann=ann,
       xlevels=xlevels,neiDstTrgs=neiDstTrgs,neiIdsTrgs=neiIdsTrgs,
       neiDstRefs=neiDstRefs,neiIdsRefs=neiIdsRefs,rfXsubsets=rfXsubsets)

   class(out)="yai"
   out
}
