DBGLM1 <-
function (y, m, groups, shrink.method = c("WEB","DEB"), contrast = c(1,2) , fdr.level = 0.05 , use.all.groups = TRUE )
 {
   if( !all(table(groups) > 1) ) stop( "Need at least 2 replicates per group" )

   yy = .DBSplitIntoGroups( y , groups )
   mm = .DBSplitIntoGroups( m , groups )
   m.offs <- sapply( mm[contrast] , rowMeans )
   Tk = mapply( FUN = function(y,m) rowSums(y)/rowSums(m) , yy[contrast] , mm[contrast] )
   K = length(unique(groups))
   z <- y/m
   zz = .DBSplitIntoGroups(z, groups)
   cols1 <- groups == unique(groups)[contrast[1]]
   cols2 <- groups == unique(groups)[contrast[2]]
   wh.cols <- c(which(cols1),which(cols2))
   n_eff <- apply(z[,wh.cols], 1, FUN=function(vec) sum(vec!=1 & vec!=0, na.rm=TRUE))
   n_eff[ n_eff<=K ] <- K+1

   if( !use.all.groups )
    {
      cols1 <- groups == unique(groups)[contrast[1]]
      cols2 <- groups == unique(groups)[contrast[2]]
      wh.cols <- c(which(cols1),which(cols2))
      y <- y[,wh.cols]
      m <- m[,wh.cols]
      groups <- as.character(groups[ wh.cols ])
      groups <- as.factor(groups)
      cols1 <- groups == unique(groups)[1]
      cols2 <- groups == unique(groups)[2]
      rows.kp <- rowSums(y)!=rowSums(m) & rowSums(m)!=0 & apply(m[,cols1,drop=FALSE],1,FUN=function(vec) sum(vec!=0)>1) & apply(m[,cols2,drop=FALSE],1,FUN=function(vec) sum(vec!=0)>1)
      shrink.method <- match.arg(shrink.method, c("WEB","DEB"))
      y <- y[rows.kp,]
      m <- m[rows.kp,]
      targets <- rownames(y)
      K = length(unique(groups))
      J = nrow(y)
    }
   else
    {
       mm = .DBSplitIntoGroups( m , groups )[ as.character(unique(groups)) ]
       mm.sums <- sapply(mm,FUN=function(mtx) apply(mtx,1,FUN=function(vec) sum(vec!=0)>1))
       rows.kp <- rowSums(y)!=rowSums(m) & rowSums(m)!=0 & rowSums(y)!=0 & apply(mm.sums,1,all)
       y <- y[rows.kp,]
       m <- m[rows.kp,]

       shrink.method <- match.arg(shrink.method, c("WEB","DEB"))
       groups <- as.factor(groups)
       targets <- rownames(y)
       K = length(unique(groups))
       J = nrow(y)
    }

   z <- y/m
   zz = .DBSplitIntoGroups(z, groups)
   neff <- apply(z, 1, FUN=function(vec) sum(vec!=1 & vec!=0, na.rm=TRUE))
   neff[ neff<=K ] <- K+1

   SNULL = .DBS(y=y,m=m,groups=NULL)
   DNULL = switch( shrink.method , "WEB" = EstimateWEBDisp(y=y,m=m,groups=NULL,neff=neff,S=SNULL),
                            "DEB" = EstimateDEBDisp(y=y,m=m,groups=NULL,neff=neff,S=SNULL) )
   S = .DBS(y=y,m=m,groups=groups)
   D = switch( shrink.method , "WEB" = EstimateWEBDisp(y=y,m=m,groups=groups,neff=neff,S=S),
                        "DEB" = EstimateDEBDisp(y=y,m=m,groups=groups,neff=neff,S=S) )
   D[D==0] <- sqrt(.Machine$double.eps)
   DNULL[DNULL==0] <- sqrt(.Machine$double.eps)
    names(D) <- targets
    names(DNULL) <- targets

   if( use.all.groups )
    {
       cols1 <- groups == unique(groups)[contrast[1]]
       cols2 <- groups == unique(groups)[contrast[2]]
       wh.cols <- c(which(cols1),which(cols2))
       y <- y[,wh.cols]
       m <- m[,wh.cols]
       groups <- as.character(groups[ wh.cols ])
       groups <- as.factor(groups)
       cols1 <- groups == unique(groups)[1]
       cols2 <- groups == unique(groups)[2]
       rows.kp2 <- rowSums(y)!=rowSums(m) & rowSums(m)!=0 & rowSums(y)!=0 & apply(m[,cols1,drop=FALSE],1,FUN=function(vec) sum(vec!=0)>1) & apply(m[,cols2,drop=FALSE],1,FUN=function(vec) sum(vec!=0)>1)
       y <- y[rows.kp2,]
       m <- m[rows.kp2,]
       K = length(unique(groups))
       J = nrow(y)
       SNULL = .DBS(y=y,m=m,groups=NULL)
       S = .DBS(y=y,m=m,groups=groups)
    }

   nj <- matrix( table(groups)[ unique(as.character(groups)) ], nrow = J, ncol = K)
   if( any(m == 0) )
    {
      wh.0 <- which(apply(m, 1, FUN=function(vec) any(vec == 0)))
      nj[wh.0,] = t(apply(m[wh.0,,drop=FALSE],1,FUN=function(vec) table(groups[which(vec!=0)])[unique(as.character(groups))]))
    }
   yy = .DBSplitIntoGroups( y , groups )
   mm = .DBSplitIntoGroups( m , groups )

   if( use.all.groups ) LRT = 2*( DNULL[rownames(y)]*SNULL - D[rownames(y)]*S + 0.5*rowSums(nj)*(log(D[rownames(y)]) - log(DNULL[rownames(y)])) )
    else LRT = 2*( DNULL*SNULL - D*S + 0.5*rowSums(nj)*(log(D) - log(DNULL)) )

   pvals = pchisq(LRT, df=K-1, lower.tail=FALSE)
   padj <- p.adjust(pvals, method="BH")
   P <- vector( length=length(rows.kp) , mode="numeric" )
    names(P) <- names(rows.kp)
   Padj <- vector( length=length(rows.kp) , mode="numeric" )
    names(Padj) <- names(rows.kp)
   model.disp <- vector( length=length(rows.kp) , mode="numeric" )
    names(model.disp) <- names(rows.kp)
   null.disp <- vector( length=length(rows.kp) , mode="numeric" )
    names(null.disp) <- names(rows.kp)
   wh.tested <- match(names(pvals),names(P))
   P[wh.tested] <- pvals
   P[-wh.tested] <- NA
   Padj[wh.tested] <- padj
   Padj[-wh.tested] <- NA
   model.disp[names(D)] <- D
   model.disp[!(names(model.disp)%in%names(D))] <- NA
   null.disp[names(DNULL)] <- DNULL
   null.disp[!(names(null.disp)%in%names(DNULL))] <- NA

   all <- cbind( Tk , P , Padj , n_eff , m.offs , model.disp , null.disp )
    colnames( all ) <- c(paste0("MLE_",unique(groups)[1]),paste0("MLE_",unique(groups)[2]),
                         "pVal","Adj.pVal","n_eff",paste0("MeanTotCt_",unique(groups)[1]),
                         paste0("MeanTotCt_",unique(groups)[2]), "Model.Disp","Null.Disp")
   wh.sig <- which( all[,"Adj.pVal"] <= fdr.level )
   sig <- all[wh.sig,,drop=FALSE]
   out <- list( "Sig" = sig , "All" = all )
   return(out)
 }
