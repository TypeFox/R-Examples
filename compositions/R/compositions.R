# (C) 2005/2008 by Gerald van den Boogaart, Greifswald
# License: GPL version 2 or newer


gsiInt <- function(x,n=NULL) {if(!is.null(n))stopifnot(length(x)==n);as.integer(x)}
gsiDouble <- function(x,n=NULL) {if(!is.null(n))stopifnot(length(x)==n);as.numeric(x)}

gsi.plain <- function(x) {
  if( is.data.frame(x) )
    unclass(data.matrix(x))
  else 
    unclass(x)
}

gsi.simshape <- function(x,oldx) {
  if(length(dim(oldx))>=2 )
    oneOrDataset(x)
  else if( length(dim(oldx)) == 1 )
    structure(c(x),dim=length(x))
  else 
    c(drop(x))
}

gsi.diagExtract <- function(x) {
  if( length(x) > 1 )
    diag(x)
  else
    c(x)
}

gsi.diagGenerate <- function(x) {
  if( length(x) > 1 )
    diag(x)
  else
    matrix(x)
}

gsi.getD  <- function(x) ncol(oneOrDataset(x))
gsi.getN  <- function(x) nrow(oneOrDataset(x))   

gsi.eq <-  function(x,y) {
  if( is.null(y) ) return(is.null(x)) # null
  if( is.finite(y) ) {           
    if( is.infinite(1/y) & (1/y)<0 )  # -0
      return(is.infinite(1/x) & (1/x)<0)
      else
        return(is.finite(x) & x==y)     # Zahlencodes 
  }
  if( is.nan(y) ) return(is.nan(x))   # NaN
  if( is.infinite(y)&y>0) return(is.infinite(x)&x>0) # Inf
  if( is.infinite(y)&y<0) return(is.infinite(x)&x<0) # -Inf
  if( is.na(y) ) return(is.na(y))                    # NA
  stop("Unkown comparison type ",y)                  # Was wurde vergessen?
}


names.acomp <- function(x) colnames(oneOrDataset(x))
names.rcomp <- names.acomp
names.aplus <- names.acomp
names.rplus <- names.acomp
names.rmult <- names.acomp
names.ccomp <- names.acomp

"names<-.acomp" <- "names<-.rcomp" <- "names<-.aplus" <- "names<-.rplus" <- "names<-.rmult" <- "names<-.ccomp" <-
  function(x,value) {
    if(is.matrix(x)) {
      colnames(x) <- value
      x
    }
    else
      NextMethod("names",x,value=value)
  }


groupparts <- function(x,...) UseMethod("groupparts",x)

groupparts.rcomp <- function(x,...,groups=list(...)) {
                                        # BDL=SZ->0, MAR->MAR, MNAR->MNAR
  x <- rmult(gsi.recodeM2C(x,clo(x),BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA))
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- rplus(x,idx)
    ss %*% rep(1,gsi.getD(ss))
  })
  rcomp(gsi.recodeC2M(erg,na=MNARvalue,nan=MARvalue))
}

groupparts.rplus <- function(x,...,groups=list(...)) {
                                        # BDL=SZ->0, MAR->MAR, MNAR->MNAR
  x <- rmult(gsi.recodeM2C(x,BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA))
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- rplus(x,idx)
    ss %*% rep(1,gsi.getD(ss))
  })
  rplus(gsi.recodeC2M(erg,na=MNARvalue,nan=MARvalue))
}

groupparts.acomp <- function(x,...,groups=list(...)) {
                                        # BDL: BDL, NA: NA, 0: 0
  x <- rmult(gsi.recodeM2C(x,clo(x),BDL=NaN,SZ=NA,MAR=Inf,MNAR=NaN))
  #SZ <- is.SZ(x)               # keep regardless of the rest NA
  #MNAR <- is.MNAR(x)|is.BDL(x) # keep if no SZ               NaN
  #MAR  <- is.MAR(x)            # keep if no SZ or MNAR are in the way Inf
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- aplus(x,idx)
    if( is.matrix(ss) )
      gsi.geometricmeanRow(ss)
    else
      gsi.geometricmean(ss)
  })
  acomp(gsi.recodeC2M(erg,na=SZvalue,nan=MNARvalue,inf=MARvalue))
  
#  SZ   <- is.na(x)&!is.na(x)   # keep regardless of the rest NA
#  MNAR <- is.nan(x)            # keep if no SZ               NaN
#  MAR  <- !is.finite(x)&!is.na(x) # keep if no SZ or MNAR are in the way Inf
}

groupparts.aplus <- function(x,...,groups=list(...)) {
  x <- rmult(gsi.recodeM2C(x,clo(x),BDL=NaN,SZ=NA,MAR=Inf,MNAR=NaN))
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- aplus(x,idx)
    if( is.matrix(ss) )
      gsi.geometricmeanRow(ss)
    else
      gsi.geometricmean(ss)
  })
  aplus(gsi.recodeC2M(erg,na=SZvalue,nan=MNARvalue,inf=MARvalue))
}

groupparts.acomp <- function(x,...,groups=list(...)) {
                                        # BDL: BDL, NA: NA, 0: 0
  x <- rmult(gsi.recodeM2C(x,clo(x),BDL=NaN,SZ=NA,MAR=Inf,MNAR=NaN))
  #SZ <- is.SZ(x)               # keep regardless of the rest NA
  #MNAR <- is.MNAR(x)|is.BDL(x) # keep if no SZ               NaN
  #MAR  <- is.MAR(x)            # keep if no SZ or MNAR are in the way Inf
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- aplus(x,idx)
    if( is.matrix(ss) )
      gsi.geometricmeanRow(ss)
    else
      gsi.geometricmean(ss)
  })
  acomp(gsi.recodeC2M(erg,na=SZvalue,nan=MNARvalue,inf=MARvalue))
  
#  SZ   <- is.na(x)&!is.na(x)   # keep regardless of the rest NA
#  MNAR <- is.nan(x)            # keep if no SZ               NaN
#  MAR  <- !is.finite(x)&!is.na(x) # keep if no SZ or MNAR are in the way Inf
}

groupparts.ccomp <- function(x,...,groups=list(...)) {
  x <- rmult(x)
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    totals(ccomp(x,idx))
  })
  ccomp(erg)
}


# groupparts(x,G1=c("Cd","S"),G2=c("Co","Ni"),G3=c("As","F"))

clo <- function(X,parts=1:NCOL(oneOrDataset(X)),total=1,
                detectionlimit=attr(X,"detectionlimit"),
                BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL,
                storelimit=!is.null(attr(X,"detectionlimit"))) {
  X <- gsi.plain(X)
  # Collect the parts
  parts  <- unique(parts)
  if( is.character(parts) ) {
    partsn <- match(parts,colnames(X))
    if( any(is.na(partsn)) )
      stop("Unknown variable name",parts[is.na(partsn)])
    parts <- partsn
  }
  nparts <- length(parts)
  Xn <- gsi.plain(oneOrDataset(X))[,parts,drop=FALSE]
  drop <- length(dim(X)) < 2
  #if( any(na.omit(c(Xn)<0)) )
  #  stop("Negative values are not valid for amounts")
  # Processing of missings
  iMAR <- if( !is.null(MAR) ) gsi.eq(Xn,MAR) else FALSE
  iMNAR<- if( !is.null(MNAR) ) gsi.eq(Xn,MNAR) else FALSE
  iSZ  <- if( !is.null(SZ) ) gsi.eq(Xn,SZ) else FALSE
  iBDL <- if( !is.null(BDL) ) gsi.eq(Xn,BDL) else FALSE
  if( is.null(detectionlimit) ) {
    if( any(iBDL) )
      Xn[iBDL]<- BDLvalue
  } else if( is.matrix(detectionlimit) ) {
    if( nrow(Xn)!=nrow(detectionlimit) | ncol(Xn)!=ncol(detectionlimit) )
      warning("Matrix of Detectionlimits does not fit x")
    Xn <- ifelse( is.finite(detectionlimit) & detectionlimit>=0,
                ifelse( (is.finite(Xn) & X <= detectionlimit)|iBDL ,
                       -detectionlimit, Xn ),
                Xn)
  } else if( length( detectionlimit) > 1 ) {
    if( ncol(Xn)!=length(detectionlimit)  )
      warning("Length of Detectionlimits does not fit x")
    detectionlimit <- outer(rep(1,nrow(Xn),detectionlimit))
    Xn <- ifelse( is.finite(detectionlimit) & detectionlimit>=0,
                ifelse( (is.finite(Xn) & X <= detectionlimit)|iBDL ,
                       -detectionlimit, Xn ),
                Xn)
  } else if( is.finite(detectionlimit) && detectionlimit > 0 ) {
    Xn <- ifelse((Xn<=detectionlimit&Xn>=0)|iBDL,-detectionlimit,Xn)
  } else
    Xn <- ifelse(iBDL,BDLvalue,Xn)
  if( any(iMAR) ) Xn[iMAR] <- MARvalue
  if( any(iMNAR)) Xn[iMNAR]<- MNARvalue
  if( any(iSZ) )  Xn[iSZ]  <- SZvalue
  # Make the sum 1 ignoring missings
  scaling <- 1
  if( length(total) > 1 || !is.na(total) ) {
    nas <- !is.NMV(Xn)&!is.BDL(Xn) # Missings are not included in closing
    bdl <- is.BDL(Xn)              # BDLs are accordingly scaled 
    naValues <- Xn[nas]
    Xn[nas]<-0
    s <- c(ifelse(bdl,0,Xn) %*% rep(1,nparts))
    scaling <-  matrix(rep(s/total,nparts),ncol=nparts)
    Xn  <- Xn / scaling
    Xn[nas] <- naValues
  }
  erg <- gsi.simshape(Xn,X)
  if( storelimit) {
    if( length(detectionlimit) == 1 )
      detectionlimit <- matrix(detectionlimit,nrow=nrow(Xn),ncol=ncol(Xn))
    detectionlimit/scaling
    detectionlimit <- gsi.simshape(detectionlimit[,parts,drop=FALSE],X)
    attr(erg,"detectionlimit") <- detectionlimit
  }
  erg
}



acomp <- function(X,parts=1:NCOL(oneOrDataset(X)),total=1,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {
  eq <- function(x,y) identical(as.numeric(x),as.numeric(y),single.NA=FALSE)
  if( is.list(X) )
    X<-data.matrix(X)
  if( is.ccomp(X) )
    X <- unclass(X)+0.5
  if( !is.null(BDL) ) {
    if( is.finite(BDL) )
      bdl <- X==BDL
    else
      bdl <- sapply(X,eq,BDL)
  } else bdl <- FALSE
  if( !is.null(MAR) ) {
    if( is.finite(MAR) )
      mar <- X==MAR
    else
      mar <- sapply(X,eq,MAR)
  } else mar <- FALSE
  if( !is.null(MNAR) ) {
    if( is.finite(MNAR) )
      mnar <- X==MNAR
    else
      mnar <- sapply(X,eq,MNAR)
  } else mnar <- FALSE
  if( !is.null(SZ) ) {
    if( is.finite(SZ) )
      sz <- X==SZ
    else
      sz <- sapply(X,eq,SZ)
  } else sz <- FALSE
  if( any( is.finite(X) & X < 0 & !(mar|mnar|bdl|sz)))
    warning("Negative values in composition are used as detection limits")
  if( !is.null(MAR) && any(mar) ) X[mar]<-BDLvalue
  if( !is.null(MNAR)&& any(mnar) ) X[mnar]<-BDLvalue
  if( !is.null(bdl) && any(bdl)) X[bdl]<-MARvalue
  if( !is.null(SZ)&&any(sz) ) X[sz]<-SZvalue
  X <-  structure(clo(X,parts,total),class="acomp")
  if( !is.null(detectionlimit) && any(X==BDLvalue) ) {
    X[sapply(X,eq,BDLvalue)]<- -detectionlimit
  }
  if( warn.na ) {
    if( any(is.SZ(X))) 
      warning("Composition has structural zeros")
    if( any(is.MAR(X) | is.MNAR(X)))
      warning("Composition has missings")
    if( any(is.BDL(X)) )
      warning("Composition has values below detection limit")
  }
  X
}



rcomp <- function(X,parts=1:NCOL(oneOrDataset(X)),total=1,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {
  X <-  structure(clo(X,parts,total,detectionlimit=detectionlimit,BDL=BDL,MAR=MAR,MNAR=MNAR,SZ=SZ),class="rcomp")
  if( warn.na ) {
    if( any(is.SZ(X))) 
      warning("Composition has structural zeros")
    if( any(is.MAR(X) | is.MNAR(X)))
      warning("Composition has missings")
    #if( any(is.BDL(X)) )
     # warning("Composition has values below detection limit")
  }
  X
}


aplus <- function(X,parts=1:NCOL(oneOrDataset(X)),total=NA,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {
  if( is.ccomp(X) )
    X <- unclass(X)+0.5
  X <- gsi.simshape(clo(X,parts,total,detectionlimit=detectionlimit,BDL=BDL,MAR=MAR,MNAR=MNAR,SZ=SZ),X)
  if( warn.na ) {
    if( any(is.SZ(X))) 
      warning("aplus has structural zeros")
    if( any(is.MAR(X) | is.MNAR(X)))
      warning("aplus has missings")
    if( any(is.BDL(X)) )
      warning("aplus has values below detection limit")
  }
  class(X) <-"aplus"
  X
}

rplus <- function(X,parts=1:NCOL(oneOrDataset(X)),total=NA,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {
  X <- gsi.simshape(clo(X,parts,total,detectionlimit=detectionlimit,BDL=BDL,MAR=MAR,MNAR=MNAR,SZ=SZ),X)
  if( warn.na ) {
    if( any(na.omit(c(X)==0)) )
      warning("rplus has structural zeros")
    if( any(is.na(c(X)) & ! is.nan(c(X))))
      warning("rplus has missings")
    if( any(is.nan(c(X))))
      warning("rplus has values below detection limit")
  }
  class(X) <-"rplus"
  X
}

ccomp <- function(X,parts=1:NCOL(oneOrDataset(X)),total=NA,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {
  X <- gsi.simshape(clo(X,parts,total,detectionlimit=detectionlimit,BDL=BDL,MAR=MAR,MNAR=MNAR,SZ=SZ),X)
  if( warn.na ) {
    if( any(na.omit(c(X)==0)) )
      warning("ccomp has structural zeros")
    if( any(is.na(c(X)) & ! is.nan(c(X))))
      warning("ccomp has missings")
    if( any(is.nan(c(X))))
      warning("ccomp has values below detection limit")
  }
  class(X) <-"ccomp"
  X
}


rmult <- function(X,parts=1:NCOL(oneOrDataset(X)),
                  orig=attr(X,"orig"),missingProjector=attr(X,"missingProjector")) {
  X <- gsi.simshape(oneOrDataset(X)[,parts,drop=FALSE],X)
  attr(X,"orig") <- orig
  attr(X,"missingProjector")<-missingProjector
  class(X) <-"rmult"
  X
}

print.rmult <- function(x,...) {
  Odata <- attr(x,"orig")
  if( ! is.null(Odata) )
    attr(x,"orig") <- missingSummary(Odata)
  mp <- attr(x,"missingProjector")
  if( ! is.null(mp) )
    attr(x,"missingProjector") <- dim(mp)
  NextMethod(x,...)
}
gsi2.invperm <- function(i,n){
  i <- unique(c(i,1:n))
  j <- numeric(length(i))
  j[i]<-1:length(i)
  j
}




rcompmargin <- function(X,d=c(1,2),name="+",pos=length(d)+1,what="data") {
  if( what=="data" ) {
    X <- rcomp(X)
    drop <- length(dim(X)) < 2
    if( mode(d)=="character" )
      d <- match(d,colnames(X))
    X <- oneOrDataset(gsi.plain(X))
    d <- unique(d)
    if( NCOL(X) <= length(d) )
      return(rcomp(X))
    else if( NCOL(X) == length(d) +1)
      return( rcomp(cbind(X[,d,drop=FALSE],X[,-d,drop=FALSE]) ))
    Xm <- gsi.recodeM2C(X[,-d,drop=FALSE],BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA)
    rest =gsi.recodeC2M(Xm %*% rep(1,NCOL(Xm)),zero=BDLvalue,nan=MARvalue,na=MNARvalue)
    tmp <- rcomp(cbind(rest=rest ,X[,d,drop=FALSE] ))
    if( !is.null(colnames(tmp)) )
      colnames(tmp)[1]<-name
    if( pos != 1 )
      tmp <- tmp[,gsi2.invperm(pos,ncol(tmp))]
    if( drop )
      tmp <- drop(tmp)
    rcomp(tmp)
  } else if( what=="var" ) {
    if( mode(d)=="character" )
      d <- match(d,colnames(X))
    Vrest <- sum(X[-d,-d])
    V     <- X[d,d,drop=FALSE]
    C     <- c(apply(X[d,-d,drop=FALSE],1,sum))
    erg   <- rbind(c(Vrest,C),cbind(C,V))
    #erg   <- erg - apply(erg,1,mean)
    #erg   <- t(t(erg)-apply(erg,2,mean))
    if( abs(sum(erg)) >1E-10)
      warning("Scaling problem in rcompmarin for variances")
    if( is.null(colnames(X)) )
      colnames(X) <- paste("var",1:ncol(X))
    colnames(erg)<-c(name,colnames(X)[d])
    row.names(erg) <-colnames(erg)
    if( pos != 1 )
      erg <- erg[gsi2.invperm(pos,ncol(erg)),gsi2.invperm(pos,ncol(erg))]
    erg
  } else
  stop("Unkown type of data in rcompmargin:",what)
}


acompmargin <- function(X,d=c(1,2),name="*",pos=length(d)+1,what="data") {
  if( what == "data" ) {
    drop <- length(dim(X)) < 2
    if( mode(d)=="character" )
      d <- match(d,colnames(X))
    X <- oneOrDataset(gsi.plain(X))
    d <- unique(d)
    if( NCOL(X) <= length(d) )
      return(X)
    else if( NCOL(X) == length(d) +1)
      return( cbind(X[,d,drop=FALSE],X[,-d,drop=FALSE]) )
    Xm <- X[,-d,drop=FALSE]
    Xm <- gsi.recodeM2C(Xm,log(ifelse(Xm>0,Xm,NA)),BDL=-Inf,MAR=NaN,MNAR=NA)
    rest <- Xm %*% rep(1/NCOL(Xm),NCOL(Xm))
    rest <- gsi.recodeC2M(rest,exp(rest),inf=BDLvalue,nan=MARvalue,na=MNARvalue)
    tmp <- acomp(cbind(rest=rest ,X[,d,drop=FALSE] ))
    if( !is.null(colnames(tmp)) )
      colnames(tmp)[1]<-name
    if( pos != 1 )
      tmp <- tmp[,gsi2.invperm(pos,ncol(tmp))]
    if( drop )
      tmp <- drop(tmp)
    acomp(tmp)
  } else if(what=="var") {
    if( mode(d)=="character" )
      d <- match(d,colnames(X))
    Vrest <- mean(X[-d,-d])
    V     <- X[d,d,drop=FALSE]
    C     <- c(apply(X[d,-d,drop=FALSE],1,mean))
    erg   <- rbind(c(Vrest,C),cbind(C,V))
    erg   <- erg - apply(erg,1,mean)
    erg   <- t(t(erg)-apply(erg,2,mean))
    if( abs(sum(erg)) >1E-10)
      warning("Scaling problem in acompmargin for variances")
    erg <- ilrvar2clr(clrvar2ilr(erg))
    if( is.null(colnames(X)) )
      colnames(X) <- paste("var",1:ncol(X))
    colnames(erg)<-c(name,colnames(X)[d])
    row.names(erg) <-colnames(erg)
    if( pos != 1 )
      erg <- erg[gsi2.invperm(pos,ncol(erg)),gsi2.invperm(pos,ncol(erg))]
    erg
  } else 
  stop("Unkown type of data in acompmargin:",what)
}


oneOrDataset <- function(W,B=NULL) {
  W <- gsi.plain(W)
  if( missing(B) || length(dim(B))!= 2 ) {
    if( length(dim(W)) == 2) {
      return( W )
    }
    else {
      tmp <- matrix(c(W),nrow=1)
      colnames(tmp) <- names(W)
      return(tmp)
    }
  } else {
    if( length(dim(W)) == 2) {
      return( W )
    }
    else {
      tmp <- matrix(c(W),nrow=NROW(B),ncol=length(W),byrow=TRUE)
      colnames(tmp)<- names(W)
      return(tmp)
    }
  }
}

gsi.geometricmean <- function(x,...) {
    exp(mean(log(c(unclass(x))),...))
}

gsi.geometricmeanRow <- function(x,...) apply(x,1,gsi.geometricmean,...)
gsi.geometricmeanCol <- function(x,...) apply(x,2,gsi.geometricmean,...)


geometricmean <- function(x,...) {
  if( any(na.omit(x==0)) )
    0
  else
    exp(mean(log(unclass(x)[is.finite(x)&x>0]),...))
}

geometricmeanRow <- function(x,...) apply(x,1,geometricmean,...)
geometricmeanCol <- function(x,...) apply(x,2,geometricmean,...)

meanCol <- function( x , ... , na.action=get(getOption("na.action"))) {
  apply(oneOrDataset(x),2,function(x,...) mean(x,na.action=na.action,...),...)
}

meanRow <- function( x , ... , na.action=get(getOption("na.action"))) {
  meanCol(t(x),...,na.action=na.action)
}


totals <- function( x , ... ) UseMethod("totals",x)

totals.acomp <- function(x,...,missing.ok=TRUE) {
  x <- oneOrDataset(x)
  if( missing.ok )
    x <- gsi.recodeM2C(x,BDL=0.0,SZ=0.0,MAR=0.0,MNAR=0.0)
  else 
    x <- gsi.recodeM2C(x,BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA)
  erg <- gsi.recodeC2M(apply(x,1,sum,...),zero=BDLvalue,nan=MARvalue,na=MNARvalue)
  erg
}

totals.aplus <- totals.acomp
totals.rcomp <- totals.acomp
totals.rplus <- totals.acomp
totals.ccomp <- totals.acomp

gsi.svdsolve <- function(a,b,...,cond=1E-10) {
  s <- svd(a,...)
  lambda1 <- s$d[1] 
  drop(s$v %*% (gsi.diagGenerate(ifelse(s$d<lambda1*cond,0,1/s$d)) %*% (t(s$u) %*% b)))
  
}

gsi.svdinverse <- function(a,...,cond=1E-10) {
  s <- svd(a,...)
  lambda1 <- s$d[1] 
  s$v %*% (gsi.diagGenerate(ifelse(s$d<lambda1*cond,0,1/s$d)) %*% t(s$u))
}


gsi.csum <- function(x) {c(rep(1,nrow(x)) %*% ifelse(is.finite(x),x,0))}
gsi.rsum <- function(x) {c(ifelse(is.finite(x),x,0)%*%rep(1,ncol(x))) }


mean.aplus<- mean.rplus <- mean.rcomp <- mean.acomp <- function( x, ... ,robust=getOption("robust")) {
  idtInv(mean(idt(x),...,robust=robust),x)
}

mean.ccomp <- function( x, ... ,robust=getOption("robust")) {
  ccomp(c(rep(1,nrow(x))%*% unclass(x)))
}

#mean.acomp <- function( x, ... ) {
#  if( has.missings(x) ) 
#    clrInv(gsi.svdsolve(sumMissingProjector(x),gsi.csum(gsi.cleanR(clr(x)))))
#  else
#    clrInv(mean(clr(x),...))
#}


#mean.rcomp <- function( x,... ) {
#  if( has.missings(x) ) 
#    cptInv(gsi.svdsolve(sumMissingProjector(x), gsi.csum(gsi.cleanR(cpt(x))) ) )
#  else
#    cptInv(meanCol(cpt(x),...))
#}

#mean.aplus <- function( x,... ) {
#  if( has.missings(x) ) 
#    iltInv(gsi.svdsolve(sumMissingProjector(x), gsi.csum(gsi.cleanR(ilt(x))) ) )
#  else
#    iltInv(meanCol(ilt(x),...))
#}

#mean.rplus <- function( x,...) {
#  if( has.missings(x) ) 
#    iitInv(gsi.svdsolve(sumMissingProjector(x), gsi.csum(gsi.cleanR(iit(x))) ) )
 # else
#    iitInv(meanCol(iit(x),...))
#}


  

mean.rmult <- function( x,...,na.action=NULL,robust=getOption("robust")) {
  if( !is.null(na.action) ) {
    x <- na.action(x)
  }
  control <- attr(robust,"control")
  if( is.logical(robust) ) robust <- if(robust) "mcd" else "pearson"
  if( has.missings(x) ) {
    if( !(is.character(robust) && robust=="pearson" ))
      warning("mean.rmult: Robust estimation currently not supported with missings")
    rmult(gsi.svdsolve(sumMissingProjector(x), gsi.csum(x,...)))
  }
  else {
    if( is.character(robust) ) {
      rmult(switch(robust,
             pearson=do.call(meanCol,c(list(x=unclass(x)),control,...)),
             mcd={
               #require("robustbase")
               if(!is.null(control)) covMcd(unclass(x),...,control=control)$center else covMcd(unclass(x),...)$center},
                   
             stop("mean.rmult: Unkown robustness type:",robust)
               ))
    } else if(is.function(robust)) {
      rmult(robust("mean",unclass(x),...,robust=robust))
    } else stop("mean.rmult: Unkown robustness type:",robust)      
  }
}


clr2ilr <- function( x , V=ilrBase(x) ) {
  gsi.simshape(gsi.recodeC2M(oneOrDataset(x),ninf=0,nan=0,na=0,inf=0) %*% V , x)
}

ilr2clr <- function( z , V=ilrBase(z=z),x=NULL ) {
  erg <- oneOrDataset(z) %*% t(V)
  if( !is.null(x) )
    colnames(erg)<-colnames(x)
  gsi.simshape( erg , z)
}


clrvar2ilr <- function( varx , V=ilrBase(D=ncol(varx)) ) {
  t(V) %*% varx %*% V
}

ilrvar2clr <- function( varz , V=ilrBase(D=ncol(varz)+1),x=NULL ) {
  erg <- V %*% varz %*% t(V)
  if( !is.null(x)) {
    colnames(erg) <- colnames(x)
    row.names(erg) <- colnames(x)
  }
  erg
}



var         <- function(x,...) UseMethod("var",x)
var.default <- function (x, y = NULL, na.rm = FALSE, use,...) stats::var(x,y,na.rm,use)

var.acomp <- function(x,y=NULL,...,
                      robust=getOption("robust"), use="all.obs",
                      giveCenter=FALSE) {
  control <- attr(robust,"control")
  if( is.logical(robust) ) robust <- if(robust) "mcd" else "pearson"
  if( has.missings(x) ||  has.missings(y) ) {
    if( !(is.character(robust) && robust=="pearson" ))
      warning("var.*: Robust estimation with losts not yet implemented")
    if(is.null(y)){
      if(use=="pairwise.complete.obs"){
        return(gsi.varwithlosts(cdt(x),giveCenter=giveCenter) )
      }else{
        tk = as.logical( gsi.geometricmeanRow( is.NMV(x) ) )
        xaux = x[tk,]
          class(xaux) = class(x)
        return(var(cdt(xaux),giveCenter=giveCenter))
      }
    }else{
      warning("Covariance with losts not yet implemented. Omitting lost values.")
        tk = as.logical( gsi.geometricmeanRow( is.NMV(cbind(x,y)) ) )
        xaux = x[tk,]
          class(xaux) = class(x)
        yaux = y[tk,]
          class(yaux) = class(y)
        return(var(unclass(cdt(xaux)),unclass(cdt(yaux)),...))
    }
  } else {
    switch(robust,
           pearson={
             erg <- var(unclass(cdt(x)),unclass(cdt(y)),...)
             if(giveCenter)
               attr(erg,"center")<-mean(x,robust=FALSE)
             erg
             },
           mcd={
             #require("robustbase")
             if(is.null(y)) {
               erg <- ilrvar2clr(if( is.null(control)) covMcd(idt(x),...)$cov else covMcd(idt(x),...,control=control)$cov,x=x)
               if(giveCenter)
                 attr(erg,"center")<-mean(x,robust=FALSE)
               erg
             } else {
               Dx<- ncol(x)
               Dy<- ncol(y)
               x1 <- idt(x)
               y1 <- idt(y)
               erg <- if( is.null(control)) covMcd(cbind(x1,y1),...) else covMcd(cbind(x1,y1),...,control=control)
               m <- erg$center
               erg <- erg$cov
               erg <- t(ilrBase(D=Dx)) %*% erg[1:(Dx-1),Dx:ncol(erg)] %*% ilrBase(D=Dy)
               row.names(erg) <- colnames(x)
               colnames(erg) <- colnames(y)
               if( giveCenter )
                 attr(erg,"center")<-idtInv(m,x)
               erg
             }},
           stop("var.*: unkown robust method",robust)
           )
  }
}
var.rcomp <- var.acomp

var.aplus  <- function(x,y=NULL,...,robust=getOption("robust"), use="all.obs",
                       giveCenter=FALSE) {
  control <- attr(robust,"control")
  if( is.logical(robust) ) robust <- if(robust) "mcd" else "pearson"
  if( has.missings(x) ||  has.missings(y) ) { 
    if( !(is.character(robust) && robust=="pearson" ))
      warning("var.*: Robust estimation with losts not yet implemented")
    if(is.null(y)){
      if(use=="pairwise.complete.obs"){
        return(gsi.varwithlosts(cdt(x),giveCenter=giveCenter) )
      }else{
        tk = as.logical( gsi.geometricmeanRow( is.NMV(x) ) )
        xaux = x[tk,]
          class(xaux) = class(x)
        return(var(cdt(xaux),giveCenter=giveCenter))
      }
    }else{
      warning("Covariance with losts not yet implemented. Omitting lost values.")
        tk = as.logical( gsi.geometricmeanRow( is.NMV(cbind(x,y)) ) )
        xaux = x[tk,]
          class(xaux) = class(x)
        yaux = y[tk,]
          class(yaux) = class(y)
        return(var(unclass(cdt(xaux)),unclass(cdt(yaux)),...))
    }
  } else {
    switch(robust,
           pearson={
             erg <- var(unclass(cdt(x)),unclass(cdt(y)),...)
             if( giveCenter) attr(erg,"center")<-mean(x,robust=FALSE)
             erg
           }
             ,
           mcd={
             #require("robustbase")
             if(is.null(y)) {
               erg <- if( is.null(control)) covMcd(idt(x),...)$cov else covMcd(idt(x),...,control=control)$cov
               if(giveCenter)
                 attr(erg,"center")<-mean(x,robust=FALSE)
               erg
             } else {
               Dx<- ncol(x)
               Dy<- ncol(y)
               x1 <- cdt(x)
               y1 <- cdt(y)
               erg <- if( is.null(control)) covMcd(cbind(x1,y1),...)$cov else covMcd(cbind(x1,y1),...,control=control)$cov
               m <- erg$center
               erg <- erg$cov
               erg <- erg[1:Dx,(Dx+1):ncol(erg)]
               row.names(erg) <- colnames(x)
               colnames(erg) <- colnames(y)
               if( giveCenter) attr(erg,"center")<-cdtInv(m,x)
               erg
             }
           },
           stop("var.*: unkown robust method",robust)
           )
  }
}
var.rplus <- var.aplus
var.rmult <- var.aplus

cov         <- function(x,y=x,...) UseMethod("cov",x)
cov.default <- function (x, y = NULL, use = "everything", method = c("pearson", 
    "kendall", "spearman"),...) stats::cov(x,y,use,method)

cov.acomp   <- var.acomp
cov.rcomp <- var.rcomp
cov.aplus <- var.aplus
cov.rplus <- var.rplus
cov.rmult <- var.rmult

cor <- function(x,y=NULL,...) UseMethod("cor",x)
cor.default <- function (x, y = NULL, use = "everything", method = c("pearson", 
    "kendall", "spearman"),...) stats::cor(x,y,use,method)




cor.acomp <- function(x,y=NULL,...,robust=getOption("robust")) {
  mat2cor <- function(x) {
    if( nrow(x) < 2 )
      return(x/x)
    sf<-diag(1/sqrt(diag(x)))
    structure( sf %*% x %*% sf ,dimnames=dimnames(x))
  }
  if( is.null(y) ) {
    mat2cor(var(x,y,...,robust=robust))
  } else {
    varX <- var(x,y=NULL,...,robust=robust)
    varY <- var(y,y=NULL,...,robust=robust)
    covXY<- var(x,y,...,robust=robust)
    sfX<-diag(1/sqrt(diag(varX)))
    sfY<-diag(1/sqrt(diag(varY)))
    structure( sfX %*% covXY %*% sfY, dimnames=list(colnames(x),colnames(y)))
  }
}

cor.rcomp <- cor.acomp 
cor.aplus <- cor.acomp
cor.rplus <- cor.acomp
cor.rmult <- cor.acomp

#  function(x,y=NULL,...) {
#  cor(unclass(x),unclass(cdt(y)),...)
#}


powerofpsdmatrix <- function(M,p,...) {
  s <- svd(M,...)
  d <- ifelse( abs(s$d)>max(abs(s$d))*1E-10, s$d^p,0)
  s$u %*% gsi.diagGenerate(d) %*% t(s$v)
}

mvar <- function(x,...) UseMethod("mvar",x)
mcov <- function(x,...) UseMethod("mcov",x)
mcor <- function(x,...) UseMethod("mcor",x)
msd  <- function(x,...) UseMethod("msd",x)

mvar.default <- function(x,y=NULL,...) {
  sum(gsi.diagExtract(var(x,y,...)))
}

mcov.default <- function(x,y=x,...) {
  sum(abs(svd(cov(idt(x),idt(y),...))$d))
}

msd.default <- function(x,y=NULL,...) {
  sqrt(mean(gsi.diagExtract(var(idt(x),y=NULL,...))))
}

mcor.default <- function(x,y,...) {
  ix <- scale(idt(x),center=TRUE,scale=FALSE)
  ix <- ix %*% powerofpsdmatrix(var(ix),-1/2)
  iy <- scale(idt(y),center=TRUE,scale=FALSE)
  iy <- iy %*% powerofpsdmatrix(var(iy),-1/2)
  mcov(ix,iy)
}



summary.acomp <- function( object,...,robust=getOption("robust") ) {
  W <- clo(gsi.plain(object))
  Wq <- apply(gsi.recodeM2C(W,BDL=NaN,SZ=NaN,MAR=NaN,MNAR=NaN),
              1,function(w) outer(w,w,"/"))
  dim(Wq)<-c(ncol(W),ncol(W),nrow(W))
  dimnames(Wq) <- list(colnames(W),colnames(W),NULL)
  vari <- if(is.null(robust) ) NULL else variation.acomp(acomp(W),robust=robust)
  narm <- function(x) x[is.finite(x)]
  structure(list(mean=if(is.null(robust)) NULL else mean(acomp(W),robust=robust),
       mean.ratio=apply(Wq,1:2,function(x) exp(mean(log(x[is.finite(x)])))),
       variation=vari,
       expsd=if( is.null(vari) ) NULL else exp(sqrt(vari)),
       invexpsd=if( is.null(vari) ) NULL else exp(-sqrt(vari)),
       min=apply(Wq,1:2,function(x) min(narm(x))),
       q1 =apply(Wq,1:2,function(x,...) quantile(narm(x),...),probs=0.25),
       med=apply(Wq,1:2,function(x,...) median(narm(x),...)),
       q3 =apply(Wq,1:2,function(x,...) quantile(narm(x),...),probs=0.75),
       max=apply(Wq,1:2,function(x,...) max(narm(x),...)),
       missingness=missingSummary(object)  
       ),class="summary.acomp")
       
}

summary.aplus <- function( object,...,digits=max(3, getOption("digits")-3),robust=NULL  ) {
  if( !missing(robust) )
    if( if(is.logical(robust)) robust else robust!="pearson" )
      warning("robustness currently not supported in summary.aplus")
  object <- ilt(object)
  erg <- sapply(data.frame(object),summary,...,digits=18)
  erg <- apply(erg,1:2,exp)
  erg <- apply(erg,1:2,signif,digits=digits)
  if( any( !is.NMV(object)) ) {
    attr(erg,"missingness")<-missingSummary(object)  
  }
  class(erg) <- c("summary.aplus",class(erg))
  erg       
}

summary.rplus <- function( object,... ,robust=NULL  ) {
  if( !missing(robust) )
    if( if(is.logical(robust)) robust else robust!="pearson" )
      warning("robustness currently not supported in summary.rplus")
  object <- iit(object)
  erg <- sapply(data.frame(object),summary,...)
  if( any( !is.NMV(object)) ) {
    attr(erg,"missingness")<-missingSummary(object)  
  }
  class(erg) <- c("summary.rplus",class(erg))
  erg       
}

summary.rmult <- function( object,...  ,robust=NULL ) {
  if( !missing(robust) )
    if( if(is.logical(robust)) robust else robust!="pearson" )
      warning("robustness currently not supported in summary.mult")
  object <- unclass(object)
  erg <- sapply(data.frame(object),summary,...)
  class(erg) <- c("summary.rmult",class(erg))
  erg       
}

summary.rcomp <- function( object,...,robust=NULL ) {
  # must support robust = NULL for no estimation with missing methods
  if( !missing(robust) )
    if( if(is.logical(robust)) robust else robust!="pearson" )
      warning("robustness currently not supported in summary.rcomp")
  object <- clo(gsi.plain(object)) 
  erg <- sapply(data.frame(object),function(x,...) summary(x[is.NMV(x)],...),...)
  if( has.missings(object) ) attr(erg,"missingness")<-missingSummary(object)
  class(erg) <- c("summary.rcomp",class(erg))
  erg       
}



vp.logboxplot <- function(x,y,...,dots=FALSE,boxes=TRUE,xlim=NULL,ylim=NULL,log=TRUE,notch=FALSE,plotMissings=TRUE,mp=~simpleMissingSubplot(missingPlotRect,
                                                                                         missingInfo,c("NM","TM",cn)),
         missingness=attr(y,"missingness")                                                                                ) {
  if(is.null(missingness))
    plotMissings <- FALSE
  fakMis <- FALSE
  if( any(is.na(x)) )  {
    fakMis<- TRUE
    levels(x) <- c(levels(x),"ERR")
    x[is.na(x)]<-"ERR"
  }
  nmv <- oneOrDataset(missingness=="NMV")
  nMis <- apply(!nmv,1,sum)
  nonmissing <- nMis==0
  lf <- length(levels(x))
  if( boxes ) {
    stats <- boxplot(split(log(ifelse(nonmissing,y,NA)),x),plot=FALSE)
    stats$stats <- exp(stats$stats)

    stats$conf  <- exp(stats$conf)
    stats$out   <- exp(stats$out)
    bxp(stats,add=TRUE,at=1:lf,width=rep(1,lf),notch=notch)
      
  }
  if( dots  ) points(x,y,...)
  if( plotMissings && !all(nmv)) {
    wM <- apply(!nmv,1,function(x) c(which(x),1)[1])
    missingTab <-  cbind(NotMissing=tapply(nonmissing,x,sum),
                     TotallyMissing=tapply(nMis>1,x,sum),
                     oneOrDataset(apply(nmv,2,function(w) tapply(nMis==1 & !w,x,sum))) 
                     )
    cn <- colnames(missingness)
    for( i in 1:length(levels(x))) {
      lev <- levels(x)[i]
      missingInfo <- missingTab[i,]
      if( sum(missingInfo[-1])>0 ) {
        usr <- par("usr")
        missingPlotRect <- c(i+0.45,i+0.5,usr[3],usr[4])
        eval(mp[[2]])
      }
    }
  }
}



#vp.boxplot <- function(x,y,...,dots=FALSE,boxes=TRUE,xlim,ylim,log,notch=FALSE,
#                       plotMissings=TRUE,
#                       mp=~simpleMissingSubplot(missingPlotRect,
#                                                missingInfo,c("NM","TM",cn)),
#                       missingness=attr(y,"missingness")) {
#    if( boxes ) boxplot(split(y,x),add=TRUE,notch=notch)
#    if( dots  ) points(x,y,...)
#}


vp.boxplot <- function(x,y,...,dots=FALSE,boxes=TRUE,xlim=NULL,ylim=NULL,log=FALSE,notch=FALSE,plotMissings=TRUE,mp=~simpleMissingSubplot(missingPlotRect,
                                                                                         missingInfo,c("NM","TM",cn)),
         missingness=attr(y,"missingness")                                                                                ) {
  if(is.null(missingness))
    plotMissings <- FALSE
  fakMis <- FALSE
  if( any(is.na(x)) )  {
    fakMis<- TRUE
    levels(x) <- c(levels(x),"ERR")
    x[is.na(x)]<-"ERR"
  }
  nmv <- oneOrDataset(missingness=="NMV")
  nMis <- apply(!nmv,1,sum)
  nonmissing <- nMis==0
  lf <- length(levels(x))
  if( boxes ) {
    stats <- boxplot(split(ifelse(nonmissing,y,NA),x),plot=FALSE)
    #stats$stats <- exp(stats$stats)

    #stats$conf  <- exp(stats$conf)
    #stats$out   <- exp(stats$out)
    bxp(stats,add=TRUE,at=1:lf,width=rep(1,lf),notch=notch)
      
  }
  if( dots  ) points(x,y,...)
  if( plotMissings && !all(nmv)) {
    wM <- apply(!nmv,1,function(x) c(which(x),1)[1])
    missingTab <-  cbind(NotMissing=tapply(nonmissing,x,sum),
                     TotallyMissing=tapply(nMis>1,x,sum),
                     oneOrDataset(apply(nmv,2,function(w) tapply(nMis==1 & !w,x,sum)) )
                     )
    cn <- colnames(missingness)
    for( i in 1:length(levels(x))) {
      lev <- levels(x)[i]
      missingInfo <- missingTab[i,]
      if( sum(missingInfo[-1])>0 ) {
        usr <- par("usr")
        missingPlotRect <- c(i+0.45,i+0.5,usr[3],usr[4])
        eval(mp[[2]])
      }
    }
  }
}



gsi.textpanel <- function(x,y,lab,...) {
  par(usr=c(0,1,0,1),xlog=FALSE,ylog=FALSE)
  text(0.5,0.5,lab,...)
}


# changed by Raimon on 2008-07-07
#   Now the verical scale of the plots has always the same length,
#   but each row of plots has its own ylim.
#   In this way, the boxplots are not too narrow, but
#   spread and location are still inter-comparable.
boxplot.acomp <- function(x,fak=NULL,...,
                          xlim=NULL,ylim=NULL,log=TRUE,
                          panel=vp.logboxplot,dots=!boxes,boxes=TRUE,
                          notch=FALSE,
                          plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                missingInfo,c("NM","TM",cn))) {
  X <- acomp(x)
  if( is.null(fak) )
    fak <- factor(rep("",nrow(X)))
  if( is.null(xlim)) {
    if( is.factor(fak) )
      xlim <- c(0,nlevels(fak)+1)
    else {
      xlim <- range(fak)
    }
  }
  if( !is.factor(fak) ) {
    boxes <- FALSE
    dots  <- TRUE
  }
  if( is.function(panel) )
    panel <- list(panel)
  ipanel <- function(x,y,...) {
    ic <- gsi.mapfrom01(log(x))
    jc <- gsi.mapfrom01(log(y))
    mis <- missingType(X[,c(ic,jc)])
    Y <- X[,c(jc)]/X[,c(ic)]
    attr(Y,"missingness") <- mis
    #a <- gsi.recodeM2Clean(unclass(X)[,gsi.mapfrom01(log(x))])
    #b <- gsi.recodeM2Clean(unclass(X)[,])
    for( thispanel in panel ) {
      thispanel(fak,Y,...,notch=notch,dots=dots,boxes=boxes,plotMissings=plotMissings,mp=mp)
    }
  }
  if( is.null(ylim) ) {
    a<-apply(x,1,function(x) {
      x <- x[is.finite(x)&x>0]
      if( length(x) > 0 ) {
        mi <- min(x)
        ma <- max(x)
        c(mi/ma,ma/mi)
      } else {
        c(1,1)
      }}
             )
    ylim <- range(a)
     nc = gsi.getD(X)
     ylims = sapply(1:nc, 
        function(i){
         aux = sapply(1:nc, function(j){
              aux = log(X[,i]/X[,j])
              aux = aux[is.finite(aux)]
              return( range(aux) )
        })
        return(c( min(aux[1,]), max(aux[2,]) ))
      })
     yrange = max(c(-1,1) %*% ylims )
     ylims = outer( c(0.5*c(1,1) %*% ylims), c(yrange*c(-1,1)/2), "+")
  }else{
     yrange = NULL
     ylims = NULL
  }

  #su <- summary.acomp(X,robust=NULL)
  #minq <- min(su$min)
  #maxq <- max(su$max)
  mm <- exp(sapply(1:NCOL(X),gsi.mapin01))
  colnames(mm) <- colnames(X)
  ipairs <- function (x, labels, panel = ipanel, ...,
                      font.main = par("font.main"),
                      cex.main = par("cex.main"), diag.panel = NULL, 
                      text.panel = textPanel,
                      label.pos = 0.5 , cex.labels = NULL, 
                      font.labels = 1, gap = 1,xlim, ylim, yrange, ylims, log) {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                     y, txt, cex = cex, font = font)
    localAxis <- function(side, xpd, bg, ...) axis(side, xpd = NA, 
                                                   ...)
    panel <- match.fun(panel)
    nc <- ncol(x)
    labels <- colnames(x)
    if (is.null(labels)) 
      labels <- paste("var", 1:nc)
    has.labs <- ! is.null(labels)
    oma <- c(4, 4, 4, 4)
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar),add=TRUE)
    for (i in 1:nc ){
     for (j in 1:nc) {
      if(!is.null(yrange)){ylim=exp(ylims[i,])}
      plot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
           type = "n", ...,xlim=xlim,ylim=ylim,log=log)
      box()
      mfg <- par("mfg")
      if (i == j) {
        if (has.labs) {
          par(usr = c(0, 1, 0, 1))
          if (is.null(cex.labels)) {
            l.wid <- strwidth(labels, "user")
            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
          }
          text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                     font = font.labels)
        }
      }
      else  
       panel(as.vector(x[, j]), as.vector(x[, i]), ...)
      if (any(par("mfg") != mfg)) 
        stop("The panel function made a new plot")
     }
    }
    invisible(NULL)
  }
  ipairs(mm,labels=colnames(X),panel=ipanel,...,log=ifelse(log,"y",""),
             ylim=ylim, yrange=yrange, xlim=xlim, ylims=ylims, text.panel=gsi.textpanel)
  
  replot(plot=match.call())  
  
}


boxplot.rcomp <- function(x,fak=NULL,...,
                         xlim=NULL,ylim=NULL,log=FALSE,panel=vp.boxplot,
                          dots=!boxes,boxes=TRUE,notch=FALSE,
                          plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                missingInfo,c("NM","TM",cn))) {
  X <- acomp(x)
  if( is.null(fak) )
    fak <- factor(rep("",nrow(X)))
  if( is.null(xlim) ) {
    if( is.factor(fak) )
      xlim <- c(0,nlevels(fak)+1)
    else {
      xlim <- range(fak)
    }
  }
  if( !is.factor(fak) ) {
    boxes <- FALSE
    dots  <- TRUE
  }
  if( is.function(panel) )
    panel <- list(panel)
  
  if( is.null(ylim) ) {
    if( log ) {
      a<-apply(x,1,function(x) {
        x <- x[is.finite(x)&x>0]
        if( length(x) > 0 ) {
          mi <- min(x)
          ma <- max(x)
          c(mi/ma,ma/mi)
        } else {
          c(1,1)
        }}
               )
      ylim <- log(range(a))
    } else ylim<-c(0,1)
  }
  ipanel <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(log(x))]
    b <- unclass(X)[,gsi.mapfrom01(log(y))]
    Y <- if(log) log(ifelse(is.NMV(b)&is.NMV(a),b/a,NA)) else ifelse(is.NMV(b)&is.NMV(a),b/(a+b),NA)
    attr(Y,"missingness")<-missingType(X[,c(gsi.mapfrom01(log(x)),gsi.mapfrom01(log(y)))])
    for( thispanel in panel ) 
      thispanel(fak,Y,
                ...,notch=notch,dots=dots,boxes=boxes,plotMissings=plotMissings,mp=mp)
  }


  #su <- summary.rcomp(X) ## Must be changed when robustness in rcomp arises
  #minq <- log(min(su[,"min"])/(min(su[,"min"])+max(su[,"max"])))
  #maxq <- log(max(su[,"max"])/(min(su[,"min"])+max(su[,"max"])))
  
  ipairs <- function (x, labels, panel = points, ..., 
                      font.main = par("font.main"),
                      cex.main = par("cex.main"), diag.panel = NULL, 
                      text.panel = textPanel,
                      label.pos = 0.5 , cex.labels = NULL, 
                      font.labels = 1, gap = 1,xlim,ylim,log="") {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                     y, txt, cex = cex, font = font)
    localAxis <- function(side, xpd, bg, ...) axis(side, xpd = NA, 
                                                   ...)
    panel <- match.fun(panel)
    nc <- ncol(x)
    #labels <- colnames(x)
    if (is.null(labels)) 
      labels <- paste("var", 1:nc)
    has.labs <- ! is.null(labels)
    oma <- c(4, 4, 4, 4)
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar),add=TRUE)
    for (i in 1:nc ) for (j in 1:nc) {
      plot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
           type = "n", ...,xlim=xlim,ylim=ylim,log=log)
      box()
      mfg <- par("mfg")
      if (i == j) {
        if (has.labs) {
          par(usr = c(0, 1, 0, 1))
          if (is.null(cex.labels)) {
            l.wid <- strwidth(labels, "user")
            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
          }
          text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                     font = font.labels)
        }
      }
      else  
        panel(as.vector(x[, j]), as.vector(x[, i]), ...)
      if (any(par("mfg") != mfg)) 
        stop("The panel function made a new plot")
    }
    invisible(NULL)
  }


  ipairs(exp(sapply(1:NCOL(X),gsi.mapin01)),labels=colnames(X),panel=ipanel,...,ylim=ylim,xlim=xlim)
  
  replot(plot=match.call())  

}

boxplot.rplus <- function(x,fak=NULL,...,ylim=NULL,log=FALSE,
                          plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                   missingInfo,
                                                   names(missingInfo))
                          ) {
  if( !is.null(fak) )
    warning("Spliting not yet implemente in boxplot.rplus")
  nmv <- oneOrDataset((is.NMV(x) ))
  xx <- ifelse(is.BDL(x),if(log) NA else 0,ifelse(nmv,x,NA))
  if( is.null(ylim) ) {
      ylim <- range(xx[nmv])
      if( !log )
        ylim[1]<-0
  }
  erg<-boxplot(as.data.frame(xx),...,ylim=ylim,log=if(log) "y" else "")
  if( plotMissings && !is.null(mp) &&
     !all(nmv)) {
    su <- missingSummary(x)
    for(i in 1:ncol(x)) {
      if(any(nmv[,i])) {
        usr <- par("usr")
        missingPlotRect <- c(i+0.45,i+0.5,usr[3],usr[4])
        cn <- colnames(x)[i]
        X<-x[,i]
        missingInfo <- su[i,]
        eval(mp[[2]])
      }
    }
  }
  replot(plot=match.call())  

  invisible(erg)
}

boxplot.aplus <- function(x,fak=NULL,...,log=TRUE,
                          plotMissings=TRUE,
                          mp=~simpleMissingSubplot(missingPlotRect,
                                                   missingInfo,
                                                   names(missingInfo))
) {
  if( !is.null(fak) )
    warning("Spliting not yet implemente in boxplot.aplus")
  
  stats <- boxplot(as.data.frame(ifelse(nmv<-is.NMV(x),ilt(x),NA)),plot=FALSE)
  delog <- function(x) {if(is.list(x)) lapply(x,delog) else exp(x)}
  stats$stats <- exp(stats$stats)
  stats$conf <- exp(stats$conf)
  stats$out <- exp(stats$out)
  erg<-bxp(stats,...,at=1:ncol(x),width=rep(1,ncol(x)),log=if(log) "y" else "")
  if( plotMissings && !is.null(mp) && !all(nmv) ) {
    su <- missingSummary(x)
    for(i in 1:ncol(x)) {
      if(any(nmv[,i])) {
        usr <- par("usr")
        missingPlotRect <- c(i+0.45,i+0.5,usr[3],usr[4])
        cn <- colnames(x)[i]
        X<-x[,i]
        missingInfo <- su[i,]
        eval(mp[[2]])
      }
    }
  }
  replot(plot=match.call())  

  invisible(erg)
}



vp.qqnorm <- function(x,y,...,alpha=NULL) {
  usr <- par("usr")
  usr[1:2] <- range(qnorm(ppoints(length(y))))
  usr[3:4] <- range(y)
  par( usr=usr )
  if( !is.null(alpha) && is.factor(x) ) 
    alpha <- alpha/nlevels(x)
  reject <- FALSE
  if( is.factor(x)) {
    for( k in split(y,x) ) {
      if( !is.null(alpha) && shapiro.test(k)$p < alpha )
        reject<-TRUE
      lines(qnorm(ppoints(length(k))),sort(k),...)
    }
  } else { 
    if( !is.null(alpha) && shapiro.test(y)$p < alpha )
        reject<-TRUE
    points(qnorm(ppoints(length(y))),sort(y),...)
  }
  qqline(y)
  if( reject )
    title(main="!",col.main="red")
    
}

qqnorm.acomp <- function(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL) {
  X <- acomp(y)
  if( !is.null(alpha) )
    alpha <- alpha/((nrow(X)*(nrow(X)-1)/2))
  if( is.function(panel) )
    panel <- list(panel)
  ipanel <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    v <- log(b/a)
    for( thispanel in panel )
      thispanel(fak,v,...,alpha=alpha)
  }
  pairs(sapply(1:NCOL(X),gsi.mapin01),labels=colnames(X),panel=ipanel,...)
  replot(plot=match.call())  

}

qqnorm.aplus <- function(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL) {
  X <- aplus(y)
  if( is.function(panel) )
    panel <- list(panel)
  if( !is.null(alpha) )
    alpha <- alpha/(nrow(X)^2)
  ipanelupper <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,log(b/a),...,alpha=alpha)
  }
  ipanellower <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,log(b*a),...,alpha=alpha)
  }
  ipaneldiag <- function(x,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    for( thispanel in panel )
      thispanel(fak,log(a),...,alpha=alpha)
  }
  itextpanel <- function(x,y,lab,...) {
    par(usr=c(0,1,0,1),xlog=FALSE,ylog=FALSE)
    text(0.1,0.9,lab,adj=c(0,1),...)
  }

  pairs(sapply(1:NCOL(X),gsi.mapin01),labels=colnames(X),lower.panel=ipanellower,upper.panel=ipanelupper,diag.panel=ipaneldiag,text.panel=itextpanel,...)
  replot(plot=match.call())  

}



qqnorm.rcomp <- function(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL) {
  X <- rcomp(y)
  if( is.function(panel) )
    panel <- list(panel)
  if( !is.null(alpha) )
    alpha <- alpha/(nrow(X)^2)
  ipanelupper <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,b-a,...,alpha=alpha)
  }
  ipanellower <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,b+a,...,alpha=alpha)
  }
  ipaneldiag <- function(x,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    for( thispanel in panel )
      thispanel(fak,a,...,alpha=alpha)
  }
  itextpanel <- function(x,y,lab,...) {
    par(usr=c(0,1,0,1),xlog=FALSE,ylog=FALSE)
    text(0.1,0.9,lab,adj=c(0,1),...)
  }

  pairs(sapply(1:NCOL(X),gsi.mapin01),labels=colnames(X),lower.panel=ipanellower,upper.panel=ipanelupper,diag.panel=ipaneldiag,text.panel=itextpanel,...)
  replot(plot=match.call())  

}

qqnorm.rplus <- function(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL) {
  X <- rplus(y)
  if( is.function(panel) )
    panel <- list(panel)
  if( !is.null(alpha) )
    alpha <- alpha/(nrow(X)^2)
  ipanelupper <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,b-a,...,alpha=alpha)
  }
  ipanellower <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,b+a,...,alpha=alpha)
  }
  ipaneldiag <- function(x,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    for( thispanel in panel )
      thispanel(fak,a,...,alpha=alpha)
  }
  itextpanel <- function(x,y,lab,...) {
    par(usr=c(0,1,0,1),xlog=FALSE,ylog=FALSE)
    text(0.1,0.9,lab,adj=c(0,1),...)
  }

  pairs(sapply(1:NCOL(X),gsi.mapin01),labels=colnames(X),lower.panel=ipanellower,upper.panel=ipanelupper,diag.panel=ipaneldiag,text.panel=itextpanel,...)
  replot(plot=match.call())  

}



gsi.drop  <-function(X,drop) if( drop ) drop(X) else X

is.acomp <- function(x) inherits(x,"acomp")

is.rcomp <- function(x) inherits(x,"rcomp")

is.aplus <- function(x) inherits(x,"aplus")

is.rplus <- function(x) inherits(x,"rplus")
   
is.rmult <- function(x) inherits(x,"rmult")

is.ccomp <- function(x) inherits(x,"ccomp")

perturbe <- function( x,y ) {
  acomp(gsi.mul(x,y))
}

perturbe.aplus <- function(x,y) {
  aplus(gsi.mul(x,y))
}



gsi.add <- function( x,y ) {
  if( length(dim(x)) == 2 )
    if( length(dim(y)) == 2 )
      unclass(x)+unclass(y)
    else
      unclass(x)+rep(c(y),rep(NROW(x),length(y)))
  else if( length(dim(y)) == 2 )
      unclass(y)+rep(c(x),rep(NROW(y),length(x)))
  else
    unclass(x)+unclass(y)
}

gsi.sub <- function( x,y ) {
 # drop <- length(dim(x)) < 2 && length(dim(y)) < 2
  if( length(dim(x)) == 2 )
    if( length(dim(y)) == 2 )
      unclass(x)-unclass(y)
    else
      unclass(x)-rep(c(y),rep(NROW(x),length(y)))
  else if( length(dim(y)) == 2 )
      unclass(y)-rep(c(x),rep(NROW(y),length(x)))
  else
    unclass(x)-unclass(y)
}

gsi.mul <- function( x,y ) {
  if( length(dim(x)) == 2 )
    if( length(dim(y)) == 2 )
      unclass(x)*unclass(y)
    else
      unclass(x)*rep(c(y),rep(NROW(x),length(y)))
  else if( length(dim(y)) == 2 )
      unclass(y)*rep(c(x),rep(NROW(y),length(x)))
  else
    unclass(x)*unclass(y)
}

gsi.div <- function( x,y ) {
  if( length(dim(x)) == 2 )
    if( length(dim(y)) == 2 )
      unclass(x)/unclass(y)
    else
      unclass(x)/rep(c(y),rep(NROW(x),length(y)))
  else if( length(dim(y)) == 2 )
      unclass(y)/rep(c(x),rep(NROW(y),length(x)))
  else
    unclass(x)/unclass(y)
}


power.acomp <- function(x,s) {
  if( is.acomp(s) || is.rcomp(s))
    stop("power.acomp is scalar multiplication only")
  if( !is.matrix(x) || nrow(x) ==1 ) {
    if( length(s)>1 )
      x <- matrix(x,byrow=TRUE,ncol=length(x),nrow=length(s))
  } else {
    if( length(s) > 1 && length(s)!= nrow(x) )
      warning("lengths do not match in power.acomp")
  }
  acomp(unclass(x)^c(s)) 
}


"+.acomp" <- function(x,y) {
  acomp(gsi.mul(x,y))
}

"-.acomp" <- function(x,y) {
  if( missing(y) )
    acomp(1/unclass(x))
  else 
    acomp(gsi.div(x,y))
}

"*.acomp" <- function(x,y) {
  if( is.acomp(x) && !is.acomp(y) )
    power.acomp(x,y)
  else if( is.acomp(y)&& !is.acomp(x) )
    power.acomp(y,x)
  else
    stop("the powertransform performed in *.acomp only operates on acomps and scalar")
}

"/.acomp" <- function(x,y) {
  if( is.acomp(x) && !is.acomp(y) )
    power.acomp(x,1/unclass(y))
  else
    stop("/.acomp only operates on acomp / numeric")
}

"+.aplus" <- function(x,y) {
    aplus(gsi.mul(x,y))
}

"-.aplus" <- function(x,y) {
  if( missing(y) )
    return(aplus(1/unclass(y)))
  else
    aplus( gsi.div(x,y) )
}

"*.aplus" <- function(x,y) {
  if( is.aplus(x)&& !is.aplus(y) )
    power.aplus(x,y)
  else if( is.aplus(y)&& !is.aplus(x) )
    power.aplus(y,x)
  else
    stop("*.aplus only operates on aplus and scalar")
}

"/.aplus" <- function(x,y) {
  if( is.aplus(x) && !is.aplus(y) )
    power.aplus(x,1/unclass(y))
  else
    stop("/.aplus only operates on aplus and scalar")
}


"+.rcomp" <- function(x,y) {
#  warning("+ is meaningless for rcomp")
  if( is.rcomp(x) )
    if( is.rcomp(y) ) {
      stop("+ is meaningless for two rcomp objects")
   } else if( is.rcomp(x)) {
      rmult(clo(x))+rmult(y)
    } else if( is.rcomp(y) ) {
      rmult(x)+rmult(clo(y))
    } else
      stop("Why are we here in +.rcomp without rcomp?")
  #rcomp(gsi.add(x,y))
}

"-.rcomp" <- function(x,y) {
  if( missing(y) )
    rmult(-unclass(x))
  else
    rmult(gsi.sub(x,y))
}

"*.rcomp" <- function(x,y) {
  if( is.rcomp(x) && is.rcomp(y) )
    rcomp(gsi.mul(x,y))
  else if( is.rcomp(x) )
    rplus(x)*y
  else if( is.rcomp(y) )
    rplus(y)*x
  else
    stop("undefined combination of arguments for *.rcomp")
}

"/.rcomp" <- function(x,y) {
  if( is.rcomp(x) && is.rcomp(y) )
    rcomp(gsi.div(x,y))
  else if( is.rcomp(x) )
    rplus(x)/y
  else
    stop("undefined combination of arguments for /.rcomp")
}

"+.rplus" <- function(x,y) {
  if( is.rplus(x) && is.rplus(y) )
    rplus(gsi.add(x,y))
  else
    rmult(gsi.add(x,y))
}

"-.rplus" <- function(x,y) {
  if( missing(y) )
    rmult(-unclass(x))
  else
    rmult(gsi.sub(x,y))
}


"*.rplus" <- function(x,y) {
  if( is.rplus(x) && is.rplus(y) )
    rplus(gsi.mul(x,y))
  else if( is.rplus(x) )
    mul.rplus(x,y)
  else if( is.rplus(y) )
    mul.rplus(y,x)
  else
    stop("undefined combination of arguments for *.rplus")
}

"/.rplus" <- function(x,y) {
  if( is.rplus(x) && is.rplus(y) )
    rplus(gsi.div(x,y))
  else if( is.rcomp(x) )
    mul.rplus(rplus(x),1/unclass(y))
  else
    stop("undefined combination of arguments for /.rcomp")
}

"+.rmult" <- function(x,y) {
  rmult(gsi.add(x,y))
}

"-.rmult" <- function(x,y) {
  if( missing(y) )
    rmult(-unclass(x))
  else
    rmult(gsi.sub(x,y))
}


"*.rmult" <- function(x,y) {
  if( is.rmult(x) && is.rmult(y) )
    rmult(gsi.mul(x,y))
  else
    rmult(unclass(x)*unclass(y))
}

"/.rmult" <- function(x,y) {
  if( is.rmult(x) && is.rmult(y) )
    rmult(gsi.div(x,y))
  else 
    rmult(unclass(x)/unclass(y))
}

"%*%" <- function(x,y) UseMethod("%*%",structure(c(),class=c(class(x),class(y))))


#gsi.internaltmp <- get("%*%",pos="package:base")
#formals(gsi.internaltmp) <- formals(get("%*%"))
#"%*%.default" <- gsi.internaltmp

"%*%.default" <- function(x,y) base::"%*%"(x,y)

"%*%.rmult" <- function(x,y) {
  if( is.rmult(y) )
    if( is.rmult(x) ) 
      c(gsi.mul(x,y) %*% rep(1,gsi.getD(x)))
    else if( is.matrix(x) ) 
      rmult(gsi.simshape(oneOrDataset(y) %*% t(x),y))
    else
      c(oneOrDataset(y) %*% x) 
  else if( is.matrix(y) )
      rmult(gsi.simshape(oneOrDataset(x) %*% y,x))
  else
      c(oneOrDataset(x) %*% y) 
  }

"%*%.acomp" <- function(x,y) {
  if( is.acomp(y) )
    if( is.acomp(x) ) 
      cdt(x) %*% cdt(y)
    else if( is.matrix(x) ) {
      if( nrow(x) == gsi.getD(y) )
        clrInv(x %*% clr(y))
      else
        ilrInv(x %*% ilr(y))
    }
    else
      stop( "%*%.acomp is only defined for special combinations I" )
  else if( is.acomp(x) ) {
    if( is.matrix(y) ) {
      if( ncol(y) == gsi.getD(x) )
        clrInv(clr(x) %*% y )
      else
        ilrInv(ilr(x) %*% y )
    }
  else
      stop( "%*%.acomp is only defined for special combinations II" )
  }
  else
      stop( "%*%.acomp is only defined for special combinations III" )
    
}

"%*%.aplus" <- function(x,y) {
  if( is.aplus(y) )
    if( is.aplus(x) ) 
      cdt(x) %*% cdt(y)
    else if( is.matrix(x) ) {
        iltInv(x %*% ilt(y))
    }
    else
      stop( "%*%.acomp is only defined for special combinations I" )
  else if( is.aplus(x) ) {
    if( is.matrix(y) ) {
        iltInv(ilt(x) %*% y )
    }
  else
      stop( "%*%.aplus is only defined for special combinations II" )
  }
  else
      stop( "%*%.aplus is only defined for special combinations III" )
    
}


convex.rcomp <- function(x,y,alpha=0.5) {
  rcomp( alpha*x + (1-alpha)*y )
}


mul.rplus <- function(x,r) {
  if( all(r>=0) )
    rplus(unclass(x)*r)
  else
    rmult(unclass(x)*r)
}

power.aplus <- function(x,r) {
  aplus(unclass(x)^r) 
}


gsi.expandrcomp <- function(x,alpha) {
  cptInv(cpt(x)*alpha)
}

endmemberCoordinates <- function(X,...) UseMethod("endmemberCoordinates")

endmemberCoordinates.default <- function(X,endmembers=diag(gsi.getD(X)),...) {
  class(endmembers) <- class(X)
  X <- oneOrDataset(idt(X))
  A <- t(unclass(idt(endmembers)))
  erg <- solve( rbind(cbind(t(A)%*%A,1),c(rep(1,ncol(A)),0)),
               rbind(t(A)%*%t(unclass(X)),1))
  erg <- rmult(t(erg[-nrow(erg),,drop=FALSE]))
  colnames(erg) <- rownames(endmembers)
  erg
}

endmemberCoordinates.acomp <- function(X,endmembers=clrInv(diag(gsi.getD(X))),...) {
  ep <- ilr(endmembers)
  rownames(ep) <- rownames(endmembers)
  endmemberCoordinates(idt(X),ep,...)
}

endmemberCoordinates.aplus <- function(X,endmembers,...) {
  ep <- ilt(endmembers)
  rownames(ep) <- rownames(endmembers)
  endmemberCoordinates(idt(X),ep,...)
}


endmemberCoordinates.rplus <- function(X,endmembers,...) {
  ep <- iit(endmembers)
  rownames(ep) <- rownames(endmembers)
  endmemberCoordinates(idt(X),ep,...)
}


endmemberCoordinatesInv <- function(K,endmembers,...) UseMethod("endmemberCoordinatesInv",endmembers)

endmemberCoordinatesInv.rmult <- function(K,endmembers,...) {
  rmult(t(t(unclass(endmembers)) %*% t(unclass(K))))
}

endmemberCoordinatesInv.acomp <- function(K,endmembers,...) {
  ilrInv(endmemberCoordinatesInv(K,ilr(endmembers)))
}

endmemberCoordinatesInv.rcomp <- function(K,endmembers,...) {
  iptInv(endmemberCoordinatesInv(K,ipt(endmembers)))
}


endmemberCoordinatesInv.aplus <- function(K,endmembers,...) {
  iltInv(endmemberCoordinatesInv(K,ilt(endmembers)))
}

endmemberCoordinatesInv.rplus <- function(K,endmembers,...) {
  iitInv(endmemberCoordinatesInv(K,iit(endmembers)))
}


formals(scale) <- c(formals(scale),alist(...= ))
formals(scale.default) <- c(formals(scale.default),alist(...= ))

scale.aplus  <- scale.acomp <- function( x,center=TRUE, scale=TRUE ,...,robust=getOption("robust")) {
  if( ! (center || scale ) ) return(x)
  va <- var(x,robust=robust,giveCenter=TRUE)
  ce <- attr(va,"center")
  if( is.logical(center) ) {
    if( center )
      x <- x-ce
  } else x <- x-center
  if( is.logical(scale) ) {
    if(scale)
      x <- (1/sqrt(mean(gsi.diagExtract(va))))*x
  } else x <- x*scale
  return(x)
#  W <- x
#  if( center ) {
#    W <- clrInv( scale(clr(W),center=center,scale=FALSE) )
#    if( scale )
#      W <- power.acomp(W,as.numeric(scale)/
#                       sqrt(sum(gsi.diagExtract(var(clr(W))))))
#  } else if( scale ) {
#    mean <- c(mean(acomp(W),robust=robust))
#    W <- perturbe(power.acomp(perturbe(W,1/mean),as.numeric(scale)/sqrt(sum(gsi.diagExtract(var(clr(W)))))),mean)
#  }
#  W
}

#scale.rcomp <- function( x,center=TRUE, scale=TRUE ) {
#  W <- x
#  if( center ) {
#    W <- cptInv( scale(cpt(W),center=center,scale=FALSE) )
#    if( scale )
#      W <- gsi.expandrcomp(W,as.numeric(scale)/sqrt(sum(gsi.diagExtract(var(cpt(W))))))
#  } else if( scale ) {
#    mean <- c(mean(rcomp(W)))
#    W <- gsi.add(mean,gsi.sub(W,mean)/sqrt(sum(gsi.diagExtract(var(cpt(W))))))
#  }
#  W
#}

#scale.aplus <- function( x,center=TRUE, scale=TRUE ) {
#  W <- x
#  if( center ) {
#    W <- iltInv( scale(ilt(W),center=center,scale=FALSE) )
#    if( scale )
#      W <- power.aplus(W,as.numeric(scale)/sqrt(sum(gsi.diagExtract(var(ilt(W))))))
#  } else if( scale ) {
#    mean <- c(mean(aplus(W)))
#    W <- perturbe.aplus(power.aplus(perturbe.aplus(W,1/mean),as.numeric(scale)/sqrt(sum(gsi.diagExtract(var(ilt(W)))))),mean)
#  }
#  W
#}

#scale.rplus <- function( x,center=TRUE, scale=TRUE ) {
#   rmult(scale(gsi.plain(x),center=center,scale=scale))
#}

scale.rcomp <- scale.rplus <- scale.rmult <- function( x,center=TRUE, scale=TRUE ,..., robust=getOption("robust")) {
  if( ! (center || scale ) ) return(x)
  var <- var(x,robust=robust,giveCenter=TRUE)
  ce <- attr(var,"center")
  if( is.logical(center) ) {
    if( center )
      x <- x-ce
  } else x <- x-center
  if( is.logical(scale) ){
    if(scale) {
      s <- gsi.diagGenerate(1/sqrt(gsi.diagExtract(var)))
      x <- s %*% x
    }
  } else if( is.matrix(s) )
     x <- s %*% x
  else if( length(s)==1)
     x <- s * x
  else
     x <- gsi.diagGenerate(s) %*% x
  return(x)
   #rmult(scale(gsi.plain(x),center=center,scale=scale))
}

normalize <- function(x,...) UseMethod("normalize",x)
normalize.default <- function(x,...) x/norm(x)

#if( !exists("norm")) norm <- function(X,...) UseMethod("norm",X)

norm.default <- function(X,...) {
  sqrt( sum(X^2) )
}

norm.acomp <- function(X,...) {
  norm.rmult(cdt(X),...)
}
norm.rcomp <- norm.acomp
norm.aplus <- norm.acomp
norm.rplus <- norm.acomp
norm.rmult <- function(X,...) {
   sqrt(X %*% X)
}

dist <- function(x,...) UseMethod("dist")
dist.default <- function(x,...) stats::dist(cdt(x),...)


scalar <- function(x,y) UseMethod("scalar")

scalar.default <- function(x,y) {
  x <- cdt(x)
  y <- cdt(y)
  tmp <- gsi.mul(oneOrDataset(x,y), oneOrDataset(y,x)) 
  c( tmp %*% rep(1,NCOL(tmp)))
}


clr <- function( x ,... ) {
  W <- oneOrDataset(x)
  nmv <- is.NMV(W)
  LOG <- unclass(log(ifelse(nmv,W,1)))
  erg <- ifelse(nmv,LOG-gsi.rsum(LOG)/gsi.rsum(nmv),0)
  #M   <- gsi.rsum(gsi.recodeM2C(W,LOG,BDL=0,SZ=0,MAR=0,MNAR=0))
  rmult(gsi.simshape(erg,x),orig=x) 
}

clrInv <- function( z,... ) {
  acomp( gsi.recodeC2M(exp(z),ninf=BDLvalue,nan=MARvalue,na=MNARvalue) )
}

ult <- function( x,... ) {
  ilt(clo(x),...)
}

ultInv <- clrInv

Kappa <- function( x, ...) {
  W <- oneOrDataset(x)
  (clr(W)-ult(W))[,1]
}

#gsi.ilrBase <- function(D) {
#  if( D==1 )
#    return(matrix(nrow=0,ncol=0))
#  tmp <- diag(D) - 1/(D)* matrix(1,ncol=D,nrow=D)
#  for( i in 1:(NCOL(tmp)-1)  ) {
#    tmp[,i] <- tmp[,i]/sqrt(sum(tmp[,i]^2))
#    rest <- (i+1):NCOL(tmp)
#    if( length(rest) != 1 ) {
#      tmp[,rest]  <-tmp[,rest,drop=FALSE] - tmp[,rep(i,length(rest)),drop=FALSE]%*%
#        gsi.diagGenerate( c(t(tmp[,i])%*%tmp[,rest,drop=FALSE] ) )
#    } 
#  }
#  tmp[,-NROW(tmp)]
#}

gsi.ilrBase <-function(D){
  if(D<=1)
    return(matrix(nrow=0,ncol=0))
  else 
    t(unclass(normalize(rmult(t(contr.helmert(n=D))))))
}




                                        # Old ilrBase
#ilrBase <- function( x=NULL , z=NULL , D = NULL ) {
#  if( missing(D) )
#    D <- if(is.null(x))
#      NCOL(oneOrDataset(z))+1
#    else
#      NCOL(oneOrDataset(x))
#  while( D > length(ilrBaseList) )
#    ilrBaseList <<- c(ilrBaseList,gsi.ilrBase(length(ilrBaseList)+1))
#  ilrBaseList[[D]]gsi.OrderIlr
#}
# ilrBaseList <- lapply(1:20,gsi.ilrBase)


ilrBase <- function(x=NULL, z=NULL, D=NULL, method="basic"){
 if (method=="basic"){
    if (missing(D))
        D <- if (is.null(x))
            NCOL(oneOrDataset(z)) + 1
        else NCOL(oneOrDataset(x))
#     while (D > length(ilrBaseList)) ilrBaseList <<- c(ilrBaseList,
#         gsi.ilrBase(length(ilrBaseList) + 1))
#     V=ilrBaseList[[D]]
  V = gsi.ilrBase(D)
 } #end if basic
 if (method=="balanced"){
    # build the merge structure (as in hclust)
    if(is.null(D)==TRUE){D=max(ncol(x),ncol(z))}
    M = c(-c(1:D),c(floor(D/2):1),c((floor(D/2)+1):(D-2)))
    M = matrix(M, byrow=TRUE, nrow=D-1, ncol=2)
    V = gsi.merge2signary(M)
    V = gsi.buildilrBase(V)
 }#end if balanced
 if (method=="optimal"){
  if(length(dim(x))<2){
   warning("method 'optimal' requires a data set")
  }
  else{
    V = gsi.optimalilrBase(x)
    V = gsi.buildilrBase(V)

  }
 } #end if optimal
 if(length(grep(pattern="PB",x=method))!=0){
   if(length(dim(x))<2){
     stop("all 'Principal Balances' methods require a data set")
   }
   gsi.PrinBal(x, method)
 } # end methods of principal balances (all having string "PB" in their name)
 return(V)
}

gsi.buildilrBase <-function(W=c(1,-1)){
  # builds an ilr base from a matrix of 1, -1 and 0 (a partition)
  if(length(W)<2){
    return(ilrBase(D=1))
  }
if(length(dim(W))==0){
  return(ilrBase(D=2))
  }
  if(length(dim(W))>0){
   W = as.matrix(W)
   nc = ncol(W)
   D = nrow(W)
   isPos = (W>0)
   isNeg = (W<0)
   nPos = matrix(1,D,D) %*% isPos
   nNeg = matrix(1,D,D) %*% isNeg
   W = (isPos * nNeg - isNeg * nPos)
   nn = sapply(1:nc,function(i){1/norm.rmult(W[,i])})
   nn = matrix(nn, ncol=ncol(W), nrow=nrow(W), byrow=TRUE)
   W = W * nn
   return(W)
  }
}
gsi.signary2ilrBase <- gsi.buildilrBase


gsi.optimalilrBase <- function(x){
  # Takes as data a binary data frame: 0=case unobserved, 1=case observed
  # performs a cluster analysis of variables on it
  # and recodes the structure to a signary matrix (pre-ilr)
  if(is.null(attr(x,"Losts"))==TRUE){
    Ones = !oneOrDataset(is.infinite(log(as.matrix(x))),x) # find zeroes
  }
  if(is.null(attr(x,"Losts"))==FALSE){
    Ones = !attr(x,"Losts")
  }
  h= hclust(dist(t(Ones)))
  V = gsi.merge2signary(h$merge)
  return(V)
}


gsi.merge2signary <- function(M){
  # takes a merge structure (as explained in hclust) and converts it to a sign matrix (encoding the ilr partition: 0=no influence, 1=numerator, -1=denominator)
 V = matrix(0,ncol=nrow(M)+1,nrow=nrow(M))
  for(i in 1:nrow(M)){
   for(j in 1:2){
    weight = (-1)^j
    k = M[i,j]
    if(k<0){
      V[i,abs(k)] = weight
    } # for singletons
    if(k>0){ # for groups
      take = as.logical(V[k,])
      V[i,take] = rep(weight,sum(take))
    }
  }
 }
return(t(V))
}




ilr    <- function( x , V=ilrBase(x),... ) {
  rmult(clr2ilr( clr(oneOrDataset(x)),V ))
}

ilrInv <- function( z, V=ilrBase(z=z),...,orig=NULL) {
  erg <- clrInv( ilr2clr(z,V) )
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}

alr <- function( x ,ivar=ncol(x),...) {
xo <- x
W <- unclass(clo(oneOrDataset(x)))
x <- gsi.recodeM2C(W,log(W),BDL=-Inf,SZ=NaN,MAR=NaN,MNAR=NA)
rmult(gsi.simshape( x[,-ivar,drop=FALSE] - c(x[,ivar]), xo))
}

alrInv <- function( z ,...,orig=NULL) {
  Z <- cbind(oneOrDataset(z),0)
  erg <- acomp(gsi.simshape( clo(gsi.recodeC2M(Z,exp(Z),
                                        ninf=BDLvalue,
                                        inf =MNARvalue,
                                        nan =MARvalue,
                                        na  =MNARvalue
                                        )) , z ))
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}


apt <- function( x ,...) {
  W <- oneOrDataset(x)
  V <- gsi.recodeM2C(W,gsi.plain(clo( W )),BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA)
  rmult(gsi.simshape(V[,-NCOL(W)],x))
}

aptInv <- function( z ,...,orig=NULL) {
  Z <- oneOrDataset(z)
  Z <- cbind(Z, 1 - gsi.recodeC2M(Z,na=0.0,nan=0.0) %*% rep(1,NCOL(Z)))
  erg <- rcomp(gsi.simshape( Z ,z ))
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}

cpt <- function( x ,...) {
  X <- oneOrDataset(x)
  rmult(clo(x) - 1/NCOL(X),orig=x)
}

cptInv <- function( z ,...) {
  if( abs(sum(z))>0.0001 )
    warning( "z not in cpt plane in cptInv")
  rcomp(z + 1/NCOL(oneOrDataset(z)))
}

ipt    <- function( x , V=ilrBase(x),...) {
  rmult(clr2ilr(cpt(x),V),orig=x)
  
}

iptInv <- function( z, V=ilrBase(z=z) ,...,orig=NULL) {
  erg<-cptInv( ilr2clr(z,V) )
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}

uciptInv <- function( z, V=ilrBase(z=z) ,...,orig=NULL) {
  tmp <- ilr2clr(z,V) + 1/(gsi.getD(z)+1)
  tmp[tmp<0]<-MNARvalue
  erg<-rcomp(tmp)
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}


ilt <- function( x ,...) {
  W <- gsi.plain(x)
  rmult(log(ifelse(is.NMV(W),W,1)),orig=x)
}

iltInv <- function( z ,...) {
  aplus(gsi.recodeC2M(z,exp(z),ninf=BDLvalue,nan=MARvalue,na=MNARvalue))
}

iit <- function( x ,...) {
  rmult( gsi.recodeM2C(x,BDL=0.0,SZ=0.0,MAR=0.0,MNAR=0.0 ) ,orig=x)
}

iitInv <- function(z,...) {
  rplus(gsi.recodeC2M(z,na=MNARvalue,nan=MARvalue))
}

idt         <- function(x,...) UseMethod("idt",x)
idt.default <- function(x,...) x
idt.acomp   <- function(x,...) ilr(x,...) 
idt.rcomp   <- function(x,...) ipt(x,...)
idt.ccomp   <- iit
idt.aplus   <- ilt 
idt.rplus   <- iit 
idt.rmult   <- function(x,...) x
idt.factor  <- function(x,...) rmult(clr2ilr(cdt(x)))


cdt         <- function(x,...) UseMethod("cdt",x)
cdt.default <- function(x,...) x
cdt.acomp   <- clr 
cdt.rcomp   <- cpt
cdt.ccomp   <- iit
cdt.aplus   <- ilt 
cdt.rplus   <- iit 
cdt.rmult   <- function(x,...) x
cdt.factor  <- function(x,...) {
  #x <- matrix(0,nrow=length(x),ncol=nlevels(x),dimnames=list(names(x),levels(x)))
  x[1:ncol(x)+unclass(x)] <- model.matrix(~-1+x)
  
  rmult(matrix(x,nrow=nrow(x),dimnames=dimnames(x)))
}

cdtInv <- function(x,orig,...) UseMethod("cdtInv",orig)
cdtInv.default <- function(x,orig,...) x
cdtInv.acomp   <- function(x,orig,...) clrInv(x,...,orig=orig)
cdtInv.rcomp   <- function(x,orig,...) cptInv(x,...,orig=orig)
cdtInv.ccomp   <- function(x,orig,...) iitInv(x,...,orig=orig)
cdtInv.aplus   <- function(x,orig,...) iltInv(x,...,orig=orig)
cdtInv.rplus   <- function(x,orig,...) iitInv(x,...,orig=orig)
cdtInv.rmult   <- function(x,orig,...) x

idtInv <- function(x,orig,...) UseMethod("idtInv",orig)
idtInv.default <- function(x,orig,...) x
idtInv.acomp   <- function(x,orig,...) ilrInv(x,...,orig=orig)
idtInv.rcomp   <- function(x,orig,...) iptInv(x,...,orig=orig)
idtInv.ccomp   <- function(x,orig,...) iitInv(x,...,orig=orig)
idtInv.aplus   <- function(x,orig,...) iltInv(x,...,orig=orig)
idtInv.rplus   <- function(x,orig,...) iitInv(x,...,orig=orig)
idtInv.rmult   <- function(x,orig,...) x


variation <- function( x, ... ) UseMethod("variation",x)

variation.acomp <- function( x,...,robust=getOption("robust") ) {
  co <-var(x,robust=robust)
  d <- NCOL(x)
  va <-gsi.diagExtract(co)
  co1 <- matrix(rep(va,each=d),ncol=d)
  co2 <- matrix(rep(va,d),ncol=d)
  -2*co+co1+co2
  
}

variation.rcomp <- function( x ,...,robust=getOption("robust")) {
  co <-var(x,robust=robust)
  d <- NCOL(x)
  va <-gsi.diagExtract(co)
  co1 <- matrix(rep(va,each=d),ncol=d)
  co2 <- matrix(rep(va,d),ncol=d)
  -2*co+co1+co2
  
}


variation.aplus <- function( x ,...,robust=getOption("robust")) {
  co <-var(x,robust=robust)
  d <- NCOL(x)
  va <-gsi.diagExtract(co)
  co1 <- matrix(rep(va,each=d),ncol=d)
  co2 <- matrix(rep(va,d),ncol=d)
  -2*co+co1+co2
  
}

variation.rmult <- function( x ,...,robust=getOption("robust")) {
  co <-var(x,robust=robust)
  d <- NCOL(x)
  va <-gsi.diagExtract(co)
  co1 <- matrix(rep(va,each=d),ncol=d)
  co2 <- matrix(rep(va,d),ncol=d)
  -2*co+co1+co2
  
}
variation.rplus <- variation.rmult


gsi.mapin01 <- function(i,min=0,max=1) {c(min,min+(max-min)/i,max)}
gsi.mapfrom01 <- function(x) {(x[3]-x[1])/(x[2]-x[1])}
gsi.mapmin <- function(x) {x[1]}
gsi.mapmax <- function(x) {x[3]}

#gsi.plots <- list()
#gsi.setPlot <- function(x) {
#  gsi.setPlot[[dev.cur()]]<-x
#}
#gsi.getPlot <- function() {
#  gsi,plots[[dev.cur()]]
#}

gsi.plots <- new.env()

gsi.setPlot <- function(x,what="XXXplots",dev=dev.cur()) {
  y<-NULL
  try(
  y <- get(what,gsi.plots),silent=TRUE
      )
  if( is.null(y))
    y <- list()
  y[[dev]]<-x
  assign(what,y,gsi.plots)
}

gsi.getPlot <- function(what="XXXplots",dev=dev.cur()) {
  y<-NULL
  try(
      y <- get(what,gsi.plots),silent=TRUE
      )
  if( is.null(y))
    y <- list()
  
  #x <- get(what,gsi.plots)
  if( dev <= length(y) )
     y[[dev]]
  else
     NULL
}

#gsi.coorInfo <- list()
#gsi.setCoorInfo <- function(...) {
#  par()
#  gsi.coorInfo[[dev.cur()]] <<- list(...)
#}
#gsi.getCoorInfo <- function() {
#  if( dev.cur() <= length(gsi.coorInfo))
#    gsi.coorInfo[[dev.cur()]]
#  else
#    NULL
#}

gsi.setCoorInfo <- function(...,all=list(...)) {
  par()
  gsi.setPlot(all,what="gsi.coorInfo")
  #gsi.coorInfo[[dev.cur()]] <<- list(...)
}
gsi.addCoorInfo <- function(...,all=list(...)) {
  coorinfo <- gsi.getPlot(what="gsi.coorInfo")
  coorinfo[names(all)]<-all
  gsi.setCoorInfo(all=coorinfo)
#  if( dev.cur() <= length(gsi.coorInfo))
#    gsi.coorInfo[[dev.cur()]]
#  else
#    NULL
}


gsi.getCoorInfo <- function() {
  gsi.getPlot(what="gsi.coorInfo")
#  if( dev.cur() <= length(gsi.coorInfo))
#    gsi.coorInfo[[dev.cur()]]
#  else
#    NULL
}


gsi.call  <- function(fkt,...) {
  if( is.character(fkt) )
    do.call(fkt,list(...))
  else
    fkt(...)
}

gsi.add2pairs <- function(x,panel,...,noplot=FALSE) {
#  if( dev.cur() <= length(gsi.plots) )
#    curplot <- gsi.plots[[dev.cur()]]
#  else
#    curplot <- NULL
  curplot <- gsi.getPlot()
  if( is.null(curplot) ) {
    panel(x[,1],x[,2],...)
  } else {
    if( !missing(panel) )
      curplot$add <- c(curplot$add,list(list(x=x,panel=panel,args=list(...))))
    gsi.setPlot(curplot)
    if(!noplot) do.call("gsi.pairs",curplot)
  }
}



##### OldVersion
#gsi.pairs <- function (x, labels, panel = points, ..., main = NULL, oma = NULL, 
#    font.main = par("font.main"), cex.main = par("cex.main"), 
#    lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
#    text.panel = textPanel, label.pos = 0.5 + has.diag/3, cex.labels = NULL, 
#    font.labels = 1, row1attop = TRUE, gap = 1,add=list(),xlim=apply(x,2,range),ylim=apply(x,2,range),log="",noplot=FALSE) 
#{
#    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
#        y, txt, cex = cex, font = font)
#    localAxis <- function(side, xpd, bg, ...) axis(side, xpd = NA, 
#        ...)
#    if (!is.matrix(x)) 
#        x <- data.matrix(x)
#    if (!is.numeric(x)) 
#        stop("non-numeric argument to pairs")
#    panel <- match.fun(panel)
#    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
#        lower.panel <- match.fun(lower.panel)
#    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
#        upper.panel <- match.fun(upper.panel)
#    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
#        diag.panel <- match.fun(diag.panel)
#    if (row1attop) {
#        tmp <- lower.panel
#        lower.panel <- upper.panel
#        upper.panel <- tmp
#        tmp <- has.lower
#        has.lower <- has.upper
#        has.upper <- tmp
#    }
#    nc <- ncol(x)
#    if (nc < 2) 
#        stop("only one column in the argument to gsi.pairs")
#    has.labs <- TRUE
#    if (missing(labels)) {
#        labels <- colnames(x)
#        if (is.null(labels)) 
#            labels <- paste("var", 1:nc)
#    }
#    else if (is.null(labels)) 
#        has.labs <- FALSE
#    if( length(dim(xlim)) < 2 ) xlim <- matrix(rep(xlim,ncol(x)),nrow=2)  
#    if( length(dim(ylim)) < 2 ) ylim <- matrix(rep(ylim,ncol(x)),nrow=2)  
#    if (is.null(oma)) {
#        oma <- c(4, 4, 4, 4)
#        if (!is.null(main)) 
#            oma[3] <- 6
#    }

#    mycall <- list(x=x,labels=labels,panel=panel,...,main=main,oma=oma,
#                   font.main=font.main,cex.main=cex.main,
#                   lower.panel=lower.panel, upper.panel=upper.panel,
#                   diag.panel=diag.panel,text.panel=text.panel,
#                   label.pos=label.pos,cex.labels=cex.labels,
#                   font.labels=font.labels,row1attop=row1attop,gap=gap,add=add,
#                   xlim=xlim,ylim=ylim,log=log)
#    if( noplot ) {
#      gsi.plots[[dev.cur()]] <<- mycall
#      return(invisible(NULL))
#    }
#    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
#    on.exit(par(opar))
#    for (i in if (row1attop) 
#        1:nc
#    else nc:1) for (j in 1:nc) {
#        plot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
#            type = "n", ...,log=ifelse(i==j,"",log),xlim=xlim[,j],ylim=ylim[,i])
#        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
#            box()
#            if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
#                localAxis(1 + 2 * row1attop, ...)
#            if (i == nc && (j%%2 || !has.upper || !has.lower)) 
#                localAxis(3 - 2 * row1attop, ...)
#            if (j == 1 && (!(i%%2) || !has.upper || !has.lower)) 
#                localAxis(2, ...)
#            if (j == nc && (i%%2 || !has.upper || !has.lower)) 
#                localAxis(4, ...)
#            mfg <- par("mfg")
#            if (i == j) {
#                if (has.diag) 
#                  diag.panel(as.vector(x[, i]))
#                if (has.labs) {
#                  par(usr = c(0, 1, 0, 1))
#                  if (is.null(cex.labels)) {
#                    l.wid <- strwidth(labels, "user")
#                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
#                  }
#                  text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
#                    font = font.labels)
#                }
#            }
#            else if (i < j) 
#                lower.panel(as.vector(x[, j]), as.vector(x[, 
#                  i]), ...)
#           else upper.panel(as.vector(x[, j]), as.vector(x[, 
#                i]), ...)
#            if( i!=j ) for( p in add ) {
#              arg <- c(list(p$panel,as.vector(p$x[,j]), as.vector(p$x[,i])),p$args)
#              do.call("gsi.call",arg)
#            }
#            if (any(par("mfg") != mfg)) 
#                stop("The panel function made a new plot")
#        }
#        else par(new = FALSE)
#    }
#    if (!is.null(main)) 
#        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
#    gsi.plots[[dev.cur()]] <<- mycall
#    invisible(NULL)
#}

noreplot <- function(expr,dev=dev.cur()) {
  curplot  <- gsi.getPlot(what="...replot",dev=dev)
  if( is.null(curplot) )
    gsi.setPlot(quote(list()),what="...replot",dev=dev)
  on.exit(gsi.setPlot(curplot,what="...replot",dev=dev),add=TRUE)
  is.null(erg <- expr)
  invisible(erg)
}

replot <- function(...,dev=dev.cur(),plot=TRUE,envir=NULL,add=FALSE) {
  if( !is.logical(plot) ) {
    if( !is.list(plot)) {
      if( !is.null(envir) )
        environment(plot)<-envir
      if( is.null(environment(plot)) )
        environment(plot)<- parent.frame(2)
    }
    if( is.logical(add) & !add ) {
      if( is.call(plot) )
        gsi.setPlot(list(main=plot),what="...replot",dev=dev)
      else 
        gsi.setPlot(plot,what="...replot",dev=dev)
      curplot <- gsi.getPlot(what="...replot",dev=dev)
    } else {
      curplot <- gsi.getPlot(what="...replot",dev=dev)
      if( is.logical(add) )
        add=length(curplot)+1
      curplot[[add]]<-plot
      gsi.setPlot(curplot,what="...replot",dev=dev)
    }
    return(invisible(curplot))
  } 
  lm <- list(...)
  if( is.logical(add) )
    if( !add ) add <- 1 else stop("Add needs to specify an additional plot") 
  curplot <- gsi.getPlot(what="...replot",dev=dev)
  if( is.null(curplot) ) {
    if( plot==FALSE && length(lm)==0 )
      return(curplot)
    stop("showOnlyPanel: The active plot does not support replot")
  }
  if(!is.null(envir))
    envir <- environment(curplot[[add]])
  if( length(lm) > 0 ) 
    curplot[[add]][names(lm)]<-lm
  if(!is.null(envir))
    environment(curplot[[add]])<-envir
  gsi.setPlot(curplot,what="...replot",dev=dev)
  if( plot ) {
    lapply(curplot,function(e) eval(e,environment(e)))
  }
  return(invisible(curplot))
}

replotable <- function(expr,add=FALSE) {
  replot(plot=substitute(expr),add=add)
  invisible(noreplot(expr))
}


#### New version.
#### New version.
# modified by Raimon on 2008-05-23:, and again on 2008-06-23
#    when margin = part, the diagonal name panel were not responding adequately
#    see lines added after "if (i == j)"
gsi.pairs <- function (x, labels, panel = points, ..., main = NULL, oma = NULL, 
    font.main = par("font.main"), cex.main = par("cex.main"), 
    lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
    text.panel = textPanel, label.pos = 0.5 + has.diag/3, cex.labels = NULL, 
    font.labels = 1, row1attop = TRUE, gap = 1,add=list(),xlim=apply(x,2,range),ylim=apply(x,2,range),log="",onlyPanel=NULL,noplot=FALSE,trimode=FALSE) 
{
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) {
      x<-0.5
      y<-0.5
      usr <- par("usr")
      if( par("xlog") )
        X <- 10^(usr[1]*(1-x)+usr[2]*x)
      else
        X <- usr[1]*(1-x)+usr[2]*x
      if( par("ylog") )
        Y <- 10^(usr[3]*(1-y)+usr[4]*y)
      else  
        Y <- usr[3]*(1-y)+usr[4]*y
      text(X,Y, txt, cex = cex, font = font)
    }
    localAxis <- function(side, xpd, bg, ...) axis(side, xpd = NA, 
        ...)
    if (!is.matrix(x)) 
        x <- data.matrix(x)
    if (!is.numeric(x)) 
        stop("non-numeric argument to pairs")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
        stop("only one column in the argument to gsi.pairs")
    has.labs <- TRUE
    if (missing(labels)) {
        labels <- colnames(x)
        if (is.null(labels)) 
            labels <- paste("var", 1:nc)
    }
    else if (is.null(labels)) 
        has.labs <- FALSE
    if( length(dim(xlim)) < 2 ) xlim <- matrix(rep(xlim,ncol(x)),nrow=2)  
    if( length(dim(ylim)) < 2 ) ylim <- matrix(rep(ylim,ncol(x)),nrow=2)  
    if (is.null(oma)) {
      if( trimode ) {
        oma <- c(3, 3, 3, 3)
        if (!is.null(main)) 
          oma[3] <- 5
      } else {
        oma <- c(4, 4, 4, 4)
        if (!is.null(main)) 
          oma[3] <- 6
      }
    }

    mycall <- list(x=x,labels=labels,panel=panel,...,main=main,oma=oma,
                   font.main=font.main,cex.main=cex.main,
                   lower.panel=lower.panel, upper.panel=upper.panel,
                   diag.panel=diag.panel,text.panel=text.panel,
                   label.pos=label.pos,cex.labels=cex.labels,
                   font.labels=font.labels,row1attop=row1attop,gap=gap,add=add,
                   xlim=xlim,ylim=ylim,log=log,onlyPanel=onlyPanel,
                   trimode=trimode)
    if( noplot ) {
      gsi.setPlot(mycall)
      #gsi.plots[[dev.cur()]] <<- mycall
      return(invisible(NULL))
    }
    
    if( is.null(onlyPanel) )
      onlyPanel <- list(1:nc,1:nc)
    if( length(onlyPanel)!=2)
        stop("onlyPanel must be a list of two")
    if( trimode ) {
      opar <- par(mfrow = c(length(onlyPanel[[1]]), length(onlyPanel[[2]])),
                  mar = rep.int(gap/2, 4)+c(1,0,0,0), oma = oma,pty="s",xaxs="i",yaxs="i")
    } else {
      opar <- par(mfrow = c(length(onlyPanel[[1]]), length(onlyPanel[[2]])),
                  mar = rep.int(gap/2, 4), oma = oma)
    }
#    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit({par(opar);gsi.addCoorInfo(d=NULL)},add=TRUE)
    for (iP in if (row1attop) 
        1:length(onlyPanel[[1]])
    else length(onlyPanel[[1]]):1) for (jP in 1:length(onlyPanel[[2]])) {
      i <- onlyPanel[[1]][iP]
      j <- onlyPanel[[2]][jP]
      gsi.addCoorInfo(d=c(i,j))
      if( trimode ){
        plot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
            type = "n", ...,log=log,xlim=c(0,1),ylim=c(0,1),usr=c(0,1,0,1))
      } else {
        plot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
            type = "n", ...,log=log,xlim=xlim[,j],ylim=ylim[,i])
      }
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            if(trimode) {
              if(i!=j) {
                #localAxis(1, ...)
              }
            } else {
              box()
              if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
                localAxis(1 + 2 * row1attop, ...)
              if (i == nc && (j%%2 || !has.upper || !has.lower)) 
                localAxis(3 - 2 * row1attop, ...)
              if (j == 1 && (!(i%%2) || !has.upper || !has.lower)) 
                localAxis(2, ...)
              if (j == nc && (i%%2 || !has.upper || !has.lower)) 
                localAxis(4, ...)
            }
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag) 
                  diag.panel(as.vector(x[, i]))
                if (has.labs) {
 #                 par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
 #                   cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                    cex.labels <- max(1, min(2, 0.9/max(l.wid)))
                  }
                  # to ensure that diagonal panels respond when margin=part
                  if(trimode){
                   mrg = gsi.getCoorInfo()$margin
                   classes = c("acomp","rcomp")
                   k = i + ifelse(mrg %in% classes, 0, ifelse( i< match(mrg,labels),0,1))
                   text.panel(0.5, label.pos, labels[k], cex = cex.labels, 
                     font = font.labels)
                  }else{
                   text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                    font = font.labels)
                  }
                }
            }
            else if (i < j) 
                lower.panel(as.vector(x[, j]), as.vector(x[, 
                  i]), ...)
            else upper.panel(as.vector(x[, j]), as.vector(x[, 
                i]), ...)
            if( i!=j ) for( p in add ) {
              arg <- c(list(p$panel,as.vector(p$x[,j]), as.vector(p$x[,i])),p$args)
              do.call("gsi.call",arg)
            }
            if (any(par("mfg") != mfg)) 
                stop("The panel function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) 
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    #gsi.plots[[dev.cur()]] <<- mycall
    gsi.setPlot(mycall)
    invisible(NULL)
}






simpleMissingSubplot <- function(loc,
                            portions,labels=NULL,
                            col=c("white","yellow","red","green","blue"),
                                 ...,border="gray60",vertical=NULL,xpd=NA) {
  x1<-loc[1]
  x2<-loc[2]
  y1<-loc[3]
  y2<-loc[4]
  opar <- par(list("xlog","ylog","xpd","pin","plt"))
  usr<-par("usr")
  #par(opar[c("xlog","ylog","xpd","pin","plt")])
  on.exit(par(opar),add=TRUE)
  par(xlog=FALSE,ylog=FALSE,xpd=xpd,par(),opar)
  if( is.null(vertical) ) {
    vertical <- opar$pin[1]*abs(x2-x1)/abs(usr[2]-usr[1]) < 
      opar$pin[2]*abs(y2-y1)/abs(usr[4]-usr[3])
    if( is.na(vertical) ) {
      print(list(opar=opar,usr=usr,args=match.call(),loc=loc))
    }
  }
  if( length(col)<length(portions) ) col <- rep(col,length(portions))
  col <- col[1:length(portions)]
  col <- col[portions>0]
  labels <- labels[1:length(portions)]
  labels<- labels[portions>0]
  portions<-portions[portions>0]
  portions <- portions/sum(portions)  
  p <- cumsum(c(0,portions))
  mp <- (p[-1]+p[-length(p)])/2
  xd <- x2-x1
  yd <- y2-y1
  if( vertical ) {
    xs <- (x2+x1)/2
    for(i in 1:length(portions))
      rect(x1,y1+yd*p[i],x2,y1+yd*p[i+1],col=col[i],border=border)
    if(!is.null(labels) ) text(xs,y1+mp*yd,labels,...,xpd=xpd)
  } else {
    ys <- (y2+y1)/2
    for(i in 1:length(portions))
      rect(x1+xd*p[i],y1,x1+xd*p[i+1],y2,col=col[i],border=border)
    if(!is.null(labels) ) text(x1+mp*xd,ys,labels,...,xpd=xpd)
  }
  invisible(mp)
}


ternaryAxis <- function(side=1:3,at=seq(0.2,0.8,by=0.2),
                        labels=if(is.list(at)) lapply(at,format) else format(at),
                        ...,
                        tick=TRUE,pos=0,
                        font.axis=par("font.axis"),
                        font.lab=par("font.lab"),
                        lty="solid",lwd=1,
                        len.tck=0.025,dist.lab=0.03,
                        dist.axis=0.03,
                        lty.tck="solid",
                        col.axis=par("col.axis"),
                        col.lab=par("col.lab"),
                        cex.axis=par("cex.axis"),
                        cex.lab=par("cex.lab"),
                        Xlab=NULL,Ylab=NULL,Zlab=NULL,small=TRUE,
                        xpd=NA,aspanel=FALSE){
  #print(match.call())
  nr <-1
  it <- function(x,Nr=nr) {
    if( length(x) == 0 )
      return(x)
    x[[((Nr-1)%%length(x))+1]]
  }
  itl <- function(x,Nr=nr) {
    if( length(x) == 0 )
      return(x)
    if( is.list(x) ) 
      x[[((Nr-1)%%length(x))+1]]
    else
      x
  }
  if( ! aspanel ) {
    axes <- match.call()
    axes$aspanel=TRUE
    environment(axes) <- parent.frame()
    replot(axes=substitute(quote(axes),list(axes=axes)))
    return(invisible(NULL))
  }
  #if (is.null(col) && length(list(...)) && !is.null(fg <- list(...)$fg)) 
  #  col <- fg
  s60 <- sin(pi/3)
  c60 <- cos(pi/3)
  s30 <- sin(pi/6)
  c30 <- cos(pi/6)
  if( small ) {
    if( !is.null(Xlab) )
      text(0,-it(dist.lab,1),Xlab,adj=c(0.5,1),font=it(font.lab,1),col=it(col.lab,1),cex=it(cex.lab,1),xpd=xpd)
    if( !is.null(Ylab) )
      text(1,-it(dist.lab,2),Ylab,adj=c(0.5,1),font=it(font.lab,2),col=it(col.lab,2),cex=it(cex.lab,2),xpd=xpd)
    if( !is.null(Zlab) )
      text(0.5,s60+it(dist.lab,3),Zlab,adj=c(0.5,0),font=it(font.lab,3),col=it(col.lab,3),cex=it(cex.lab,3),xpd=xpd)
  } else {
    if( !is.null(Xlab) )
      text(-c30*it(dist.lab,1),-s30*it(dist.lab,1),
           Xlab,adj=c(1,1),font=it(font.lab,1),col=it(col.lab,1),cex=it(cex.lab,1),xpd=xpd)
    if( !is.null(Ylab) )
      text(1+c30*it(dist.lab,2),
           -s30*it(dist.lab,2),Ylab,adj=c(0,1),font=it(font.lab,2),col=it(col.lab,2),cex=it(cex.lab,2),xpd=xpd)
    if( !is.null(Zlab) )
      text(c60,s60+it(dist.lab,3),Zlab,
           adj=c(0.5,0),font=it(font.lab,3),col=it(col.lab,3),cex=it(cex.lab,3),xpd=xpd)
  }
  if(length(side)>0) for( nr in  1:length(side) ) {
    b <- it(pos)
    a <- 1-b
    if( is.list(at) ) {
      att<-at[[nr]]
    } else {
      att<-at
    }
    switch(as.character(side[nr]),
           "0"={},
           "1"={
             if( !is.na(it(lty)) )
               segments(0*a+c60*b,
                        0*a+s60*b,
                        1*a+c60*b,
                        0*a+s60*b,lwd=it(lwd),lty=it(lty),
                        col=it(col.axis),
                        xpd=xpd)
             if( !is.na(it(lty.tck)) && it(tick))
               segments(att*a+c60*b,
                        0*a+s60*b,
                        att*a+c60*b,
                 0*a+s60*b-it(len.tck),
                        lwd=it(lwd),lty=it(lty.tck),
                        col=it(col.axis),
                        xpd=xpd)
             if( length(itl(labels))>0 )
               text(att*a+c60*b,
                    0*a+s60*b-it(dist.axis),
                    as.graphicsAnnot(itl(labels)),adj=c(0.5,1),
                    font=it(font.axis),col=it(col.axis),
                    cex=it(cex.axis),xpd=xpd)
           },
           "2"={
             if( !is.na(it(lty)) )
               segments(1*a,
                        0*a,
                        c60*a,
                        s60*a,lwd=it(lwd),lty=it(lty),
                        col=it(col.axis),
                        xpd=xpd)
             if( !is.na(it(lty.tck)) && it(tick))
               segments((1+att*(c60-1))*a,
                        (0+att*(s60-0))*a,
                        (1+att*(c60-1))*a+it(len.tck)*c30,
                        (0+att*(s60-0))*a+it(len.tck)*s30,
                        lwd=it(lwd),lty=it(lty),
                        col=it(col.axis)
                        ,xpd=xpd)
             if( length(itl(labels))>0 )
               text((1+att*(c60-1))*a+it(dist.axis)*c30,
                    (0+att*(s60-0))*a+it(dist.axis)*s30,
                    adj=c(0,0),
                    as.graphicsAnnot(itl(labels)),
                    font=it(font.axis),
                    col=it(col.axis),
                    cex=it(cex.axis),
                    xpd=xpd)
           },
           "3"={
             if( !is.na(it(lty)) )
               segments(c60*a+1*b,
                        s60*a+0*b,
                        0*a+1*b,
                        0*a+0*b,lwd=it(lwd),lty=it(lty),
                        col=it(col.axis),
                        xpd=xpd)
             if( !is.na(lty.tck) && it(tick))
               segments((c60+att*-c60)*a+1*b,
                        (s60+att*-s60)*a+0*b,
                        (c60+att*-c60)*a+1*b+it(len.tck)*-c30,
                        (s60+att*-s60)*a+0*b+it(len.tck)*s30,
                        lwd=it(lwd),lty=it(lty.tck),
                        col=it(col.axis),
                        xpd=xpd)
             if( length(itl(labels))>0 )
               text((c60+att*-c60)*a+1*b+it(len.tck)*-c30,
                    (s60+att*-s60)*a+0*b+it(len.tck)*s30,
                    as.graphicsAnnot(itl(labels)),
                    adj=c(1,0),
                    font=it(font.axis),
                    col=it(col.axis),
                    cex=it(cex.axis),
                    xpd=xpd)
           },
           "-1"={
             if( !is.na(it(lty)) )
               segments(1*a+c60*b,
                        0*a+s60*b,
                        (1-c30*s60)*a+(c60-c30*s60)*b,
                        (0-s30*s60)*a+(s60-s30*s60)*b,lwd=it(lwd),lty=it(lty),
                        col=it(col.axis),
                        xpd=xpd)
             if( !is.na(it(lty.tck)) && it(tick))
               segments(1*a+c60*b-att*c30*s60,
                        0*a+s60*b-att*s30*s60,
                        1*a+c60*b-att*c30*s60-c60*it(len.tck),
                        0*a+s60*b-att*s30*s60+s60*it(len.tck),
                        lwd=it(lwd),lty=it(lty.tck),
                        col=it(col.axis),
                        xpd=xpd)
             if( length(itl(labels))>0 )
               text(1*a+c60*b-att*c30*s60-c60*it(dist.axis),
                    0*a+s60*b-att*s30*s60+s60*it(dist.axis),
                    as.graphicsAnnot(itl(labels)),adj=c(1,0),
                    font=it(font.axis),col=it(col.axis),
                    cex=it(cex.axis),xpd=xpd)
           },
           "-2"={
             if( !is.na(it(lty)) )
               segments(c60*a+0*b,
                        s60*a+0*b,
                        (c60+c30*s60)*a+(0+c30*s60)*b,
                        (s60-s30*s60)*a+(0-s30*s60)*b,lwd=it(lwd),lty=it(lty),
                        col=it(col.axis),
                        xpd=xpd)
             if( !is.na(it(lty.tck)) && it(tick))
               segments(c60*a+0*b+att*c30*s60,
                        s60*a+0*b-att*s30*s60,
                        c60*a+0*b+att*c30*s60+c60*it(len.tck),
                        s60*a+0*b-att*s30*s60+s60*it(len.tck),
                        lwd=it(lwd),lty=it(lty.tck),
                        col=it(col.axis),
                        xpd=xpd)
             if( length(itl(labels))>0 )
               text(c60*a+0*b+att*c30*s60+c30*it(dist.axis),
                    s60*a+0*b-att*s30*s60+s30*it(dist.axis),
                    as.graphicsAnnot(itl(labels)),adj=c(0,0),
                    font=it(font.axis),col=it(col.axis),
                    cex=it(cex.axis),xpd=xpd)
           },
           "-3"={
             if( !is.na(it(lty)) )
               segments(0*a+1*b,
                        0*a+0*b,
                        (0-0*s60)*a+(1-0*s60)*b,
                        (0+1*s60)*a+(0+1*s60)*b,lwd=it(lwd),lty=it(lty),
                        col=it(col.axis),
                        xpd=xpd)
             if( !is.na(it(lty.tck)) && it(tick))
               segments(0*a+1*b-att*0*s60,
                        0*a+0*b+att*1*s60,
                        0*a+1*b-att*0*s60-1*it(len.tck),
                        0*a+0*b+att*1*s60+0*it(len.tck),
                        lwd=it(lwd),lty=it(lty.tck),
                        col=it(col.axis),
                        xpd=xpd)
             if( length(itl(labels))>0 )
               text(0*a+1*b-att*0*s60-1*it(dist.axis),
                    0*a+0*b+att*1*s60+0*it(dist.axis),
                    as.graphicsAnnot(itl(labels)),adj=c(1,0.5),
                    font=it(font.axis),col=it(col.axis),
                    cex=it(cex.axis),xpd=xpd)
           },

           warning("Unkown axis side",side[nr])
           )
  }
}

## Paradigma:
#
# 1) Prepare the parameters
# 2) (if !add and !aspanel) Set Up Coordinate system
# 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
# 3a) Set up gsi.pairs
# 3a.1) Store original data in environment 
# 3a.2) Define an infkt
# 3a.3) Call gsi.add2pairs or gsi.pairs
# 3b) Plot directly
# 3b.1) if( newPlot ) gsi.setPlot(NULL)
# 3b.2) Transform data
# 3b.3) prepare plotting
# 3b.4) Eventually create the plot
# 3b.5) plot coordinate System
# 3b.6) Analyse Missings
# 3b.7) Plot Nonmissings
# 3b.8) Plot Missings
# 4) Postprocessing: create dependent subplotting
# 5) if( ! aspanel ) replot(plot=match.call(),add=add)
# 7) return(invisible(replot(plot=FALSE))) 
plot.acomp <- function(x,...,labels=names(x),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),margin="acomp",add=FALSE,triangle=!add,col=par("col"),axes=FALSE,plotMissings=TRUE,lenMissingTck=0.05,colMissingTck="red",mp=~simpleMissingSubplot(c(0,1,0.95,1),missingInfo,c("NM","TM",cn)),robust=getOption("robust"))
 {
   ## Prepare the parameters
   va <- NULL
   col <- unclass(col)
   if( is.null(colMissingTck) ) colMissingTck<-col
   if( is.null(labels) ) labels <- paste("var",1:gsi.getD(x))
   #colnames(x) <- labels
   newPlot <- ! aspanel && !add
   ## oX <- X
   ## Setting up the coordinate system
   if( newPlot ) {
     D  <- gsi.getD(x)
     ce <- acomp(structure(rep(1,D),names=names(x)))
     ms <- 1
     va <- NULL
     if( is.logical(center) && is.logical(scale) ) {
       if( center || scale ) {
         if( is.null(va) )
           va <- var(x,robust=robust,giveCenter=TRUE)
         if( center )
           ce <- attr(va,"center")
         if( scale ) {
           ms <- (1/sqrt(mean(gsi.diagExtract(va))))
         }
       }
     } else {
       if( !is.logical(center) )
         ce <- center
       if( !is.logical(scale) )
         ms <- scale
     }
     if( gsi.getD(x) > 3 ){
       nc <- gsi.getD(x) - if( margin %in% c("acomp","rcomp") ) 0 else 1
     }else{
       nc <- 1}
     gsi.setCoorInfo(mean=ce,
                     scale=ms,
   #                  var=va,
                     geo="acomp",
                     margin=margin,nc=nc)
     ## Set up panels
   }
   ## Setting up a gsi.pairs or plotting directly?
#   if( newPlot && gsi.getD(x) > 3 ) {
   if( !aspanel && gsi.getD(x) > 3 ) {
     ## Set up gsi.pairs
     X <- x # Store the dataset for later use in panel.
     infkt <- function(x,y,...) {
       plot.acomp(x=X,aspanel=TRUE,labels=labels,center=center,scale=scale,col=col,plotMissings=plotMissings,lenMissingTck=lenMissingTck,colMissingTck=colMissingTck,mp=mp,robust=robust,axes=axes,...)
     }
# ##    if( margin=="rcomp" )
# ##      infkt <- function(x,y,...) {
# ##        plot.acomp(rcompmargin(X,d=c(gsi.mapfrom01(x),gsi.mapfrom01(y)),pos=3),...,aspanel=TRUE,center=center,scale=scale,col=col,plotMissings=plotMissings,lenMissingTck=lenMissingTck,colMissingTck=colMissingTck,mp=mp,robust=robust,axes=axes)
# ##      }
# ##    else if(margin=="acomp") {
# ##      infkt <- function(x,y,...) {
# ##        plot.acomp(acompmargin(X,d=c(gsi.mapfrom01(x),gsi.mapfrom01(y)),pos=3),...,aspanel=TRUE,center=center,scale=scale,col=col,plotMissings=plotMissings,lenMissingTck=lenMissingTck,colMissingTck=colMissingTck,mp=mp,robust=robust,axes=axes)
# ##      }
# ##      
# ##    } else {
# ##      if( !is.numeric(margin))
# ##        margin <- match(margin,colnames(X))
# ##      fest <- X[,margin,drop=FALSE]
# ##      X    <- X[,-margin]
# ##      infkt <- function(x,y,...) {
# ##        plot.acomp(acomp(cbind(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],fest)),...,aspanel=TRUE,center=center,scale=scale,col=col,plotMissings=plotMissings,lenMissingTck=lenMissingTck,colMissingTck=colMissingTck,mp=mp,robust=robust,axes=axes)
# ##      }
# ##    }
     nc <- gsi.getCoorInfo()$nc
     if( add )
       noreplot(gsi.add2pairs(sapply(1:nc,gsi.mapin01),infkt,...))
     else
       noreplot(gsi.pairs(sapply(1:nc,gsi.mapin01),labels=labels,panel=infkt,...,trimode=TRUE))
   } else {
     ## Dont set up gsi.pairs, but plot directly
     if( newPlot )
       gsi.setPlot(NULL) ## not a gsi.pairs plot
     ## From here on we have a single ternary diagram
     ## and a coorInfo set up
     ## in this else part we are actually plotting
     
     ## Transform into local system
#     coorInfo <- gsi.getCoorInfo()
     Y <- oneOrDataset(gsi.pltTrafo(x,what="data"))
     cn<- gsi.pltTrafo(labels,what="names")
     if( is.null(cn) ) {
       cn <- c("x","y","z")
     }
     ## Prepare plotting
     s60 <- sin(pi/3)
     c60 <- cos(pi/3)
     s30 <- sin(pi/6)
     c30 <- cos(pi/6)
     
     ## Create the plot
     if( newPlot ) {
       usr <- par(list("pty")); on.exit(par(usr),add=TRUE)
       par( pty="s" )
       plot(x=c(0,c60,1,0),y=c(0,s60,0,0),
            xlim=c(0,1),ylim=c(0,1),type="n",xlab="",ylab="",
            axes=FALSE)
     }
     ## Plot triangle, axes and annotations
     if( triangle ) {
       lines(x=c(0,c60,1,0),y=c(0,s60,0,0))
     }
     if( is.logical(axes) ) {
       if( axes ) 
         ternaryAxis(1:3,Xlab=cn[1],Ylab=cn[2],Zlab=cn[3],aspanel=TRUE,small=aspanel)
       else if( !add )
         ternaryAxis(0,Xlab=cn[1],Ylab=cn[2],Zlab=cn[3],aspanel=TRUE,small=aspanel)
     } else {
       if( is.null(axes$Xlab) ) axes$Xlab<- cn[1]
       if( is.null(axes$Ylab) ) axes$Ylab<- cn[2]
       if( is.null(axes$Zlab) ) axes$Zlab<- cn[3]
       if( is.null(axes$small) )
         axes$small <- aspanel
       if( is.call(axes) ) {
         env <- environment(axes)
         if( is.null(env) )
           eval(axes)
         else
           eval(axes,env)
       } else {
         axes$aspanel<-TRUE
         do.call("ternaryAxis",as.list(axes))
       }
     }
     ## Analyse Missingness
     nmv <- is.NMV(Y)
     bdl <- is.BDL(Y)
     nMis <- apply(if( plotMissings ) !nmv else !(nmv | bdl),1,sum)
     nonmissing   <- nMis == 0
     Y<-oneOrDataset(clo(structure(c(ifelse(nmv,Y,0)),dim=dim(Y))))
     ## plot  Nonmissings
     x1 <- Y[,2]+Y[,3]*c60
     y1 <- Y[,3]*s60
     points(ifelse(nonmissing,x1,NA),ifelse(nonmissing,y1,NA),...,col=col)
     ## plot Missings
     if( plotMissings && !all(nonmissing) ) {
       ## Missing ticks
       totallyMissing   <- nMis > 1
       partiallyMissing <- nMis == 1
       wM <- apply(!nmv,1,function(x) c(which(x),0)[1])
       if( lenMissingTck != 0 && any(partiallyMissing) ) {
         xD <- -(x1-c(NA,0,1,c60)[wM+1])*lenMissingTck
         yD <- -(y1-c(NA,0,0,s60)[wM+1])*lenMissingTck
         segments(x1,y1,x1+xD,y1+yD,col=colMissingTck,...,xpd=TRUE)
       }
       ## Missing-Panel-Plot
       missingInfo <-  c(NotMissing=sum(nonmissing),
                         TotallyMissing=sum(totallyMissing),
                         Missing1=sum(partiallyMissing&!nmv[,1]),
                         Missing2=sum(partiallyMissing&!nmv[,2]),
                         Missing3=sum(partiallyMissing&!nmv[,3]))
       eval(mp[[2]])
     }
     ## Id 
     if( id ) {
       if( is.null(idlabs) )
         nam <- names(x)
       idlabs <- apply(x,1,function(k) paste(names(x),k,sep="=",col="\n"))
                                        #paste(cn[1],"=",round(Y[,1],2),",\n",
                                        #              cn[2],"=",round(Y[,2],2),",\n",
                                        #              cn[3],"=",round(Y[,3],2))
       return( identify(x1,y1,idlabs,col=idcol,xpd=NA))
     }
   }
### After the actual plot:
   ## Create PCA
  if( pca && ! aspanel) {
       if( is.null(va) ) va <- var(acomp(x),robust=robust,giveCenter=TRUE)
       pca.d <- acomp(princomp(acomp(x),robust=robust,covmat=va)$Loadings[1,])
       pca.c <- attr(va,"center")
       noreplot(straight.acomp(pca.c,pca.d,col=col.pca))
   }
   ## Set up replotting information:
   if( !aspanel )
     replot(plot=match.call(),add=add)
   return(invisible(replot(plot=FALSE)))
}

plot.ccomp <- function(x,...) {
  x<-unclass(x)
  x<-rcomp(structure(abs(x+runif(length(x),-0.3,0.3)),dim=dim(x)))
  plot(x,...)
}


# modified by Raimon in May 2008
# more modifications by Gerald
### plot.rcomp
### Differences
###  a) Scaling equals Acomp centering, no mean shift. *
###  b) Uses rcomp geometry in finding the centering mean *
###      but not in doing it !!!
###  c) Principal Components are in rcomp geometry
###  d) default margin is acomp
###  e) BDL is not treated as nonmissing if missings are not plotted

plot.rcomp <- function(x,...,labels=names(x),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),margin="rcomp",add=FALSE,triangle=!add,col=par("col"),axes=FALSE,plotMissings=TRUE,lenMissingTck=0.05,colMissingTck="red",mp=~simpleMissingSubplot(c(0,1,0.95,1),missingInfo,c("NM","TM",cn)),robust=getOption("robust"))
{
# 1) Prepare the parameters
  col <-unclass(col)
  if( is.null(colMissingTck) ) colMissingTck<-col
  if( !is.logical(center) || !is.logical(scale) || center || scale )
    warning("Scaling and centering meaningless for rcomp-compositions");
# 2) (if !add and !aspanel) Set Up Coordinate system
   newPlot <- ! aspanel && !add
   ## oX <- X
   ## Setting up the coordinate system
   if( newPlot ) {
     D  <- gsi.getD(x)
     ce <- rcomp(structure(rep(1,D),names=names(x)))
     ms <- 1
     #va <- NULL
     #if( is.logical(center) && is.logical(scale) ) {
     #  if( center || scale ) {
     #    if( is.null(va) )
     #      va <- var(x,robust=robust,giveCenter=TRUE)
     #    if( scale )
     #      ce <- attr(va,"center")
     #    if( center ) {
     #      warning("Centering not supported for rcomp-compositions")
     #      ms <- 1 # (1/sqrt(mean(gsi.diagExtract(va))))
     #    }
     #  }
     #} else {
     #  if( !is.logical(center) )
     #    ce <- center ## Both are the same thing
     #  if( !is.logical(scale) )
     #    ce <- scale  ## Both are the same thing
     #}
     if( gsi.getD(x) > 3){
       nc <- gsi.getD(x) - if( margin %in% c("acomp","rcomp") ) 0 else 1
     }else{
       nc <- 1
     }
     gsi.setCoorInfo(mean=ce,scale=ms,geo="rcomp", margin=margin,nc=nc)
   }
# 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
  #X <- oneOrDataset(x)
  #oX<-X
  if( !aspanel && gsi.getD(x) > 3 ) {
# 3a) Set up gsi.pairs
# 3a.1) Store original data in environment
    X<-x
# 3a.2) Define an infkt
    infkt <- function(x,y,...) {
      plot.rcomp(x=X,aspanel=TRUE,labels=labels,center=center,scale=scale,col=col,plotMissings=plotMissings,lenMissingTck=lenMissingTck,colMissingTck=colMissingTck,mp=mp,robust=robust,axes=axes,...)
     }
# 3a.3) Call gsi.add2pairs or gsi.pairs
#    if( margin=="rcomp" )
#      infkt <- function(x,y,...) {
#        plot.rcomp(rcompmargin(X,d=c(gsi.mapfrom01(x),gsi.mapfrom01(y)),pos=3),...,aspanel=TRUE,col=col,plotMissings=plotMissings,lenMissingTck=lenMissingTck,colMissingTck=colMissingTck,mp=mp,robust=robust,axes=axes)
#      }
#  else if(margin=="acomp") {
#    infkt <- function(x,y,...) {
#      plot.rcomp(acompmargin(X,d=c(gsi.mapfrom01(x),gsi.mapfrom01(y)),pos=3),...,aspanel=TRUE,col=col,plotMissings=plotMissings,lenMissingTck=lenMissingTck,colMissingTck=colMissingTck,mp=mp,robust=robust,axes=axes)
#    }
#    
#  } else {
#    if( !is.numeric(margin))
#      margin <- match(margin,colnames(X))
#    fest <- X[,margin,drop=FALSE]
#    X    <- X[,-margin]
#    infkt <- function(x,y,...) {
#                                        # plot.rcomp(acomp(cbind(X[,c(gsi.mapfrom01(y),gsi.mapfrom01(x))],fest)),...,aspanel=TRUE)
#      plot.rcomp(rcomp(cbind(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],fest)),...,aspanel=TRUE, col=col,plotMissings=plotMissings,lenMissingTck=lenMissingTck,colMissingTck=colMissingTck,mp=mp,robust=robust,axes=axes)
#    }
#  }
    nc <- gsi.getCoorInfo()$nc 
     if( add ){
       noreplot(gsi.add2pairs(sapply(1:nc,gsi.mapin01),infkt,...))
     }else{
       noreplot(gsi.pairs(sapply(1:nc,gsi.mapin01),labels=labels,panel=infkt,...,trimode=TRUE))
     }
#    if( add )
#      gsi.add2pairs(sapply(1:nc,gsi.mapin01),infkt,...)
#    else
#      gsi.pairs(sapply(1:NCOL(X),gsi.mapin01),labels=labels,panel=infkt,...,
#              trimode=TRUE)
  } else {
### 3b) Plot directly
### 3b.1) if( newPlot ) gsi.setPlot(NULL)
    if( newPlot ) gsi.setPlot(NULL)
### 3b.2) Transform data
    Y <- oneOrDataset(gsi.pltTrafo(x,what="data"))
    cn<- gsi.pltTrafo(labels,what="names")
    if( is.null(cn) ) {
      cn <- c("x","y","z")
    }
    
### 3b.3) prepare plotting
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
    s30 <- sin(pi/6)
    c30 <- cos(pi/6)
### 3b.4) Eventually create the plot
    if( newPlot ) {
      usr <- par(list("pty")); on.exit(par(usr),add=TRUE)
      par( pty="s" )
      plot(x=c(0,c60,1,0),y=c(0,s60,0,0),
           xlim=c(0,1),   #c(min(c(0,x)),max(c(1,x))),
           ylim=c(0,1),   #c(min(c(0,y)),max(c(1,y))),
           type="n",xlab="",ylab="",
           axes=FALSE)
     # gsi.setPlot(NULL)
    }
  
### 3b.5) plot coordinate System
    if( triangle ) {
      segments(x0=c(0,1,c60),y0=c(0,0,s60),x1=c(1,c60,0),y1=c(0,s60,0))
    }
    if( is.logical(axes) ) {
      if( axes ) 
        ternaryAxis(1:3,Xlab=cn[1],Ylab=cn[2],Zlab=cn[3],aspanel=TRUE)
      else
        ternaryAxis(0,Xlab=cn[1],Ylab=cn[2],Zlab=cn[3],aspanel=TRUE)
    } else {
      if( is.null(axes$Xlab) )  axes$Xlab<- cn[1]
      if( is.null(axes$Ylab) ) axes$Ylab<- cn[2]
      if( is.null(axes$Zlab) ) axes$Zlab<- cn[3]
      if( is.null(axes$small) )
        axes$small <- aspanel
      if( is.call(axes) )
        eval(axes)
      else {
        axes$aspanel<-TRUE
        do.call("ternaryAxis",as.list(axes))
      }
    }
    
### 3b.6) Analyse Missings
 
#    X <- rcomp(oneOrDataset(X),c(1,2,3))
#    Y <- X
#    ce <- acomp(c(1,1,1))
#    ms <- 1
#    va <- 1 ## Vorsicht!!!!!
#    gsi.setCoorInfo(mean=ce,scale=ms)
    nmv <- is.NMV(Y) | (is.finite(Y) & Y==0)
    bdl <- is.BDL(Y)
    nMis <- apply(if( plotMissings ) !nmv else !(nmv | is.BDL(Y)),1,sum)
    nonmissing   <- nMis == 0
    Y<-oneOrDataset(clo(structure(c(ifelse(nmv,Y,0)),dim=dim(Y))))
  #  names(Y) <- cn
### 3b.7) Plot Nonmissings
    x1 <- Y[,2]+Y[,3]*c60
    y1 <- Y[,3]*s60
    points(ifelse(nonmissing,x1,NA),ifelse(nonmissing,y1,NA),...,col=col)
### 3b.8) Plot Missings
    if( plotMissings && !all(nonmissing) ) {
      ## Missing ticks
      totallyMissing   <- nMis > 1
      partiallyMissing <- nMis == 1
      wM <- apply(!nmv,1,function(x) c(which(x),0)[1])
      if( lenMissingTck != 0 && any(partiallyMissing) ) {
        xD <- -c(NA,c30,-c30,0)[wM+1]*lenMissingTck*s60
        yD <- -c(NA,s30,s30,-1)[wM+1]*lenMissingTck*s60
        segments(x1,y1,x1+xD,y1+yD,col=colMissingTck,...,xpd=TRUE)
      }
      ## equal area plot
      missingInfo <-  c(NotMissing=sum(nonmissing),
                        TotallyMissing=sum(totallyMissing),
                        Missing1=sum(partiallyMissing&!nmv[,1]),
                        Missing2=sum(partiallyMissing&!nmv[,2]),
                        Missing3=sum(partiallyMissing&!nmv[,3]))
      eval(mp[[2]])
    }
    ##### !!! Special activity
    if( id ) {
      if( is.null(idlabs) )
        idlabs <- paste(cn[1],"=",round(Y[,1],2),",\n",
                        cn[2],"=",round(Y[,2],2),",\n",
                        cn[3],"=",round(Y[,3],2))
      if( !aspanel ) replot(plot=match.call(),add=add)  
      return( identify(x1,y1,idlabs,col=idcol,xpd=NA))
    }
  }
## 4) Postprocessing: create dependent subplotting
  if( pca && ! aspanel) {
    va <- var(rcomp(x),robust=robust,giveCenter=TRUE)
    pca.d <- rcomp(princomp(rcomp(x),covmat=va,robust=robust)$Loadings[1,])
    pca.c <- attr(va,"center")
    straight.rcomp(pca.c,pca.d,col=col.pca)
  }
# 5) if( ! aspanel ) replot(plot=match.call(),add=add)
  if( !aspanel ) replot(plot=match.call(),add=add)  
# 7) return(invisible(replot(plot=FALSE)))   
  return(invisible(replot(plot=FALSE)))
}



isoPortionLines <- function(...) {
  UseMethod("isoPortionLines",gsi.getCoorInfo()$mean)
}

isoProportionLines <- function(...) {
  UseMethod("isoProportionLines",gsi.getCoorInfo()$mean)
}

####  simplified for pure adding (with panels):
### 1) Prepare the parameters
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment 
### 3a.2) Define an infkt
### 3a.3) gsi.add2pairs
### 3b) Plot directly
### 3b.1) --
### 3b.2) Transform data
### 3b.3) prepare plotting
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings 
### 3b.7) Plot Nonmissings
### 3b.8) Plot Missings    # only with data-type plotting
### 4) Postprocessing: create dependent subplotting
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
### 7) return(invisible(replot(plot=FALSE))) 

isoPortionLines.acomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="") {
  coor <- gsi.getCoorInfo()
  if( coor$scale != 1 )
    stop("Scaling not implemented  in isoPortionLines.acomp")
  for(k in parts) {
    isoCollaps <- function(kkw,k) {
      s <- sum(kkw[-k])
      rcomp(switch(k,c(kkw[k],0,s),c(s,kkw[k],0),c(0,s,kkw[k])))
    }
    dir   <- rep(0,3)
    dir[k]<- 0
    dir[-k]<-c(1,-1)
    dirt = rmult(rcomp(dir/unclass(coor$mean))) # directions must be also perturbed!
    diru = rmult(rcomp(dir)) # directions must be also per
    for(p in at) {
      start <- rep((1-p)/2,3)
      start[k] <- p
      tx <- rep(0,3)    # to p,center=TRUElace the text, it is easy to...
      tx[k] <- p         # ... take points on the sides of the diagram
      tx[c(3,1,2)[k]] <- 1-p
      if( p>0 && p<1 ) {
        #kwt <- clo(start/unclass(coor$mean)) # perturbation must be positive
        kwu <- clo(start)#/unclass(coor$mean)) # perturbation must be positive
        tx <- clo(tx/unclass(coor$mean))
        if(lines) noreplot(straight.rcomp(kwu,diru,...))
        kwt <- isoCollaps(tx,k)
        if( labs ){
            text(kwt[2]+cos(60*pi/180)*kwt[3],kwt[3]*sin(60*pi/180),
            paste(p*total,unit[(k-1)%%length(unit)+1]),pos=c(2,1,4)[k],...,xpd=TRUE)
           }
      }
    }
  }
  replot(plot=match.call(),add=TRUE)  

  invisible(NULL)
}

#isoPortionLines.acomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="") {
#  coor <- gsi.getCoorInfo()
#  if( coor$scale != 1 )
#    stop("Scaling not implemented  in isoPortionLines.acomp")
#  for(k in parts) {
#    isoCollaps <- function(kw,k) {
#      s <- sum(kw[-k])
#      rcomp(switch(k,c(kw[k],0,s),c(s,kw[k],0),c(0,s,kw[k])))
#    }
#    dir   <- rep(0,3)
#    dir[k]<- 0
#    dir[-k]<-c(1,-1)
#    for(p in at) {
#      start <- rep((1-p)/2,3)
#      start[k] <- p
#      if( p>0 && p<1 ) {
#        kw<-rcomp(acomp(start)-coor$mean)
#        if(lines) straight(kw,rmult(dir),...)
#        kw <- isoCollaps(kw,k)
#        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
#             paste(p*total,unit[(k-1)%%length(unit)+1]),...,pos=c(2,1,4)[k],xpd=TRUE)
#      }
#    }
#  }
#  invisible(NULL)
#}

isoProportionLines.acomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,labs=TRUE,lines=TRUE) {
  coor <- gsi.getCoorInfo()
  for(k in parts) {
    dir   <- rep(0.25,3)
    dir[k]<- 0.5
    dir = acomp(dir) # * coor$scale # not needed
    for(p in at) {
      if( p>0 && p<1) {
        start <- acomp(switch(k,c(1,1-p,p),c(p,1,1-p),c(1-p,p,1)))
#        start <- acomp(switch(k,c(1E-17,1-p,p),c(p,1E-17,1-p),c(1-p,p,1E-17)))
        kw<-(start+coor$mean)*(coor$scale)  # perturbation must be positive
        if(lines) noreplot(straight(kw,dir,...))
        kw[k]<- 1E-17
        kw <- acomp(kw)
        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
             paste(p),...,pos=c(4,2,1)[k],xpd=TRUE)
      }
    }
  }
  replot(plot=match.call(),add=TRUE)  
  invisible(NULL)
}

#isoProportionLines.acomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,labs=TRUE,lines=TRUE) {
#  coor <- gsi.getCoorInfo()
#  for(k in parts) {
#    dir   <- rep(0.25,3)
#    dir[k]<- 0.5
#    for(p in at) {
#      if( p>0 && p<1) {
#        start <- acomp(switch(k,c(1,1-p,p),c(p,1,1-p),c(1-p,p,1)))
#        kw<-(start-coor$mean)*coor$scale
#        if(lines) straight(kw,acomp(dir),...)
#        kw[k]<- 1E-17
#        kw <- acomp(kw)
#        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
#             paste(p),...,pos=c(4,2,1)[k],xpd=TRUE)
#      }
#    }
#  }
#  invisible(NULL)
#}


isoPortionLines.rcomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="") {
  coor <- gsi.getCoorInfo()
  if( coor$scale != 1 || norm(acomp(coor$mean))>0.001 )
    stop("Scaling and centering not implemented  in isoPortionLines.rcomp")
  for(k in parts) {
    isoCollaps <- function(kw,k) {
      s <- sum(kw[-k])
      rcomp(switch(k,c(kw[k],0,s),c(s,kw[k],0),c(0,s,kw[k])))
    }
    dir   <- rep(0,3)
    dir[k]<- 0
    dir[-k]<-c(1,-1)
    for(p in at) {
      start <- rep((1-p)/2,3)
      start[k] <- p
      if( p>0 && p<1 ) {
        try({
        kw <- rcomp(start) # isoCollaps(rcomp(start)-coor$mean)
        if(lines ) noreplot(straight(kw,rmult(dir),...))
        kw <- isoCollaps(kw,k)
        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
             paste(p*total,unit[(k-1)%%length(unit)+1]),...,pos=c(2,1,4)[k],xpd=TRUE)
      },silent=FALSE)
      }
    }
  }
  replot(plot=match.call(),add=TRUE)  
  invisible(NULL)
}


isoProportionLines.rcomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,labs=TRUE,lines=TRUE) {
  coor <- gsi.getCoorInfo()
  if( coor$scale != 1 || norm(acomp(coor$mean))>0.001)
    stop("Scaling not implemented  in isoPortionLines.rcomp")
  for(k in parts) {
    dir   <- rep(0.25,3)
    dir[k]<- 0.5
    for(p in at) {
      if( p>0 && p<1) {
        start <- acomp(switch(k,c(1,1-p,p),c(p,1,1-p),c(1-p,p,1)))
        kw<-start
        if(lines) noreplot(straight(kw,acomp(dir),...))
        kw[k]<- 1E-17
        kw <- acomp(kw)
        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
             paste(p),...,pos=c(4,2,1)[k],xpd=TRUE)
      }
    }
  }
  replot(plot=match.call(),add=TRUE)  
  invisible(NULL)
}

#plot.aplus <- function (x, ..., labels = colnames(x), cn = colnames(x), aspanel = FALSE,
# The X not x is here to apply that to the dataset if only a single point is given.
# Is there any situation in which this is a problem?
plot.aplus <- function (x, ..., labels = colnames(X), cn = colnames(X), aspanel = FALSE,
    id = FALSE, idlabs = NULL, idcol = 2, center = FALSE, scale = FALSE,
    pca = FALSE, col.pca = par("col"), add = FALSE, logscale = TRUE,
    xlim = NULL, ylim = xlim, col = par("col"), plotMissings = TRUE,
    lenMissingTck = 0.05, colMissingTck = "red", mp = ~simpleMissingSubplot(missingPlotRect,
        missingInfo, c("NM", "TM", cn)), robust = getOption("robust")){
    col <- unclass(col)
    if (is.null(colMissingTck))
        colMissingTck <- col
    if (!aspanel && (center || scale))
        x <- scale(x, center = center, scale = scale, robust = robust)
    if (!aspanel && !add) {
        if (is.null(xlim))
            xlim <- if (logscale)
                apply(ifelse(is.NMV(x), x, NA), 2, function(x) {
                  erg <- range(x, na.rm = TRUE)
                  if (erg[1] == erg[2])
                    erg <- erg * c(1/1.1, 1.1)
                  erg
                })
            else apply(ifelse(is.NMV(x), x, 0), 2, function(x) c(0,
                max(x)))
        if (is.null(ylim))
            ylim <- xlim
    }
    X <- oneOrDataset(x)
    oX <- X
    if (NCOL(X) > 2) {
        infkt <- function(x, y, ...) {
            plot.aplus(X[, c(x[1], y[1]), drop = FALSE], ...,
                aspanel = TRUE, center = center, scale = scale,
                logscale = logscale, add = add, col = col, plotMissings = plotMissings,
                lenMissingTck = lenMissingTck, colMissingTck = colMissingTck,
                mp = mp, robust = robust)
        }
        usr <- par(list("xlog", "ylog"))
        on.exit(par(usr), add = TRUE)
        if (add)
            gsi.add2pairs(matrix(1:NCOL(X), nrow = 1), panel = infkt,
                ...)
        else {
            gsi.pairs(matrix(1:NCOL(X), nrow = 1), labels = labels,
                panel = infkt, ..., log = ifelse(logscale, "xy",
                  ""), xlim = xlim, ylim = ylim)
        }
    }
    else {
        if (is.null(cn)) {
            cn <- c("x", "y")
        }
        if (!add && !aspanel) {
            plot(x = c(1), y = c(1), xlim = xlim[, 1], ylim = ylim[,
                2], type = "n", log = ifelse(logscale, "xy",
                ""), xlab = cn[1], ylab = cn[2])
            gsi.setPlot(NULL)
        }
        nmv <- is.NMV(X)
        nMis <- apply(if (plotMissings || logscale)
            !nmv
        else !(nmv | is.BDL(x)), 1, sum)
        nonmissing <- nMis == 0
        Y <- oneOrDataset(structure(c(ifelse(nmv, X, 0)), dim = dim(x)))
        x1 <- ifelse(nmv[, 1], Y[, 1], if (par("xlog"))
            10^par("usr")[1]
        else par("usr")[1])
        y1 <- ifelse(nmv[, 2], Y[, 2], if (par("ylog"))
            10^par("usr")[3]
        else par("usr")[3])
        points(ifelse(nonmissing, x1, NA), ifelse(nonmissing,
            y1, NA), ..., col = col)
        if (plotMissings && !all(nmv)) {
            opar <- par(list("xlog", "ylog", "pin", "plt"))
            usr <- par("usr")
            try({
                if (lenMissingTck != 0) {
                  yH <- nmv[, 2] & !nmv[, 1]
                  if (any(yH)) {
                    par(xlog = FALSE)
                    segments(usr[1], ifelse(yH, y1, NA), usr[1] -
                      (usr[1] - usr[2]) * lenMissingTck, ifelse(yH,
                      y1, NA), ..., col = colMissingTck, xpd = TRUE)
                    par(xlog = opar$xlog)
                  }
                  xH <- nmv[, 1] & !nmv[, 2]
                  if (any(xH)) {
                    par(ylog = FALSE)
                    segments(ifelse(xH, x1, NA), usr[3], ifelse(xH,
                      x1, NA), usr[3] - (usr[3] - usr[4]) * lenMissingTck,
                      ..., col = colMissingTck, xpd = TRUE)
                    par(ylog = opar$ylog)
                  }
                }
                if (!is.null(mp)) {
                  missingPlotRect <- c(usr[2], usr[2] + (usr[2] -
                    usr[1]) * (1 - opar$plt[2])/(opar$plt[2] -
                    opar$plt[1]), usr[3], usr[4])
                  missingInfo <- c(NotMissing = sum(nonmissing),
                    TotallyMissing = sum(nMis == 2), Missing1 = sum(nMis ==
                      1 & !nmv[, 1]), Missing2 = sum(nMis ==
                      1 & !nmv[, 2]))
                  eval(mp[[2]])
                }
            }, silent = FALSE)
            par(opar)
        }
        if (id) {
            if (is.null(idlabs))
                idlabs <- paste(cn[1], "=", round(X[, 1], 2),
                  ",\n", cn[2], "=", round(X[, 2], 2))
            if (!aspanel)
                replot(plot = match.call(), add = add)
            return(identify(x1, y1, idlabs, col = idcol, xpd = NA))
        }
    }
    if (pca && !aspanel) {
        pca.d <- iltInv(princomp(ilt(oX), robust = robust)$loadings[,
            1])
        pca.c <- mean(aplus(oX), robust = robust)
        straight.aplus(pca.c, pca.d, col = col.pca)
    }
    if (!aspanel)
        replot(plot = match.call(), add = add)
    return(invisible(NULL))
}


plot.rplus <- function (x, ..., labels = colnames(X), cn = colnames(X), aspanel = FALSE,
    id = FALSE, idlabs = NULL, idcol = 2, center = FALSE, scale = FALSE,
    pca = FALSE, col.pca = par("col"), add = FALSE, logscale = FALSE,
    xlim = NULL, ylim = xlim, col = par("col"), plotMissings = TRUE,
    lenMissingTck = 0.05, colMissingTck = "red", mp = ~simpleMissingSubplot(missingPlotRect,
        missingInfo, c("NM", "TM", cn)), robust = getOption("robust"))
{
    col <- unclass(col)
    if (is.null(colMissingTck))
        colMissingTck <- col
    if (!aspanel && !add) {
        if (is.null(xlim))
            xlim <- if (logscale)
                apply(ifelse(is.NMV(x), x, NA), 2, function(x) {
                  erg <- range(x, na.rm = TRUE)
                  if (erg[1] == erg[2])
                    erg <- erg * c(1/1.1, 1.1)
                })
            else apply(ifelse(is.NMV(x), x, 0), 2, function(x) c(0,
                max(x)))
        if (is.null(ylim))
            ylim <- xlim
    }
    if (scale)
        warning("Scaling has no graphical effect in rplus-amounts")
    X <- oneOrDataset(x)
    oX <- X
    if (NCOL(X) > 2) {
        if (!is.matrix(xlim))
            xlim <- replicate(ncol(oneOrDataset(x)), xlim)
        if (!is.matrix(xlim))
            ylim <- replicate(ncol(oneOrDataset(x)), xlim)
        infkt <- function(x, y, ...) {
            plot.rplus(X[, c(x[1], y[1]), drop = FALSE], ...,
                aspanel = TRUE, center = center, scale = scale,
                logscale = logscale, plotMissings = plotMissings,
                lenMissingTck = lenMissingTck, colMissingTck = colMissingTck,
                mp = mp, col = col)
        }
        if (add)
            gsi.add2pairs(matrix(1:NCOL(X), nrow = 1), infkt,
                ...)
        else {
            gsi.pairs(matrix(1:NCOL(X), nrow = 1), labels = labels,
                panel = infkt, ..., xlim = xlim, ylim = ylim,
                log = if (logscale)
                  "xy"
                else "")
        }
    }
    else {
        if (is.null(cn)) {
            cn <- c("x", "y")
        }
        if (!add && !aspanel) {
            plot(x = c(1), y = c(1), xlim = xlim[, 1], ylim = ylim[,
                2], type = "n", log = ifelse(logscale, "xy",
                ""), xlab = cn[1], ylab = cn[2])
            gsi.setPlot(NULL)
        }
        nmv <- is.NMV(X) | (is.finite(X) & X == 0)
        nMis <- apply(if (plotMissings || logscale)
            !nmv
        else !(nmv | is.BDL(x)), 1, sum)
        nonmissing <- nMis == 0
        x1 <- ifelse(nmv[, 1], X[, 1], if (par("xlog"))
            10^par("usr")[1]
        else par("usr")[1])
        y1 <- ifelse(nmv[, 2], X[, 2], if (par("ylog"))
            10^par("usr")[3]
        else par("usr")[3])
        points(ifelse(nonmissing, x1, NA), ifelse(nonmissing,
            y1, NA), ..., col = col)
        if (plotMissings && !all(nmv)) {
            opar <- par(list("xlog", "ylog", "pin", "plt"))
            usr <- par("usr")
            try({
                if (lenMissingTck != 0) {
                  yH <- nmv[, 2] & !nmv[, 1]
                  if (any(yH)) {
                    par(xlog = FALSE)
                    segments(usr[1], ifelse(yH, y1, NA), usr[1] -
                      (usr[1] - usr[2]) * lenMissingTck, ifelse(yH,
                      y1, NA), ..., col = colMissingTck, xpd = TRUE)
                    par(xlog = opar$xlog)
                  }
                  xH <- nmv[, 1] & !nmv[, 2]
                  if (any(xH)) {
                    par(ylog = FALSE)
                    segments(ifelse(xH, x1, NA), usr[3], ifelse(xH,
                      x1, NA), usr[3] - (usr[3] - usr[4]) * lenMissingTck,
                      ..., col = colMissingTck, xpd = TRUE)
                    par(ylog = opar$ylog)
                  }
                }
                if (!is.null(mp)) {
                  missingPlotRect <- c(usr[2], usr[2] + (usr[2] -
                    usr[1]) * (1 - opar$plt[2])/(opar$plt[2] -
                    opar$plt[1]), usr[3], usr[4])
                  missingInfo <- c(NotMissing = sum(nonmissing),
                    TotallyMissing = sum(nMis == 2), Missing1 = sum(nMis ==
                      1 & !nmv[, 1]), Missing2 = sum(nMis ==
                      1 & !nmv[, 2]))
                  eval(mp[[2]])
                }
            }, silent = FALSE)
            par(opar)
        }
        if (id) {
            if (is.null(idlabs))
                idlabs <- paste(cn[1], "=", round(X[, 1], 2),
                  ",\n", cn[2], "=", round(X[, 2], 2))
            if (!aspanel)
                replot(plot = match.call(), add = add)
            return(identify(x1, y1, idlabs, col = idcol, xpd = NA))
        }
    }
    if (pca && !aspanel) {
        pca.d <- princomp(iit(oX), robust = robust)$loadings[,
            1]
        pca.c <- mean(rplus(oX), robust = robust)
        straight.rplus(pca.c, pca.d, col = col.pca)
    }
    if (!aspanel)
        replot(plot = match.call(), add = add)
    return(invisible(NULL))
}


plot.rmult <- function(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),add=FALSE,logscale=FALSE,col=par("col"),robust=getOption("robust")) {
X <- oneOrDataset(x)
oX <- X
if( NCOL(X) > 2 ) {
    infkt <- function(x,y,...) {
      plot.rmult(X[,c(x[1],y[1]),drop=FALSE],...,aspanel=TRUE,center=center,scale=scale,pca=pca,col.pca=col.pca,logscale=logscale,col=col)
    }
    if( add )
      gsi.add2pairs(matrix(1:NCOL(X),nrow=1),infkt,...)
    else
      gsi.pairs(matrix(1:NCOL(X),nrow=1),labels=labels,panel=infkt,...,xlim=apply(x,2,range),ylim=apply(x,2,range),log=if(logscale) "xy" else "")
  } else {
    if( is.null(cn) ) {
      cn <- c("x","y")
    }
    x <- X[,1]
    y <- X[,2]
    if( aspanel && ! add ) {
      usr <- par(list("usr")); on.exit(par(usr),add=TRUE)
      #if( logscale )
       # par( xlog=TRUE,ylog=TRUE,usr=c(log10(min(x)),log10(max(x)),log10(min(y)),log10(max(y))))
      #else
      #  par( usr=c(min(x),max(x),min(y),max(y)) )
                                        #axis(1)
                                        #axis(2)
    } else {
      if( !add ) {
        plot(x=c(1),y=c(1),
             xlim=range(x),ylim=range(y),type="n",
             log=ifelse(logscale,"xy",""),xlab=cn[1],ylab=cn[2])
        #gsi.plots[[dev.cur()]]<<-NULL
        gsi.setPlot(NULL)
      }
    }
    points(x,y,...,col=col)
    if( id ) {
      if( is.null(idlabs) )
        idlabs <- paste(cn[1],"=",round(X[,1],2),",\n",
                        cn[2],"=",round(X[,2],2))
      if( !aspanel ) replot(plot=match.call(),add=add)  
      return( identify(x,y,idlabs,col=idcol,xpd=NA))
    }
  }
if( pca && ! aspanel ) {
  pca.d <- iitInv(princomp(iit(oX),robust=robust)$loadings[,1])
  pca.c <- mean(oX,robust=robust)
  straight.rmult(pca.c,pca.d,col=col.pca)
}
    if( !aspanel ) replot(plot=match.call(),add=add)  

return( invisible(NULL))
}



#### Plot Function Paradigma
####  simplified for pure adding (with panels):
### 1) Prepare the parameters
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment 
### 3a.2) Define an infkt
### 3a.3) gsi.add2pairs
### 3b) Plot directly
### 3b.1) --
### 3b.2) Transform data
### 3b.3) prepare plotting
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings 
### 3b.7) Plot Nonmissings
### 3b.8) Plot Missings    # only with data-type plotting
### 4) Postprocessing: create dependent subplotting
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
### 7) return(invisible(replot(plot=FALSE))) 

### Called in three modes:
### a) User Mode: Add lines to plot
### b) internal panel mode: Draws lines in global coordinates
### c) external panel mode: Draws lines in local coordinates

lines.acomp <- function(x,...,steps=30,aspanel=FALSE) {
### 1) Prepare the parameters
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
  if( !aspanel && gsi.getD(x) > 3 ) {
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment 
    X <- oneOrDataset(x)
### 3a.2) Define an infkt
    infkt <- function(x,y,...) {
      lines.acomp(x=X,...,steps=steps,aspanel=TRUE)
    }
### 3a.3) gsi.add2pairs
    nc <- gsi.pltTrafo(x,"mfrow")[1]
    gsi.add2pairs(sapply(1:nc[1],gsi.mapin01),infkt,...)
  } else {
### 3b) Plot directly
### 3b.1) --
### 3b.2) Transform data
    Xt <- oneOrDataset(gsi.pltTrafo(structure(acomp(x),trafoed=attr(x,"trafoed")),"data"))
### 3b.3) prepare plotting
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings
    nmv <-is.NMV(Xt)
    NNN <- apply(nmv,1,all)
    if(!all(NNN))
       Xt[!NNN,]<-NA
### --- 
### 3b.7) Plot Nonmissings
    nx <- NROW(Xt)
    Xe <- Xt[-1,,drop=FALSE]
    Xs <- Xt[-nx,,drop=FALSE]
    nx <- NROW(Xs)
    l   <- rep((0:steps)/steps,nx)
    i   <- rep(1:nx,each=steps+1)
    XP  <- unclass(ilrInv((1-l)*ilr(Xs[i,,drop=FALSE]) +
                           l*ilr(Xe[i,,drop=FALSE])))
    x1 <- XP[,2]+XP[,3]*c60
    y1 <- XP[,3]*s60
    lines(x1,y1,...)

### 3b.8) Plot Missings    # only with data-type plotting
### --- is geometry type plotting    
  }
### 4) Postprocessing: create dependent subplotting
### ---
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
### 7) return(invisible(replot(plot=FALSE))) 
return(invisible(replot(plot=FALSE)))
}


lines.rcomp <- function(x,...,steps=30,aspanel=FALSE) {
#### Plot Function Paradigma
####  simplified for pure adding (with panels):
### 1) Prepare the parameters
  X <- oneOrDataset(clo(x))
  trafoed <- attr(x,"trafoed")
  if(is.null(trafoed)){trafoed=FALSE}
  if( (!aspanel||trafoed) && any( !is.finite(X)|x<0) ) {
    # Remove out of box problems
    X[is.BDL(x)]<-0.00000000001
    Xa <- X[-nrow(X),,drop=FALSE]
    Xb <- X[-1,,drop=FALSE]
    Xd <- Xb-Xa
    nm <- apply(is.NMV(Xa),1,all)&apply(is.NMV(Xb),1,all) 
    XaOk <-nm & apply(Xa>=0,1,all)
    XbOk <-nm & apply(Xb>=0,1,all)
    nm <- nm | (!XaOk & !XbOk)
    if( any(ra <- !XaOk & nm ) ) {
      delta <- apply(Xb[ra,,drop=FALSE]/-Xd[ra,,drop=FALSE],2,function(x) min(x[x>0]))
      Xa[ra,] <- Xb[ra,]-Xd[ra,]*delta
    }
    if( any(rb <- !XbOk & nm ) ) {
      delta <- apply(Xa[rb,,drop=FALSE]/Xd[rb,,drop=FALSE],2,function(x) min(x[x>0]))
      Xb[rb,] <- Xb[rb,]-Xd[rb,]*delta
    }
    Xa[!nm,]<-NA
    Xb[!nm,]<-NA
    LNK <- apply(Xa[-1,]==Xb[-nrow(Xb),],1,all)
    PAD <- !is.finite(LNK)
    LNK <- ifelse(is.finite(LNK),LNK,FALSE)
    padLNK <- apply(!is.finite(Xa[-1,])&!is.finite(Xb[-nrow(Xb),]),1,all)
    Redundant <- LNK | padLNK
    # Problem: Not LNK and NOT PAD would require extra pad
    gp <- rbind(c(TRUE,!Redundant),TRUE,
                c(!LNK & apply(is.finite(Xb[-nrow(Xb),])&is.finite(Xb[-nrow(Xb),]),1,all) ,FALSE))
    X<-t(structure(structure(t(cbind(Xa,Xb,Xb*NA)),dimnames=NULL),dim=c(ncol(Xb),3*nrow(Xb))))[c(gp),,drop=FALSE]
  }
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
  if( !aspanel && gsi.getD(X) > 3 ) {
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment
    #X<-X
### 3a.2) Define an infkt
    infkt <- function(x,y,...) {
      lines.rcomp(x=X,...,steps=steps,aspanel=TRUE)
    }
### 3a.3) gsi.add2pairs
    nc <- gsi.pltTrafo(X,"mfrow")
    gsi.add2pairs(sapply(1:nc[1],gsi.mapin01),infkt,...)
  } else {
### 3b) Plot directly
### 3b.1) --
### 3b.2) Transform data
#    Xt <- oneOrDataset(gsi.pltTrafo(structure(rcomp(X),trafoed=attr(X,"trafoed")),"data"))
    Xt <- oneOrDataset(gsi.pltTrafo(structure(rcomp(X),trafoed=trafoed),"data"))
#    Xt <- gsi.pltTrafo(structure(rcomp(X),trafoed=trafoed),"data")
### 3b.3) prepare plotting
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings
###   --- predone    
### 3b.7) Plot Nonmissings
    Xe <- Xt[-1,,drop=FALSE]
    Xs <- Xt[-nrow(X),,drop=FALSE]
    l   <- rep(c((0:steps)/steps),NROW(Xs))
    i   <- rep(1:NROW(Xs),each=steps+1)
    XP  <- unclass(convex.rcomp(Xs[i,,drop=FALSE],Xe[i,,drop=FALSE],l))
    x <- XP[,2]+XP[,3]*c60
    y <- XP[,3]*s60
    lines(x,y,...)
### 3b.8) Plot Missings    # only with data-type plotting
### ----
  }
### 4) Postprocessing: create dependent subplotting
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
### 7) return(invisible(replot(plot=FALSE))) 
  return(invisible(replot(plot=FALSE)))
}




lines.aplus <- function(x,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      lines.aplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...,aspanel=TRUE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    Y <- X[-1,,drop=FALSE]
    X <- X[-nrow(X),,drop=FALSE]
    l   <- rep((0:steps)/steps,NROW(X))
    i   <- rep(1:NROW(X),each=steps+1)
    XP  <- unclass(iltInv((1-l)*ilt(X[i,,drop=FALSE]) + l*ilt(Y[i,,drop=FALSE])))
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
}

lines.rplus <- function(x,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x)
  if( !aspanel && any( !is.finite(X)|x<0) ) {
    # Remove out of box problems
    X[is.BDL(x)]<-0.00000000001
    Xa <- X[-nrow(X),,drop=FALSE]
    Xb <- X[-1,,drop=FALSE]
    Xd <- Xb-Xa
    nm <- apply(is.NMV(Xa),1,all)&apply(is.NMV(Xb),1,all) 
    XaOk <-nm & apply(Xa>=0,1,all)
    XbOk <-nm & apply(Xb>=0,1,all)
    nm <- nm | (!XaOk & !XbOk)
    if( any(ra <- !XaOk & nm ) ) {
      delta <- apply(Xb[ra,,drop=FALSE]/-Xd[ra,,drop=FALSE],2,function(x) min(x[x>0]))
      Xa[ra,] <- Xb[ra,]-Xd[ra,]*delta
    }
    if( any(rb <- !XbOk & nm ) ) {
      delta <- apply(Xa[rb,,drop=FALSE]/Xd[rb,,drop=FALSE],2,function(x) min(x[x>0]))
      Xb[rb,] <- Xb[rb,]-Xd[rb,]*delta
    }
    Xa[!nm,]<-NA
    Xb[!nm,]<-NA
    LNK <- apply(Xa[-1,]==Xb[-nrow(Xb),],1,all)
    PAD <- !is.finite(LNK)
    LNK <- ifelse(is.finite(LNK),LNK,FALSE)
    padLNK <- apply(!is.finite(Xa[-1,])&!is.finite(Xb[-nrow(Xb),]),1,all)
    Redundant <- LNK | padLNK
    # Problem: Not LNK and NOT PAD would require extra pad
    gp <- rbind(c(TRUE,!Redundant),TRUE,
                c(!LNK & apply(is.finite(Xb[-nrow(Xb),])&is.finite(Xb[-nrow(Xb),]),1,all) ,FALSE))
    X<-t(structure(t(cbind(Xa,Xb,Xb*NA)),dim=c(ncol(Xb),3*nrow(Xb))))[c(gp),,drop=FALSE]
  }
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      lines.rplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...,aspanel=TRUE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    Y <- X[-1,,drop=FALSE]
    X <- X[-nrow(X),,drop=FALSE]
    l   <- rep((0:steps)/steps,NROW(X))
    i   <- rep(1:NROW(X),each=steps+1)
    XP  <- unclass(iitInv((1-l)*iit(X[i,,drop=FALSE]) + l*iit(Y[i,,drop=FALSE])))
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
}

lines.rmult <- function(x,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      lines.rmult(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...,aspanel=TRUE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    Y <- X[-1,,drop=FALSE]
    X <- X[-nrow(X),,drop=FALSE]
    l   <- rep((0:steps)/steps,NROW(X))
    i   <- rep(1:NROW(X),each=steps+1)
    XP  <- (1-l)*unclass(X)[i,,drop=FALSE] + l*unclass(Y)[i,,drop=FALSE]
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
}

#segments <- function(x0,...) UseMethod("segments",x0)
#segments.default <- graphics::segments

segments.acomp <- function(x0,y,...,steps=30,aspanel=FALSE) {
#### Plot Function Paradigma
####  simplified for pure adding (with panels):
### 1) Prepare the parameters
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
  if( !aspanel && gsi.getD(x0) > 3 ) {
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment
    x0Store<-x0
    yStore <-y
### 3a.2) Define an infkt
    infkt <- function(x,y,...) {
      segments.acomp(x0=x0Store,y=yStore,...,steps=steps,aspanel=TRUE)
    }
### 3a.3) gsi.add2pairs
    nc <- gsi.getCoorInfo()$nc 
    gsi.add2pairs(sapply(1:nc,gsi.mapin01),infkt,...) 
   } else {
### 3b) Plot directly
### 3b.1) --
### 3b.2) Transform data
    Xs <- oneOrDataset(gsi.pltTrafo(x0,"data"))
    Xe <- oneOrDataset(gsi.pltTrafo(y,"data"))
### 3b.3) prepare plotting
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings
    ## Unclear
    Xsm <- !apply(is.NMV(Xs),1,all)
    Xem <- !apply(is.NMV(Xe),1,all)
    if(any(Xsm)) Xs[Xsm,]<-NA
    if(any(Xem)) Xe[Xem,]<-NA
### 3b.7) Plot Nonmissings
    l   <- rep(c((0:steps)/steps,NA),NROW(Xs))
    i   <- rep(1:NROW(Xs),each=steps+2)
    XP  <- unclass(ilrInv((1-l)*ilr(Xs[i,]) + l*ilr(Xe[i,])))
    #print(XP)
    x <- XP[,2]+XP[,3]*c60
    y <- XP[,3]*s60
    lines(x,y,...)
### 3b.8) Plot Missings    # only with data-type plotting
    ## geometry type plotting
  }
### 4) Postprocessing: create dependent subplotting
  ## --
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
  if( ! aspanel ) replot(plot=match.call(),add=TRUE)
### 7) return(invisible(replot(plot=FALSE))) 
  return(invisible(replot(plot=FALSE)))
}


segments.rcomp <- function(x0,y,...,steps=30,aspanel=FALSE) {
### 1) Prepare the parameters
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
  if( !aspanel && gsi.getD(x0) > 3 ) {
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment
    x0Store<-x0
    yStore <-y
### 3a.2) Define an infkt
    infkt <- function(x,y,...) {
      segments.rcomp(
                     x0=x0Store,y=yStore,
                                 ...,steps=steps,aspanel=TRUE)
    }
### 3a.3) gsi.add2pairs
    nc <- gsi.getCoorInfo()$nc
    gsi.add2pairs(sapply(1:nc,gsi.mapin01),infkt,...)
  } else {
### 3b) Plot directly
### 3b.1) --
### 3b.2) Transform data
    Xs <- oneOrDataset(gsi.pltTrafo(x0,"data"))
    Xe <- oneOrDataset(gsi.pltTrafo(y,"data"))
### 3b.3) prepare plotting
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings
    ## remove those data points with BDL (in the future, maybe a flag allows to zeroreplace them)
    Xsm <- !apply(!is.BDL(Xs),1,all)
    Xem <- !apply(!is.BDL(Xe),1,all)
    if(any(Xsm)) Xs[Xsm,]<-0
    if(any(Xem)) Xe[Xem,]<-0
    ## remove those data points with Missings
    Xsm <- !apply(is.NMV(Xs)|is.BDL(Xs),1,all)
    Xem <- !apply(is.NMV(Xe)|is.BDL(Xe),1,all)
    if(any(Xsm)) Xs[Xsm,]<-NA
    if(any(Xem)) Xe[Xem,]<-NA
### 3b.7) Plot Nonmissings
    l   <- rep(c((0:steps)/steps,NA),NROW(Xs))
    i   <- rep(1:NROW(Xs),each=steps+2)
    XP  <- unclass(convex.rcomp(Xs[i,],Xe[i,],l))
    x <- XP[,2]+XP[,3]*c60
    y <- XP[,3]*s60
    lines(x,y,...)
### 3b.8) Plot Missings    # only with data-type plotting
    ## -- geometry type plot
  }
### 4) Postprocessing: create dependent subplotting
  ##--
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
### 7) return(invisible(replot(plot=FALSE))) 
  return(invisible(replot(plot=FALSE)))
}





segments.aplus <- function(x0,y,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x0,y)
  Y <- oneOrDataset(y,x0)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      segments.aplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],
                     Y[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...,aspanel=TRUE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    l   <- rep(c((0:steps)/steps,NA),NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- unclass(iltInv((1-l)*ilt(X[i,]) + l*ilt(Y[i,])))
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
}

segments.rplus <- function(x0,y,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x0,y)
  Y <- oneOrDataset(y,x0)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      segments.rplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],
                     Y[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...,aspanel=TRUE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    l   <- rep(c((0:steps)/steps,NA),NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- unclass(iitInv((1-l)*iit(X[i,]) + l*iit(Y[i,])))
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  

}

segments.rmult <- function(x0,y,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x0,y)
  Y <- oneOrDataset(y,x0)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      segments.rmult(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],
                     Y[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...,aspanel=TRUE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    l   <- rep(c((0:steps)/steps,NA),NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- (1-l)*unclass(X)[i,] + l*unclass(Y)[i,]
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  

}



gsi.closespread <- function(spread) {
  if(length(dim(spread))>3) {
    return(apply(spread,1,gsi.closespread))
  }
  d <- nrow(spread)
  Pmat <- diag(d)-1/d
  row.names(Pmat)<-row.names(spread)
  Pmat %*% spread %*% t(Pmat)
}

gsi.spreadToIsoSpace <- function(spread) {
  if(length(dim(spread))>3) {
    return(apply(spread,1,gsi.closespread))
  }
  d <- nrow(spread)
  V <- ilrBase(D=d)
  t(V) %*% spread %*% V
}

ellipses <- function(mean,...) UseMethod("ellipses",mean)



#ellipses.rcomp <- function(mean,var,r=1,...,steps=360,thinRatio=1/100) {
#mean <- oneOrDataset(ipt(mean))
#sp <- var
#w  <- seq(0,2*pi,length.out=steps)
#for(i in 1:nrow(mean)) {
#    if( length(dim(var))==3 )
#      sp<-var[i,,]
#    isp <- gsi.spreadToIsoSpace(sp)
#    mi <- mean[i,]
#    eisp <- eigen(isp,TRUE)
#    MEV <- eisp$values[1]
#    NE  <- max(sum(eisp$values>MEV*thinRatio),2)
#    X<-NULL
#    for(r1 in 1:(NE-1))
#      for(r2 in (r1+1):NE) {
#        X <- rbind(X,t(mi+ t(sqrt(eisp$values[1])*r*cos(w) %o% eisp$vectors[,1] +
#             sqrt(eisp$values[2])*r*sin(w) %o% eisp$vectors[,2])),NA)
#      }
#    lines.rcomp(uciptInv(X),...)
#  }
#}

#ellipses.rplus <- function(mean,var,r=1,...,steps=360,thinRatio=1/100) {
#mean <- oneOrDataset(iit(mean))
#sp <- var
#w  <- seq(0,2*pi,length.out=steps)
#for(i in 1:nrow(mean)) {
#    if( length(dim(var))==3 )
#      sp<-var[i,,]
##    isp <- gsi.spreadToIsoSpace(sp)
#    mi <- mean[i,]
#    eisp <- eigen(sp,TRUE)
#    MEV <- eisp$values[1]
#    NE  <- max(sum(eisp$values>MEV*thinRatio),2)
#    X<-NULL
#    for(r1 in 1:(NE-1))
#      for(r2 in (r1+1):NE) {
#        X <- rbind(X,t(mi+ t(sqrt(eisp$values[1])*r*cos(w) %o% eisp$vectors[,1] +
#                     sqrt(eisp$values[2])*r*sin(w) %o% eisp$vectors[,2])),NA)
#      }
#    lines.rmult(X,...)
#  }
#}



#ellipses.rmult <- function(mean,var,r=1,...,steps=360,thinRatio=1/100) {
#mean <- oneOrDataset(mean)
#sp <- var
#w  <- seq(0,2*pi,length.out=steps)
#for(i in 1:nrow(mean)) {
#  if( length(dim(var))==3 )
#    sp<-var[i,,]
#  mi <- mean[i,]
#  eisp <- eigen(sp,TRUE)
#  MEV <- eisp$values[1]
#  NE  <- max(sum(eisp$values>MEV*thinRatio),2)
#  X<-NULL
#  for(r1 in 1:(NE-1))
#    for(r2 in (r1+1):NE) {
#      X <- rbind(X,t(mi+ t(sqrt(eisp$values[r1])*r*cos(w) %o% eisp$vectors[,r1] +
#                           sqrt(eisp$values[r2])*r*sin(w) %o% eisp$vectors[,r2])),NA)
#    }
#  lines.rmult(X,...)
#}
#}
# var given as idt!!! --->>> This gives some failures: var given better as clr:
gsi.DrawCompEllipses = function(mean,var,r,steps=72,...) {
  # mean and variance in clr/cpt
  ei   <- eigen(clrvar2ilr(var),symmetric=TRUE)
  w <- seq(0,2*pi,length.out=steps+1)
  sw<-sin(w)
  cw<-cos(w)
  if( min(ei$values) / max(ei$values) < -1E-8) {
    warning("Non positive Semidefinite Matrix used in Ellipses")
    print(list(problem="Non positive Semidefinite Matrix used in Ellipses",var=var,eigen=ei))
  }
  rs <- sqrt(abs(ei$values))*r
 # Loop over ellipse centers
  meFull <- oneOrDataset(idt(mean))
  for(k in 1:nrow(meFull) ) {
    # me   <- structure(meFull[k,],class="rmult")
    me   <- meFull[k,]
    X <- idtInv(cbind(me[1]+rs[1]*ei$vectors[1,1]*sw+rs[2]*ei$vectors[1,2]*cw,
                      me[2]+rs[1]*ei$vectors[2,1]*sw+rs[2]*ei$vectors[2,2]*cw
                      ), orig=mean)
    noreplot(lines(structure(X,trafoed=TRUE),...,aspanel=TRUE))
  }
}




gsi.ellipsesCompPanel <- function(i,j,margin,mean,var,r=1,...,steps=72) {
  va <- gsi.pltTrafo(var, what="var", geo=class(mean))
  me <- gsi.pltTrafo(mean, what="data", geo=class(mean))
  gsi.DrawCompEllipses(structure(me,class=class(mean)),va,r,steps=steps,...)
}


ellipses.rcomp <- ellipses.acomp <- function(mean,var,r=1,...,steps=72,thinRatio=NULL,aspanel=FALSE) {
### 1) Prepare the parameters
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
  if( !aspanel && gsi.getD(mean) > 3 ) {
    if( is.null(thinRatio) ) {
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment 
### 3a.2) Define an infkt
      infkt <- function(x,y,...) {
      gsi.ellipsesCompPanel(mean,var,r,...)
    }
### 3a.3) gsi.add2pairs
      nc <- gsi.getCoorInfo()$nc 
      gsi.add2pairs(sapply(1:nc,gsi.mapin01),infkt,...,mean=mean,var=var,r=r,steps=steps,thinRatio=thinRatio)
    }
  } else {
### 3b) Plot directly
    if( is.null(thinRatio) ) {
      gsi.DrawCompEllipses(mean,var,r=r,steps=steps,...)
### Includes:
### 3b.1) --
### 3b.2) Transform data
### 3b.3) prepare plotting
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings 
### 3b.7) Plot Nonmissings
### 3b.8) Plot Missings    # only with data-type plotting
    }
  }
### 4) Postprocessing: create dependent subplotting
  if( !is.null(thinRatio) ) {
    meFull <- oneOrDataset(idt(mean))
    sp <- var
    w  <- seq(0,2*pi,length.out=steps)
    for(i in 1:nrow(meFull)) {
      if( length(dim(var))==3 )
        sp<-var[i,,]
      isp <- gsi.spreadToIsoSpace(sp)
      mi <- meFull[i,]
      eisp <- eigen(isp,TRUE)
      MEV <- eisp$values[1]
      NE  <- max(sum(eisp$values>MEV*thinRatio),2)
      X<-mi
#      X<-mi
      for(r1 in 1:(NE-1))
        for(r2 in (r1+1):NE) {
#           X <- rbind(X,t(mi + t(sqrt(eisp$values[r1])*r*cos(w) %o% eisp$vectors[,r1] +
#                                 sqrt(eisp$values[r2])*r*sin(w) %o% eisp$vectors[,r2])),NA)
            X <- rbind(X,
                         t(mi + t(sqrt(eisp$values[r1])*r*cos(w) %o% eisp$vectors[,r1] +
                           sqrt(eisp$values[r2])*r*sin(w) %o% eisp$vectors[,r2])), mi)
        }
      if( class(mean) == "rcomp")
        noreplot(lines.rcomp(uciptInv(X),...))
      else
        noreplot(lines(ilrInv(X),...))
    }
  }
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
  replot(plot=match.call(),add=TRUE)  
### 7) return(invisible(replot(plot=FALSE))) 
  return(invisible(replot(plot=FALSE)))
}



gsi.ellipsesRealPanel <- function(i,j,mean,var,r=1,...,steps=72) {
  # Preparation
  va   <- var[c(i,j),c(i,j)]
  ei   <- eigen(va)
  w <- seq(0,2*pi,length.out=steps+1)
  # Loop over ellipse centers
  meFull <- unclass(oneOrDataset(cdt(mean)))
  for(k in 1:nrow(meFull) ) {
    me   <- meFull[k,c(i,j)]
    rs <- sqrt(ei$values)*r
    X <- cdtInv(cbind(me[1]+rs[1]*ei$vectors[1,1]*sin(w)+rs[2]*ei$vectors[1,2]*cos(w),
                      me[2]+rs[1]*ei$vectors[2,1]*sin(w)+rs[2]*ei$vectors[2,2]*cos(w)
                      ),mean)
    if( class(mean) == "rplus" ) {
      noreplot(lines(structure(X,class="rplus"),...))
    } else noreplot(lines(X,...))
  }
}



ellipses.rmult<-ellipses.rplus<-ellipses.aplus <- function(mean,var,r=1,...,steps=72,thinRatio=NULL) {
    if( is.null(thinRatio) ) {
      if( gsi.getD(mean) > 2 ) {
        # Plot BigEllipse via panel
        infkt <- function(x,y,...) {
          gsi.ellipsesRealPanel(gsi.mapfrom01(x),gsi.mapfrom01(y),...)
        }
        gsi.add2pairs(sapply(1:gsi.getD(mean),gsi.mapin01),infkt,...,mean=mean,var=var,r=r,steps=steps)
      } else {
        # Plot BigEllipse directly
        gsi.ellipsesRealPanel(1,2,mean=mean,var=var,r=r,steps=steps,...)
      }
    } else {
      # Plot SmallEllipses directly
      meFull <- oneOrDataset(idt(mean))
      sp <- var
      w  <- seq(0,2*pi,length.out=steps)
      for(i in 1:nrow(meFull)) {
        if( length(dim(var))==3 )
          sp<-var[i,,]
                                        # isp <- gsi.spreadToIsoSpace(sp)
        mi <-meFull[i,]
###############
        eisp <- eigen(sp,TRUE)
        MEV <- eisp$values[1]
        NE  <- max(sum(eisp$values>MEV*thinRatio),2)
        X<-NULL
        for(r1 in 1:(NE-1))
          for(r2 in (r1+1):NE) {
            X <- rbind(t(mi + t(sqrt(eisp$values[r1])*r*cos(w) %o% eisp$vectors[,r1] +
                                  sqrt(eisp$values[r2])*r*sin(w) %o% eisp$vectors[,r2])))
         #   X <- rbind(X,t(mi + t(sqrt(eisp$values[r1])*r*cos(w) %o% eisp$vectors[,r1] +
         #                    sqrt(eisp$values[r2])*r*sin(w) %o% eisp$vectors[,r2])),NA)
            noreplot(lines(idtInv(X,mean),...))
          }
         # noreplot(lines(idtInv(X,mean),...))
      }
    }
    replot(plot=match.call(),add=TRUE)  

}




straight  <- function(x,...) UseMethod("straight",x)

#
# straight.acomp deals with partially missing directions by using
# a early transformation scheme 
#

#
# straight.rcomp deals with incompatible projections by using a
# late projection scheme
#

straight.acomp <- function(x,d,...,steps=30,aspanel=FALSE) {
### 1) Prepare the parameters
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
  if( !aspanel && gsi.getD(x) > 3 ) {
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment
    X<-x
    D<-d
### 3a.2) Define an infkt
    infkt <- function(x,y,...) {
      straight.acomp(x=X,d=D,
                     ...,steps=steps,aspanel=TRUE)
    }
### 3a.3) gsi.add2pairs
    nc <- gsi.getCoorInfo()$nc 

    gsi.add2pairs(sapply(1:nc,gsi.mapin01),infkt,...)
  } else {
### 3b) Plot directly
### 3b.1) --
### 3b.2) Transform data
    X <- oneOrDataset(gsi.pltTrafo(x,"data"),d)
    D <- oneOrDataset(gsi.pltTrafo(d,"direction"),x)
### 3b.3) prepare plotting
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings
    Xm <- !apply(is.NMV(X),1,all)
    Dm <- !apply(is.NMV(D),1,all)
    if(any(Xm)) X[Xm,]<-NA
    if(any(Dm)) D[Dm,]<-NA
    
### 3b.7) Plot Nonmissings
    D <- normalize(acomp(D)) 
    X <- perturbe(X,power.acomp(D,-scalar(acomp(X),acomp(D)))) 
    l   <- rep(2*c((0:steps)/steps,NA)-1,NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- acomp(clrInv(clr(X[i,]) + (2*l)^3*clr(D[i,])))
    x <- XP[,2]+XP[,3]*c60
    y <- XP[,3]*s60
    lines(x,y,...)
### 3b.8) Plot Missings    # only with data-type plotting
  }
### 4) Postprocessing: create dependent subplotting
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
### 7) return(invisible(replot(plot=FALSE))) 
  return(invisible(replot(plot=FALSE)))
}

straight.aplus <- function(x,d,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x,d)
  d <- oneOrDataset(d,x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      straight.aplus(
                     aplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                     aplus(d[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                                 ...,steps=steps,aspanel=TRUE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    s <- log(d[,2])/log(d[,1])
#  if( par("xlog") && par("ylog") )
#    abline(exp(log(x[,2])-log(x[,1])/s),s,untf=FALSE,...)
#  else {
    l   <- rep(2*c((0:steps)/steps,NA)-1,NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    r   <- par("usr")
    if( ! par("xlog") ) r[1:2] <- log(c(r[2]/100,r[2]))
    if( ! par("ylog") ) r[3:4] <- log(c(r[4]/100,r[4]))
    r   <- abs(r[2]-r[1])+abs(r[4]-r[3])
    XP  <- aplus(X[i,]) + r*l*normalize(aplus(d[i,]))
    noreplot(lines.rmult(XP,...)) # lines.aplus(XP,...) produced double lines
                                        #    warning("straight.aplus not yet implemented in nonlog coordinates");
                                        #  }
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  

}


straight.rcomp <- function(x,d,...,steps=30,aspanel=FALSE) {
### 1) Prepare the parameters
### 2) -
### 3) if( !aspanel && ncol>3 ) Set up gsi.pairs else plot directly
  if( !aspanel && gsi.getD(x) > 3 ) {
### 3a) Set up gsi.pairs
### 3a.1) Store original data in environment
    X<-x
    D<-d
### 3a.2) Define an infkt
    infkt <- function(x,y,...) {
      straight.rcomp(x=X,d=D,...,steps=steps,aspanel=TRUE)
    }
### 3a.3) gsi.add2pairs
    nc <- gsi.getCoorInfo()$nc 
    gsi.add2pairs(sapply(1:nc,gsi.mapin01),infkt,...)
  } else { 
### 3b) Plot directly
    X <- oneOrDataset(x,d)
    d <- oneOrDataset(d,x)
    l1 <- apply(-X/d,1,function(x) {max(x[x<=0])})*0.99999 
    l2 <- apply(-X/d,1,function(x) {min(x[x>=0])})*0.99999
    X1 <- rcomp(rmult(X)+(l1*rmult(d))) 
    X2 <- rcomp(rmult(X)+(l2*rmult(d)))
    segments.rcomp(X1,X2,...,aspanel=TRUE) ## noreplot < aspanel
### Includes:
### 3b.1) --
### 3b.2) Transform data
### 3b.3) prepare plotting
### 3b.4) --
### 3b.5) -- 
### 3b.6) Analyse Missings 
### 3b.7) Plot Nonmissings
### 3b.8) Plot Missings    # only with data-type plotting
  }
### 4) Postprocessing: create dependent subplotting
### 5) if( ! aspanel ) replot(plot=match.call(),add=TRUE)
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  
### 7) return(invisible(replot(plot=FALSE))) 
  return(invisible(replot(plot=FALSE)))
}




#straight.rplus <- function(x,d,...,steps=30) {
#  x <- oneOrDataset(x,d)
#  d <- oneOrDataset(d,x)
#  s <- d[,2]/d[,1]
#  abline(x[,2]-x[,1]/s,s,untf=TRUE)
#  l1 <- apply(-x/d,1,function(x) {max(c(-100,x[x<=0]))}) 
#  l2 <- apply(-x/d,1,function(x) {min(c(100,x[x>=0]))})
#  X1 <- rcomp(gsi.add(x,gsi.mul(l1,d))) 
#  X2 <- rcomp(gsi.add(x,gsi.mul(l2,d)))
#  segments.rplus(X1,X2,...)
#}

straight.rplus <- function(x,d,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x,d)
  d <- oneOrDataset(d,x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      straight.rplus(
                     rplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                     rmult(d[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                                 ...,steps=steps,aspanel=FALSE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {  
    #s <- log(d[,2])/log(d[,1])
                                        #  if( par("xlog") && par("ylog") )
                                        #    abline(exp(log(x[,2])-log(x[,1])/s),s,untf=FALSE,...)
                                        #  else {
    l   <- rep(2*c((0:steps)/steps,NA)-1,NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    r   <- par("usr")
    if( par("xlog") ) r[1:2] <- 10^(r[1:2])
    if( par("ylog") ) r[3:4] <- 10^(r[3:4])
    r   <- abs(r[2]-r[1])+abs(r[4]-r[3])
    XP  <- rmult(X[i,]) + r*l*normalize(rmult(d[i,]))
    noreplot(lines.rmult(XP,...))
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  

}

straight.rmult <- function(x,d,...,steps=30,aspanel=FALSE) {
  X <- oneOrDataset(x,d)
  d <- oneOrDataset(d,x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      straight.rmult(
                     aplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                     aplus(d[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                                 ...,steps=steps,aspanel=TRUE)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {  
    #s <- log(d[,2])/log(d[,1])
                                        #  if( par("xlog") && par("ylog") )
                                        #    abline(exp(log(x[,2])-log(x[,1])/s),s,untf=FALSE,...)
#  else {
    l   <- rep(2*c((0:steps)/steps,NA)-1,NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    r   <- par("usr")
    if( par("xlog") ) r[1:2] <- 10^(r[1:2])
    if( par("ylog") ) r[3:4] <- 10^(r[3:4])
    r   <- abs(r[2]-r[1])+abs(r[4]-r[3])
    XP  <- rmult(X[i,]) + r*l*normalize(rmult(d[i,]))
    noreplot(lines.rmult(XP,...))
                                        #    warning("straight.aplus not yet implemented in nonlog coordinates");
                                        #  }
  }
  if( !aspanel ) replot(plot=match.call(),add=TRUE)  

}


#straight.panel.acomp <- function(x,d,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    x <- margin(x,what)
#    d <- margin(d,what)
#    straight.acomp(x,d,...,steps=steps)
#  }
#}

#straight.panel.rcomp <- function(x,d,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    x <- margin(x,what)
#    d <- margin(d,what)
#    straight.rcomp(x,d,...,steps=steps)
#  }
#}



rDirichlet.acomp <- function(n,alpha) {
  acomp(sapply(alpha,rgamma,n=n))
}

rDirichlet.rcomp <- function(n,alpha) {
  rcomp(sapply(alpha,rgamma,n=n))
}

rpois.ccomp <- function(n,p,lambda) {
  if( missing(p) ) {
    L <- lambda
    if( missing(n) )
      n <- nrow(rplus(lambda))
  } else {
    L <- rplus(rcomp(p))*lambda
    if( missing(n) )
      n <- max(nrow(p),length(lambda))
  }
  if( length(dim(L)) == 2 ) 
    L <- L[ ((1:n)-1) %% nrow(L) +1,]
  else
    L <- t(structure(rep(L,n),dim=c(length(L),n)))
  ccomp(structure(rpois(length(L),c(L)),dim=dim(L)))
}

rmultinom.ccomp <- function(n,p,N) {
  if( missing(n) )
    n <- max(nrow(p),length(N))
  erg <- matrix(0,ncol=gsi.getD(p),nrow=n)
  if( length(dim(p))> 1) {
    for(i in 1:n) {
      erg[i,]<-rmultinom(1,N[(i-1)%%length(N)+1],p[i,])
    }
  } else {
    for(i in 1:n) {
      erg[i,]<-rmultinom(1,N[(i-1)%%length(N)+1],p)
    }
  }
  ccomp(erg)
}

runif.acomp <- function(n,D) rDirichlet.acomp(n,rep(1,D))
runif.rcomp <- function(n,D) rDirichlet.rcomp(n,rep(1,D))

rnorm.aplus <- function(n,mean,var) {
  D <- NCOL(oneOrDataset(mean))
  perturbe.aplus(iltInv(matrix(rnorm(n*length(mean)),ncol=D) %*% chol(var)),
                mean)
}

dnorm.aplus <- function(x,mean,var) {
  x <- aplus(x)
  mean <- aplus(mean)
  w <- ilt(x-mean)
  D <- ncol(oneOrDataset(x))
  if( length(dim(w)) == 2 ) 
    u <- c(rep(1,ncol(w))%*%t((solve(var,t(w)))*t(w)))
  else
    u <- sum(solve(var,w)*w)
  exp(-u/2)/sqrt(2^D*pi^D*det(var))
}

dnorm.rmult <- function(x,mean,var) {
  w <- rmult(x)-rmult(mean)
  D <- gsi.getD(x)
  if( length(dim(w)) == 2 ) 
    u <- c(rep(1,ncol(w))%*%t((solve(var,t(w)))*t(w)))
  else
    u <- sum(solve(var,w)*w)
  exp(-u/2)/sqrt(2^D*pi^D*det(var))
}




#rnorm.acomp <- function(n,mean,var) {
#  D <- NCOL(oneOrDataset(mean))
#  perturbe(ilrInv(matrix(rnorm(n*(D-1)),ncol=D-1) %*%
#                         chol(clrvar2ilr(var))),
#                mean)
#}

rnorm.acomp <-function (n, mean, var){
    D <- gsi.getD(mean)-1
    erg <- perturbe(ilrInv(matrix(rnorm(n*D), ncol = D) %*% chol(clrvar2ilr(var))), mean)
    colnames(erg) <- colnames(oneOrDataset(mean))
    erg
}




dnorm.acomp <- function(x,mean,var) {
  x <- acomp(x)
  mean <- acomp(mean)
  w <- ilr(x-mean)
  D <- ncol(oneOrDataset(x))
  if( length(dim(w)) == 2 ) 
    u <- c(rep(1,D-1)%*%t((solve(clrvar2ilr(var),t(w)))*t(w)))
  else
    u <- sum(solve(clrvar2ilr(var),w)*w)
  exp(-u/2)/sqrt(2*pi*det(clrvar2ilr(var)))
}


rnorm.ccomp <- function(n,mean,var,lambda) {
  rpois.ccomp(n,rnorm.acomp(n,mean=mean,var=var),lambda)

}

rlnorm.rplus <- function(n,meanlog,varlog) {
  D <- NCOL(oneOrDataset(meanlog))
  erg<-rplus(perturbe.aplus(exp(matrix(rnorm(n*length(meanlog)),ncol=D) %*% chol(varlog)),
                exp(meanlog)))
  colnames(erg) <- colnames(oneOrDataset(mean))
  erg
}

dlnorm.rplus <- function(x,meanlog,varlog) {
  xx <- oneOrDataset(x)
  w <- ilt(x)-meanlog
  if( length(dim(w)) == 2 ) {
    u <- c(rep(1,ncol(w))%*%t((solve(varlog,t(w)))*t(w)))
    v <- c(exp(log(xx) %*% rep(1,ncol(xx)))) 
  }
  else {
    u <- solve(varlog,w)%*%w
    v <- prod(x)
  }
  exp(-u/2)/sqrt(2*pi*det(varlog))/v
}


rnorm.rplus <- function(n,mean,var) {
  D <- ncol(var)
  erg <- rplus(pmax(rmult(matrix(rnorm(n*D),ncol=D) %*% chol(var))+rmult(mean),0))
  colnames(erg) <- colnames(oneOrDataset(mean))
  erg
}

rnorm.rmult <- function(n,mean,var) {
  D <- ncol(var)
  erg<-rmult(matrix(rnorm(n*D),ncol=D) %*% chol(var))+rmult(mean)
  colnames(erg) <- colnames(oneOrDataset(mean))
  erg
}

rnorm.rcomp <- function(n,mean,var) {
  D <- ncol(var)
  erg<-rcomp(pmax(ilr2clr(matrix(rnorm(n*D),ncol=D-1) %*% chol(clrvar2ilr(var)))+rplus(rcomp(mean)),0))
  colnames(erg) <- colnames(oneOrDataset(mean))
  erg
}



gsi.plotmargin <- function(X,d,margin,what="data") {
  X <- oneOrDataset(X)
  if( margin=="rcomp" )
    rcompmargin(X,d=d,pos=3)
  else if( margin=="acomp" )
    acompmargin(X,d=d,pos=3)
  else {
    if( ! is.numeric(margin))
      margin <- match(margin,colnames(X))
    fest <- X[,margin,drop=FALSE]
    X    <- X[,-margin,drop=FALSE]
    acomp(cbind(X[,d,drop=FALSE],fest))
  }
}


gsi.pltTrafo <- function(X,what="data",coorinfo=gsi.getCoorInfo(),geo="acomp",...) {
  # Types: data, var, mfrow, names, mfrow, direction
  trafoed <- attr(X,"trafoed")
  if( !is.null(trafoed) && trafoed )
    return(X)
  d <- coorinfo$d             # what panel
  margin  <- coorinfo$margin   # margin type
  D       <- coorinfo$D        # Dimension from which to transform
  mean    <- coorinfo$mean     # mean 
  scale   <- coorinfo$scale    # scaling
  coorgeo <- coorinfo$geo      # The geometry of scaling
  if( coorgeo == "acomp")  {
    zeroCenter <- max(coorinfo$mean)==min(coorinfo$mean)
    noScale    <- scale==1
    noTrafo    <- noScale&& zeroCenter
  } else {noTrafo<-noScale<-zeroCenter<-TRUE}
  if( !(geo %in% c("acomp","rcomp") ))
    stop("gsi.pltTrafo: Unkown data geometry",geo)
  if( !(coorgeo %in% c("acomp","rcomp") ))
    stop("gsi.pltTrafo: Unkown plot geometry",coorgeo)
  if( !(what %in% c("data","var","direction","mfrow","names")))
    stop("gsi.pltTrafo: Unkown transformation type ",what)
    
  #if( class(geo)!= "acomp" ) return(UseMethod("gsi.pltTrafo",geo))
  ### Scaling and centering
  if( what== "mfrow") {
    if( coorgeo %in% c("acomp","rcomp") )
      return(rep(gsi.getD(X),2))
    else
      return(rep(gsi.getD(X)-1,2))
  }
  if( coorgeo == "acomp" ) {
    if( what == "data" ) {
      Xn <- scale*(acomp(X)-acomp(mean))
    } else if( what =="direction") {
      if( geo == "acomp") { 
        dir <- scale*acomp(X)
      } else if( geo =="rcomp" ) {
        if( ! noTrafo )
          warning("gsi.pltTrafo: rcomp direction in acomp geometry ")
        dir <- X
      }
    } else if( what=="var" ){
      if( geo=="acomp" )
        var <- scale^2*X
      else if( geo=="rcomp" ) {
        if( !noTrafo )
          warning("gsi.pltTrafo: ",geo," variance in acomp geometry ")
        var <- X
      } else stop("Unkown geometry",geo)
    } else if(what=="names") {
      nam <- X
    } else if(what=="info") {
      orig <- X
    } else stop("gsi.pltTrafo: unkown request")
  } else if( coorgeo == "rcomp" ) {
    if( what == "data" )
      Xn <- X # No centering or scaling !!!
    else if(what=="direction") {
      if( geo == "acomp") { 
        dir <- X
      } else if( geo =="rcomp" ) {
        dir <- X
      } else stop("unkown geometry of data")
    } else if( what == "var" ) {
      if( geo == "acomp" )
        var <- X
      else if( geo=="rcomp" ) {
        var <- X
      } else stop("unkown geometry of data")
    } else if(what=="names") {
      nam <- X
    } else stop("gsi.pltTrafo: unkown request")
  } else stop("gsi.pltTrafo: unkown geometry in plot")
  ####  Marginalisation
  # is marginalisation necessary? (are we in a panel?)
  if( is.null(d) ) {
    if( what == "data" ) {
      return(Xn)
    } else if( what == "var" ){
      return(var)
    } else if( what == "names"){
      return(nam)
    } else if( what=="direction") {
      return(dir)
    } else stop("gsi.pltTrafo: Unkown request")
  }
  # acomp margin
  if( margin == "acomp" ) {
    if( what == "data" ) {
      return(acompmargin(Xn,d=d,pos=3,what="data"))
    } else if(what == "direction") {
      if( geo == "acomp" )
        return(acompmargin(dir,d=d,pos=3,what="data"))
      else if( geo=="rcomp" ) {
        warning("Incompatible geometries in plot")
        return(unclass(dir)*NA)
      }
    } else if( what == "var" ) {
      if( geo!="acomp" )
        warning("pltTrafo: Can not correctly marginalise a rcomp-variance with acomp-margin")
      return(acompmargin(var,d=d,pos=3,what="var"))
    } else if(what =="names") {
      return( c(nam[d],"*") )
    } else stop("unkown transformation type")
  # rcomp margin
  } else if( margin== "rcomp" ) {
    if( what == "data" ) {
      return(rcompmargin(Xn,d=d,pos=3,what="data"))
    } else if( what=="direction") {
      if( geo=="rcomp" )
        return(cpt(rcompmargin(rcomp(unclass(dir)+1),d=d,pos=3)))
      else {
        warning("gsi.pltTrafo: Can not marginalise acomp direction with rcomp-marginals")
        return(rcompmargin(rcomp(dir*NA),d=d,pos=3))
      }
    } else if( what == "var" ) {
      if( geo!="rcomp" )
        warning("pltTrafo: Can not correctly marginalise an acomp-variance with rcomp-margins")
      return(rcompmargin(var,d=d,pos=3,what="var"))
    } else if(what =="names") {
      return( c(nam[d],"*") )
    } else stop("unkown transformation type")
  # single variable margin
  } else { 
    if( is.character(margin))
      margin <- match(margin,names(mean))
    if( what == "data" ) {
      Xn <- oneOrDataset(Xn)
      return( structure(clo(cbind(Xn[,-margin,drop=FALSE][,d,drop=FALSE],
                                  Xn[,margin])),class=class(X)))
    } else if(what=="direction") {
#      Xn <- oneOrDataset(Xn)
      Xn <- oneOrDataset(dir)
      return( cbind(Xn[,-margin,drop=FALSE][,d,drop=FALSE],Xn[,margin]))      
    } else if( what =="var" ) {
      if( geo=="rcomp")
        warning("Can not correctly marginalise a rcomp-variances to subcompositions")
      varmar = rbind(cbind(var[-margin,-margin][d,d],var[-margin,margin][d]),
                   cbind(t(var[margin,-margin][d]),var[margin,margin]))
      return(varmar)
    } else if(what =="names") {
      return( c(nam[-margin][d],nam[margin]) )
    } 
  }
}




gsi.isSingleRow <- function(X) {
  return( NROW(X) == 1 || NCOL(X) ==1 )
}




barplot.acomp <- function(height,...,legend.text=TRUE,beside=FALSE,total=1,plotMissings=TRUE,missingColor="red",missingPortion=0.01) {
  nmv <- is.NMV(oneOrDataset(height))
  if( is.null(total) ) {
    X <- oneOrDataset(gsi.plain(ifelse(nmv,height,0)))
  } else {
    X <- oneOrDataset(gsi.plain(rcomp(ifelse(nmv,height,0),total=total)))
  }
  if( plotMissings && missingPortion>0 ) {
    X <- oneOrDataset(ifelse(nmv,X,missingPortion))
  }
  if( !beside && gsi.isSingleRow(X) ) {
    erg <- barplot(t(rbind(X,0)),c(1,0),...,legend.text=legend.text,
            beside=beside)
    ht <- t(X)
    ergDelta <- (erg[2]-erg[1])
    erg <- erg[1]
  }
  else {
    erg <- barplot(ht<-t(X),rep(1,ncol(nmv)),...,legend.text=legend.text,beside=beside)
    ergDelta <- (erg[2]-erg[1])/2
    erg 
  }
  ergDelta<-0.5
  if( plotMissings && any(!nmv) ) {
    if( is.matrix(erg) ) {
      rect(erg[!nmv]-ergDelta,
           0,
           erg[!nmv]+ergDelta,
           t(ht)[!nmv],
           col=missingColor
           )
    } else {
      ergMat <- matrix(erg,nrow=nrow(X),ncol=ncol(X))
      htMat  <- t(apply(rbind(0,ht),2,cumsum))
                                        #recover()
      rect(ergMat[!nmv]-ergDelta,
           htMat[cbind(!nmv,FALSE)],
           ergMat[!nmv]+ergDelta,
           htMat[cbind(FALSE,!nmv)],
           col=missingColor
           )
    }                                  #segments(ergMat[!nmv]-ergDelta,
    #      htMat[!nmv],
    #      ergMat[!nmv]+ergDelta,
    #      htMat[!nmv],
    #      lwd=3,
    #      col=missingColor,pch="-")
  }
  replot(plot=match.call(),add=FALSE)  

  invisible(erg)
}
barplot.rcomp <- barplot.acomp
barplot.ccomp <- barplot.acomp
#barplot.rcomp <- function(height,...,legend.text=TRUE,beside=FALSE,total=1) {
#  X <- height
#  if( gsi.isSingleRow(X) )
#    barplot(t(rbind(gsi.plain(rcomp(X,total=total)),0)),c(1,0),...,legend.text=legend.text,beside=beside)
#  else
#    barplot(gsi.plain(t(rcomp(X,total=total))),...,legend.text=legend.text,beside=beside);
#}
barplot.aplus <- barplot.acomp
#formals(barplot.aplus)$beside <- TRUE
formals(barplot.aplus)[c("beside","total")]  <- list(beside=TRUE,total=NULL)

#barplot.aplus <- function(height,...,legend.text=TRUE,beside=TRUE) {
#  X <- height
#  if( gsi.isSingleRow(X) )
#    barplot(t(rbind(gsi.plain(aplus(X)),0)),c(1,0),...,legend.text=legend.text,
#            beside=beside)
#  else
#    barplot(gsi.plain(t(aplus(X))),...,legend.text=legend.text,beside=beside);
#}

barplot.rplus <- barplot.aplus
#barplot.rplus <- function(height,...,legend.text=TRUE,beside=TRUE) {
#  X <- height
#  if( gsi.isSingleRow(X) )
#    barplot(t(rbind(gsi.plain(rplus(X)),0)),c(1,0),...,
#            legend.text=legend.text,beside=beside)
#  else
#    barplot(gsi.plain(t(rplus(X))),...,legend.text=legend.text,beside=beside);
#}



split.acomp <- function(x,f,drop=FALSE,...) {
  cls <- class(x)
  lapply(split(1:NROW(x),f,drop=drop,...),function(i) structure(x[i,],class=cls))
}
split.rcomp <- split.acomp
split.aplus <- split.acomp
split.rplus <- split.acomp
split.rmult <- split.acomp
split.ccomp <- split.acomp

as.data.frame.acomp <- function(x,...) as.data.frame.matrix(unclass(x))
as.data.frame.rcomp <- function(x,...) as.data.frame.matrix(unclass(x))
as.data.frame.aplus <- function(x,...) as.data.frame.matrix(unclass(x))
as.data.frame.rplus <- function(x,...) as.data.frame.matrix(unclass(x))
as.data.frame.rmult <- function(x,...) as.data.frame.matrix(unclass(x))

gsi.addclass <- function(x,cls) {
  class(x) <- c(cls,attr(x,"class"))
  x
}


gsi <- new.env(hash=TRUE,emptyenv())




princomp.acomp <- function(x,...,scores=TRUE,center=attr(covmat,"center"),
                           covmat=var(x,robust=robust,giveCenter=TRUE),robust=getOption("robust")) {
  cl <- match.call()
  D <- gsi.getD(x)
  tmp <- princomp(clr(x),...,center=clr(center),covmat=covmat,scores=scores)
  tmp$sdev        <- tmp$sdev[-D]
  tmp$loadings    <- structure(tmp$loadings[,-D],class="loadings")
  tmp$Center      <- clrInv(tmp$center)
  tmp$Loadings  <- acomp(clrInv(t(tmp$loadings)),total=D)
  tmp$DownLoadings<- acomp(clrInv(t(-tmp$loadings)),total=D)
  tmp$call <- cl
  gsi.addclass(tmp,"princomp.acomp")
}

print.princomp.acomp <- function(x,...) {
  NextMethod("print",x,...)
  cat("Mean (compositional):\n")
  print(x$Center)
  cat("+Loadings (compositional):\n")
  print(x$Loadings)
  cat("-Loadings (compositional):\n")
  print(x$DownLoadings)
  invisible(x)
}

plot.princomp.acomp <- function(x,y=NULL,...,npcs=min(10,length(x$sdev)),
                                        type=c("screeplot","variance",
                                          "biplot","loadings","relative"),
                                main=NULL,scale.sdev=1) {
  if( missing(main) ) main <- deparse(substitute(x))
  type <- match.arg(type)
  if( type=="biplot" )
    biplot(x,...,main=main)
  else if( type=="loadings" ) {
    if( is.na(scale.sdev) )
      scl<-1
    else
      scl<-scale.sdev*x$sdev
    barplot(acomp((x$Loadings*scl)[1:npcs,]),...,
            main=main,total=gsi.getD(x$Loadings))
  }
  else if(type=="relative") {
    tmp <- relativeLoadings(x,scale.sdev=scale.sdev)[,1:npcs]
    tmp <- barplot(t(tmp),...,main=main,beside=TRUE,legend=TRUE)
    abline(h=1)
    invisible(tmp)
  } else  {
  if( type=="screeplot" ) type <- "lines"
  if( type=="variance" ) type <- "barplot"
  screeplot(x,...,npcs=npcs, main=main,type=type)
  }
}

predict.princomp.acomp <- function(object,newdata,...) {
  NextMethod("predict",object,newdata=clr(newdata),...)
}

                                
#panel.princomp.acomp <- function(x,choice,t,...){
#  straight.panel.acomp(x$Center,x$Loadings)
#}


princomp.rcomp <- function(x,...,scores=TRUE,center=attr(covmat,"center"),
                           covmat=var(x,robust=robust,giveCenter=TRUE),robust=getOption("robust")) {
  cl <- match.call()
  D <- gsi.getD(x)
  tmp <- princomp(cpt(x),...,scores=scores,covmat=covmat,center=cpt(center))
  tmp$sdev        <- tmp$sdev[-D]
  tmp$loadings    <- structure(tmp$loadings[,-D],class="loadings")
  tmp$Center      <- cptInv(tmp$center,x)
  tmp$Loadings    <- cptInv(t(tmp$loadings),x)
  tmp$call <- cl
  gsi.addclass(tmp,"princomp.rcomp")
}


print.princomp.rcomp <- function(x,...) {
  NextMethod("print",x,...)
}

plot.princomp.rcomp <- function(x,y=NULL,...,npcs=min(10,length(x$sdev)),
                                        type=c("screeplot","variance",
                                          "biplot","loadings","relative"),
                                main=NULL,scale.sdev=1) {
  if( missing(main) ) main <- deparse(substitute(x))
  type <- match.arg(type)
  if( type=="biplot" )
    erg<-biplot(x,...,main=main)
  else if( type=="loadings" ) {
    if( is.na(scale.sdev) )
      scl<-1
    else
      scl<-scale.sdev*x$sdev
    erg<-barplot((rmult(t(x$loadings))*scl)[1:npcs,],...,
            main=main)
  }
  else if(type=="relative") {
    tmp <- relativeLoadings(x,scale.sdev=scale.sdev)[,1:npcs]
    tmp <- barplot(t(tmp),...,main=main,beside=TRUE,legend=TRUE)
    erg<-invisible(tmp)
  } else  {
  if( type=="screeplot" ) type <- "lines"
  if( type=="variance" ) type <- "barplot"
  erg<-screeplot(x,...,npcs=npcs, main=main,type=type)
  }
  replot(plot=match.call(),add=FALSE)  
  invisible(erg)
}

predict.princomp.rcomp <- function(object,newdata,...) {
  NextMethod("predict",object,newdata=cpt(newdata),...)
}

princomp.aplus <- function(x,...,scores=TRUE,center=attr(covmat,"center"),
                           covmat=var(x,robust=robust,giveCenter=TRUE),robust=getOption("robust")) {
  cl <- match.call()
  D <- gsi.getD(x)
  tmp <- princomp(ilt(x),...,scores=scores,covmat=covmat,center=ilt(center))
  tmp$Center      <- center
  tmp$Loadings  <- iltInv(t(tmp$loadings))
  tmp$DownLoadings<- iltInv(t(-tmp$loadings))
  tmp$call <- cl
  gsi.addclass(tmp,"princomp.aplus")
}

print.princomp.aplus <- function(x,...) {
  NextMethod("print",x,...)
  cat("Mean (compositional):\n")
  print(x$Center)
  cat("+Loadings (compositional):\n")
  print(x$Loadings)
  cat("-Loadings (compositional):\n")
  print(x$DownLoadings)
  invisible(x)
}

plot.princomp.aplus <- function(x,y=NULL,...,npcs=min(10,length(x$sdev)),
                                        type=c("screeplot","variance",
                                          "biplot","loadings","relative"),
                                main=NULL,scale.sdev=1) {
  if( missing(main) ) main <- deparse(substitute(x))
  type <- match.arg(type)
  if( type=="biplot" )
    erg<-biplot(x,...,main=main)
  else if( type=="loadings" ) {
        if( is.na(scale.sdev) )
      scl<-1
    else
      scl<-scale.sdev*x$sdev
    erg<-barplot(aplus((x$Loadings*scl)[1:npcs,]),...,
            main=main)
  }
  else if(type=="relative") {
    tmp <- relativeLoadings(x,scale.sdev=scale.sdev)[,1:npcs]
    tmp <- barplot(t(tmp),...,main=main,beside=TRUE,legend=TRUE)
    abline(h=1)
    erg<-invisible(tmp)
  } else  {
  if( type=="screeplot" ) type <- "lines"
  if( type=="variance" ) type <- "barplot"
  erg<-screeplot(x,...,npcs=npcs, main=main,type=type)
  }
  replot(plot=match.call(),add=FALSE)  
  invisible(erg)
}

predict.princomp.aplus <- function(object,newdata,...) {
  NextMethod("predict",object,newdata=ilt(newdata),...)
}

princomp.rplus <- function(x,...,scores=TRUE,center=attr(covmat,"center"),
                           covmat=var(x,robust=robust,giveCenter=TRUE),robust=getOption("robust")) {
  cl <- match.call()
  robust<-robust
  tmp <- princomp(iit(x),...,scores=scores,covmat=covmat,center=iit(center))
  tmp$Center      <- iitInv(center)
  tmp$call <- cl
  tmp$Loadings <- rmult(t(tmp$loadings))
  gsi.addclass(tmp,"princomp.rplus")
}

print.princomp.rplus <- function(x,...) {
  NextMethod("print",x,...)
  cat("Mean:\n")
  print(x$Center)
  cat("Loadings:\n")
  print(x$Loadings)
  invisible(x)
}

plot.princomp.rplus <- function(x,y=NULL,...,npcs=min(10,length(x$sdev)),
                                        type=c("screeplot","variance",
                                          "biplot","loadings","relative"),
                                main=NULL,scale.sdev=1) {
  if( missing(main) ) main <- deparse(substitute(x))
  type <- match.arg(type)
  if( type=="biplot" )
    erg<-biplot(x,...,main=main)
  else if( type=="loadings" ) {
    if( is.na(scale.sdev) )
      scl<-1
    else
      scl<-scale.sdev*x$sdev
    erg<-barplot((rmult(t(x$loadings))*scl)[1:npcs,],...,
            main=main)
  }
  else if(type=="relative") {
    tmp <- relativeLoadings(x,scale.sdev=scale.sdev)[,1:npcs]
    tmp <- barplot(t(tmp),...,main=main,beside=TRUE,legend=TRUE)
    erg<-invisible(tmp)
  } else  {
  if( type=="screeplot" ) type <- "lines"
  if( type=="variance" ) type <- "barplot"
  erg<-screeplot(x,...,npcs=npcs, main=main,type=type)
  }
  replot(plot=match.call(),add=FALSE)  
  invisible(erg)
}

predict.princomp.rplus <- function(object,newdata,...) {
  NextMethod("predict",object,newdata=iit(newdata),...)
}


princomp.rmult <- function(x,cor=FALSE,scores=TRUE,
                           covmat=var(rmult(x[subset,]),robust=robust,giveCenter=TRUE),center=attr(covmat,"center"),  subset = rep(TRUE, nrow(x)),...,robust=getOption("robust")) {
Nprincomp <-
    function(x, cor = FALSE, scores = TRUE, covmat = NULL,center=center,
             subset = rep(TRUE, nrow(as.matrix(x))), ...)
{
    cl <- match.call()
    cl[[1]] <- as.name("princomp")
    #if(!missing(x) && !missing(covmat))
    #    warning("both 'x' and 'covmat' were supplied: 'x' will be ignored")
    z <- as.matrix(unclass(x))[subset, , drop = FALSE]
    cv <- covmat
    n.obs <- nrow(x) # NA is maybe more appropriate
    cen <- center
    dn <- dim(z)
    if(dn[1] < dn[2])
      stop("'princomp' can only be used with more units than variables")
    if(!is.numeric(cv)) stop("PCA applies only to numerical variables")
    if (cor) {
        sds <- sqrt(diag(cv))
        if(any(sds == 0))
            stop("cannot use cor=TRUE with a constant variable")
        cv <- cv/(sds %o% sds)
    }
    edc <- eigen(cv, symmetric = TRUE)
    ev <- edc$values
    if (any(neg <- ev < 0)) { # S-PLUS sets all := 0
        ## 9 * : on Solaris found case where 5.59 was needed (MM)
        if (any(ev[neg] < - 9 * .Machine$double.eps * ev[1]))
            stop("covariance matrix is not non-negative definite")
        else
            ev[neg] <- 0
    }
    cn <- paste("Comp.", 1:ncol(cv), sep = "")
    names(ev) <- cn
    dimnames(edc$vectors) <- if(missing(x))
        list(dimnames(cv)[[2]], cn) else list(dimnames(x)[[2]], cn)
    sdev <- sqrt(ev)
    sc <- if (cor) sds else rep(1, ncol(cv))
    names(sc) <- colnames(cv)
    scr <- #if (scores && !missing(x) && !is.null(cen))
      (unclass((rmult(z)-rmult(cen))/sc))%*%edc$vectors
    #  scale(z, center = cen, scale = sc) %*% edc$vectors
    if (is.null(cen)) cen <- rep(NA_real_, nrow(cv))
    edc <- list(sdev = sdev,
                loadings = structure(edc$vectors, class="loadings"),
                center = c(unclass(cen)), scale = sc, n.obs = n.obs,
                scores = scr, call = cl)
    class(edc) <- "princomp"
    edc
}
Nprincomp(unclass(x),cor=cor,scores=scores,covmat=covmat,center=unclass(center),subset=subset,...)
}

gsi.pairrelativeMatrix <- function(names) {
  n  <- length(names)
  ii <- rep(1:n,n)
  jj <- rep(1:n,each=n)
  jgi <- jj>ii
  ii <- ii[jgi]
  jj <- jj[jgi]
  N <- length(ii)
  erg <- rep(0,N*n)
  erg[1:N+N*(ii-1)]<-  1
  erg[1:N+N*(jj-1)]<- -1
  erg <- matrix(erg,nrow=N,ncol=n)
  colnames(erg)  <- names
  row.names(erg) <- paste(names[ii],names[jj],sep="/") 
  erg
}

relativeLoadings <- function(x,...) UseMethod("relativeLoadings",x)
relativeLoadings.princomp.acomp <- function(x,...,log=FALSE,scale.sdev=NA,cutoff=0.1) {
  if( is.na(scale.sdev) ) 
    scale <- rep(1,length(x$sdev))
  else 
    scale <- scale.sdev*x$sdev
  bl <- gsi.pairrelativeMatrix(row.names(x$loadings)) %*% (unclass(x$loadings) %*% diag(scale))
  colnames(bl) <- colnames(x$loadings)
  if( !log )
    bl <- exp(bl)
  structure(bl,class="relativeLoadings.princomp.acomp",cutoff=cutoff,scale=scale,log=log)
}

relativeLoadings.princomp.aplus <- function(x,...,log=FALSE,scale.sdev=NA,cutoff=0.1) {
  if( is.na(scale.sdev) ) 
    scale <- rep(1,length(x$sdev))
  else 
    scale <- scale.sdev*x$sdev
  bl <- gsi.pairrelativeMatrix(row.names(x$loadings)) %*% (unclass(x$loadings) %*% diag(scale))
  colnames(bl) <- colnames(x$loadings)
  if( !log )
    bl <- exp(bl)
  structure(bl,class="relativeLoadings.princomp.aplus",cutoff=cutoff,scale=scale,log=log)
}

relativeLoadings.princomp.rcomp <- function(x,...,scale.sdev=NA,cutoff=0.1) {
  if( is.na(scale.sdev) ) 
    scale <- rep(1,length(x$sdev))
  else 
    scale <- scale.sdev*x$sdev
  bl <- gsi.pairrelativeMatrix(row.names(x$loadings)) %*% (unclass(x$loadings) %*% diag(scale))
  colnames(bl) <- colnames(x$loadings)
  structure(bl,class="relativeLoadings.princomp.rcomp",cutoff=cutoff,scale=scale)
}

relativeLoadings.princomp.rplus <- function(x,...,scale.sdev=NA,cutoff=0.1) {
  if( is.na(scale.sdev) ) 
    scale <- rep(1,length(x$sdev))
  else 
    scale <- scale.sdev*x$sdev
  bl <- gsi.pairrelativeMatrix(row.names(x$loadings)) %*% (unclass(x$loadings) %*% diag(scale))
  colnames(bl) <- colnames(x$loadings)
  structure(bl,class="relativeLoadings.princomp.rplus",cutoff=cutoff,scale=scale)
}


print.relativeLoadings.princomp.acomp <- function(x,...,
                                                 cutoff=attr(x,"cutoff"),
                                                 digits=2
                                                 ) {
  bt <- format(x,digits=digits)
  if( !attr(x,"log") ) { 
    if( !is.na(cutoff) )
      bt[abs(log(x))<cutoff] <- ""
  } else {
    if( !is.na(cutoff) )
      bt[abs(x)<cutoff] <- ""
  }
  print(bt,quote=FALSE)
  x
}

print.relativeLoadings.princomp.aplus <- print.relativeLoadings.princomp.acomp

print.relativeLoadings.princomp.rcomp <- function(x,...,
                                                 cutoff=attr(x,"cutoff"),
                                                 digits=2
                                                 ) {
  bt <- format(x,digits=digits)
  bt[abs(x)<cutoff] <- ""
  print(bt,quote=FALSE)
  x
}

print.relativeLoadings.princomp.rplus <- function(x,...,
                                                 cutoff=attr(x,"cutoff"),
                                                 digits=2
                                                 ) {
  bt <- format(x,digits=digits)
  bt[abs(x)<cutoff] <- ""
  print(bt,quote=FALSE)
  x
}


plot.relativeLoadings.princomp.acomp<- function(x,...) {
  barplot(t(x),...,beside=TRUE)
  if( !attr(x,"log") )
    abline(h=1)
  replot(plot=match.call(),add=FALSE)  
  invisible(x)
}
plot.relativeLoadings.princomp.aplus<- function(x,...) {
  barplot(t(x),...,beside=TRUE)
  if( !attr(x,"log") )
    abline(h=1)
  replot(plot=match.call(),add=FALSE)  
  invisible(x)
}

plot.relativeLoadings.princomp.rcomp<- function(x,...) {
  barplot(t(x),...,beside=TRUE)
  replot(plot=match.call(),add=FALSE)  
  invisible(x)
}

plot.relativeLoadings.princomp.rplus<- function(x,...) {
  barplot(t(x),...,beside=TRUE)
  replot(plot=match.call(),add=FALSE)  
  invisible(x)
}



read.geoeas <- function (file){
#read title
print("reading title:...")
title <- read.table(file, nrows=1, quote="", colClasses="character", sep="\t")
print("reading title: OK")

#read number of variables
print("reading number of variables:...")
nvars <- read.table(file, skip=1, nrows=1, sep="\t")
nvars <- as.integer(nvars)
#read labels of the variables
print("reading variables:...")
labels <- scan(file, skip=2, nmax=nvars, nlines= nvars, sep="\t", quote="", what="character")

#read data matrix
print("reading dataset:...")
dades <- read.table(file, skip=2+nvars)
print("reading dataset: OK")
#relate variables with their names
colnames(dades)=as.vector(labels)
#add title as an attribute
attr(dades,"title") <- as.character(title)
return(dades)
}

read.geoEAS <- function(file){ read.geoeas(file) }



#read.standard <- function (file){
##read title
#print("reading title:...")
#title <- read.table(file, nrows=1, quote="", colClasses="character", sep="\t")
#print("reading title: OK")

#read data matrix
#print("reading dataset:...")
#dades <- read.table(file, skip=1, header=T)
#print("reading dataset: OK")#

#add title as an attribute
#attr(dades,"title") <- as.character(title)
#return(dades)
#}

 #  =========================================================================#
# Based on the the tedrahedron plot from MixeR
# http://vlado.fmf.uni-lj.si/pub/MixeR/
# Vladimir Batagelj & Matevz Bren
# plots mix object 'm' into tetrahedron and exports it in kinemage format
#   clu - partition -> color of points
#   vec - vector of values -> size of points
#   col - color of points if clu=NULL
#   scale - relative size of points
#   king - FALSE (for Mage); TRUE (for King)
kingTetrahedron <- function(X,parts=1:4,file="tmptetrahedron.kin",clu=NULL,vec=NULL,king=TRUE,scale=0.2,col=1,title="Compositional Tetrahedron"){
  m <- list()
  m$mat <- data.frame(clo(X,parts=parts))
  m$tit <- title
  kinColors <- c('white', 'red', 'blue', 'yellow',
    'green', 'cyan', 'magenta', 'pink', 'lime',
    'orange', 'peach', 'gold', 'purple', 'sea',
    'gray', 'sky', 'brown', 'lilac', 'hotpink',
    'yellowtint', 'pinktint', 'peachtint',
    'greentint', 'bluetint', 'lilactint',
    'deadwhite', 'deadblack', 'invisible')
  warn <- ""
  if (king) warn <- "\n*** works only with KiNG viewer: http://kinemage.biochem.duke.edu/"
  head <- paste("@text\n",
    file," / ",date(),
    "\nDataset: ", m$tit,warn,
"\nLayout obtained using composition
based on  MixeR
http://vlado.fmf.uni-lj.si/pub/MixeR/
Vladimir Batagelj & Matevz Bren
Institute of Mathematics, Physics and Mechanics
Ljubljana, Slovenia
@kinemage 1
@caption\n", m$tit,
"\n@fontsizeword 18
@zoom 1.00
@thinline
@zclipoff
@group {Complete} dominant animate movieview = 1
@spherelist 0 radius= 0.20\n",sep="")
  foot <-
"@vectorlist {}  color=  blue
P 9.500, 0.500, 9.500
0.500, 9.500, 9.500
P 9.500, 0.500, 9.500
0.500, 0.500, 0.500
P 9.500, 0.500, 9.500
9.500, 9.500, 0.500
P 0.500, 9.500, 9.500
0.500, 0.500, 0.500
P 0.500, 9.500, 9.500
9.500, 9.500, 0.500
P 0.500, 0.500, 0.500
9.500, 9.500, 0.500\n"
  kin <- file(file,"w")
  n <- nrow(m$mat)
  if (is.null(clu)) {clu <- rep(col,n)}
  clu <- c(clu,0,0,0,0)
  if (is.null(vec)) {vec <- rep(1,n)}
  vec <- c(vec,1,1,1,1)
  size <- scale/max(vec)
  rn <- c(dimnames(m$mat)[[1]],"A","B","C","D")
  a <- 10/9+c(m$mat[,1],1,0,0,0)
  b <- c(m$mat[,2],0,1,0,0)
  c <- c(m$mat[,3],0,0,1,0)
  d <- c(m$mat[,4],0,0,0,1)
  x <- (a - b - c + d)*0.45
  y <- (a - b + c - d)*0.45
  z <- (a + b - c - d)*0.45
  cat(head,file=kin)
  for (i in 1:length(rn)) {
    color <- "white"
    if (clu[i]>0) color <- kinColors[2+(clu[i]-1)%%18]
    cat("{",rn[i],"} ", color," ",file=kin,sep="")
    if (king) cat(" r=", vec[i]*size," ",file=kin,sep=" ")
    cat(10*x[i],10*(1-y[i]),10*z[i],"\n",file=kin,sep=" ")
  }
  cat(foot,file=kin)
  close(kin)
}


gsi.getBalStruct <- function(descr,names,allowMinus=FALSE,allowOne=FALSE) {
    if( !is.call(descr))
      return(structure(list(),all=list(as.character(descr))))
    command <- as.character(descr[[1]])
                                        # next two lines added!
    if( command == "~" )
      return( Recall(descr[[2]]) )
    if( command == "(" )
      return( Recall(descr[[2]]) )
    if( allowOne && command == "1" )
      return( structure(list(),all=list(as.character(descr))) )
    else if( (allowMinus && command == "-") || command == "/" || command =="*" || command =="+" || command==":") {
      a <- Recall(descr[[2]])
      b <- Recall(descr[[3]])
      aall <- unlist(c(attr(a,"all"),attr(b,"all")))
      if( command == "/" || (command == "-" && allowMinus) ) {
        erg <- c( list(list(attr(a,"all"),attr(b,"all"))),
                 c(a),c(b)
                 )
      } else {
        erg <- c( c(a),c(b) )
        
      }
      attr(erg,"all") <- aall
      erg
    } else {
      structure(list(),all=list(as.character(command)))
    }
  }


balanceBase <- function(X,...) UseMethod("balanceBase",X)

balanceBase.acomp <- function(X,expr,...) {
  if( !is.character(X) )
    X <- colnames(oneOrDataset(X))
  bs <- gsi.getBalStruct(expr)
  nam <- sapply(bs,function(k) {
    paste(paste(unlist(k[[1]]),sep="",collapse=""),
          paste(unlist(k[[2]]),sep="",collapse=""),sep="/")
  })
  bas <- sapply(bs,function(k) {
    x <- rep(0,length(X))
    x[match(k[[1]],X)]<-1/length(k[[1]])
    x[match(k[[2]],X)]<- - 1/length(k[[2]])
    x <- x/sqrt(sum(x^2))
    x
  })
  colnames(bas)<-nam
  rownames(bas)<-X
  bas
}

balanceBase.aplus <- function(X,expr,...) {
  if( !is.character(X) )
    X <- colnames(oneOrDataset(X))
  bs <- gsi.getBalStruct(expr,allowOne=TRUE)
  nam <- sapply(bs,function(k) {
    paste(paste(unlist(k[[1]]),sep="",collapse=""),
          paste(unlist(k[[2]]),sep="",collapse=""),sep="/")
  })
  bas <- sapply(bs,function(k) {
    x <- rep(0,length(X))
    if(length(k[[1]])>0) x[match(k[[1]],X)]<-1/length(k[[1]])
    if(length(k[[2]])>0) x[match(k[[2]],X)]<- - 1/length(k[[2]])
    x <- x/sqrt(sum(x^2))
    x
  })
  colnames(bas)<-nam
  rownames(bas)<-X
  bas
}


balanceBase.rcomp <- function(X,expr,...) {
  if( !is.character(X) )
    X <- colnames(oneOrDataset(X))
  bs <- gsi.getBalStruct(expr,allowMinus=TRUE)
  nam <- sapply(bs,function(k) {
    paste(paste(unlist(k[[1]]),sep="",collapse=""),
          paste(unlist(k[[2]]),sep="",collapse=""),sep="-")
  })
  bas <- sapply(bs,function(k) {
    x <- rep(0,length(X))
    x[match(k[[1]],X)]<-1
    x[match(k[[2]],X)]<- -1
    x
  })
  colnames(bas)<-nam
  rownames(bas)<-X
  bas
}

balanceBase.rplus <- function(X,expr,...) {
  if( !is.character(X) )
    X <- colnames(oneOrDataset(X))
  bs <- gsi.getBalStruct(expr,allowMinus=TRUE,allowOne=TRUE)
  nam <- sapply(bs,function(k) {
    paste(paste(unlist(k[[1]]),sep="",collapse=""),
          paste(unlist(k[[2]]),sep="",collapse=""),sep="-")
  })
  bas <- sapply(bs,function(k) {
    x <- rep(0,length(X))
    if(length(k[[1]])>0) x[match(k[[1]],X)]<-1
    if(length(k[[2]])>0) x[match(k[[2]],X)]<- -1
    x
  })
  colnames(bas)<-nam
  rownames(bas)<-X
  bas
}


balance <- function(X,...) UseMethod("balance",X)

balance.acomp <- function(X,expr,...) {
  ilr(X,balanceBase(X,expr))
}

balance.aplus <- function(X,expr,...) {
  clr2ilr(ilt(X),balanceBase(X,expr))
}

balance.rcomp <- function(X,expr,...) {
  clr2ilr(iit(rcomp(X)),balanceBase(X,expr))
}

balance.rplus <- function(X,expr,...) {
  clr2ilr(ilt(X),balanceBase(X,expr))
}

balance01 <- function(X,...) UseMethod("balance01",X)

balance01.acomp <- function(X,expr,...) {
  ilr(X,balanceBase(X,expr))
}


balance01.acomp <- function(X,expr,...) {
bal01 <- function(X) {
  A <- exp(unclass(oneOrDataset(X))*sqrt(2))
  A / (1+A)
}
bal01(balance(X,expr))
}


balance01.rcomp <- function(X,expr,...) {
  bb <- balanceBase(X,expr)
  p <- clr2ilr(iit(rcomp(X)),ifelse(bb>0,1,0))
  t <- clr2ilr(iit(rcomp(X)),ifelse(bb!=0,1,0))
  p/t
}
  
print.acomp <- function(x,...,replace0=TRUE) {
pad <- function(x,n) {
  maxN <- max(n)
  padding <- substring(rep(paste(rep(" ",maxN),collapse="",sep=""),max(length(n),length(x))),1,pmax(n-nchar(x),0))
  paste(padding,x,sep="")
}

subser <- function(erg,x) {
  at <- attributes(x)
  nc <- nchar(erg)
  if( any(nc<4) )
    erg[nc<4]<-pad(erg[nc<4],4)  
  nc <- nchar(erg)
  erg<-gsub(" NaN$"," MAR",erg)
  erg<-gsub("  NA$","MNAR",erg)
  erg<-gsub("-Inf$","  SZ",erg)
  zero <- is.finite(x)&x==0
  if( any(zero) )
    erg[zero]<-pad("BDL",nchar(erg[zero]))
  erg<-gsub("^ *-","<",erg)
  erg<-gsub("Inf$","ERR",erg)
  attributes(erg) <- at
  erg
}

if( is.matrix(x) )
     erg<- apply(x,2,format,...)
  else 
     erg<-format(x,...)
erg <- subser(erg,x)
if( !is.null(attr(x,"Losts"))) {
  L <- oneOrDataset(attr(x,"Losts"))
   attr(erg,"Losts")<-noquote(c("..."=paste("Containing ", sum(L)," lost and replaced values!")))
 }

                             
                                        #  dm <-dim(erg)
                                        #  dim(erg)<-dm
                                        #  dimnames(erg)<-dimnames(x)
quote <- FALSE
if( !is.null(list(...)$quote ))
  quote <- list(...)$quote
print.default(x=erg,...,quote=FALSE)
invisible(x)
}



print.aplus<- print.acomp
print.rcomp<- print.acomp
formals(print.rcomp)$replace0 <- FALSE
print.rplus<- print.rcomp

#"(.acomp" <- function(x,...) { 
#  if(is.matrix(x))
#     return(acomp(x[...]))
#  else
#     return(acomp(oneOrDataset(x)[...]))
#  acomp(oneOrDataset(x)[...])
#}

AitchisonDistributionIntegrals <- function(
         theta=alpha+sigma %*% clr(mu),
         beta=-1/2*gsi.svdinverse(sigma),
         alpha=mean(theta),
         mu=clrInv(c(sigma%*%(theta-alpha))),
         sigma=-1/2*gsi.svdinverse(beta),
         grid=30,
         mode=3) {
if( !xor( missing(theta) , (missing(mu) || missing(alpha)) ) )
  stop("Specify either theta or mu and alpha")
if( !xor( missing(beta) , missing(sigma) ) )
  stop("Specify either beta or sigma")
D <- length(theta)
if( nrow(beta) == D-1 )
  beta <- ilrvar2clr(beta)
if( any(abs(beta-t(beta))>1E-10) ) {
  warning("AitchisonDistributionIntegrals: beta was not symmetric 1")
  print(beta)
}
gsiInt(dim(beta),2)
stopifnot( length(dim(beta))==2,ncol(beta)==D,nrow(beta)==D)
erg <- .C(gsiAitchisonDistributionIntegral,
          D   =gsiInt(D,1),
          grid=gsiInt(grid,1),
          mode=gsiInt(mode,1),
          theta=gsiDouble(theta,D),
          beta =gsiDouble(beta,D*D),
          expKappa =numeric(1),
          loggxMean=numeric(1),
          clrMean  =numeric(D),
          clrVar   =numeric(D*D)
          )
erg$beta <- matrix(erg$beta,nrow=D)
erg$SqIntegral <- matrix(erg$clrVar,nrow=D,ncol=D)
erg$alpha=alpha
erg$mu=mu
erg$sigma=sigma
erg$clrSqExpectation <- matrix(erg$clrVar,nrow=D,ncol=D)
dim(erg$clrVar) <- c(D,D)
  return(erg[c("theta","beta","alpha","mu","sigma",if( mode>=0 ) c("expKappa","loggxMean") else c(),if(mode >= 1) "clrMean" else c(),if(mode==2) "clrSqExpectation" else if(mode>=3) "clrVar" else c())])
}



dAitchison <- function(x,theta=alpha+sigma %*% clr(mu),beta=-1/2*gsi.svdinverse(sigma),alpha=mean(theta),mu=clrInv(c(sigma%*%(theta-alpha))),sigma=-1/2*gsi.svdinverse(beta),grid=30,
realdensity=FALSE,expKappa=AitchisonDistributionIntegrals(theta,beta,grid=grid,mode=1)$expKappa) {
if( !xor( missing(theta) , (missing(mu) || missing(alpha)) ) )
  stop("Specify either theta or mu and alpha")
if( !xor( missing(beta) , missing(sigma) ) )
  stop("Specify either beta or sigma")
D <- length(theta)
if( any(abs(beta-t(beta))>1E-10) ) {
  warning("dAitchison: beta was not symmetric 1")
  print(beta)
}
if( nrow(beta) == D-1 )
  beta <- ilrvar2clr(beta)
if( any(abs(beta-t(beta))>1E-10) ) {
  warning("dAtichison: beta was not symmetric 2")
  print(beta)
}
stopifnot(gsi.getD(x)==D)
clrx <- clr(x)
cf <- if(realdensity) 1 else 0
exp((clrx%*%beta)%*%clrx+ult(x)%*%rmult(theta-cf))/expKappa
}


rAitchison <- function(n,theta=alpha+sigma %*% clr(mu),beta=-1/2*gsi.svdinverse(sigma),alpha=mean(theta),mu=clrInv(c(sigma%*%(theta-alpha))),sigma=-1/2*gsi.svdinverse(beta),withfit=FALSE) {
#withfit=FALSE # in the future, this should allow to simulate more efficiently by playing with the decomposition Ait = Normal x Dirichlet
if( !xor(missing(theta),missing(mu) || missing(alpha)) )
  stop("Specify either theta or mu and alpha")
if( !xor(missing(beta),missing(sigma)) )
  stop("Specify either beta or sigma")
if( !missing(theta) ) nam <- names(theta) else if( !missing(mu) ) nam<- names(mu)
D <- length(theta)
if( nrow(sigma) == D-1 )
  sigma <- ilrvar2clr(sigma)
  else
  sigma <- ilrvar2clr(clrvar2ilr(sigma))
# Prepare Sigma
if( any(abs(sigma-t(sigma))>1E-10) )
  warning("rAtichison: sigma was not symmetric")
SVD <- svd(sigma)
if( any(SVD$d < -1E-8 ) )
  warning("rAitchison currently only works correctly with positive semidefinit sigmas. Results wrong!!!")
sqrtSigma <- with(SVD,u%*% gsi.diagGenerate(sqrt(d)) %*% t(v))
if( withfit ) {
     # find the best decomposition Ait = Normal x Dirichlet; where both have the same mode
	bestfit <- gsiFindSolutionWithDerivative(
                           function(alpha) exp(alpha)+sigma%*%alpha-theta,
                           function(alpha) diag(exp(alpha))+sigma,
                           c(theta),
                           iter=20)
	if( attr(bestfit,"status")!="ok" ) {
	warning("Problems in finding optimal simulation algorith, using fallback")
	}
	# Compute best fitter
	mu <- bestfit - sum(bestfit)
	SimTheta <- exp(bestfit)
	stopifnot(abs(sum(SimTheta)-sum(theta))<1E-6)
} else {
    # use as normal(0,sigma) and as dirichlet(theta)
	SimTheta <- theta
	mu <- theta*0
	# SimTheta <- rep(alpha, length(theta))
	# mu <- mu
}
if( ! all(SimTheta>0) ) {
  if( all(theta>0) ) {
	SimTheta <- theta
	mu <- theta*0
        warning("rAitchison: withfit ignored");
      } else {
        stop("rAitchison: This implementation only works with positive theta")
      }
}
# Compute with rejection sampling
erg <- .C(gsirandomClr1Aitchison,
   D=gsiInt(D,1),
   n=gsiInt(n,1),
   erg=numeric(D*n),
   theta=gsiDouble(SimTheta,D),
   mu=gsiDouble(mu,D),
   sqrtSigma=gsiDouble(sqrtSigma,D*D)
   )
# Format result from CLR
res <- clrInv(matrix(erg$erg,nrow=n))
names(res) <- nam
res
}


gsiFindSolutionWithDerivative <- function(f,Der,start,iter=30) {
  nstart <- start
  try({
    it <- 0
    for( i in 1:iter) {
      y  <- f(c(nstart))
      ny = norm(c(y))
      if( i == 1) firstnorm = ny
      Div  <- Der(c(nstart))
      nstart <- nstart-solve(Div,y)
      it <- i
      if( ny <1E-14 && ny/firstnorm<1E-6 )
        break
    }
    return(structure(nstart,value=y,status=if(ny/firstnorm<1E-6 || ny <1E-15) "ok" else "not converged",iterations=it))
  },silent=FALSE)
  return(structure(start,value=f(c(start)),status="failed",iterations=it))
}

R2 <- function(object,...) UseMethod("R2",object)
R2.lm <- function(object,...,adjust=TRUE,ref=0) {
  if( !is.numeric(ref))
    ref<-Recall(ref,...,adjust=adjust)
  pr <- predict(object)
  re <- resid(object)
  n  <- nrow(re)
    if (is.null(n)) n <- length(re) # consider the case of one-dimensional response
  y  <- pr+re
  erg <- if( adjust ) {
    dfres<-object$df.residual
    1-(mvar(re)*(n-1)/dfres)/mvar(y)
  } else
    1-mvar(re)/mvar(y)
  (erg-ref)/(1-ref)
}
R2.default <- function(object,...,ref=0) {
  if( !is.numeric(ref))
    ref<-Recall(ref,...)
  pr <- predict(object)
  re <- resid(object)
  n  <- nrow(re)
    if (is.null(n)) n <- length(re)  # consider the case of one-dimensional response
  y  <- pr+re
  erg <- 1-mvar(re)/mvar(y)
  (erg-ref)/(1-ref)
}
var.mlm <- function(x,...) {r<-unclass(resid(x));(t(r)%*%r)/x$df.residual} 
var.lm <- function(x,...) {r<-unclass(resid(x));sum(r^2)/x$df.residual} 


vcovAcomp <- function(object,...){
  co <- coef(object)
  aperm(structure(vcov(object,...),
        dim=c(dim(co),dim(co)),
        dimnames=c(dimnames(co),dimnames(co))),
        c(2,4,1,3))
}
qHotellingsTsq <- function(p,n,m){
  qf(p,n,m-n+1)*(n*m)/(m-n+1)
}
pHotellingsTsq <- function(q,n,m){
  pf(q/((n*m)/(m-n+1)),n,m-n+1)
}


ConfRadius <- function(model,prob=1-alpha,alpha) {
  sqrt(qHotellingsTsq(prob,ncol(coef(model)),model$df.residual))
}

