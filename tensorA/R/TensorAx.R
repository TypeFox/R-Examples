
#.First.lib <- function(lib,pkg) {
#  library.dynam("tensorA",pkg,lib)
#}

#.Last.lib <- function(libpath) {
#  library.dynam.unload("tensorA",libpath)
#}

gsi.cat <- function(...,l=list(...)) {
  n <- names(l)
  if( length(l)>1 || !is.null(names(l)) ) {
    cat("{")
    for( i in 1:length(l) ) {
      if( !is.null(n) ) cat(n[i],"=",sep="")
      gsi.cat(l=l[[i]])
      cat(" ")
    }
    cat("}")
  } else if( length(l)>0 )
    try(cat(l[[1]]))
}

if( FALSE ) {
  gsi.debug <- function(...) {try(gsi.cat(...))}
  gsi.debugn<- function(...) {try(gsi.cat(...))}
  gsi.debugf<- function(...) {print(match.call("",call=sys.call(sys.parent())))}
  gsi.debugr<- function(X) {
    cat(as.character(as.list(sys.call(sys.parent()))[[1]]),"\n");
    gsi.cat(dim(X),"\n");
    X
  }
} else {
  gsi.debug <- function(...) {list(...)}
  gsi.debugn<- function(...) {}
  gsi.debugf <- function(...) {}
  gsi.debugr <- function(X) {X}
}

# Some nice conveniance functions for matrices not worth a package
#Aidx<- function(n,m) {rep(1:n,m)}
#Bidx<- function(n,m) {rep(1:m,each=n)}

# Modify the formal arguments of a function
gsi.setarg <-function(fun,...) {
  env <- environment(fun)
  P <- formals(fun)
  L <- list(...)
  P[names(L)]<-L
  formals(fun)<-P
  environment(fun)<-env
  fun
}
#
# The index described by what is restructured to a tensor of the structure
# given by dims.    
#

to.tensor <- function(X,...) UseMethod("to.tensor",X)

  
to.tensor.default  <- function(X,dims=NULL,ndimnames=NULL,what=1,addIndex=FALSE,...){
                                        # Make what the first unstructured index 
  gsi.debugn("to.tensor$dims=",dims)
  if( !missing(what) ) 
    X <- untensor(X,what)
                                        # Process prior structure
  if( is.null(d<- dim(X)) )
    dim(X) <- d <- c("I"=length(X))
  if( is.null( odimnames <- dimnames(X) ) ) 
    odimnames <- rep(list(NULL),length(d))
  else
    attr(X,"dimnames") <- NULL
  if( missing(dims) ) {
    if( is.null(ndimnames) ) {
      if( is.null(names(dim(X))) )
        if( is.null(names(dimnames(X))) )
          names(dim(X))<-gsi.stdnames(length(dim(X)),"I")
        else
          names(dim(X))<-names(dimnames(X))
      names(odimnames) <- names(dim(X))
      attr(X,"dimnames") <- odimnames
      if( length(dim(X))==2 ) 
        class(X) <- c("tensor","matrix")
      else
        class(X) <- "tensor"
      return(gsi.debugr(X))
    } else {
      dims <- sapply(ndimnames,length) 
      names(dims) <- names(ndimnames)
    } 
  } else if( is.list(dims) ) {
    ndimnames <- dims
    dims <- sapply(ndimnames,length)
    names(dims) <- names(ndimnames)
  }
  
                                        # Process new structure
                                        #   Process dimension names
  if( !is.null(dims) && length(dims) == 0 && length(X) <2 )
    return( c(X) )
    
  newnames <- names(dims)
  if( is.null(newnames) ) {
    if( ! is.null(ndimnames) ) 
      newnames <- names(ndimnames)
    if( is.null(newnames))
      newnames <- gsi.stdnames(length(dims),"I",avoid=names(dim(X)[-1]))
    names(dims) <- newnames
  }
  
  pl <- prod(dims)
                                        # Repeat elements in case needed
  if( d[1]==1 ){
    X <- rep(c(X),each=pl)
    d[1]<-pl
  }
                                        # Check
  if( is.null(ndimnames) )
    ndimnames <- rep(list(NULL),length(dims))
  if( length(ndimnames) != length(dims) )
    stop("ndimnames and dims don't match: ",length(ndimnames)," ",length(dims))
                                        # Define new structure
  if( addIndex==FALSE  ) {
    if( d[1] != prod(dims))
      stop("to.tensor: Not right number of elements (1)",d[1]," ",dims)
    dim(X) <- c(dims,d[-1])
    dimnames(X) <- c(ndimnames,odimnames[-1])
  } else {
    if( d[1] %% prod(dims))
      stop("to.tensor: Not right number of elements (2)",d[1]," ",dims)
    dim(X) <- c(dims,"@"=d[1]/pl,d[-1])
    if( is.character(addIndex) )
      names(dim(X))[length(dims)+1] <- addIndex
    dimnames(X) <- c(ndimnames,list(NULL),odimnames[-1])
  }
                                        # Return
  names(dimnames(X))<-names(dim(X))
  if(  length(dim(X)) == 2 )
    class(X) <-c("tensor","matrix")
  else
    class(X) <- c("tensor")
  if( length(names(X))!=length(unique(names(X))))
    warning("Tensor with duplicated names generated: ", paste(names(X),col=" "))
  gsi.debugr(X)
}

ftable.tensor <- function(x,...) {
  n  <- length(dim(x))
  dn <- dimnames(x)
  nams <- names(dim(x))
  if( is.null(nams) )
    nams <- letters[1:n]
  if( is.null(dn) )
    dn <- rep(list(NULL),n)
  dimnames(x) <- lapply(1:n,function(i) {
    if( is.null(dn[[i]]) )
      gsi.stdnames(dim(x)[i],nams[i])
    else
      dn[[i]]
  })
  class(x) <- "table"
  ftable(x,...)
}

names.tensor <- function(x) {
  return(names(dim(x)))
}

"names<-.tensor" <- function(x,value) {
  if( length(value) != length(unique(value)))
    warning("Duplicated names assigned to tensor ",paste(value,col=" "))
  dn <- dimnames(x)
  names(dim(x)) <- value   # gültig
  names(dn) <- value
  dimnames(x) <- dn
  x
}

"dim<-.tensor" <- function(x,value) {
  if( length(value) > 0 )
    NextMethod()
  else if( length(x) > 1 )
    stop("Not fitting dimension")
  else
    NextMethod(value=NULL)
}

"dimnames<-.tensor" <- function(x,value) {
  if( is.null(value) )
    value <- rep(list(NULL),length(dim(x)))
  names(value) <- names(dim(x))
  if( is.null(dim(x)) )
    return(NULL)
  NextMethod(x,value)
}

dimnames.tensor <- function(x) {
  dn <- NextMethod("dimnames",x)
  if( is.null(dn) )
    return( structure(rep(list(NULL),length(dim(x))),names=names(x)) )
  else
    structure(dn,names=names(x))
}


# Euklidische Norm des Tensors
norm <- function(X,...) UseMethod("norm")

norm.tensor <- function(X,i=NULL,...,by=NULL) {
  gsi.debugn("to.tensor$dimX=",dim(X),"i=",i," by=",by)
  if( missing(i) )
    if(is.null(by))
      return(sqrt(sum(abs(X)^2)))
    else
      i <- (1:level.tensor(X))[-toPos.tensor(X,by)]
  i <- toPos.tensor(X,i)
  sqrt(margin.tensor(abs(X)^2,i))
}

opnorm <- function(X,...) UseMethod("opnorm")

opnorm.tensor <- function(X,i=NULL,...,by=NULL) {
  svd.tensor(X,i=i,by=by,name="opnormdim")$d[[opnormdim=1]]
}

margin.tensor <- function(X,i=NULL,by=NULL) {
  i <- toPos.tensor(X,i,by=by)
  A <- gsi.matrify(X,i)
  to.tensor( c(rep(1,nrow(A)) %*% A),dim(X)[-i], dimnames(X)[-i])  
}

diagmul.tensor <- function(X,i=names(D),D,j=i,by=NULL) {
  gsi.debugn("diagmul.tensor$dim=",dim(X)," i=",i," j=",j," by"=by)
  if( !is.tensor(X) ) warning("X must be tensor in diagmul.tensor")
  if( !is.tensor(D) ) warning("D must be tensor in diagmul.tensor")
  if( is.tensor(i) || is.tensor(j) )
    stop("Tensor provided as index in diagmul.tensor")
  if( is.null(by) )  byx<-byy<-c() else {
    byx <- toPos.tensor(X,by,missing.ok=TRUE)
    byy <- toPos.tensor(D,by,missing.ok=TRUE)
    by  <- by[!is.na(byx)&!is.na(byy)]
  }
  ndims <- paste(":!",1:length(c(i,by)),sep="")
  odims <- names(X)[c(toPos.tensor(X,i),toPos.tensor(X,by))]
  names(X)[c(toPos.tensor(X,i),toPos.tensor(X,by))] <- ndims
  names(D)[c(toPos.tensor(D,j),toPos.tensor(D,by))] <- ndims
  erg <- mul.tensor(X,c(),D,c(),by=ndims)
  #print(dim(erg))
  names(erg)[toPos.tensor(erg,ndims)]<-odims
  erg
}

#diagmul.tensor <- function(X,i=names(D),D,j=i,by=NULL) {
#  if( !is.tensor(X) ) X <- as.tensor(X)
#  if( !is.tensor(D) ) D <- as.tensor(D)
#  if( length(by) > 0 )
#    xby <- (!is.na(toPos.tensor(X,by,missing.ok=TRUE)) &
#            !is.na(toPos.tensor(D,by,missing.ok=TRUE)))
#  else
#    xby<-numeric(0)
#  j <- toPos.tensor(D,c(j,by[xby]))
#  i <- toPos.tensor(X,c(i,by[xby]))
#  D <- reorder.tensor(D,j)
#  if( length(i)!=length(dim(D)) || any(dim(X)[i]!=dim(D)[j]) )
#    stop("diagmul.tensor dimension mismatch",length(i)!=length(j),dim(X)[i]!=dim(D)[j])
#  d  <- dim(X)
#  dn <- dimnames(X)
#  gsi.unmatrify(gsi.matrify(X,i)*c(D),d,i,dn)
#}

is.tensor <- function(X) return( "tensor" %in% class(X) )

# Erzeugt alle Tensorpositionen als Matrix
# in Speicherreihenfolge aus einer Dimension
pos.tensor <- function(d) {
  if( is.list(d) )
    d <- sapply(d,length)
  for(i in 1:length(d)){
    if(i>1)
      E <- cbind(E[rep(1:dim(E)[1],d[i]),],
                 rep(1:d[i],rep(dim(E)[1],d[i]))) 
		else
                  E <- cbind(1:d[i])
  }
  colnames(E) <- names(d)
  E
}

# Sortiert die Tensorindices so um, daß 
# die in i genannten in dieser Reihenfolge zuerst kommen

reorder.tensor <- function(x,i=NULL,...,by=NULL)	{
                                        #cat(i)
  d <- dim(x)
  i <- toPos.tensor(x,i,by=by)
  i <- unique(c(i,1:length(d)))
  odimnames <- dimnames(x)
  ndim <- d[i]
  ndx <- c(1)
  weights<- gsi.weights(d)
  for(j in 1:length(d))
    ndx <- rep(ndx,ndim[j])+
      rep((1:ndim[j]-1)*weights[i[j]],
				  rep(length(ndx),ndim[j]))
  to.tensor(x[ndx],ndim,ndimnames=odimnames[i])
}

"^.tensor" <- function(x,y) {
  if( is.character(y) && is.tensor(x) ) {
    if( substr(y[1],1,1) == "$" ) {
      y[1] <- substr(y[1],2,10000)
      if( length( grep("[.]",y) ) > 0  )
        i <- unlist(strsplit(y,"[.]"))
      else 
        i <- unlist(strsplit(y,""))
    } else
    i <- unlist(strsplit(y,"[.]"))
    names(x)[1:length(i)]<-i
    x
  } else if( length(y) == 0 ) {
    x
  } else NextMethod()
}

"$.tensor" <- renamefirst.tensor <- function(x,y) {
  if( is.character(y) && is.tensor(x) ) {
    if( length(grep("[.]",y))>0 ) 
      i <- unlist(strsplit(y,"[.]"))
    else
      i <- unlist(strsplit(y,""))
    names(x)[1:length(i)]<-ifelse(is.covariate(names(x)[1:length(i)]),
                                  i,contraname(i))
    x
  } else if( length(y) == 0 ) {
    x
  } else NextMethod()
}


"|.tensor" <- function(x,y) {
  if( is.character(y) && is.tensor(x) ) {
    if( substr(y[1],1,1) == "$" ) {
      y[1] <- substr(y[1],2,10000)
      if( length( grep("[.]",y) ) > 0  )
        i <- toPos.tensor(x,unlist(strsplit(y,"[.]")))
      else 
        i <- toPos.tensor(x,unlist(strsplit(y,"")))
    } else
    i <- toPos.tensor(x,unlist(strsplit(y,"[.]")))    
    reorder.tensor(x,i)
  } else if( length(y) == 0 ) {
    x
  } else NextMethod()
}



# Gibt die zur Umsortierung im Speicher nötige Indexreihen
# - folge um die in i genannten Indices nach vorn zu holen

reorder.tidx <- function(x,i,...){
  d <- x
  i <- unique(c(i,1:length(d)))
  weights<- gsi.weights(d)
  ndim <- d[i]
  ndx <- c(1)
  for(j in 1:length(d))
    ndx <- rep(ndx,ndim[j])+
      rep((1:ndim[j]-1)*weights[i[j]],
          rep(length(ndx),ndim[j]))
  ndx
}

# Macht aus einer Reihe nach vorne zu holender Indices
# die komplette umsortierung, so daß alle in ihrer neuen 
# Reihenfolge auftreten

gsi.fullreorder <- function(d,first=NULL,last=NULL){
  if( length(last)>0 )
    unique(c(first,rev(unique(rev(c(1:length(d),last))))))
  else
    unique(c(first,1:length(d)))
}

#
# Wie gsi.fullreorder, nur daß die indices i nach hinten
# geschoben werden
#
#gsi.fullreorder.anti <- function(d,i){
#  i <- unique(i)
#  li<-length(i)
#  ld<-length(d)
#  if(li==ld) return(gsi.fullreorder(d,i))
#  unique(c(i,1:ld))[c((li+1):ld,1:li)]
#}

gsi.weightedndx <- function(d,w){
  ndx <- c(1)
  for(j in 1:length(d))
    {
      ndx <- rep(ndx,d[j])+
        rep( ((1:d[j]-1)*w[j]),rep(length(ndx),d[j]) )
    }
  ndx
}

# erzeugt den multiplicator für jeden Arrayindex

gsi.weights   <- function(d){
  cumprod(c(1,d[-length(d)]))
}

# Gibt den Vektor der übrigen indices zurück
gsi.rest <- function(d,i){
  if( length(i)==0 )
    return(1:length(d))
  unique(c(i,1:length(d)))[-(1:length(i))]
}

# invertiert eine Permutation

gsi.invperm <- function(i){
  if( length(i) == 0 )
    return(1)
  i <- unique(c(i,1:max(i)))
  j <- numeric(length(i))
  j[i]<-1:length(i)
  j
}

gsi.matrify <- function(X,i){
  matrix(reorder.tensor(X,i),nrow=prod(dim(X)[i]))
}

gsi.unmatrify <- function(X,d,i,dn=NULL){
  i <- gsi.fullreorder(d,i)
  ii<- gsi.invperm(i)
  dim(X) <- d[i]
  reorder.tensor(X,ii)
}

gsi.lefts <-function(X)	{
  n <- length(dim(X))
  if( n== 2 )
    return(1)
  else
    c(1,3:((n/2)+1))
}

gsi.rights <- function(X) {
  n <-length(dim(X))
  if( n== 2 )
    return(2)
  else
    (1:n)[-c(1,3:(n/2+1))]
}

gsi.without <- function(l,i) if(length(i)>0) l[-i] else l



#
# Collabiert zwei tensoren über die genannten Indices
# i für X und Y für j
#

mul.tensor <- function(X,i=c(),Y,j=i,by=NULL){
  gsi.debugn("mul.tensor$dimX=",dim(X)," i=",i," dimY=",dim(Y)," j=",j," by"=by)
  if( ! is.tensor(X) || ! is.tensor(Y) ) warning("Nontensors multiplied with mul.tensor")
  X <- as.tensor(X)
  Y <- as.tensor(Y)
  if( is.tensor(i) || is.tensor(j) )
    stop("Tensor provided as index in mul.tensor")

  if(length(i)==0 && is.null(by)) { # "Outer Product" 
      if( is.null(dim(X))  ) dim(X) <- length(X)
		if( is.null(dim(Y)) ) dim(Y) <- length(Y)
      E <- rep(c(X),length(Y))*rep(c(Y),rep(length(X),length(Y)))
      return( gsi.debugr(to.tensor(c(E),c(dim(X),dim(Y)),ndimnames=c(dimnames(X),dimnames(Y)))))
    }
  #if( missing(j))
  #  j <- i
  j <- toPos.tensor(Y,j)
  i <- toPos.tensor(X,i)
  if( ! is.null(by) ) {
    # parallization is only done when present in both and not used 
    byi <- toPos.tensor(X,by,missing.ok=TRUE)
    byj <- toPos.tensor(Y,by,missing.ok=TRUE)
    parallel <- !is.na(byi) & !is.na(byj) & !(byi %in% i) & !(byj %in% j)
    byi <- byi[parallel]
    byj <- byj[parallel]
  } else {
    byi <- byj <- c()
  }
  dx <- dim(X)
  dy <- dim(Y)
  if(length(i)!=length(j) || any(dx[i]!=dy[j]) )
    stop("mul.tensor: i incompatible to j")
  if(any(dx[byi]!=dy[byj]) )
    stop("mul.tensor: by incompatible between tensors")
  rix <- gsi.fullreorder(dx,last=c(i,byi))
  riy <- gsi.fullreorder(dy,j,last=byj)
  inner <- prod(dx[i])
  para  <- prod(dy[byj])
  outerx <- prod(gsi.without(dx,c(i,byi)))
  outery <- prod(gsi.without(dy,c(j,byj)))
  dime  <- c(gsi.without(dx,c(i,byi)),gsi.without(dy,c(j,byj)),dx[byi])
  dnX <- dimnames(X)
  dnY <- dimnames(Y)
  ndim  <- c(gsi.without(dnX,c(i,byi)),gsi.without(dnY,c(j,byj)),dnX[byi])
  if( is.null(dnX) )
    dnX <- rep(list(NULL),length(dim(X)))
  if( is.null(dnY) )
    dnY <- rep(list(NULL),length(dim(Y)))
  xtidx <- reorder.tidx(dx,rix)
  ytidx <- reorder.tidx(dy,riy)
  dimx <- c(outerx,inner,para)
  if( is.complex(X) || is.complex(Y) ) {
    E <- .C(tensoraCmulhelper,
            dimx=as.integer(c(outerx,inner,para)),
            as.integer(c(inner,outery,para)),
            as.integer(c(outerx,outery,para)),
            as.complex(X[xtidx]),
            as.complex(Y[ytidx]),
            erg=complex(outerx*outery*para),
            NAOK=TRUE,DUP=FALSE
            )$erg
  } else {
    E <- .C(tensoramulhelper,
            dimx=as.integer(c(outerx,inner,para)),
            as.integer(c(inner,outery,para)),
            as.integer(c(outerx,outery,para)),
            as.numeric(X[xtidx]),
            as.numeric(Y[ytidx]),
            erg=numeric(outerx*outery*para),
            NAOK=TRUE,DUP=FALSE
            )$erg
       #E <- matrix(X[xtidx],ncol=inner)%*%
      #  matrix(Y[ytidx],nrow=inner)
  }
  return( gsi.debugr(to.tensor(c(E),dime,ndimnames=ndim)))
  
} 
##
rep.tensor <- function(x,times,pos=1,name="i",...) {
  gsi.debugn("rep.tensor$dimX=",dim(x)," times=",times," pos=",pos," name=",name)
  if( length(times) > 1 || is.na(name)) {
    pos <- toPos.tensor(x,pos)
    x <- reorder.tensor(x,pos)
    m <- prod(dim(x)[-1])
    ndimnames <- dimnames(x)
    ndimnames[[1]] <- rep(ndimnames[[1]],times)
    tmp <- to.tensor(rep(c(unclass(x)),rep(times,m)),
                     c(gsi.namednumber(name,sum(times)),dim(x)[-1]),
                     ndimnames=dimnames(x))
  } else {
    tmp <- to.tensor(rep(c(unclass(x)),each=times),
                     c(gsi.namednumber(name,times),dim(x)),
                     ndimnames=c(gsi.namedlist(name,NULL),dimnames(x)))
  }
  if( !missing(pos) )
    tmp <- reorder.tensor(tmp,gsi.invperm(pos))
  gsi.debugr(tmp)
}


#
# collabiert Tensor (spur bilden) über die Indices i,j
# (auch listen sind erlaubt,
#  dann werden mehere paare collabiert)

trace.tensor <- function(X,i,j) {
    gsi.debugn("trace.tensor$dimX=",X," i=",i," dimY=",dim(X)," j=",j," by"=by)

  i <- toPos.tensor(X,i)
  j <- toPos.tensor(X,j)
  
  il <- length(i)
  jl <- length(j)
  d  <- dim(X)
  if( il!=jl )
    stop("trace.tensor needs pairs to collaps")
  if( length(i)==0 )
    return(X)
  if( is.null(d) )
    stop("trace.tensor only works on tensors")		
  if( length( unique(c(i,j)) )!=2*il ) 
    stop("trace.tensor only allows different indices")
  if( any( d[i]!=d[j] )) 
    stop("trace.tensor non conformable indices used")
  if( 2*il == 1 )
    return(sum(X))
  odimnames <- dimnames(X)
  weight <- gsi.weights(d)
  rest   <- gsi.rest(d,c(i,j))
  nweight <- c(weight[rest],weight[i]+weight[j])
  dime <- d[rest]
  dimZ <- c(dime,d[i])
  rep(1,prod(d[j]))
  collapsdim<- prod(d[j])
  E <- matrix(
              X[gsi.weightedndx(d[c(rest,i)],nweight)],
              ncol=collapsdim)%*%rep(1,collapsdim)
  gsi.debugr(to.tensor(c(E),dime,ndimnames=odimnames[rest]))  
}



delta.tensor <-function(d,mark="'",dn=NULL,by=NULL) {
  if( is.list(d) ) {
    ndimnames <- d
    d <- sapply(ndimnames,length)
  } else {
    ndimnames <- rep(list(NULL),length(d))
    names(ndimnames) <- names(d)
  }
  if( !is.null(by) ) {
    by <- toPos.tensor(,by,mnames=names(d))
    d2 <- d[by]
    d  <- d[-by]
    ndimnames2 <-ndimnames[by] 
    ndimnames <- ndimnames[-by] 
  } else {
    d2 <- c()
  }
  X <- diag(prod(d))
  if( is.null(by) )
    gsi.debugr(to.tensor(c(X),c(d,mark(d,mark)),ndimnames=c(ndimnames,ndimnames)))
  else
    gsi.debugr(mul.tensor(to.tensor(c(X),c(d,mark(d,mark)),ndimnames=c(ndimnames,ndimnames)),c(),one.tensor(d2,ndimnames2)))
    
}

diag.tensor <- function(X,mark="'",dn=NULL,by=NULL) {
#  if( length(by) > 0 ) {
  by <- toPos.tensor(X,by,missing.ok=TRUE)
  by <- by[!is.na(by)]
  nby <- gsi.without(1:level.tensor(X),toPos.tensor(X,by))
  dx <- dim(X)
  if( is.null(dn) ) dn <- dimnames(X)
  XX <- reorder.tensor(X,c(nby,by))
  pnb <- prod(dx[nby])
    pb  <- prod(dx[by])
  if( is.complex(X) )
    tt <- complex(pnb^2*pb)
  else
    tt <- numeric(pnb^2*pb)
  tt[ rep( (1:pnb)*(1+pnb)-pnb , pb) + pnb^2*(rep((1:pb)-1,each=pnb)) ]<-c(XX)
  gsi.debugr(to.tensor(tt,c(dx[nby],mark(dx[nby],mark=mark),dx[by]),
            c(dn[nby],dn[nby],dn[by])))
  
#  } else gsi.debugr(to.tensor(c(diag(c(X))),c(dim(X),mark(dim(X),mark)),if(is.null(dn)) c(dimnames(X),dimnames(mark(X,mark))) else c(dn,dn) ))
}


tripledelta.tensor <- function(d,mark1="'",mark2="*",dn=NULL) {
  p <- prod(d)
  tt <- rep(0,p^3)
  tt[1+(1+p+p^2)*(0:(p-1))]<-1
  gsi.debugr(to.tensor(tt,c(d,mark(d,mark1),mark(d,mark2)),if(is.null(dn)) NULL else c(dn,dn,dn)))
  }	

one.tensor <- function(d=NULL,dn=NULL) {
  if( length( d )==0 )
    return(1)
  gsi.debugr(to.tensor(rep(1,prod(d)),d,dn))
}

mark <- function(X,mark,...) UseMethod("mark")

mark.tensor <- function(X,mark="'",i=1:level.tensor(X),...,by=NULL) {
  i <- toPos.tensor(X,i)
  if( ! is.null(by) ) {
    not <- toPos.tensor(X,by)
    i <- i[! (i %in% by[!is.na(by)])]
  }
  nam <- names(X)
  nam[i] <- paste(nam[i],mark,sep="")
  names(X) <- nam
  gsi.debugr(X)
}

mark.numeric <- function(X,mark="'",i=1:length(X),...,by=NULL) {
  nam <- names(X)
  if( is.character(i) ) {
    oi<-i
    i <- charmatch(i,nam)
    if( any(is.na(i)) || any(i==0) )
      stop("mark.numeric: No mach found for ", oi[is.na(oi)|i==0])
  }
  if(! is.null(by) ) {
    if( is.character(by) )
      by <- charmatch(by,nam)
    i <- i[!(i %in% by[!is.na(by)])]
  }
  nam[i] <- paste(nam[i],mark,sep="")
  names(X) <- nam
  X
}

mark.character <- function(X,mark="'",i=1:length(X),...,by=NULL) {
  nam <- X
  if( is.character(i) ) {
    oi<-i
    i <- charmatch(i,nam)
    if( any(is.na(i)) || any(i==0) )
      stop("mark.character: No mach found for ", oi[is.na(oi)|i==0])
  }
  if( !is.null(by) ) {
    if( is.character(by) )
      by <- charmatch(by,nam)
    i <- i[!(i %in% by[!is.na(by)])]
  }
  nam[i] <- paste(nam[i],mark,sep="")
  nam
}


inv.tensor <- function(X,i,...,allowSingular=FALSE,eps=1E-10,by=NULL) {
  gsi.debug("inv.tensor$dimX=",dim(X)," i=",i)
  if( ! is.tensor(X) )
    warning("Not a tensor given")
  if( is.tensor(i) )
    stop("Tensor given as index")
  X <- as.tensor(X)
  i  <- toPos.tensor(X,i)
  by <- toPos.tensor(X,by,missing.ok=TRUE)
  by <- by[!is.na(by)]
  j <- gsi.without(1:length(dim(X)),c(i,by))
  dx <-dim(X)
  dn <-dimnames(X)
  XX <- reorder.tensor(X,c(i,j,by))
  XX <- c(unclass(XX))
  dim(XX) <- c(prod(dx[i]),prod(dx[j]),prod(dx[by]))
  if( allowSingular ) {
    XX <- unlist(lapply(as.list(1:dim(XX)[3]),
                        function(i) {
                          udv <- svd(XX[,,i])
                          if( abs(udv$d[1])^2 == 0 )
                            rank <- 0
                          else
                            rank <- sum( abs(udv$d[1])*eps < abs(udv$d) )
                          # cat("The rank is ",rank,"\n")
                          if( rank == 0 )
                            return(matrix(0,nrow=dim(XX)[1],ncol=dim(XX)[2]))
                            Conj(udv$u[,1:rank,drop=FALSE] %*%
                              ( 1/udv$d[1:rank] * t(Conj(udv$v[,1:rank,drop=FALSE])) ))
                        }
                        ))
  } else {
    XX <- unlist(lapply(as.list(1:dim(XX)[3]),
                        function(i) {
                          t(solve(XX[,,i]))
                        }
                        ))
  }
  gsi.debugr(to.tensor(c(XX),dx[c(i,j,by)],dn[c(i,j,by)]))
}



solve.tensor <- function(a,b,i,j=i,...,allowSingular=FALSE,eps=1E-10,by=NULL) {
  gsi.debug("mul.tensor$dima=",dim(a)," i=",i," dimb=",dim(b)," j=",j," by"=by)
  if( !is.tensor(a) || ! is.tensor(b) )
    warning("Nontensors provided as tensors in solve.tensor")
  if( is.tensor(i) || is.tensor(j) )
    stop("Tensor provided as index in solve.tensor")
  X <- as.tensor(a)
  b <- as.tensor(b)
  i  <- toPos.tensor(X,i)
  j  <- toPos.tensor(b,j)
  byi <-toPos.tensor(X,by,missing.ok=TRUE)
  byj<- toPos.tensor(b,by,missing.ok=TRUE)
  byNob <- by[!is.na(byi) & is.na(byj)]
  byBoth<- by[!is.na(byi)&!is.na(byj)]
  byi <- toPos.tensor(X,byBoth)
  byj <- toPos.tensor(b,byBoth)
  byi2<- toPos.tensor(X,byNob)
  k <- gsi.without(1:length(dim(X)),c(i,byi,byi2))
  l <- gsi.without(1:length(dim(b)),c(j,byj))
  dx <-dim(X)
  db <-dim(b)
  dnx <-dimnames(X)
  dnb <-dimnames(b)

  XX <- reorder.tensor(X,c(i,k,byi,byi2))
  bb <- reorder.tensor(b,c(j,l,byj))
  XX <- c(unclass(XX))
  bb <- c(unclass(bb))
  if( length(i)!=length(j) || any(dx[i]!=db[j] ) )
    stop("Not fitting dimensions in solve.tensor da=",paste(dx[i]),
         " db=",paste(db[j]))
  if( ! allowSingular && ( length(k)!=length(i) || prod(dx[i])!=prod(dx[k])) )
    stop("Not square tensor in solve.tensor da=",paste(dx[i]),
         " db=",paste(dx[k]))
    
  dim(XX) <- c(prod(dx[i]),prod(dx[k]),prod(dx[byi]),prod(dx[byi2]))
  dim(bb) <- c(prod(db[j]),prod(db[l]),prod(db[byj]))
  if( allowSingular ) {
    EE <- unlist(lapply(as.list(1:dim(XX)[3]),
                        function(i) {
                          unlist(lapply(as.list(1:dim(XX)[4]),function(j) {
                          udv <- svd(XX[,,i,j])
                          if( abs(udv$d[1])^2 == 0 )
                            rank <- 0
                          else
                            rank <- sum( abs(udv$d[1]*eps) < abs(udv$d) )
                          #cat("Eps is ",eps,"\n")
                          #cat("The rank is ",rank,"\n")
                          # ud t(C(v)) x = b
                          # x = v 1/d t(C(u))b
                          udv$v[,1:rank] %*%
                            ( 1/udv$d[1:rank] * (t(Conj(udv$u[,1:rank])) %*%
                                                   bb[,,i]
                                                   )) 
                        }))}
                        ))
  } else {
    EE <- unlist(lapply(as.list(1:dim(XX)[3]),
                        function(i) {
                          unlist(lapply(as.list(1:dim(XX)[4]),function(j) {
                          solve(XX[,,i,j],bb[,,i])
                        }))}
                        ))
  }
  d <- c(dx[c(k)],db[c(l,byj)],dx[byi2])
  pn <- mapply(function(x,y) if(is.null(x)) y else x,dnb[byj],dnx[byi])
  dn <- c(dnx[k],dnb[l],pn,dnx[byi2])
  gsi.debugr(to.tensor(c(EE),d,dn))
}
  

#solve.tensor <- function(a,b,i,j=i,...,by=c()){
#  if( length(by) > 0 ) {
#    i  <- toPos.tensor(X,i)
#    j  <- toPos.tensor(b,j)
#    byi <- toPos.tensor(X,by,missing.ok=TRUE)
#    byj <- toPos.tensor(X,by,missing.ok=TRUE)
#    
#    
#  }
#  X <- to.tensor(a)
#  b <- to.tensor(b)
#  i <- toPos.tensor(X,i)
#  j <- toPos.tensor(X,j)
#  
#  dx <- dim(X)
#  db <- dim(b)
#  ri <- gsi.rest(dx,i)
#  rj <- gsi.rest(db,j)
#  X<-gsi.matrify(X,i)
#  b<-gsi.matrify(b,j)
#  x<-solve(X,b)
#  x <- to.tensor(c(x),c(dx[ri],db[rj]),c(dimnames(X)[ri],dimnames(b)[rj]))
#  x	
#}

chol.tensor <- function(X,i,j,...,name="lambda") {
  gsi.debug("chol.tensor$dimX=",dim(X)," i=",i," j=",j)
  if( !is.tensor(X) )
    warning("Nontensors provided as tensors in chol.tensor")
  if( is.tensor(i) || is.tensor(j) )
    stop("Tensor provided as index in chol.tensor")

  X <- as.tensor(X)
  i  <- toPos.tensor(X,i)
  j  <- toPos.tensor(X,j)
  by <- gsi.without((1:length(dim(X))),c(i,j))
  if( length(i)!=length(j) || any(dim(X)[i]!=dim(X)[j]) )
    stop("No symmetry in chol")
  dx <-dim(X)
  dn <-dimnames(X)
  XX <- reorder.tensor(X,c(i,j,by))
  XX <- c(unclass(XX))
  dim(XX) <- c(prod(dx[i]),prod(dx[j]),prod(dx[by]))
  ld <- gsi.namednumber(name,prod(dx[j]))
  XX <- unlist(lapply(as.list(1:dim(XX)[3]),
                      function(i,...) {
                        chol(XX[,,i],...)
                      }
                      ,...))
  return(gsi.debugr(to.tensor(XX,c(ld,dx[c(i,by)]),c(list(NULL),dn[c(i,by)]))))
}


#chol.tensor <- function(X,i,j,...,name="lambda") {
#  X <- as.tensor(X)
#  i <- toPos.tensor(X,i)
#  j <- toPos.tensor(X,j)
#  d <- dim(X)
#  dn <- dimnames(X)
#  if( !is.null(dn) )
#    dn <- c(gsi.namedlist(name,NULL),dn[i])
#  if( length(i)!=length(j) ||
#     length(unique(c(i,j)))!=length(d) )
#    stop("chol.tensor: i or j wrong")
#  if( any(d[i]!=d[j]) )
#    stop("chol.tensor: tensor not square")
#  ii <- gsi.invperm(c(i,j))
#  X <- reorder.tensor(X,c(i,j))
#  X <- gsi.matrify(X,1:length(i))
#  print(svd(X))
#  X <- chol(X)
#  X <- to.tensor(c(X),c(dim(X)[1],d[i]),dn)
#  X
#}

level.tensor <- function(X,...) {
  if( is.null(dim(X)) )
    if( length(X) != 1 )
      return(1)
    else
      return(0)
  return(length(dim(X)))
}


svd.tensor <- function(X,i,j=NULL,...,name="lambda",by=NULL) {
  gsi.debug("svd.tensor$dimX=",dim(X)," i=",i," j=",j," by"=by)
  if( !is.tensor(X) )
    warning("Nontensors provided as tensors in svd.tensor")
  if( is.tensor(i) || is.tensor(j) )
    stop("Tensor provided as index in svd.tensor")

  X <- as.tensor(X)
  i <- toPos.tensor(X,i)
  by <- toPos.tensor(X,by,missing.ok=TRUE)
  by <- by[!is.na(by)]
  if( is.null(j) )
    j <- gsi.without((1:length(dim(X))),c(i,by))
  else
    j <- toPos.tensor(X,j)
  if( missing(i) )
    i <- gsi.without((1:length(dim(X))),c(j,by))
  by <- gsi.without((1:length(dim(X))),c(i,j))
  dx <- dim(X)
  dn<-  dimnames(X)
  di <- dim(X)[i]
  din<- dimnames(X)[i]
  dj <- dim(X)[j]
  djn<- dimnames(X)[j]
  XX <- unclass(reorder.tensor(X,c(i,j,by)))
  dim(XX) <- c(prod(dx[i]),prod(dx[j]),prod(dx[by]))
  EE <- lapply(1:(dim(XX)[3]),function(i) {
    svd(XX[,,i])
  },...)
  ld <- gsi.namednumber(name,min(prod(dx[i]),prod(dx[j])))
  u <- to.tensor(unlist(lapply(EE,function(x) x$u)),
                 c(dx[i],ld,dx[by]),
                 c(dn[i],list(NULL),dn[by]))
  d <- to.tensor(unlist(lapply(EE,function(x) x$d)),
                 c(ld,dx[by]),
                 c(list(NULL),dn[by]))
  v <- to.tensor(unlist(lapply(EE,function(x) x$v)),
                 c(dx[j],ld,dx[by]),
                 c(dn[j],list(NULL),dn[by]))
  list(u=gsi.debugr(u),
       d=gsi.debugr(d),
       v=gsi.debugr(v))
}

power.tensor <- function(X,i,j,p=0.5,by=NULL)	{
  gsi.debug("power.tensor$dimX=",dim(X)," i=",i," j=",j," by"=by)
  if( !is.tensor(X) )
    warning("Nontensors provided as tensors in power.tensor")
  if( is.tensor(i) || is.tensor(j) )
    stop("Tensor provided as index in power.tensor")

  X <- as.tensor(X)
  udv <- svd.tensor(X,i,j,by=by,name="<<<")
  gsi.debugr(einstein.tensor(udv$u,diag=udv$d^p,udv$v,only="<<<",by=by))
}



#power.tensor <- function(X,i,j,p=0.5)	{
#  X <- as.tensor(X)
#  i <- toPos.tensor(X,i)
#  j <- toPos.tensor(X,j)
#  d <- dim(X)
#	if( length(i)!=length(j) ||
#           length(unique(c(i,j)))!=length(d) )
#		stop("root.tensor: i or j wrong")
#  if( any(d[i]!=d[j]) )
#    stop("root.tensor: tensor not square")
#  ii <- gsi.invperm(c(i,j))
 # udv <- svd.tensor(X,i,j)
 # k <- level.tensor(udv$u)
 # X <- mul.tensor(mul.tensor(diag(udv$d^p),1,udv$u,k),1,udv$v,k) 
#  reorder.tensor(X,ii)
#}

to.matrix.tensor <- function(X,i,j,by=NULL) {
  if( !is.tensor(X) )
    warning("Nontensors provided as tensors in to.matrix.tensor")
  if( (!missing(i) && is.tensor(i)) || (!missing(j) && is.tensor(j)) )
      stop("Tensor provided as index in to.matrix.tensor")

  X <- as.tensor(X)
  d <- dim(X)
  dn <- dimnames(X)
  by = toPos.tensor(X,by,missing.ok=TRUE)
  by <- by[!is.na(by)]
  if( missing(i) ) {
    j <- toPos.tensor(X,j)
    i <- gsi.without((1:length(dim(X))),c(j,by))
  }
  i <- toPos.tensor(X,i)
  if( missing(j) )
    j <- gsi.without((1:length(dim(X))),c(i,by))
  j <- toPos.tensor(X,j)
  structure(unclass(c(reorder.tensor(X,c(i,j,by)))),
            dim=c(i=prod(d[i]),j=prod(d[j]),d[by]),
            dimnames=c(list(NULL),list(NULL),dn[by])
            )
}

gsi.checkduplicate <- function(x) {
  value <- names(x)
  if( length(value) > length(unique(value)) ) {
    warning("Duplicated names in tensor",paste(value,col=" "))
    TRUE
  } else
  FALSE
}

gsi.stdnames <- function(k,prefix="I",avoid=NULL) {
  if( k==0 )
    character(0)
  else {
    tmp <- paste(prefix,1:k,sep="")
    if( !is.null(avoid) ) {
      while( any( a <- tmp %in% avoid ) ) {
        na <- sum(a)
        tmp <- c(tmp[!a],paste(prefix,(k+1):(k+na),sep=""))
        k <- k+na 
      }
    }
    tmp
  }
}



gsi.namedlist <- function(nam,...) {
  tmp <- list(...)
  names(tmp) <- nam
  tmp
}

gsi.namednumber <- function(nam,...) {
  tmp <- c(...)
  names(tmp) <- nam
  tmp
}



untensor <- function(X,i=NULL,name=NULL,pos=1,by=NULL){
  if( is.list(i) ) {
    if( !is.null(name))
      names(i) <- name
    i <- lapply(i,function(i) names(X)[toPos.tensor(X,i)])
    for(a in 1:length(i)) 
      X <- untensor(X,i=i[[a]],names(i)[[a]],pos=1)
    return(X)
  }
  i <- toPos.tensor(X,i,by=by)
  d <- dim(X)
  r <- gsi.rest(d,i)
  X<- reorder.tensor(X,i)
  odimnames <- dimnames(X)
  k <- prod(d[i])
  if( is.null(name))
    names(k) <- gsi.stdnames(1,"I",avoid=names(d[-i]))
  else
    names(k) <- name
  dimnames(X) <- NULL
  dim(X) <- c(k,d[r])
  dimnames(X) <- c(structure(list(gsi.untensornames(odimnames[i],d[i])),names=name),odimnames[-(1:length(i))])
  if( pos!=1 ) X <- reorder.tensor(X,gsi.invperm(pos))
  gsi.debugr(X)
} 


gsi.untensornames <- function(X,d=sapply(X,length)) {
  if( is.null(X) )
    return(NULL)
  if( length(X) == 0 )
    return(NULL)
  wrong <- sapply(X,is.null)
  if( all(wrong) )
    return(NULL)
  if( any(wrong ) )
    X[wrong]<- lapply(d[wrong],seq)
  s<-""
  for(i in 1:length(X) ) {
    s <- outer(s,X[[i]],paste,sep="")
  }
  return(c(s))
  
}

as.tensor <- function(X,...) UseMethod("as.tensor")
as.tensor.default <- function(X,...,dims=NULL) {
#	if( !is.numeric(X))
#		stop("as.tensor.default: Only numeric tensors\
#			   supported")
  if( is.null(dims) )
    gsi.debugr(to.tensor(X))
  else
    gsi.debugr(to.tensor(c(X),dims))
} 

as.tensor.tensor <- function(X,...) {X}

#renorm.rows <- function(X) {
#  X/c(sqrt(X^2%*%rep(1,dim(X)[2])))
#}

#renorm.tensor<-function(X,i){
#  if( missing(i) ) 
#    return(to.tensor(c(X/norm.tensor(X)),dim(X),dimnames(X)))
#  i <- toPos.tensor(X,i)
 # d <- dim(X)
 # inorm <- 1/sqrt(mul.tensor(X^2,i,one.tensor(d[i])))
 # diagmul.tensor(X,gsi.rest(dim(X),i),inorm,1:level.tensor(inorm))
#}

slice.tensor <-function(X,i,what,drop=FALSE) {
  i <- toPos.tensor(X,i)
  if( length( i ) > 1 ) {
    stop("multiple slice not supported yet")
  }
  d <- dim(X)
  dn <- dimnames(X)
  dimnames(X) <- NULL
  if( is.character(what) ) {
    if(!is.null(dn) && !is.null(dn[[i]]))
      what <- match(what,dn[[i]])
    else
      stop("slice.tensor: Missing names with named subscript")
  }
  if( any(what>d[i]))
    stop( "slice.tensor: subscript out of bound: ",what)
  if( any(what<1) ) {
    if( all(what<0) )
      what <- 1:(d[i])[what]
    else
      stop( "slice.tensor: mixed positiv and negativ subscribts" )
  }
  w <- gsi.weights(d)
  nd <- d
  nd[i]<-length(what)
  ndx <- c(1)
  for(j in 1:length(nd))
    if( j!=i )
      {
        ndx <- rep(ndx,nd[j])+
          rep( ((1:nd[j]-1)*w[j]),rep(length(ndx),nd[j]) )
      }
    else
      {
        ndx <- rep(ndx,nd[j])+
          rep( ((what[1:nd[j]]-1)*w[j]),
              rep(length(ndx),nd[j]) )
      } 
  X <- X[ndx]
  if( !is.null(dn) && !is.null(dn[[i]]) )
    dn[[i]] <- dn[[i]][what]
  if( drop && length(what) == 1 ) {
    nd <- nd[-i]
    if( ! is.null(dn) )
      dn <- dn[-i]
  }
  dim(X)<-nd
  names(dn) <- names(nd)
  dimnames(X) <- dn
  gsi.debugr(X)
}

"[[.tensor" <- function(X,...,drop=TRUE) {
  namedargs <- list(...)
  if( is.null(names(namedargs)) )
    return(NextMethod("[[",X))
  for(n in names(namedargs)) {
    if( is.call(namedargs[[n]])) {
      k <- match(n,names(X))
      if( is.na(k) )
        stop("noexisting dimension in [[.tensor")
      names(X)[k]<- as.character(namedargs[[n]][[2]])
    } else
      X <- slice.tensor(X,n,namedargs[[n]],drop=drop)
  }
  gsi.debugr(X)
}

"[.tensor" <- function(X,...,drop=TRUE) {
  # as.tensor.default(NextMethod("[",X))
  dimnames(X) <- dimnames(X)
  r <- NextMethod("[",X)
  names(dim(r))<-names(dimnames(r))
  as.tensor.default(r)
}


"[[<-.tensor" <- function(X,...,value) {
  U <- to.tensor(1:length(X),dim(X),dimnames(X))
  U <- U[[...]]
  if( !is.null(dim(value)) && (length(dim(value))!=length(dim(U)) || any(dim(U)!=dim(value))))
      warning("non matching arrays in [<-.tensor",paste(dim(U),col=","),":",paste(dim(value),col=","))
  X[c(U)]<- c(value)
  gsi.debugr(X)
}


undrop.tensor <- function(A,name,pos=1) {
  dn <- c(dimnames(A),"!intern"=list(NULL))
  dm <- c(dim(A),"!intern"=1)
  names(dn)[length(names(dn))]<-name
  names(dm)[length(names(dm))]<-name
  attr(A,"dimnames") <- NULL
  dim(A) <- dm
  dimnames(A) <- dn
  if( !missing(pos) ) {
    A <- reorder.tensor(A,2:pos)
  }
  A
}

#combineCF.tensor <- function(A,i,B,j) {
#  i <- toPos.tensor(A,i)
#  j <- toPos.tensor(B,j)
#  dA <- dim(A)
#  dB <- dim(B)
#  if( length(i) !=length(j) || any(dA[i]!=dB[j] ))
#    stop("combineCF: Dimensions don't match", paste(da[i]),":",paste(dB[j]))
#  A <- gsi.matrify(A,i)
#  B <- gsi.matrify(B,j)
#  nB <-dim(B)[2]
#  O <- matrix(0,nrow=nB,ncol=nB)
#  C <- rbind(cbind(A,B),cbind(t(B),O))
#  C
#}


bind.tensor <- function(A,dA=NULL,B,dB=dA) {
  if( is.null(A) )
    return(B)
  if( is.null(B) )
    return(A)
  if( is.null(dA) ) {
    A <- undrop.tensor(A,"i")
    dA <- length(dim(A))
  }
  dA <- toPos.tensor(A,dA)
  if( is.null(dB) ) {
    B <- undrop.tensor(B,"i")
    dB <- length(dim(B))
  }
  dB <- toPos.tensor(B,dB)
  matching <- match(names(A)[-dA],names(B)[-dB])
  if( !all(is.na(matching))) {
    if(any(is.na(matching))) {
      stop("wrong match in bind.tensor")
    }
    B  <- reorder.tensor(B,c(dB,matching))
    dB <- 1
  }
  if( length(dB) != 1
  || 1!=length(dA) || length(dim(B))!=length(dim(A)) ||
     any(dim(B)[-dB]!=dim(A)[-dA])
     ) {
    stop("bind.tensor dimensions don't match")
  }
  dAo <- dA
  dBo <- dB
  A <- reorder.tensor(A,(1:length(dim(A)))[-dA])
  B <- reorder.tensor(B,(1:length(dim(B)))[-dB])
  dA <- length(dim(A))
  dB <- length(dim(A))
  dn <- dimnames(A)
  dnB <- dimnames(dB)[[dB]]
  if(! is.null(dn[[dA]]) || ! is.null(dnB)) {
    if( is.null(dnB))
      dnB <- gsi.stdnames(dim(B)[dB],"B")
    if( is.null(dn[[dA]]))
      dn[[dA]] <- gsi.stdnames(dim(A)[dA],"A")
    dn[[dA]] <- c(dn[[dA]],dnB)
  }
  dm <- dim(A)
  dm[dA] <- dm[dA]+dim(B)[dB]
  erg <- c(A,B)
  dim(erg) <- dm
  dimnames(erg) <- dn
  reorder.tensor(erg,c(gsi.vonbis(1,dAo-1),length(dA)))
}

gsi.vonbis <- function(a,b) {
  if( b>=a )
    a:b
  else
    c()
}

                       
toPos.tensor <- function(M,l=NULL,mnames=names(dim(M)),by=NULL,...,both=FALSE,missing.ok=FALSE) {
  if( is.null(l) && ! is.null(by) ) 
    return( (1:length(mnames))[-toPos.tensor(M,by,mnames,both=both)])
  if( length(l) == 0 )
    return(numeric(0))
  if( is.name(l) )
    l <- as.character(l)
  if( is.character(l) ) {
    if( both )
      e <- charmatch(as.covariate(l),as.covariate(mnames))
    else
      e <- charmatch(l,mnames)
    if( any(is.na(e)) & !missing.ok ) 
      stop("No match found for ",paste(l[is.na(e)],col=","))
    if( any(e[!is.na(e)]==0) )
      stop("Match not unique for ", paste(l[e==0],col=",") )
    return(e);
  } else if( is.numeric(l) ) {
    return(l)
  } else {
    n <- sapply(l,is.character)
    if( any( n ))
      l[n]<-toPos.tensor(M,l=as.character(l[n]),names=mnames)
    return(l);
  }
}

einstein.tensor <- function(...,only=NULL,by=NULL) {
  ts <- list(...)
  if( length(ts) == 0 )
    return(NULL)
  if( length(ts) < 2 )
    return(ts[[1]])
  tmp <- NULL
  nams <- names(ts)
  if( is.null(nams) ) {
    nam <- rep("",length(ts))
  }
    
  for(k in 1:length(ts)) {
    ten <- ts[[k]]
    nam <- nams[k]
    if( !is.null(nam) && nam=="diag" ) {
      if( !is.tensor(ten) )
        warning("diagmul with nontensor in einstein.tensor")
      tmp <- diagmul.tensor(tmp,names(ten),ten,names(ten))
    } else if( is.character(ten) ) {
      olds <- match(nam,names(tmp))
      if( is.na(olds) )
        stop("Unknown dimension ",ten," in einstein.tensor")
      if( ten %in% names(tmp)  ) {
        tmp <- trace.tensor(tmp,nam,ten)
      } else {
        names(tmp)[olds]<-ten
      }
    } else if( is.null(tmp) ) {
      tmp <- ten
    } else if( is.null(ten)) {
      tmp <- tmp
    } else if( length(ten) == 1 && is.null(dim(ten))) {
      tmp <- tmp*ten
    } else if( length(tmp) == 1 && is.null(dim(tmp))) {
      tmp <- tmp*ten
    } else {
      n1 <- names(tmp)
      n2 <- names(ten)
      if( is.null(n1) || is.null(n2)) {
        tmp <- mul.tensor(tmp,c(),ten,c())
      } else {
        jm <- match(n1,n2)
        i <- n1[!is.na(jm)]
        if( !missing(only) )
          i <- i[i %in% only]
        if( !is.null(by) )
          i <- i[! (i %in% by)]
        tmp<- mul.tensor(tmp,i,ten,i,by=by)
      }
    }
  }
  gsi.debugr(tmp)
}

"%e%" <- function(x,y) UseMethod("%e%")
"%e%.tensor" <- function(x,y) einstein.tensor(x,y)
"%r%" <- function(x,y) UseMethod("%r%")
"%r%.tensor" <- function(x,y) riemann.tensor(x,y)

#"%+%" <- function(x,y) UseMethod("%+%")
"+.tensor" <- function(x,y) add.tensor(x,y)
#"-" <- function(x,y) UseMethod("%-%")
"-.tensor" <- function(x,y) {
  if( missing(y) ) {
    oc <- class(x)
    structure(-unclass(x),class=oc)
  } else add.tensor(x,y,"-")
}

"*.tensor" <- function(x,y) add.tensor(x,y,"*")
"/.tensor" <- function(x,y) add.tensor(x,y,"/")


add.tensor <- function(X,Y,op="+",only=NULL) {
  if( !is.tensor(X)  )
    if( length(X) == 1 )
      return(gsi.debugr(to.tensor(c(do.call(op,list(unclass(X),unclass(Y)))),
                                  dim(Y),dimnames(Y))))
    else stop("Tensor operation with nontensor")
  if( !is.tensor(Y))
    if( length(Y) == 1)
    return(gsi.debugr(to.tensor(c(do.call(op,list(unclass(X),unclass(Y)))),
                                dim(X),dimnames(X))))
    else stop("Tensor operation with nontensors")
  wk <- names(X) %in% names(Y)
  if( !is.null(only) )
    wk <- wk & (names(X) %in% only)
  nams <- names(X)[wk]
  #if( length(nams) < 1 )
  #  warning("No match in add.tensor",names(X),":",names(Y),":",only)
  i <- toPos.tensor(X,nams)
  j <- toPos.tensor(Y,nams)
  if( level.tensor(Y)>length(nams)  ) {
    Xb <- mul.tensor(X,c(),one.tensor(gsi.without(dim(Y),j),gsi.without(dimnames(Y),j)),c())
#    dimnames(Xb)[(length(dim(X))+1) : length(dim(Xb))]<-dimnames(Y)[-j] 
  } else
    Xb <- X
  if( level.tensor(X)>length(nams)  ) 
    Yb <- mul.tensor(Y,c(),one.tensor(gsi.without(dim(X),i),
                                      gsi.without(dimnames(X),i)),c())
  #to.tensor( unclass(one.tensor(dim(X)[-i])) %o% unclass(Y) )
  else
    Yb <- Y
  Yb <- reorder.tensor(Yb,names(Xb))
  gsi.debugr(to.tensor(c(do.call(op,list(unclass(Xb),unclass(Yb)))),dim(Xb),dimnames(Xb)))
}




#secondnames <- function(n,tag="'") {
#  paste(n,tag,sep="")
#}

#firstnames  <- function(n,tag) {
#  substr(n,1,nchar(n)-nchar(tag))
#}

contraname   <- function(x) ifelse(substr(x,1,1)=="^",substr(x,2,100000),paste("^",x,sep=""))

is.covariate     <- function(x,...) UseMethod("is.covariate")
is.covariate.tensor <- function(x,...) {is.covariate(names(dim(x)))}
is.covariate.numeric <- function(x,...) {is.covariate(names(x))}
is.covariate.character <- function(x,...) {substr(x,1,1)!="^"}
as.covariate <- function(x,...) UseMethod("as.covariate")
as.covariate.character <- function(x,...) ifelse(substr(x,1,1)=="^",substr(x,2,100000),x)

is.contravariate <- function(x,...) UseMethod("is.contravariate")
is.contravariate.tensor <- function(x,...) {is.contravariate(names(x))}
is.contravariate.numeric <- function(x,...) {is.contravariate(names(x))}
is.contravariate.character <- function(x,...) {substr(x,1,1)=="^"}
as.contravariate <- function(x,...) UseMethod("as.contravariate")
as.contravariate.character<- function(x,...) ifelse(substr(x,1,1)=="^",x,paste("^",x,sep=""))

drag.tensor <- function(x,g,d) {
  gsi.debug("drag.tensor$dimx=",dim(x)," d=",d," dimg=",dim(g))
  cg <- is.covariate(g)
  if( (any(cg ) && ! all(cg)) || level.tensor(g)!= 2 )
    stop("g must be either covariate or contravariate")
  gcov <- gcon <- g
  if( all(cg) )
    gcon <- inv.tensor(g,1)
  else
    gcov <- inv.tensor(g,1)
  for(i in d) {
    na <- names(x)[toPos.tensor(x,i,both=TRUE)]
    xko <- is.covariate(na)
    if( xko ) {
      names(gcon) <- c(na,contraname(na))
      x <- mul.tensor(x,na,gcon,na)
    } else {
      names(gcov) <- c(na,contraname(na))
      x <- mul.tensor(x,na,gcov,na)
    }
  }
  gsi.debugr(x)
}




riemann.tensor <- function(...,only=NULL,by=NULL) {
  ts <- list(...)
  if( length(ts) == 0 )
    return(NULL)
  if( length(ts) < 2 )
    return(ts[[1]])
  tmp <- NULL
  nams <- names(ts)
  if( is.null(nams) ) {
    nam <- rep("",length(ts))
  }
    
  for(k in 1:length(ts)) {
    ten <- ts[[k]]
    nam <- nams[k]
    if( !is.null(nam) && nam=="diag" ) {
      tmp <- diagmul.tensor(tmp,contraname(names(ten)),ten,names(ten))
    } else if( identical(is.character(ten),TRUE) ) {
      olds <- match(nam,names(tmp))
      news <- match(contraname(ten),names(tmp))
      if( is.na(olds) )
        if( is.na(news) )
          stop("Unknown dimension ",nam,"and",contraname(ten),
               " in riemann.tensor")
        else  
          names(tmp)[news]<-contraname(nam)
      else if( is.na(news) )
        names(tmp)[olds]<-ten
      else
        tmp <- trace.tensor(tmp,nam,contraname(ten))
    } else if( is.null(tmp) ) {
      tmp <- ten
    } else if( is.null(ten)) {
      tmp <- tmp
    } else if( length(ten) == 1 && is.null(dim(ten))) {
      tmp <- tmp*ten
    } else if( length(tmp) == 1 && is.null(dim(tmp))) {
      tmp <- tmp*ten
    } else {
      n1 <- names(tmp)
      n2 <- names(ten)
      if( is.null(n1) || is.null(n2)) {
        tmp <- mul.tensor(tmp,c(),ten,c())
      } else {
        jm <- match(contraname(n1),n2)
        i <- n1[!is.na(jm)]
        if( !missing(only) )
          i <- i[i %in% only]
        if( !is.null(by) )
          i<- i[! (i%in% by) & ! (contraname(i) %in% by)]
        tmp<- mul.tensor(tmp,i,ten,contraname(i))
      }
    }
  }
  gsi.debugr(tmp)
}



mean.tensor <- function(x,along,...,na.rm=FALSE) {
  if( !na.rm && !all(is.finite(x)) ) 
    stop("Missings in mean.tensor")
  N <- to.tensor(as.numeric(is.finite(x)),dim(x))
  along <- toPos.tensor(x,along)
  one <-one.tensor(dim(x)[along])
  x[!is.finite(x)]<- 0
  S <- mul.tensor( x  ,along , one , 1:level.tensor(one) ) 
  n <- mul.tensor( N  ,along , one , 1:level.tensor(one) )
  S/n # to.tensor(c(S)/c(n),dim(S),dimnames(S))
}

var.tensor <- function(x,y=NULL,...,along,by=NULL,na.rm=FALSE,mark="'") {
  if( is.null(y) ) {
    if( !na.rm && !all(is.finite(x)) ) 
      stop("Missings in bar.tensor")
    Nx <- to.tensor(as.numeric(is.finite(x)),dim(x))
    by <- toPos.tensor(x,by)
    by <- by[!is.na(by)]
    along <- toPos.tensor(x,along)
    one <-one.tensor(dim(x)[along])
    x[!is.finite(x)]<- 0
    S <- mul.tensor( x  ,along , one , 1:level.tensor(one) ) 
    n <- mul.tensor( Nx ,along , one , 1:level.tensor(one) )
    #M <-one.tensor(dim(x)[along])
    #S <- mul.tensor( x,along , N , 1:level.tensor(M) )
    #d <- gsi.without(1:level.tensor(X),c(along,by))
    x <- x - S/n # to.tensor(c(S)/c(n),dim(S),dimnames(S))
    y <- mark(x,mark,by=c(along,by))
    Ny<- mark(Nx,mark,by=c(along,by))
    S2 <- mul.tensor( x,along, y,along,by=by)
    N2 <- mul.tensor(Nx,along,Ny,along,by=by)
    S2 / (N2-1) #to.tensor(c(S2)/(c(N2)-1),dim(S2),dimnames(S2))
  } else {
    if( !na.rm && !all(is.finite(x)) ) 
      stop("Missings in bar.tensor")
    Nx <- to.tensor(as.numeric(is.finite(x)),dim(x))
    Ny <- to.tensor(as.numeric(is.finite(y)),dim(y))
    byx <- toPos.tensor(x,by)
    byy <- toPos.tensor(y,by)
    by <- by[!is.na(byx)&!is.na(byy)]
    byx <- toPos.tensor(x,by)
    byy <- toPos.tensor(y,by)
    alongx <- toPos.tensor(x,along)
    alongy <- toPos.tensor(y,along)
    if( any(dim(x)[alongx]!=dim(y)[alongy]) )
      stop("Data dimensions dont match in covariance of tensors")
    onex <-one.tensor(dim(x)[alongx])
    oney <-one.tensor(dim(y)[alongy])
    x[!is.finite(x)]<- 0
    y[!is.finite(y)]<- 0
    Sx <- mul.tensor( x  ,along , onex , 1:level.tensor(onex) ) 
    Sy <- mul.tensor( y  ,along , oney , 1:level.tensor(oney) ) 
    nx <- mul.tensor( Nx  ,along , onex , 1:level.tensor(onex) )
    ny <- mul.tensor( Ny  ,along , oney , 1:level.tensor(oney) )
    MeanX <- to.tensor(c(Sx)/c(nx),dim(Sx),dimnames(Sx))
    MeanY <- to.tensor(c(Sy)/c(ny),dim(Sy),dimnames(Sy))
    x <- x-MeanX # x <- x %-% MeanX
    y <- y-MeanY # y %-% MeanY
    S2 <- mul.tensor( x,along, y,along,by=by)
    N2 <- mul.tensor(Nx,along,Ny,along,by=by)
    S2/(N2-1)#to.tensor(c(S2)/(c(N2)-1),dim(S2),dimnames(S2))
  }
}


