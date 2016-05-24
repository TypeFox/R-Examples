testang.ltraj <- function(x, which=c("absolute","relative"),
                          nrep=999,
                          alter=c("two-sided","less","greater"))
{
    if (!inherits(x,"ltraj"))
        stop("x should be of class ltraj")
    if ((!is.regular(x))&attr(x,"typeII"))
        stop("x should be regular or of type I")

    whi.ang <- match.arg(which)

    foo <- function(df,whi.ang=whi.ang,
                    nrepet=nrep,
                    alter=match.arg(alter)){
        if(whi.ang=="absolute"){
            y <- df[,"abs.angle"]
        }
        else if(whi.ang=="relative"){
            y <- df[,"rel.angle"]
        }

        y <- .autourna(y)
        res <- .C("testindepangl", sim=as.double(rep(0,nrepet+1)),
                  as.double(y[[1]]), as.integer(length(y[[1]])),
                  as.integer(y[[2]]),as.integer(y[[3]]),
                  as.integer(length(y[[2]])),as.integer(nrepet),
                  PACKAGE = "adehabitatLT")
        return(as.randtest(obs=res$sim[1],sim=res$sim[-1],alter=alter))
    }
    return(lapply(x,foo,whi.ang=whi.ang,nrepet=nrep,alter=match.arg(alter)))
}

testdist.ltraj <- function(x,nrep=999,
                           alter=c("two-sided","less","greater"))
{
    if (!inherits(x,"ltraj"))
        stop("x should be of class ltraj")
    if ((!is.regular(x))&attr(x,"typeII"))
        stop("x should be regular or of type I")

    foo <- function(df,nrepet=nrep,alter=match.arg(alter)){
        y <- .autourna(df[,"dist"])
        res <- .C("testindepdist", sim=as.double(rep(0,nrepet+1)),
                  as.double(y[[1]]), as.integer(length(y[[1]])),
                  as.integer(y[[2]]), as.integer(y[[3]]),
                  as.integer(length(y[[2]])), as.integer(nrepet),
                  PACKAGE = "adehabitatLT")

        return(as.randtest(obs=res$sim[1],sim=res$sim[-1],alter=alter))
    }
    return(lapply(x,foo,alter=match.arg(alter),nrepet=nrep))
}

indmove.detail <- function(x, detail=c("dx","dy"),
                           nrep=999, alter=c("two-sided","less","greater"))
{
    if (!inherits(x,"ltraj"))
        stop("x should be of class ltraj")
    if ((!is.regular(x))&attr(x,"typeII"))
        stop("x should be regular or of type I")

    detail <- match.arg(detail)
    foo <- function(df,nrepet=nrep,alter=match.arg(alter)){
        y <- .autourna(df[,detail])
        res <- .C("testindepdist",sim=as.double(rep(0,nrepet+1)),
                  as.double(y[[1]]), as.integer(length(y[[1]])),
                  as.integer(y[[2]]),as.integer(y[[3]]),
                  as.integer(length(y[[2]])),as.integer(nrepet),
                  PACKAGE = "adehabitatLT")

        return(as.randtest(obs=res$sim[1],sim=res$sim[-1],alter=alter))
    }
    return(lapply(x,foo,nrepet=nrep,alter=match.arg(alter)))
}


.autourna <- function(x){
  w<-is.na(x)
  w2<-which(w)
  w3<-which(!is.na(x))
  f1 <- function(i){
    ii <- NA
    indx <- which(w3>w2[i])
    if(length(indx)>0) ii <- min(w3[indx])
    return(ii)
  }
  f2 <- function(i){
    ii <- NA
    indx <- which(w3<w2[i])
    if(length(indx)>0) ii <-  max(w3[indx])
    return(ii)
  }

  ind.sup<-min(w3)
  if(length(w2)>0) ind.sup <- c(ind.sup,unique(sapply(1:length(w2),f1)))
  ind.inf  <- max(w3)
  if(length(w2)>0) ind.inf <- c(unique(sapply(1:length(w2),f2)),ind.inf)
  res <- cbind(ind.sup,ind.inf)
  res <- na.omit(res)
  new.sup <- res[,1]-sapply(res[,1],function(x) sum(w2<x))
  new.inf <- res[,2]-sapply(res[,2],function(x) sum(w2<x))

  return(list(na.omit(x),new.sup,new.inf))

}
