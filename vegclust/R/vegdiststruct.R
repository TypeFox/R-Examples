vegdiststruct<-function(x, type="profile", method="bray", transform=NULL, classWidths=NULL) {
  method= match.arg(method,c("bray","ruzicka","kulczynski","ochiai", "canberra","relman"))
  type= match.arg(type,c("profile","abundance","volume"))
  if(inherits(x,"stratifiedvegdata")) {
    x = CAP(x, transform=transform)
  } else if(inherits(x,"CAP")) {
    if(!is.null(transform)) x = lapply(x, FUN=transform)
    class(x)<-c("list","CAP")
  } else{
    stop("Wrong data type for 'x'. Please use a 'stratifiedvegdata' object, or a 'CAP' object.")
  }
  Y = CAP2matrix(x,type=type, classWidths=classWidths)
  n = nrow(Y)
  res = matrix(0,n,n)
  Ysums = rowSums(Y)
  if(method=="bray") {
    for(i in 2:n) {
      for(j in 1:(n-1)) {
        A = sum(pmin(Y[i,],Y[j,]))
        res[i,j] = 1-(2*A/(Ysums[i]+Ysums[j]))
      }
    }
  } else if(method=="relman") {
    return(dist(decostand(Y,method="total"), method="manhattan"))
  } else if(method=="ruzicka") {
    for(i in 2:n) {
      for(j in 1:(n-1)) {
        A = sum(pmin(Y[i,],Y[j,]))
        res[i,j] = 1-(A/(Ysums[i]+Ysums[j]-A))
      }
    }
  } else if(method=="kulczynski") {
    for(i in 2:n) {
      for(j in 1:(n-1)) {
        A = sum(pmin(Y[i,],Y[j,]))
        res[i,j] = 1-0.5*((A/Ysums[i])+(A/Ysums[j]))
      }
    }
  } else if(method=="ochiai") {
    for(i in 2:n) {
      for(j in 1:(n-1)) {
        A = sum(pmin(Y[i,],Y[j,]))
        res[i,j] = 1-(A/sqrt(Ysums[i]*Ysums[j]))
      }
    }
  } else if(method=="canberra") {
    p = nrow(x[[1]]) #Number of species
    s = ncol(x[[1]]) #Number of strata
    print(p)
    print(s)
    for(i in 2:n) {
      for(j in 1:(n-1)) {
        num = colSums(matrix((abs(Y[i,]-Y[j,])), nrow=s, ncol=p))
        den = colSums(matrix((Y[i,]+Y[j,]), nrow=s, ncol=p))
        sel = (den>0)
        pp = sum(sel)
        num = num[sel]
        den = den[sel]
        res[i,j] = sum(num/den)/pp
      }
    }
  }
  rownames(res) = colnames(res) = names(x)
  return(as.dist(res))
}