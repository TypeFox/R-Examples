## Input: Variance matrix (or list of them). Output: Collapsed versions
varCollapse <- function(x,nl){
  if(is.list(x)){
    colx <- replicate(length(x),matrix(0,2*nl+1,2*nl+1),simplify=FALSE)
    names(colx) <- names(x)
  }
  else colx <- matrix(0,2*nl+1,2*nl+1)
  for(i in 0:nl){
    for(j in 0:(nl-i)){
      for(l in 0:nl){
        for(k in 0:(nl-l)){
          if(is.list(x)){
            for(n in 1:length(x)){
              colx[[n]][2*i+j+1,2*l+k+1] <- colx[[n]][2*i+j+1,2*l+k+1]+x[[n]][m2v(i=i+1,j=j+1,L=nl),m2v(i=l+1,j=k+1,L=nl)]
            }
          }
          else colx[2*i+j+1,2*l+k+1] <- colx[2*i+j+1,2*l+k+1]+x[m2v(i=i+1,j=j+1,L=nl),m2v(i=l+1,j=k+1,L=nl)]
        }
      }
    }
  }
  colx
}

## Updates alpha (the potens of the ps) which is used by Sa
updateAlpha <- function(a){
  aa <- replicate(length(a)-1,a[-length(a)],simplify=FALSE)
  for(i in 1:length(aa)) aa[[i]][i] <- aa[[i]][i]+a[length(a)]
  aa
}

## counts how many copies there is of each alpha-vector and returns the
## unique vectors together with their counts
permSummary <- function(x){
  y <- unlist(lapply(x,function(z) paste(sort(z),collapse="")))
  y <- table(y)
  z <- lapply(strsplit(names(y),""),as.numeric)
  list(type=z,terms=as.numeric(y))
}

## generates a vector of length x with a 1 in the k'th position
ek <- function(x,k) (0+(seq(along=x)==k))

## computes the sum over sum[] P(AAAAA) including theta correction
Sab <- function(a=rep(0,length(b)),b,t,p){
  if(all(t==0)) return(Sa(a+b,p=p))
  if(length(b)==0){
    ## cat(paste("Sa(a=c(",paste(a,collapse=","),"),p=p)*(1/(1-t))\n\n",sep=""))
    return(Sa(a,p)*(1/(1-t)))
  }
  else if(any(b==0)) return(Sab(a,b[b!=0],t,p))
  else{
    if(b[length(b)]==1){
      ## cat(paste("(1-",t,")*Sab(a=c(",paste(a+ek(a,length(b)),collapse=","),"),b=c(",paste(b[-length(b)],collapse=","),"),t=",t,",p=dkfreq[[1]]):\n",sep=""))
      return((1-t)*Sab(a+ek(a,length(b)),b[-length(b)],t,p)) #return((1-t)*Sab(a+ek(a,length(b)),b-ek(b,length(b)),t,p))
    }
    else{
      ## cat(paste((b[length(b)]-1),"*",t,"*Sab(a=c(",paste(a,collapse=","),"),b=c(",paste(b-ek(b,length(b)),collapse=","),"),t=",t,"p=dkfreq[[1]]) + (1-",t,")*Sab(a=c(",paste(a+ek(a,length(b)),collapse=","),"),b=c(",paste(b-ek(b,length(b)),collapse=","),"),t=",t,",p=dkfreq[[1]]):\n",sep=""))
      return((b[length(b)]-1)*t*Sab(a,b-ek(b,length(b)),t,p) + (1-t)*Sab(a+ek(a,length(b)),b-ek(b,length(b)),t,p))
    }
  }
}

## computes sum[i,...,j]^{different} P(A[i]...A[j])
Sa <- function(a,p){
  if(length(a)==1) return(sum((p^a)))
  else{
    perm <- permSummary(updateAlpha(a))
    return(Sa(a[length(a)],p=p)*Sa(a[-length(a)],p=p) - sum(perm$terms*unlist(lapply(perm$type,Sa,p=p))))
  }
}

## computes sum[i,...,j]^{different} P(A[i]...A[j])
Sa.list <- function(a,p){
  if(length(a)==1) return(unlist(lapply(p,function(x,y) sum(x^y),y=a)))
  else{
    perm <- permSummary(updateAlpha(a))
    perm.mat <- do.call("rbind",perm$type)
    return(Sa.list(a[length(a)],p=p)*Sa.list(a[-length(a)],p=p) - rowSums( apply(perm.mat,1,Sa.list,p=p)*rep(perm$terms,each=length(p)) ))
  }
}


## computes the sum over sum[] P(AAAAA) including theta correction
Sab.list <- function(a=rep(0,length(b)),b,t,p){
  if(all(t==0)) return(Sa.list(a+b,p=p))
  if(length(b)==0) return(Sa.list(a,p)*(1/(1-t)))
  else if(any(b==0)) return(Sab.list(a,b[b!=0],t,p))
  else{
    if(b[length(b)]==1) return((1-t)*Sab.list(a+ek(a,length(b)),b[-length(b)],t,p)) #return((1-t)*Sab(a+ek(a,length(b)),b-ek(b,length(b)),t,p))
    else return((b[length(b)]-1)*t*Sab.list(a,b-ek(b,length(b)),t,p) + (1-t)*Sab.list(a+ek(a,length(b)),b-ek(b,length(b)),t,p))
  }
}

Sab.listmat <- function(a=rep(0,length(b)),b,t,p){
  if(all(t==0)) return(matrix(Sa.list(a+b,p=p),length(t),length(p),byrow=TRUE,dimnames=list(theta=t,locus=names(p))))
  if(length(b)==0) return(matrix(Sa.list(a,p),length(t),length(p),byrow=TRUE,dimnames=list(theta=t,locus=names(p)))*(1/(1-t)))
  else if(any(b==0)) return(matrix(Sab.listmat(a,b[b!=0],t,p),length(t),length(p),byrow=FALSE,dimnames=list(theta=t,locus=names(p))))
  else{
    if(b[length(b)]==1) return(matrix(Sab.listmat(a+ek(a,length(b)),b[-length(b)],t,p),length(t),length(p),byrow=FALSE,dimnames=list(theta=t,locus=names(p)))*(1-t)) 
    else return((b[length(b)]-1)*matrix(Sab.listmat(a,b-ek(b,length(b)),t,p),length(t),length(p),byrow=FALSE,dimnames=list(theta=t,locus=names(p)))*t +
                matrix(Sab.listmat(a+ek(a,length(b)),b-ek(b,length(b)),t,p),length(t),length(p),byrow=FALSE,dimnames=list(theta=t,locus=names(p)))*(1-t))
  }
}

EE <- function(p,t){ ## computes both PsVar3 and PsVar4
  d <- unlist(lapply(t,function(x) prod(1+1:4*x)))
  D <- unlist(lapply(t,function(x) prod(1+1:6*x)))
### PsVar3
  pp <- list("111111"=Sab(b=rep(1,6),t=t,p=p),"11112"=Sab(b=c(rep(1,4),2),t=t,p=p),"1113"=Sab(b=c(1,1,1,3),t=t,p=p),
             "1122"=Sab(b=c(1,1,2,2),t=t,p=p),"114"=Sab(b=c(1,1,4),t=t,p=p),"123"=Sab(b=c(1,2,3),t=t,p=p),"33"=Sab(b=c(3,3),t=t,p=p),
             "15"=Sab(b=c(1,5),t=t,p=p),"222"=Sab(b=c(2,2,2),t=t,p=p),"24"=Sab(b=c(2,4),t=t,p=p),"6"=Sab(b=6,t=t,p=p))
### PsVar4
  PP <- list("11111111"=Sab(b=rep(1,8),t=t,p=p),"1111112"=Sab(b=c(rep(1,6),2),t=t,p=p),"111113"=Sab(b=c(rep(1,5),3),t=t,p=p),
             "111122"=Sab(b=c(rep(1,4),2,2),t=t,p=p),"11114"=Sab(b=c(rep(1,4),4),t=t,p=p),"1115"=Sab(b=c(rep(1,3),5),t=t,p=p),
             "116"=Sab(b=c(rep(1,2),6),t=t,p=p),"17"=Sab(b=c(1,7),t=t,p=p),"11123"=Sab(b=c(1,1,1,2,3),t=t,p=p),"11222"=Sab(b=c(1,1,2,2,2),t=t,p=p),
             "1124"=Sab(b=c(1,1,2,4),t=t,p=p),"1133"=Sab(b=c(1,1,3,3),t=t,p=p),"125"=Sab(b=c(1,2,5),t=t,p=p),
             "1223"=Sab(b=c(1,2,2,3),t=t,p=p),"134"=Sab(b=c(1,3,4),t=t,p=p),"2222"=Sab(b=rep(2,4),t=t,p=p),"224"=Sab(b=c(2,2,4),t=t,p=p),
             "233"=Sab(b=c(2,3,3),t=t,p=p),"26"=Sab(b=c(2,6),t=t,p=p),"35"=Sab(b=c(3,5),t=t,p=p),"44"=Sab(b=c(4,4),t=t,p=p),"8"=Sab(b=8,t=t,p=p))
### PsVar3  
  ## P_{0/0,0/0} = P_{0,0}:
  p11 <- pp$"111111" + 7*pp$"11112" + 4*pp$"1113" + 9*pp$"1122" + pp$"114"+ 4*pp$"123" + 3*pp$"222" + pp$"24"
  ## P_{0/1,0/0} = P_{1,0}: 
  p21 <- 4*pp$"11112" + 4*pp$"1113" + 12*pp$"1122" + 12*pp$"123" + 2*pp$"33"
  ## P_{1/0,0/0} = P_{2,0}: 
  p31 <- 2*pp$"1122" + pp$"114" + 2*pp$"222" + pp$"24"
  ## P_{0/0,0/1} = P_{0,1}: = P_{1,0} due to symmetry
  p12 <- p21
  ## P_{0/0,1/0} = P_{0,2}: = P_{2,0} due to symmetry
  p13 <- p31
  ## P_{0/1,0/1} = P_{1,1}: 
  p22 <- 8*pp$"1113" + 8*pp$"1122" + 12*pp$"114" + 16*pp$"123" + 2*pp$"15" + 8*pp$"222" + 4*pp$"24" + 2*pp$"33"
  ## P_{1/0,0/1} = P_{2,1}: 
  p32 <- 8*pp$"123" + 2*pp$"15" + 4*pp$"24"
  ## P_{0/1,1/0} = P_{1,2}: = P_{2,1} due to symmetry
  p23 <- p32
  ## P_{1/0,1/0} = P_{2,2}: 
  p33 <- 4*pp$"33" + pp$"6"
### PsVar4
  # P_{0/0,0/0} = P_{0,0}: 
  P11 <- PP$"11111111" + 20*PP$"1111112" + 16*PP$"111113" + 110*PP$"111122" + 4*PP$"11114" + 128*PP$"11123" + 164*PP$"11222" +
    24*PP$"1124" + 40*PP$"1133" + 144*PP$"1223" + 16*PP$"134" + 33*PP$"2222" + 12*PP$"224" + 24*PP$"233" + 2*PP$"44"      
  # P_{0/1,0/0} = P_{1,0}: 
  P21 <- 4*PP$"1111112" + 20*PP$"111113" + 40*PP$"111122" + 24*PP$"11114" + 152*PP$"11123" + 8*PP$"1115" + 84*PP$"11222" +
    104*PP$"1124" + 40*PP$"1133" + 196*PP$"1223" + 24*PP$"125" + 32*PP$"134" + 16*PP$"2222" + 32*PP$"224" + 48*PP$"233" + 8*PP$"35"     
  # P_{1/0,0/0} = P_{2,0}: 
  P31 <- 2*PP$"111122" + PP$"11114" + 16*PP$"11123" + 4*PP$"1115" + 4*PP$"11222" + 10*PP$"1124" + 24*PP$"1133" + 2*PP$"116" +
    16*PP$"1223" + 4*PP$"125" + 16*PP$"134" + 2*PP$"2222" + 9*PP$"224" + 8*PP$"233" + 2*PP$"26" + 4*PP$"44"
  # P_{0/0,0/1} = P_{0,1}: = P_{1,0} due to symmetry
  P12 <- P21
  # P_{0/0,1/0} = P_{0,2}: = P_{2,0} due to symmetry
  P13 <- P31
  # P_{0/1,0/1} = P_{1,1}: 
  P22 <- 16*PP$"111122" + 16*PP$"11114" + 96*PP$"11123" + 32*PP$"1115" + 64*PP$"11222" + 128*PP$"1124" + 112*PP$"1133" + 16*PP$"116" +
    192*PP$"1223" + 64*PP$"125" + 96*PP$"134" + 32*PP$"2222" + 96*PP$"224" + 80*PP$"233" + 16*PP$"26" + 16*PP$"44"
  # P_{1/0,0/1} = P_{2,1}: 
  P32 <- 8*PP$"11222" + 20*PP$"1124" + 4*PP$"116" + 40*PP$"1223" + 24*PP$"125" + 36*PP$"134" + 4*PP$"17" + 32*PP$"233" + 20*PP$"35"  
  # P_{0/1,1/0} = P_{1,2}: = P_{2,1} DUE TO SYMMETRY
  P23 <- P32
  # P_{1/0,1/0} = P_{2,2}: 
  P33 <- 4*PP$"2222" + 20*PP$"224" + 8*PP$"26" + 9*PP$"44" + PP$"8"
  Pij <- as.list(rep(0,length(t)))
  if(length(t)>1){
    for(i in 1:length(t)){
      Pij[[i]] <- list("E3"=matrix(c(p11[i],p12[i],p13[i],p21[i],p22[i],p23[i],p31[i],p32[i],p33[i]),3,3,byrow=TRUE)/d[i],
                       "E4"=matrix(c(P11[i],P12[i],P13[i],P21[i],P22[i],P23[i],P31[i],P32[i],P33[i]),3,3,byrow=TRUE)/D[i])
    }
  }
  else Pij <- list("E3"=matrix(c(p11,p12,p13,p21,p22,p23,p31,p32,p33),3,3,byrow=TRUE)/d,
                   "E4"=matrix(c(P11,P12,P13,P21,P22,P23,P31,P32,P33),3,3,byrow=TRUE)/D)
  Pij
}

EE.list <- function(p,t){ ## computes both PsVar3 and PsVar4
  d <- unlist(lapply(t,function(x) prod(1+1:4*x)))
  D <- unlist(lapply(t,function(x) prod(1+1:6*x)))
### PsVar3
  pp <- list("111111"=Sab.listmat(b=rep(1,6),t=t,p=p),"11112"=Sab.listmat(b=c(rep(1,4),2),t=t,p=p),"1113"=Sab.listmat(b=c(1,1,1,3),t=t,p=p),
             "1122"=Sab.listmat(b=c(1,1,2,2),t=t,p=p),"114"=Sab.listmat(b=c(1,1,4),t=t,p=p),"123"=Sab.listmat(b=c(1,2,3),t=t,p=p),"33"=Sab.listmat(b=c(3,3),t=t,p=p),
             "15"=Sab.listmat(b=c(1,5),t=t,p=p),"222"=Sab.listmat(b=c(2,2,2),t=t,p=p),"24"=Sab.listmat(b=c(2,4),t=t,p=p),"6"=Sab.listmat(b=6,t=t,p=p))
### PsVar4
  PP <- list("11111111"=Sab.listmat(b=rep(1,8),t=t,p=p),"1111112"=Sab.listmat(b=c(rep(1,6),2),t=t,p=p),"111113"=Sab.listmat(b=c(rep(1,5),3),t=t,p=p),
             "111122"=Sab.listmat(b=c(rep(1,4),2,2),t=t,p=p),"11114"=Sab.listmat(b=c(rep(1,4),4),t=t,p=p),"1115"=Sab.listmat(b=c(rep(1,3),5),t=t,p=p),
             "116"=Sab.listmat(b=c(rep(1,2),6),t=t,p=p),"17"=Sab.listmat(b=c(1,7),t=t,p=p),"11123"=Sab.listmat(b=c(1,1,1,2,3),t=t,p=p),"11222"=Sab.listmat(b=c(1,1,2,2,2),t=t,p=p),
             "1124"=Sab.listmat(b=c(1,1,2,4),t=t,p=p),"1133"=Sab.listmat(b=c(1,1,3,3),t=t,p=p),"125"=Sab.listmat(b=c(1,2,5),t=t,p=p),
             "1223"=Sab.listmat(b=c(1,2,2,3),t=t,p=p),"134"=Sab.listmat(b=c(1,3,4),t=t,p=p),"2222"=Sab.listmat(b=rep(2,4),t=t,p=p),"224"=Sab.listmat(b=c(2,2,4),t=t,p=p),
             "233"=Sab.listmat(b=c(2,3,3),t=t,p=p),"26"=Sab.listmat(b=c(2,6),t=t,p=p),"35"=Sab.listmat(b=c(3,5),t=t,p=p),"44"=Sab.listmat(b=c(4,4),t=t,p=p),"8"=Sab.listmat(b=8,t=t,p=p))
### PsVar3  
  ## P_{0/0,0/0} = P_{0,0}:
  p11 <- pp$"111111" + 7*pp$"11112" + 4*pp$"1113" + 9*pp$"1122" + pp$"114"+ 4*pp$"123" + 3*pp$"222" + pp$"24"
  ## P_{0/1,0/0} = P_{1,0}: 
  p21 <- 4*pp$"11112" + 4*pp$"1113" + 12*pp$"1122" + 12*pp$"123" + 2*pp$"33"
  ## P_{1/0,0/0} = P_{2,0}: 
  p31 <- 2*pp$"1122" + pp$"114" + 2*pp$"222" + pp$"24"
  ## P_{0/0,0/1} = P_{0,1}: = P_{1,0} due to symmetry
  p12 <- p21
  ## P_{0/0,1/0} = P_{0,2}: = P_{2,0} due to symmetry
  p13 <- p31
  ## P_{0/1,0/1} = P_{1,1}: 
  p22 <- 8*pp$"1113" + 8*pp$"1122" + 12*pp$"114" + 16*pp$"123" + 2*pp$"15" + 8*pp$"222" + 4*pp$"24" + 2*pp$"33"
  ## P_{1/0,0/1} = P_{2,1}: 
  p32 <- 8*pp$"123" + 2*pp$"15" + 4*pp$"24"
  ## P_{0/1,1/0} = P_{1,2}: = P_{2,1} due to symmetry
  p23 <- p32
  ## P_{1/0,1/0} = P_{2,2}: 
  p33 <- 4*pp$"33" + pp$"6"
### PsVar4
  # P_{0/0,0/0} = P_{0,0}: 
  P11 <- PP$"11111111" + 20*PP$"1111112" + 16*PP$"111113" + 110*PP$"111122" + 4*PP$"11114" + 128*PP$"11123" + 164*PP$"11222" +
    24*PP$"1124" + 40*PP$"1133" + 144*PP$"1223" + 16*PP$"134" + 33*PP$"2222" + 12*PP$"224" + 24*PP$"233" + 2*PP$"44"      
  # P_{0/1,0/0} = P_{1,0}: 
  P21 <- 4*PP$"1111112" + 20*PP$"111113" + 40*PP$"111122" + 24*PP$"11114" + 152*PP$"11123" + 8*PP$"1115" + 84*PP$"11222" +
    104*PP$"1124" + 40*PP$"1133" + 196*PP$"1223" + 24*PP$"125" + 32*PP$"134" + 16*PP$"2222" + 32*PP$"224" + 48*PP$"233" + 8*PP$"35"     
  # P_{1/0,0/0} = P_{2,0}: 
  P31 <- 2*PP$"111122" + PP$"11114" + 16*PP$"11123" + 4*PP$"1115" + 4*PP$"11222" + 10*PP$"1124" + 24*PP$"1133" + 2*PP$"116" +
    16*PP$"1223" + 4*PP$"125" + 16*PP$"134" + 2*PP$"2222" + 9*PP$"224" + 8*PP$"233" + 2*PP$"26" + 4*PP$"44"
  # P_{0/0,0/1} = P_{0,1}: = P_{1,0} due to symmetry
  P12 <- P21
  # P_{0/0,1/0} = P_{0,2}: = P_{2,0} due to symmetry
  P13 <- P31
  # P_{0/1,0/1} = P_{1,1}: 
  P22 <- 16*PP$"111122" + 16*PP$"11114" + 96*PP$"11123" + 32*PP$"1115" + 64*PP$"11222" + 128*PP$"1124" + 112*PP$"1133" + 16*PP$"116" +
    192*PP$"1223" + 64*PP$"125" + 96*PP$"134" + 32*PP$"2222" + 96*PP$"224" + 80*PP$"233" + 16*PP$"26" + 16*PP$"44"
  # P_{1/0,0/1} = P_{2,1}: 
  P32 <- 8*PP$"11222" + 20*PP$"1124" + 4*PP$"116" + 40*PP$"1223" + 24*PP$"125" + 36*PP$"134" + 4*PP$"17" + 32*PP$"233" + 20*PP$"35"  
  # P_{0/1,1/0} = P_{1,2}: = P_{2,1} DUE TO SYMMETRY
  P23 <- P32
  # P_{1/0,1/0} = P_{2,2}: 
  P33 <- 4*PP$"2222" + 20*PP$"224" + 8*PP$"26" + 9*PP$"44" + PP$"8"
  if(length(t)>1){
    Pij <- replicate(length(t),as.list(rep(0,length(p))),simplify=FALSE)
    for(i in 1:length(t)){
      for(j in 1:length(p))
        Pij[[i]][[j]] <- list("E3"=matrix(c(p11[i,j],p12[i,j],p13[i,j],p21[i,j],p22[i,j],p23[i,j],p31[i,j],p32[i,j],p33[i,j]),3,3,byrow=TRUE)/d[i],
                              "E4"=matrix(c(P11[i,j],P12[i,j],P13[i,j],P21[i,j],P22[i,j],P23[i,j],P31[i,j],P32[i,j],P33[i,j]),3,3,byrow=TRUE)/D[i])
      names(Pij[[i]]) <- names(p)
    }
    names(Pij) <- paste(t)
  }
  else{
    Pij <- as.list(rep(0,length(p)))
    for(j in 1:length(p))
      Pij[[j]] <- list("E3"=matrix(c(p11[j],p12[j],p13[j],p21[j],p22[j],p23[j],p31[j],p32[j],p33[j]),3,3,byrow=TRUE)/d,
                       "E4"=matrix(c(P11[j],P12[j],P13[j],P21[j],P22[j],P23[j],P31[j],P32[j],P33[j]),3,3,byrow=TRUE)/D)
    names(Pij) <- names(p)
  }
  Pij
}

## maps from M-matrix to vector
m2v <- function(i,j,L) ((i-1)*(L+1) - (i-1)*(i-2)/2 + j) ## minus one from each index (i,j) because input format is i+1 and j+1

## Makes the recursion over loci. I.e. builds up the expected values by adding one locus at the time
rareVar <- function(qu){ #(q,u){
  q <- lapply(qu,function(z) z[[1]])
  u <- lapply(qu,function(z) z[[2]])
  S <- length(q)
  M <- replicate(S,matrix(0,(S+1)*(S+2)/2,(S+1)*(S+2)/2),simplify=FALSE)
  N <- replicate(S,matrix(0,(S+1)*(S+2)/2,(S+1)*(S+2)/2),simplify=FALSE)
  locus1a <- m2v(rep(c(1,1,2),3),rep(c(1,2,1),3),S)
  locus1b <- m2v(rep(c(1,2),c(6,3)),rep(c(1,2,1),each=3),S)
  q1 <- as.numeric(q[[1]])
  u1 <- as.numeric(u[[1]])
  for(i in 1:9){ ## 1:9 comes from the length of locus1a/b
    M[[1]][locus1a[i],locus1b[i]] <- q1[i]
    N[[1]][locus1a[i],locus1b[i]] <- u1[i]
  }
  for(s in 2:S){
    for(m in 1:(s+1)){
      for(mm in 1:(s+1)){
        for(p in 1:(s-m+2)){
          for(pp in 1:(s-mm+2)){
            if((m==1 & p==1) & (mm==1 & pp==1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)]
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)]
            }          
            else if((m==1 & p==1) & (mm==1 & pp>1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][1,2]*M[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][1,2]*N[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)])
            }
            else if((m==1 & p==1) & (mm>1 & pp==1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][1,3]*M[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][1,3]*N[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)])
            }
            else if((m==1 & p==1) & (mm>1 & pp>1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][1,2]*M[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + q[[s]][1,3]*M[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][1,2]*N[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + u[[s]][1,3]*N[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)])
            }
            else if((m==1 & p>1) & (mm==1 & pp==1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][2,1]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][2,1]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)])
            }
            else if((m==1 & p>1) & (mm==1 & pp>1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][2,1]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] +
                                                  q[[s]][1,2]*M[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + q[[s]][2,2]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp-1,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][2,1]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] +
                                                  u[[s]][1,2]*N[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + u[[s]][2,2]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp-1,S)])
            }
            else if((m==1 & p>1) & (mm>1 & pp==1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][1,3]*M[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] +
                                                  q[[s]][2,1]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + q[[s]][2,3]*M[[s-1]][m2v(m,p-1,S),m2v(mm-1,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][1,3]*N[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] +
                                                  u[[s]][2,1]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + u[[s]][2,3]*N[[s-1]][m2v(m,p-1,S),m2v(mm-1,pp,S)])
            }
            else if((m==1 & p>1) & (mm>1 & pp>1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][2,1]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + q[[s]][1,2]*M[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] +
                                                  q[[s]][1,3]*M[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + q[[s]][2,2]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp-1,S)] + q[[s]][2,3]*M[[s-1]][m2v(m,p-1,S),m2v(mm-1,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][2,1]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + u[[s]][1,2]*N[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] +
                                                  u[[s]][1,3]*N[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + u[[s]][2,2]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp-1,S)] + u[[s]][2,3]*N[[s-1]][m2v(m,p-1,S),m2v(mm-1,pp,S)])
            }
            else if((m>1 & p==1) & (mm==1 & pp==1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][3,1]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][3,1]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)])
            }
            else if((m>1 & p==1) & (mm==1 & pp>1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][3,1]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  q[[s]][1,2]*M[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + q[[s]][3,2]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp-1,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][3,1]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  u[[s]][1,2]*N[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + u[[s]][3,2]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp-1,S)])
            }
            else if((m>1 & p==1) & (mm>1 & pp==1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][3,1]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  q[[s]][1,3]*M[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + q[[s]][3,3]*M[[s-1]][m2v(m-1,p,S),m2v(mm-1,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][3,1]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  u[[s]][1,3]*N[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + u[[s]][3,3]*N[[s-1]][m2v(m-1,p,S),m2v(mm-1,pp,S)])
            }
            else if((m>1 & p==1) & (mm>1 & pp>1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][3,1]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] + q[[s]][1,2]*M[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] +
                                                  q[[s]][1,3]*M[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + q[[s]][3,2]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp-1,S)] + q[[s]][3,3]*M[[s-1]][m2v(m-1,p,S),m2v(mm-1,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][3,1]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] + u[[s]][1,2]*N[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] +
                                                  u[[s]][1,3]*N[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + u[[s]][3,2]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp-1,S)] + u[[s]][3,3]*N[[s-1]][m2v(m-1,p,S),m2v(mm-1,pp,S)])
            }
            else if((m>1 & p>1) & (mm==1 & pp==1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][2,1]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + q[[s]][3,1]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][2,1]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + u[[s]][3,1]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)])
            }
            else if((m>1 & p>1) & (mm==1 & pp>1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][2,1]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + q[[s]][3,1]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  q[[s]][1,2]*M[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + q[[s]][2,2]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp-1,S)] + q[[s]][3,2]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp-1,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][2,1]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + u[[s]][3,1]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  u[[s]][1,2]*N[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + u[[s]][2,2]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp-1,S)] + u[[s]][3,2]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp-1,S)])
            }
            else if((m>1 & p>1) & (mm>1 & pp==1)){
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][2,1]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + q[[s]][3,1]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  q[[s]][1,3]*M[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + q[[s]][2,3]*M[[s-1]][m2v(m,p-1,S),m2v(mm-1,pp,S)] + q[[s]][3,3]*M[[s-1]][m2v(m-1,p,S),m2v(mm-1,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][2,1]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + u[[s]][3,1]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  u[[s]][1,3]*N[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + u[[s]][2,3]*N[[s-1]][m2v(m,p-1,S),m2v(mm-1,pp,S)] + u[[s]][3,3]*N[[s-1]][m2v(m-1,p,S),m2v(mm-1,pp,S)])
            }
            # if((m>1 & p>1) & (mm>1 & pp>1)){} ... SAME AS ELSE:
            else{
              M[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (q[[s]][1,1]*M[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + q[[s]][2,1]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + q[[s]][3,1]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  q[[s]][1,2]*M[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + q[[s]][1,3]*M[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + q[[s]][2,2]*M[[s-1]][m2v(m,p-1,S),m2v(mm,pp-1,S)] +
                                                  q[[s]][3,2]*M[[s-1]][m2v(m-1,p,S),m2v(mm,pp-1,S)] + q[[s]][2,3]*M[[s-1]][m2v(m,p-1,S),m2v(mm-1,pp,S)] + q[[s]][3,3]*M[[s-1]][m2v(m-1,p,S),m2v(mm-1,pp,S)])
              #
              N[[s]][m2v(m,p,S),m2v(mm,pp,S)] <- (u[[s]][1,1]*N[[s-1]][m2v(m,p,S),m2v(mm,pp,S)] + u[[s]][2,1]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp,S)] + u[[s]][3,1]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp,S)] +
                                                  u[[s]][1,2]*N[[s-1]][m2v(m,p,S),m2v(mm,pp-1,S)] + u[[s]][1,3]*N[[s-1]][m2v(m,p,S),m2v(mm-1,pp,S)] + u[[s]][2,2]*N[[s-1]][m2v(m,p-1,S),m2v(mm,pp-1,S)] +
                                                  u[[s]][3,2]*N[[s-1]][m2v(m-1,p,S),m2v(mm,pp-1,S)] + u[[s]][2,3]*N[[s-1]][m2v(m,p-1,S),m2v(mm-1,pp,S)] + u[[s]][3,3]*N[[s-1]][m2v(m-1,p,S),m2v(mm-1,pp,S)])
            }
          }
        }
      }
    }
  }
  x <- M[[S]]
  y <- N[[S]]
  dimnames(x) <- dimnames(y) <- list(dbCats(S,vector=TRUE),dbCats(S,vector=TRUE))
  list(E3=x,E4=y)
}

dbVariance <-  function(probs,theta=0,n=1,collapse=FALSE){
  EEs <- EE.list(probs,t=theta)
  ##  Es <- rareVar(q=lapply(q,PsVar3,t=theta),u=lapply(q,PsVar4,t=theta))
  if(length(theta)>1){
    V <- as.list(rep(0,length(theta)))
    for(t in 1:length(theta)){
      E <- dbExpect(probs,theta=theta[t],vector=TRUE)
      Es <- rareVar(EEs[[t]])
      V1 <- diag(E) - E%*%t(E)
      V2 <- Es$E3 - E%*%t(E)
      V3 <- Es$E4 - E%*%t(E)
      if(n==1) V[[t]] <- list(V1=V1,V2=V2,V3=V3)
      else V[[t]] <- choose(n,2)*V1 + 6*choose(n,3)*V2 + 6*choose(n,4)*V3
      if(collapse) V[[t]] <- varCollapse(V[[t]],nl=length(probs))
    }
    names(V) <- paste(theta)
  }
  else{
    E <- dbExpect(probs,theta=theta,vector=TRUE)
    Es <- rareVar(EEs)
    V1 <- diag(E) - E%*%t(E)
    V2 <- Es$E3 - E%*%t(E)
    V3 <- Es$E4 - E%*%t(E)
    if(n==1) V <- list(V1=V1,V2=V2,V3=V3)
    else V <- choose(n,2)*V1 + 6*choose(n,3)*V2 + 6*choose(n,4)*V3
    if(collapse) V <- varCollapse(V,nl=length(probs))
  }
  V
}

