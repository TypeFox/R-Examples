#Calculate Z-statistic
ztest<-function(n1,n2,p1,p2){
  return((p1-p2)/sqrt((p1*(1-p1)/n1)+(p2*(1-p2)/n2)))
}
#Calculate t-statistic
ttest<-function(x,y){
  return((mean(x)-mean(y))/sqrt((var(x)/length(x))+(var(y)/length(y))))
}
#Construct the set Omega (Hung, 2000)
Omega2.fun<-function(a,b){
  binary.repr<-function(x,digits){
    if(missing(digits)){
      mx<-max(x)
      digits<-if(mx>0) 1+floor(log(mx,base=2))
      else 1
    }
    ans<-0:(digits-1)
    lx<-length(x)
    bin<-(rep(x,rep(digits,lx))%/%2^ans)%%2
    dim(bin)<-c(digits,lx)
    return(bin)
  }
  in.Omega<-function(x,a,b){
    if(a==1 | b==1) return(TRUE)
    else{
      for(i in 1:(a-1)){for(j in 1:(b-1)){for(l in (i+1):a){
        pi<-c(x[i,j],x[i,b],x[l,j],x[l,b])
        if(all(pi==c(1,0,0,1)) | all(pi==c(0,1,1,0))) return(FALSE)
      }}}
    }
    return(TRUE)
  }
  if(b>1){
    Omega.reduce<-Omega2.fun(a,b-1)
    Omega<-array(dim=c(a,b,dim(Omega.reduce)[3]*2^a))
    Omega.add<-binary.repr(0:(2^a-1))
    for(i in 1:(2^a)){
      for(j in 1:dim(Omega.reduce)[3]){
        Omega[,,(i-1)*dim(Omega.reduce)[3]+j]<-cbind(Omega.reduce[,,j],Omega.add[,i])
      }
    }
    ab.in<-apply(Omega,3,in.Omega,a,b)
    return(array(Omega[,,ab.in],c(a,b,sum(ab.in))))
  }
  else return(array(binary.repr(0:(2^a-1)),c(a,1,2^a)))
}
#Check if data are binary
setGeneric("is.binary",function(x){standardGeneric("is.binary")})
setMethod("is.binary","numeric",function(x){
  y<-as.numeric(levels(as.factor(x)))
  if(length(y)>2) return(FALSE)
  return(all(y==c(0,1)) || all(y==0) || all(y==1))
})
setMethod("is.binary","list",function(x){
  for(i in 1:length(x)){
    y<-as.numeric(levels(as.factor(x[[i]])))
    if(length(y)>2) return(FALSE)
    if(!(all(y==c(0,1)) || all(y==0) || all(y==1))) return(FALSE)
  }
  return(TRUE)
})
setMethod("is.binary","ANY",function(x){return(FALSE)})
#Extract data from 'carpet' object
carpetdaten<-function(C){
  data<-numeric(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){
    data<-c(data,C@data[[(a*(C@D[2]+1))+(b+1)]])
  }}
  data
}
#Extract data from 'cube' object
cubedaten<-function(C){
  data<-numeric(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){for(c in 0:C@D[3]){
    data<-c(data,C@data[[(a*(C@D[2]+1)*(C@D[3]+1))+(b*(C@D[3]+1))+c+1]])
  }}}
  data
}
#Construct contrast matrix for comparisons
#of combinations to the marginals
carpetmatrix<-function(C){
  indexplus<-indexminus<-numeric(0)
  ctrnames<-grpnames<-character(0)
  cmatrix<-matrix(0,nrow=2*C@D[1]*C@D[2],ncol=(C@D[1]+1)*(C@D[2]+1))
  for(i in 1:C@D[1]){for(j in 1:C@D[2]){
    indexplus=c(indexplus,rep(i*(C@D[2]+1)+j+1,2))
    indexminus=c(indexminus,(i*(C@D[2]+1)+1),j+1)
    ctrnames=c(ctrnames,paste(i,j,"-",i,"0",sep=""),paste(i,j,"-0",j,sep=""))
  }}
  for(l in 1:(2*C@D[1]*C@D[2])){
    cmatrix[l,indexplus[l]]=1
    cmatrix[l,indexminus[l]]=-1
  }
  for(i in 0:C@D[1]){for(j in 0:C@D[2]){
    grpnames[i*(C@D[2]+1)+j+1]=paste(i,j,sep="")
  }}
  dimnames(cmatrix)=list(ctrnames,grpnames)
  cmatrix
}
#Construct character vector to represent combination group labels in a data frame
carpetkombi<-function(C){
  kombi<-character(0)
  cb<-function(a,b){(a*(C@D[2]+1))+b+1}
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){
    kombi<-c(kombi,rep(paste("(",a,",",b,")",sep=""),length(C@data[[cb(a,b)]])))
  }}
  return(as.factor(kombi))
}
#Construct characters to represent combination group labels
kombi<-function(D){
  kombi<-character(0)
  if(length(D)==2){
    for(a in 1:D[1]){for(b in 1:D[2]){kombi<-c(kombi,paste("(",a,",",b,")",sep=""))}}
  }
  if(length(D)==3){
    for(a in 1:D[1]){for(b in 1:D[2]){for(c in 1:D[3]){kombi<-c(kombi,paste("(",a,",",b,",",c,")",sep=""))}}}
  }
  return(kombi)
}
#Truncate of extend strings to a fixed length
fixdigit<-function(x,dig){
  long<-paste(x,"              ",sep="")
  substr(long,1,dig)
}
