require(tensorA)
if(FALSE) {
  # Commands for testing
  debugger()
  rm(list=objects())
  options(error=dump.frames)
  detach("package:tensorA")
  dyn.unload("/home/boogaart/R/tensorA/tests/../../tensorA.Rcheck/tensorA/libs/tensorA.so")
  library(tensorA,lib.loc="../../tensorA.Rcheck")
}

set.seed(23)


summary.tensor <- function(x,...) {
  n <- level.tensor(x)
  d <- dim(x)
  dm <- pmin(d,1)
  dm2<- pmin(d,2)
  print(dim(x))
  if( !all(sapply(dimnames(x),is.null))) print(dimnames(x))
  if( prod(dim(x)) < 10 )
    print(x)
  else if( n == 1 )
    print(x[1:min(length(x),10)])
  else if( n == 2 )
    print(x[1:dm[1],1:dm[2]])
  else if( n == 3 )
    print(x[1:dm[1],1:dm[2],1:dm2[3],drop=FALSE])
  else if( n == 4 )
    print(x[1:(dm2[1]),1:dm[1],1,1,drop=FALSE])
  else if( n == 5 )
    print(x[1:(dm2[1]),1:dm[1],1,1,1,drop=FALSE])
  else if( n == 6 )
    print(x[1:(dm2[1]),1:dm[1],1,1,1,1,drop=FALSE])
  else
    print(x[1:min(length(x),10)])
}

checker <- function(x,y) UseMethod("checker")

checker.default <- function(x,y) {
    cat("Wrong type\n",deparse(match.call()),"\n---------------y=\n")
    print(y)
    cat("\n------------x=\n")
    print(x)
    stop("Unkown type")
}


checker.tensor <- function(x,y) {
  if( !cmp(x,y) ) {
    cat("Misfit\n",deparse(match.call()),"\n------------------y=\n")
    print(y)
    cat("\n-----------------------x=\n")
    print(x)
#    print(summary(x))
#    print(summary(y))
    stop("Missfit");
  } else {
    x
  }
}

print.tensor <- function(x){
  print.default(unclass(x))
  print(dim(x))
}

cmp.tensor <- function(x,y) {
  if(!is.null(names(y))) {
    mat <- match(names(y),names(x))
    if( any(is.na(mat)) ) {
      print(names(y))
      print(names(x))
      return(FALSE)
    }
    else
      x <- reorder.tensor(x,mat)
  }
  return( length(x)==length(y) && !any(is.na(x)) && 
         sum(abs(c(x)-c(y))^2)<1E-10 &&
         length(dim(x))==length(dim(y)) &&
         !any(is.na(dim(x))) && !any(is.na(dim(y))) &&         
         all(dim(x)==dim(y)))
}
  
cmp <- function(x,y) UseMethod("cmp")

cmp.character <- function(x,y) {
  return( length(x)==length(y) && all(x==y) )
}


cmp.numeric <- function(x,y){
  if( length(x) != length(y) )
    return(FALSE)
  if(!is.null(names(y))) {
    if( is.null(names(y)) || !all( names(x) %in% names(y) ) )
      return(FALSE)
    if( ! identical(any( x[match(names(y),names(x))]!=y ),FALSE) )
      return(FALSE)
    return(TRUE)
  } else {
    return(all(x==y))
  }
}

cmp.default <- function(x,y) identical(x,y)

checker.numeric <- function(x,y) {
  if( !cmp(x,y) ) {
    cat("Misfit\n",deparse(match.call()),"\n-----------------y=\n")
    print(y)
    cat("\n-------------------------x=\n")
    print(x)
    stop("Missfit");
  } else print(x)
}


checker.character <- function(x,y) {
  if( !cmp(x,y) ) {
    cat("Misfit\n",deparse(match.call()),"\n-----------------y=\n")
    print(y)
    cat("\n----------------------------y=\n")
    print(x)
    stop("Missfit");
  } else print(x)
}

checker.list <- function(x,y) {
  if( length(x) != length(y) || identical(!all(mapply(cmp,x,y)),TRUE) ) {
    cat("Misfit\n",deparse(match.call()),"\n---------------y=\n")
    print(y)
    cat("\n----------------------x=\n")
    print(x)
    stop("Missfit");
  } else print(x)    
}

nn <- function(...,nn=list(...)) {
  lapply(1:length(nn),function(i) paste(names(nn)[i],1:nn[[i]],sep=""))
 }



# to.tensor
I4 <- diag(4)
names(dim(I4))<-c("Q","z")
checker(to.tensor(I4,c(a=2,b=2)),
        to.tensor(c(diag(4)),c(a=2,b=2,z=4)))

dimnames(I4) <- nn(Q=4,z=4)
checker(to.tensor(I4,c(a=2,b=2),nn(a=2,b=2)),
        to.tensor(c(diag(4)),c(a=2,b=2,z=4),nn(a=2,b=2,z=4)))


#

KT5 <- to.tensor(rnorm(30),c(a=3,b=2,c=5),nn(a=3,b=2,c=5))
dim(KT5[[c=2]])
dim(KT5[,,1:2])
summary(KT5)

KT5 <- to.tensor(rnorm(30),c(a=3,b=10,c=1),nn(a=3,b=10,c=1))
drop(KT5)

checker(as.tensor(KT5),KT5)
checker(as.tensor.default(KT5),KT5)
checker(to.tensor(KT5),KT5)

# 
I <- to.tensor(diag(3),c(a=3,b=3),what=1:2)
checker(to.tensor(I),I)
checker(as.tensor(I),I)
checker(inv.tensor(I,"a","b"),I)

R1  <- matrix(rnorm(9),nrow=3)
R1i <- solve(R1)
R2 <- to.tensor(R1,c(a=3,b=3),what=1:2)
R2i <- to.tensor(R1i,c(b=3,a=3),what=1:2)

checker(inv.tensor(R2,"a","b"),R2i)
checker(inv.tensor(R2,"a","b",allowSingular=TRUE),R2i)

checker(inv.tensor(rep(R2,4,1,"K"),"a","b",by="K"),rep(R2i,4,1,"K"))
checker(inv.tensor(rep(R2,4,1,"K"),"a","b",by="K",allowSingular=TRUE),rep(R2i,4,3,"K"))

R3 <- to.tensor(rnorm(15),c(a=3,z=5))

checker(mul.tensor(R2i,"b",mul.tensor(R2,"a",R3)),R3)

checker(solve.tensor(R2i,R3[[z=1]],"a"),mul.tensor(R2,"a",R3[[z=1]]))
checker(solve.tensor(R2i,R3,"a"),mul.tensor(R2,"a",R3))

checker(solve.tensor(R2i,R3[[z=1]],"a",allowSingular=TRUE),mul.tensor(R2,"a",R3[[z=1]]))
checker(solve.tensor(R2i,R3,"a",allowSingular=T),mul.tensor(R2,"a",R3))

checker(solve.tensor(rep(R2i,4,1,"K"),R3[[z=1]],"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(rep(R2i,4,1,"K"),rep(R3[[z=1]],4,1,"K"),"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(R2i,rep(R3[[z=1]],4,1,"K"),"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))

checker(solve.tensor(rep(R2i,4,1,"K"),R3[[z=1]],"a",by="K",allowSingular=TRUE),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(rep(R2i,4,1,"K"),rep(R3[[z=1]],4,1,"K"),"a",by="K",allowSingular=TRUE),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(R2i,rep(R3[[z=1]],4,1,"K"),"a",by="K",allowSingular=TRUE),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))



checker(solve.tensor(R2i,R3,"a"),mul.tensor(R2,"a",R3))

checker(solve.tensor(R2i,R3[[z=1]],"a",allowSingular=TRUE),mul.tensor(R2,"a",R3[[z=1]]))
checker(solve.tensor(R2i,R3,"a",allowSingular=TRUE),mul.tensor(R2,"a",R3))



summary(I)


A <- to.tensor(c(diag(3)),c(a=3,b=3))

A <- to.tensor(c(1,1,1,0,1,1,0,0,1),c(a=3,b=3))
checker(mul.tensor(A,"b",A,"a"),to.tensor(c(A%*%A),c(a=3,b=3)))

A <- to.tensor(c(1,1,1,1,0,1,1,1,0,0,1,1),c(a=4,b=3))
checker(mul.tensor(A,"b",A[[a=~c]],"b"),to.tensor(c(A%*%t(A)),c(a=4,c=4)))

A <- to.tensor(rnorm(15),c(a=5,b=3))
checker(mul.tensor(A,"b",A[[a=~c]],"b"),to.tensor(c(A%*%t(A)),dim=c(a=5,c=5)))

A <- to.tensor(rnorm(5*3),c(a=5,b=3))
B <- to.tensor(rnorm(3*17),c(a=3,b=17))
checker(mul.tensor(A,"b",B,"a"),to.tensor(c(A%*%B),c(a=5,b=17)))

A <- to.tensor(c(1,1,1,0,1,1,0,0,1),c(a=3,b=3))

# 
A <- to.tensor(c(1,1,1,0,1,1,0,0,1),c(a=3,b=3))

solve.tensor(mul.tensor(I,"b",I,"a"),I,"a","a")
checker(solve.tensor(mul.tensor(I,"b",I,"a"),I,"a","b"),I)
mul.tensor(A,"b",I,"a")


checker(solve.tensor(mul.tensor(A,"b",I[[a=~c]],"b"),A,"a","a",allowSingular=TRUE),I[[a=~c]])

B <- A[[b=~z]]
checker(solve.tensor(B,mul.tensor(A,"a",B,"a"),"z","z",allowSingular=TRUE),A)
checker(solve.tensor(B,mul.tensor(A,"a",B,"a"),"z","z",allowSingular=TRUE),
        solve.tensor(B,mul.tensor(A,"a",B,"a"),"z","z"))



solve.tensor(mul.tensor(A,"b",A,"a"),A,"a","a")
checker(solve.tensor(A[[b=~c]],mul.tensor(A,"b",A,"a"),"a","a"),structure(A,dim=c(c=3,b=3)))

A <-  to.tensor(c(1,1,1,1,0,1,1,1,0,0,1,1),c(a=4,b=3))
A
solve.tensor(A[[b=~c]],mul.tensor(A,"a",A[[b=~c]],"a"),"c",allowSingular=TRUE)


A <- to.tensor(rnorm(100),c(a=4,b=5,c=5))
An <- to.tensor(rnorm(100),c(a=4,b=5,c=5),nn(a=4,b=5,c=5))
B <- to.tensor(rnorm(100),c(d=4,e=5,f=5))
An <- to.tensor(rnorm(100),c(a=4,b=5,c=5),nn(a=4,b=5,c=5))

mt <- mul.tensor(A,2,B,2)


checker(einstein.tensor(A,b="e",B),mt)


G <- to.tensor(rnorm(20*20*3),c(a=4,b=5,a1=4,b1=5,I=3))
mul.tensor(A,c("a","b"),G,c("a","b"))



checker(einstein.tensor(A,G),einstein.tensor(G,A))

checker(einstein.tensor(A,G),einstein.tensor(G,A))

checker(einstein.tensor(einstein.tensor(A,G),inv.tensor(G,c("a","b"),by="I"),by=c("I")),rep(A,dim(G)["I"],1,"I"))

einstein.tensor(A,G)

## chol

A <- to.tensor(rnorm(15),c(a=3,b=5))
AAt <- einstein.tensor(A,mark(A,i="a"))
ch <- chol.tensor(AAt,"a","a'",name="lambda")
#names(ch)[1]<-"lambda"
checker(einstein.tensor(ch,mark(ch,i="a")),AAt)

A <- to.tensor(rnorm(30),c(a=3,b=5,c=2))
AAt <- einstein.tensor(A,mark(A,i="a"),by="c")
ch <- chol.tensor(AAt,"a","a'",name="lambda")
checker(einstein.tensor(ch,mark(ch,i="a"),by="c"),AAt)

ftable(A)

# norm

A <- to.tensor(c(1,1,1,1,0,1,1,1,0,0,1,1),c(a=4,b=3))
checker(norm.tensor(A),sqrt(9))
checker(norm.tensor(A,c(1,2)),sqrt(9))
checker(norm.tensor(A,"b"),to.tensor(sqrt(c(1,2,3,3)),c(a=4)))
checker(norm.tensor(A,"a"),to.tensor(sqrt(c(4,3,2)),c(b=3)))
checker(norm.tensor(A,by="a"),to.tensor(sqrt(c(1,2,3,3)),c(a=4)))
checker(norm.tensor(A,by="b"),to.tensor(sqrt(c(4,3,2)),c(b=3)))

# opnorm

A <- to.tensor(c(1,0,0,0,1,0,0,1),c(a=2,b=2,s=2))
checker(opnorm(A,"a",by="s"),to.tensor(c(1,1),c(s=2)))


# margin


A <- to.tensor(rnorm(30),c(a=3,b=2,c=5))
checker( margin.tensor(A,c("a","c")),einstein.tensor(A,one.tensor(c(a=3,c=5))))
checker( margin.tensor(A,by=c("a","c")),einstein.tensor(A,one.tensor(c(b=2))))
checker(one.tensor(c(a=3,c=5)),to.tensor(rep(1,15),c(a=3,c=5)))

# diagmul

A <- to.tensor(rnorm(30),c(a=3,b=2,c=5))
B <- to.tensor(rnorm(6),c(a=3,b=2))

checker(einstein.tensor(A,diag.tensor(B,mark="m")),
        diagmul.tensor(A,B,i=c("a","b"))[[a=~am,b=~bm]])

# is.tensor
if( !identical(c(
is.tensor(FALSE),
is.tensor(TRUE),
is.tensor(matrix(1:3)),
is.tensor(as.tensor(matrix(1:3))),
is.tensor(A)),c(FALSE,FALSE,FALSE,TRUE,TRUE)))
  stop("Fehler")

#

A <- to.tensor(rnorm(30),c(a=3,b=2,c=5))
checker(reorder.tensor(reorder.tensor(reorder.tensor(A,c("c","b","a")),c(1,3,2)),c("a","b","c")),A)




mul.tensor(A,c(),B,c(),by=c("a","b"))
#### complex

A <- to.tensor( c(1+1i,1,0,2-13i) , c(a=2,b=2) )
B <- to.tensor( c(1+1i,1,0,2-13i) , c(a=2,b=2) )
A
Am <- matrix(c(A),nrow=nrow(A))
Bm <- matrix(c(B),nrow=nrow(B))
checker(mul.tensor(A,"b",B,"a"),to.tensor(c(Am%*%Bm),c(nrow(A),ncol(B))))


R1  <- matrix(rnorm(9)+rnorm(9)*1i,nrow=3)
R1i <- solve(R1)
R2 <- to.tensor(R1,c(a=3,b=3),what=1:2)
R2i <- to.tensor(R1i,c(b=3,a=3),what=1:2)

checker(inv.tensor(R2,"a","b"),R2i)
checker(inv.tensor(R2,"a","b",allowSingular=TRUE),R2i)

checker(inv.tensor(rep(R2,4,1,"K"),"a","b",by="K"),rep(R2i,4,1,"K"))
checker(inv.tensor(rep(R2,4,1,"K"),"a","b",by="K",allowSingular=TRUE),rep(R2i,4,3,"K"))

R3 <- to.tensor(rnorm(15),c(a=3,z=5))

checker(mul.tensor(R2i,"b",mul.tensor(R2,"a",R3)),R3)



checker(solve.tensor(R2i,R3[[z=1]],"a"),mul.tensor(R2,"a",R3[[z=1]]))
checker(solve.tensor(R2i,R3,"a"),mul.tensor(R2,"a",R3))

checker(solve.tensor(R2i,R3[[z=1]],"a",allowSingular=TRUE),mul.tensor(R2,"a",R3[[z=1]]))
checker(solve.tensor(R2i,R3,"a",allowSingular=T),mul.tensor(R2,"a",R3))

checker(solve.tensor(rep(R2i,4,1,"K"),R3[[z=1]],"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(rep(R2i,4,1,"K"),rep(R3[[z=1]],4,1,"K"),"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(R2i,rep(R3[[z=1]],4,1,"K"),"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))

checker(solve.tensor(rep(R2i,4,1,"K"),R3[[z=1]],"a",by="K",allowSingular=T),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(rep(R2i,4,1,"K"),rep(R3[[z=1]],4,1,"K"),"a",by="K",allowSingular=T),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(R2i,rep(R3[[z=1]],4,1,"K"),"a",by="K",allowSingular=T),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))

# trace.tensor

A <- to.tensor(rep(1,16),c(a=4,b=4))
checker(trace.tensor(A,"a","b"),4)

A <- to.tensor(rep(1,16*3),c(a=4,c=3,b=4))
checker(trace.tensor(A,"a","b"),to.tensor(rep(4,3),c(c=3)))

A <- to.tensor(1:(2*2*3*3*5),c(a=2,b=3,c=2,d=3,e=5))
erg <- sapply(1:5,function(e) sum(diag(matrix(A[[e=e]],nrow=6))))
checker(trace.tensor(reorder(A,c(3,4,1,5,2)),c("a","b"),c("c","d")),to.tensor(c(erg),c(e=5)))

# delta.tensor

checker(delta.tensor(c(a=2,b=3)),to.tensor(c(diag(6)),c(a=2,b=3,"a'"=2,"b'"=3)))
# diag.tensor

A <- to.tensor(1:6,c(a=2,b=3))
checker(diag.tensor(A),to.tensor(c(diag(c(A))),c(a=2,b=3,"a'"=2,"b'"=3)))
A <- to.tensor(1:6,c(a=2,b=3))
checker(diag.tensor(A,by="b"),to.tensor(c(1,0,0,2,3,0,0,4,5,0,0,6),c(a=2,"a'"=2,b=3)))

# tripledelta.tensor
checker(tripledelta.tensor(c(a=2)),to.tensor(c(1,0,0,0,0,0,0,1),c(a=2,"a'"=2,"a*"=2)))

checker(tripledelta.tensor(c(a=2,b=4)),
        einstein.tensor(tripledelta.tensor(c(a=2)),tripledelta.tensor(c(b=4)) )
        )


# one.tensor
checker( one.tensor(c(a=3,b=4)), to.tensor(rep(1,12),c(a=3,b=4)))


checker(level.tensor(A),length(dim(A)))


# svd.tensor

A <- to.tensor(rnorm(120),c(a=2,b=2,c=5,d=3,e=2))

SVD <- svd.tensor(A,c("a","d"),c("b","c"),by="e")
dim(SVD$v)
# Kompositionseigenschaft
checker(einstein.tensor(SVD$v,diag=SVD$d,SVD$u,by="e"),A)
# Orthogonalitaet:
checker( SVD$v %e% SVD$v[[lambda=~"lambda'"]],2*delta.tensor(c(lambda=6)))
checker( SVD$u %e% SVD$u[[lambda=~"lambda'"]],2*delta.tensor(c(lambda=6)))
checker( SVD$u %e% mark(SVD$u,"'",c("a","d")),2*delta.tensor(c(a=2,d=3)))


# power.tensor
A <- to.tensor(rnorm(120),c(a=2,b=2,c=5,d=3,e=2))
AAt <- A %e% mark(A,"'",c("a","b"))

checker(power.tensor(AAt,c("a","b"),c("a'","b'"),-1),inv.tensor(AAt,c("a","b")))
checker(power.tensor(AAt,c("a","b"),c("a'","b'"),2),
        mul.tensor(AAt,c("a","b"),AAt,c("a'","b'")))

checker(power.tensor(power.tensor(AAt,c("a","b"),c("a'","b'"),1/pi),
                     c("a","b"),c("a'","b'"),pi),AAt)


AAt <- einstein.tensor(A , mark(A,"'",c("a","b")),by="e")

checker(power.tensor(AAt,c("a","b"),c("a'","b'"),-1,by="e"),
        inv.tensor(AAt,c("a","b"),by="e"))
checker(power.tensor(AAt,c("a","b"),c("a'","b'"),2,by="e"),
        mul.tensor(AAt,c("a","b"),AAt,c("a'","b'"),by="e"))

checker(power.tensor(power.tensor(AAt,c("a","b"),c("a'","b'"),1/pi,by="e"),
                     c("a","b"),c("a'","b'"),pi,by="e"),AAt)


# to.matrix.tensor                               # 
A <- reorder.tensor(to.tensor(1:30,c(a=2,b=3,c=5)),c("c","a","b"))

checker(to.matrix.tensor(A,"a",c("b","c")),matrix(1:30,nrow=2))

checker(to.matrix.tensor(A,c("a","b"),c("c")),matrix(1:30,nrow=6))

checker(to.matrix.tensor(A,c("a","b"),by=c("c")),structure(1:30,dim=c(6,1,5)))
checker(to.matrix.tensor(A,c("a"),by=c("c")),structure(1:30,dim=c(2,3,5)))

# untensor

A <- reorder.tensor(to.tensor(1:30,c(a=2,b=3,c=5)),c("c","a","b"))

checker(untensor(A,c("a","b"),pos=1),to.tensor(1:30,c(I1=6,c=5)))
checker(untensor(A,c("a","b"),"new",pos=2),reorder(to.tensor(1:30,c(new=6,c=5)),2:1))
checker(untensor(A,list(u=c("a","b"),v=c(c="c"))),to.tensor(1:30,c(u=6,v=5)))


# as.tensor

checker(as.tensor(diag(5)),to.tensor(c(diag(5)),c(I1=5,I2=5)))

# Slice tensor

A <- reorder.tensor(to.tensor(1:30,c(a=2,b=3,c=5)),c("c","a","b"))
checker(slice.tensor(A,"c",1:2),to.tensor(1:12,c(a=2,b=3,c=2)))
checker(slice.tensor(A,"c",1,drop=TRUE),to.tensor(1:6,c(a=2,b=3)))
checker(slice.tensor(A,"c",1,drop=FALSE),to.tensor(1:6,c(a=2,b=3,c=1)))

# Indexing with [[]]

A <- reorder.tensor(to.tensor(1:30,c(a=2,b=3,c=5)),c("c","a","b"))
checker(A[[b=2]],slice.tensor(A,"b",2,drop=TRUE))
checker(A[[b=2:3]],slice.tensor(A,"b",2:3))
checker(A[[b=2:3,c=3:4]],slice.tensor(slice.tensor(A,"b",2:3),"c",3:4))
checker(A[[b=~q]],to.tensor(1:30,c(a=2,q=3,c=5)))

#undrop.tensor

checker(undrop.tensor(slice.tensor(A,"c",2,drop=TRUE),"c"),
        slice.tensor(A,"c",2,drop=FALSE))

# bind.tensor

A <- to.tensor(1:6,c(a=2,b=3))
checker( bind.tensor(A,"a",A), to.tensor(c(1,2,1,2,3,4,3,4,5,6,5,6),c(a=4,b=3)))
checker( bind.tensor(A,"b",A), to.tensor(c(1:6,1:6),c(a=2,b=6)))




# einstein.tensor
A <- to.tensor(1:6,c(a=2,b=3))
checker( einstein.tensor(A,A), mul.tensor(A,c("a","b"),A))
checker( einstein.tensor(A,A,by="b"), mul.tensor(A,c("a"),A,by="b"))
checker( einstein.tensor(A,diag=A,by="b"), diagmul.tensor(A,c("a"),A,by="b"))


# adding
A <- to.tensor(1:30,c(a=2,b=3,c=5))

checker(A + reorder(A,c("b","a","c")),2*A)
checker(A - reorder(A,c("b","a","c")),to.tensor(rep(0,30),c(a=2,b=3,c=5)))

B <- to.tensor(rnorm(6),c(b=3,a=2))

checker( A + B , rep.tensor(B,5,name="c") + A ) 
checker( B + A , rep.tensor(B,5,name="c") + A ) 
checker( A + (B - A),  rep.tensor(B,5,name="c") ) 
checker( B + (A - B),  A ) 

C <- to.tensor(1:42,c(a=2,b=3,d=7))

checker( A+C , rep.tensor(A,7,name="d") + rep.tensor(C,5,name="c")) 

# Einstein 

checker( C %e% A , einstein.tensor(C,A))

# Riemann + drag
checker( C %r% mark(A), mul.tensor(C,c(),mark(A),c()))

A <- to.tensor(1:16,c(a=2,b=2,c=2,d=2))
gij <- to.tensor(c(1,0.5,0.5,1),c(i=2,j=2))

ginv <- c(solve(matrix(c(1,0.5,0.5,1),2)))
gij1 <- to.tensor(ginv,c("^a"=2,a=2))
gij2 <- to.tensor(ginv ,c("^b"=2,b=2))
is.covariate(gij)

checker( drag.tensor(gij,gij,c("i")) , to.tensor(c(1,0,0,1),c("^i"=2,"j"=2)))
checker( drag.tensor(gij,gij,c("j")) , to.tensor(c(1,0,0,1),c("i"=2,"^j"=2)))
checker( drag.tensor(gij,gij,c("i","j")), to.tensor(ginv,c("^i"=2,"^j"=2)))
checker( drag.tensor(gij,gij,c("i","j"))[["^i"=~a,"^j"=~b]], riemann.tensor(gij[[i=~a,j=~b]],gij1,gij2))

checker( drag.tensor(A,gij,c("a","b")), einstein.tensor(A,gij1,gij2) )
###########################################################################
### Names






##########################################################################

# 
I <- to.tensor(diag(3),c(a=3,b=3),nn(a=3,b=3),what=1:2)
checker(to.tensor(I),I)
checker(as.tensor(I),I)
checker(inv.tensor(I,"a","b"),I)

R1  <- matrix(rnorm(9),nrow=3)
R1i <- solve(R1)
R2 <- to.tensor(R1,c(a=3,b=3),nn(a=3,b=3),what=1:2)
R2i <- to.tensor(R1i,c(b=3,a=3),nn(b=3,a=3),what=1:2)

checker(inv.tensor(R2,"a","b"),R2i)
checker(inv.tensor(R2,"a","b",allowSingular=TRUE),R2i)


R3 <- to.tensor(rnorm(15),c(a=3,z=5),nn(a=3,z=5))

checker(mul.tensor(R2i,"b",mul.tensor(R2,"a",R3)),R3)

checker(solve.tensor(R2i,R3[[z=1]],"a"),mul.tensor(R2,"a",R3[[z=1]]))
checker(solve.tensor(R2i,R3,"a"),mul.tensor(R2,"a",R3))

checker(solve.tensor(R2i,R3[[z=1]],"a",allowSingular=TRUE),mul.tensor(R2,"a",R3[[z=1]]))
checker(solve.tensor(R2i,R3,"a",allowSingular=TRUE),mul.tensor(R2,"a",R3))

checker(solve.tensor(rep(R2i,4,1,"K"),R3[[z=1]],"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(rep(R2i,4,1,"K"),rep(R3[[z=1]],4,1,"K"),"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(R2i,rep(R3[[z=1]],4,1,"K"),"a",by="K"),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))

checker(solve.tensor(rep(R2i,4,1,"K"),R3[[z=1]],"a",by="K",allowSingular=TRUE),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(rep(R2i,4,1,"K"),rep(R3[[z=1]],4,1,"K"),"a",by="K",allowSingular=TRUE),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))
checker(solve.tensor(R2i,rep(R3[[z=1]],4,1,"K"),"a",by="K",allowSingular=TRUE),rep(mul.tensor(R2,"a",R3[[z=1]]),4,1,"K"))



checker(solve.tensor(R2i,R3,"a"),mul.tensor(R2,"a",R3))

checker(solve.tensor(R2i,R3[[z=1]],"a",allowSingular=TRUE),mul.tensor(R2,"a",R3[[z=1]]))
checker(solve.tensor(R2i,R3,"a",allowSingular=TRUE),mul.tensor(R2,"a",R3))



summary(I)


A <- to.tensor(c(diag(3)),c(a=3,b=3),nn(a=3,b=3))

A <- to.tensor(c(1,1,1,0,1,1,0,0,1),c(a=3,b=3),nn(a=3,b=3))
checker(mul.tensor(A,"b",A,"a"),to.tensor(c(A%*%A),c(a=3,b=3),nn(a=3,b=3)))

A <- to.tensor(c(1,1,1,1,0,1,1,1,0,0,1,1),c(a=4,b=3),nn(a=4,b=3))
checker(mul.tensor(A,"b",A[[a=~c]],"b"),
        to.tensor(c(A%*%t(A)),dim=c(a=4,c=4),nn(a=4,a=4)))

A <- to.tensor(rnorm(15),c(a=5,b=3),nn(a=5,b=3))
checker(mul.tensor(A,"b",A[[a=~c]],"b"),
        to.tensor(c(A%*%t(A)),dim=c(a=5,c=5),nn(a=5,a=5)))

# 
A <- to.tensor(c(1,1,1,0,1,1,0,0,1),c(a=3,b=3))

solve.tensor(mul.tensor(I,"b",I,"a"),I,"a","a")
checker(solve.tensor(mul.tensor(I,"b",I,"a"),I,"a","b"),I)
mul.tensor(A,"b",I,"a")

# + - * /
A <- to.tensor(1:4,c(i=4))
B <- to.tensor(1:4,c(j=4))

checker((A+B-A-B),0*(A %e% B))
checker((A-B-A+B),0*(A %e% B))
checker((A+B-B-A),0*(A %e% B))
checker((A-B+B-A),0*(A %e% B))
checker((A-A+B-B),0*(A %e% B))
checker((A*B/A/B),one.tensor(c(dim(A),dim(B))))
checker((A/B/A*B),one.tensor(c(dim(A),dim(B))))
checker((A*B/B/A),one.tensor(c(dim(A),dim(B))))
checker((A/B*B/A),one.tensor(c(dim(A),dim(B))))
checker((A/A*B/B),one.tensor(c(dim(A),dim(B))))

# $,|,^

A <- to.tensor(1:6,c(a=1,b=2,c=3))
checker( names(A$ijk) , c("i","j","k") )
checker( names(A$i.j.k) , c("i","j","k") )
checker( A$ijk^c("a","b","c") , A )
checker( A$ijk^"a.b.c" , A )
checker( A$ijk^"$abc" , A )
checker( names(A|"$bca"), c("b","c","a"))
checker( names(A|"b.c.a"), c("b","c","a"))
checker( names(A|"$b.c.a"), c("b","c","a"))
checker( names(A|c("$b","c.a")), c("b","c","a"))
checker( names(A|c("$b","c","a")), c("b","c","a"))
checker( names(A|c("$","b","c","a")), c("b","c","a"))
checker( names(A|c("b.c.a")), c("b","c","a"))


# slice.tensor<-
#A <- to.tensor(1:24,c(a=2,b=3,c=4))
#B <- to.tensor(25:48,c(a=2,b=3,c=4))
#slice.tensor(A,"c",1:2)<- slice.tensor(B,"c",1:2)
#checker(slice.tensor(A,"c",1:2),slice.tensor(B,"c",1:2))
#slice.tensor(A,"b",2)<- slice.tensor(B,"b",2)
#checker(slice.tensor(A,"b",2),slice.tensor(B,"b",2))
