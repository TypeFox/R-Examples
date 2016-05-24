require(tensorA)

if(FALSE) {
  # Commands for testing
  debugger()
  options(error=dump.frames)
  detach("package:tensorA")
  library(tensorA,lib.loc="../../tensorA.Rcheck")
}

#Aidx(3,4) # c(1,2,3,1,2,3,1,2,3,1,2,3)
#Bidx(3,4) # c(1,1,1,2,2,2,3,3,3,4,4,4)

# gsi.eps
#gsi.eps  # 1E-10
# gs.setarg
#tmp <- function(a=2,b=3) {a*b}
#if( gsi.setarg(tmp,b=5)()!=10 )
#  stop("Fehler gs.setarg")
  
# to.tensor
to.tensor(c(1,2,3))
dim(to.tensor(c(1,2,3)))

set.seed(23)

A <- to.tensor(1:20,c(U=2,V=2,W=5))
A
dim(A)
names(A)
dimnames(A)

ftable(to.tensor(A))
ftable(to.tensor(c(A),dim(A)))
ftable(to.tensor(c(A),dim(A),dimnames(A)))
ftable(to.tensor(A,dim(A),what=1:3))
ftable(to.tensor(A,dim(A)[1:2],dimnames(A)[1:2],1:2))
ftable(to.tensor(A,dim(A)[1],dimnames(A)[1],1,addIndex=TRUE))


Anamed <- A
#dimnames(Anamed)[["U"]]<- gsi.stdnames(2,"u")
#dimnames(Anamed)[["V"]]<- gsi.stdnames(2,"v")
#dimnames(Anamed)[["U"]]<- gsi.stdnames(dim(A)["w"],"w")
ftable(Anamed)
ftable(to.tensor(Anamed))
ftable(to.tensor(c(Anamed),dim(Anamed)))
ftable(to.tensor(c(Anamed),dim(Anamed),dimnames(Anamed)))
ftable(to.tensor(Anamed,dim(Anamed),what=1:3))
ftable(to.tensor(Anamed,dim(Anamed)[1:2],dimnames(Anamed)[1:2],1:2))
ftable(to.tensor(Anamed,dim(Anamed)[1],dimnames(Anamed)[1],1,addIndex=TRUE))
ftable(to.tensor(Anamed,dim(Anamed)[1],dimnames(Anamed)[1],1,addIndex="I"))


B <- to.tensor(1:30,list(U=c("a","b","c"),V=c("B1","B2"),W=1:5))
B
dim(B)
names(B)
dimnames(B)

B2 <- to.tensor(B,c(x=1,y=3,z=1))
ftable(B2)


C <- to.tensor(1:20,c(A=4,B=5))
C
C <- to.tensor(C,c(A1=2,A2=2))
C
D <- C
C <- to.tensor(C,c(A1=1,A2=4),what=c("A2","A1"))
C

C <- to.tensor(C,c(A1=2),what=c("A2","A1"),addIndex="Q")
C

Cnamed <- to.tensor(1:20,c(A=4,B=5))
dimnames(Cnamed)[["B"]]<-LETTERS[1:5]
dimnames(Cnamed)[["A"]]<-letters[1:4]
Cnamed
Cnamed <- to.tensor(Cnamed,c(A1=2,A2=2))
dimnames(Cnamed)
Dnamed <- Cnamed
Cnamed <- to.tensor(Cnamed,c(A1=1,A2=4),what=c("A2","A1"))
ftable(Cnamed)

Cnamed <- to.tensor(Cnamed,c(A1=2),what=c("A2","A1"),addIndex="Q")
ftable(Cnamed)



# names
names( C )
names( C ) <- c("A1","A2","A3")
C
names(C)

names( Cnamed )
names( Cnamed ) <- c("A1","A2","A3")
ftable(Cnamed)
names(Cnamed)


# norm
norm.tensor(C) - sqrt(sum((1:20)^2))
norm.tensor(Cnamed) - sqrt(sum((1:20)^2))

# margin.tensor
margin.tensor(C,"A1")
margin.tensor(Cnamed,"A1")

# diagmul.tensor
norm.tensor( diagmul.tensor(C,D=C)-C^2 ) # 0
norm.tensor( diagmul.tensor(Cnamed,D=Cnamed)-Cnamed^2 ) # 0

Cr <- rep.tensor(C,10,1,name="K")
#norm(diagmul.tensor(Cr,names(C),C,names(C),by="K") - diagmul.tensor(Cr,names(C),D=Cr,names(C),by="K")) # 0

# pos.tensor
pos.tensor(dim(C))

# reorder.tensor
reorder.tensor(C,c(2))
reorder.tensor(C,c("A3","A2"))
ftable(reorder.tensor(Cnamed,c(2)))
ftable(reorder.tensor(Cnamed,c("A3","A2")))

# mul.tensor
AA <- mul.tensor(A,1,mark(A),1)
AA
dim(AA)

AA <- mul.tensor(Anamed,1,mark(Anamed),1)
ftable(AA)
dim(AA)

A <- to.tensor(as.complex(1:20),c(U=2,V=2,W=5))
AA <- mul.tensor(A,1,mark(A),1)


#


# rep.tensor
rep(A,3,name="Q",pos=4)
rep(Anamed,3,name="Q",pos=4)
ftable(rep(Anamed,3,name="Q",pos=4))

# trace.tensor

trace.tensor(D,"A1","A2")
dim(trace.tensor(D,"A1","A2"))

trace.tensor(Dnamed,"A1","A2")
dim(trace.tensor(Dnamed,"A1","A2"))


#delta.tensor

delta.tensor(c(a=2))
delta.tensor(c(a=2,b=3))
dim(delta.tensor(c(a=2,b=4),"."))

ftable(delta.tensor(c(a=2,b=2)))
ftable(delta.tensor(c(a=2,b=3)))

# tripledelta.tensor

tripledelta.tensor(c(a=2))
ftable(tripledelta.tensor(c(a=2)))
dim(tripledelta.tensor(c(a=2),"1","2"))


# mark.tensor
A
mark(A,"*")
names(mark(A,"*"))
mark(A,"*","U")
names(mark(A,"*","U"))
names(mark(A,"*",2:3))

Anamed
mark(Anamed,"*")
names(mark(Anamed,"*"))
mark(Anamed,"*","U")
names(mark(Anamed,"*","U"))
names(mark(Anamed,"*",2:3))


# mark.numeric
mark(dim(A),"*")
mark(dim(A),"*","U")
mark(dim(A),"*",2:3)

# mark.character
mark(names(A),"*")
mark(names(A),"*","U")
mark(names(A),"*",2:3)


# inv.tensor
E <- to.tensor(rnorm(16),c(a=2,b=2,A=2,B=2))
Ei <- inv.tensor(E,c("A","B"))
E
dim(Ei)
names(Ei)
dim(mul.tensor(E,c("A","B"),mark(Ei,c("a","b"))))
dim(delta.tensor(c(a=2,b=2)))
norm.tensor(mul.tensor(E,c("A","B"),mark(Ei,c("a","b")))-delta.tensor(c(a=2,b=2)))

AA <- to.tensor(rnorm(20),c(a=2,b=2,C=5))
inv.tensor(AA,"a",by="C")



E <- to.tensor(rnorm(16),c(a=2,b=2,A=2,B=2),list(c("a","b"),NULL,c("A1","A2"),c("B1","B2")))
Ei <- inv.tensor(E,c("A","B"))
E
dim(Ei)
names(Ei)
norm.tensor(mul.tensor(E,c("A","B"),mark(Ei,c("a","b")))-delta.tensor(c(a=2,b=2)))
ftable(round(mul.tensor(E,c("A","B"),mark(Ei,c("a","b")))),3)

# solve.tensor 

X1 <- to.tensor(
                c( 0 , 10 , 20 ,30 , 2,2,2,2, 45,21,34,5, 67,0,0,0 ),
                c(a=2,b=2,"a'"=2,"b'"=2) )
b1 <- to.tensor(c( 0,10,20,30),c(a=2,b=2))
solve.tensor(X1,b1,c("a","b"))
norm.tensor(solve.tensor(X1,b1,c("a","b"))  - to.tensor(c( 1,0,0,0 ),c(a=2,b=2))) #0

dimnames(X1)[["a"]]<- 1:2
dimnames(X1)[["b"]]<- 1:2
dimnames(X1)[["a'"]]<- 1:2
dimnames(b1)[["a"]] <-c("M","K")

solve.tensor(X1,b1,c("a","b"))
norm.tensor(solve.tensor(X1,b1,c("a","b"))  - to.tensor(c( 1,0,0,0 ),c(a=2,b=2))) #0


# chol.tensor
Xs <- einstein.tensor(X1,a="a*",b="b*",X1)
dim(Xs)
norm.tensor( Xs - mul.tensor(X1,c("a'","b'"),mark(X1,"*",c("a","b")),c("a'","b'"))) #0
RXs <- chol.tensor(Xs,c("a","b"),c("a*","b*"))
norm.tensor(mul.tensor(RXs,1,mark(RXs),1) - Xs ) # 0

dimnames(X1)[["a"]]<-c("q1","q2")
dimnames(X1)[["a'"]]<-c("p1","p2")

Xs <- einstein.tensor(X1,a="a*",b="b*",X1)
dimnames(Xs)
dimnames(Xs)["a*"]<-list(NULL)
dim(Xs)
norm.tensor( Xs - mul.tensor(X1,c("a'","b'"),mark(X1,"*",c("a","b")),c("a'","b'"))) #0
RXs <- chol.tensor(Xs,c("a","b"),c("a*","b*"))
norm.tensor(mul.tensor(RXs,1,mark(RXs),1) - Xs ) # 0


# level.tensor
dim(X1)
level.tensor(X1) # 4

# svd.tensor
SXs <- svd.tensor(X1,c("a","b"))
SXs
DXs <- to.tensor(c(diag(SXs$d)),c(lambda=length(SXs$d),"lambda'"=length(SXs$d)))
SXs2<-mul.tensor(mul.tensor(SXs$u,"lambda",DXs,"lambda'"),"lambda",SXs$v) #
SXs2
norm.tensor(SXs2-X1) #0

# power.tensor
Xs
rXs <- power.tensor(Xs,c(1,2),c(3,4),p=0.5)
Xs2 <- mul.tensor(rXs,c("a","b"),rXs,c("a*","b*"))
Xs2
norm.tensor((Xs-Xs2)) # 0

# to.matrix.tensor
to.matrix.tensor(X1,j=c("a'","b'"))
to.matrix.tensor(X1,i=c("a","b")  )

# gsi.stdnames
#gsi.stdnames(7,avoid=gsi.stdnames(5))
#gsi.stdnames(3,"B")
#gsi.stdnames(3,"B",avoid=gsi.stdnames(5))


# gsi.namedlist
#gsi.namedlist(c("Holla","Hoppla"),1,3)
# gsi.lefts/gsi.rights
#gsi.lefts(X1)
#gsi.rights(X1)

# untensor
tmp <- untensor(X1,c("a","b"),"A",1)
X1 <- to.tensor(
                c( 0 , 10 , 20 ,30 , 2,2,2,2, 45,21,34,5, 67,0,0,0 ),
                c(a=2,b=2,"a'"=2,"b'"=2) )
tmp <- untensor(X1,c("a","b"),"A",1)
tmp
untensor(tmp,c("a'","b'"),"B",2)

# gsi.untensornames
untensor(A,c(1,2))
untensor(Anamed,c(1,2))

# as.tensor
as.tensor(X1)
as.tensor(c(1,2,3))
dim(as.tensor(c(1,2,3)))

# renorm.matrix / renorm.tensor
#renorm.rows(X1[,,1,1])
#renorm.tensor(X1[,,1,1])
#norm.tensor(renorm.tensor(X1,c(1,2)),c(1,2)) # 1 1 1 1

# slice.tensor
slice.tensor(B,"U","a")
slice.tensor(B,"U","a",drop=TRUE)
slice.tensor(B,"U",c("a","c"))

# [[]].tensor
B[[U="a"]]
B[[U="a",drop=FALSE]]
old <- B[[U=c("a","c")]]
B[1:2,,]
B[[1]]
B[[U=c("a","c")]]<-7
B
B[[U=c("a","c")]]<- old

# undrop.tensor
undrop.tensor(A,name="Q",2)

# combineCF.tensor
# combineCF.tensor(Xs,c("a","b"),Xs,c("a","b"))

# bind.tensor
bind.tensor(A,1,einstein.tensor(A,V="q","U"="V",q="U"),2)
bind.tensor(A,1,A,1)
bind.tensor(Anamed,1,einstein.tensor(Anamed,V="q","U"="V",q="U"),2)
bind.tensor(Anamed,1,Anamed,1)

# toPos.tensor
#toPos.tensor(A,c("A","U","I")) # error
toPos.tensor(A,rev(names(A)))

# einstein.tensor
einstein.tensor(A,U="U'",B)
einstein.tensor(A,U="U'",mark(B,"k"))
einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk")
einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk",1/10)
einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk",diag=to.tensor(c(1,1/10,1/100),c(Uk=3)))

ftable(einstein.tensor(A,U="U'",B))
ftable(einstein.tensor(A,U="U'",mark(B,"k")))
ftable(einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk"))
ftable(einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk",1/10))
ftable(einstein.tensor(A,U="U'",mark(B,"k"),V="Vk",W="Wk",diag=to.tensor(c(1,1/10,1/100),c(Uk=3))))

# %e% [[]]
dim(A[[U=~M]])
A[[U=~M]] %e% B
A[[U=~M,V=~"L"]] %e% B


# firstnames/secondnames
#firstnames(secondnames(names(A),"2"),"2")
contraname(names(A))

dim(A)
dim(B)

add.tensor(A,A)/2-A
norm.tensor(add.tensor(A,A)/2-A)

tmp <- add.tensor( A,A[[U=~K]] ,op="-")
dim(tmp)
ftable(tmp)

A %e% A
norm.tensor(reorder(A,c(2,3,1)) - A)

# Dragging 

g <- to.tensor(c(1,2,0,1),c(i=2,j=2))
A <- to.tensor(rnorm(8),c(a=2,b=2,c=2))
A2 <- drag.tensor(A,g,c("b","c"))
A2
names(A2)
as.covariate(names(A2))
as.contravariate(names(A2))
is.covariate(A2)
is.contravariate(A2)
riemann.tensor(A2,g)


#  #######################

 d1 <- c(a=2,b=2)
 d2 <- c("a'"=2,"b'"=2)
 m <- to.tensor(1:4,d1)             
 V <- delta.tensor(d1)+one.tensor(c(d1,d2))
 V
 X <- (power.tensor(V,c("a","b"),c("a'","b'"),p=1/2)  %e%
      to.tensor(rnorm(1000*2*2),c(i=1000,d2))) + m
 # The mean
 mean.tensor(X,along="i")
 # Full tensorial covariance:
 var.tensor(X,along="i")
 # Variance of the slices  X[[b=1]] and X[[b=2]] :
 var.tensor(X,along="i",by="b")
 # Covariance of the slices X[[b=1]] and X[[b=2]] :
 var.tensor(X[[b=1]],X[[a=~"a'",b=2]],along="i")


##  #########################

Var <- function(x,along) {
  one <- one.tensor(dim(x)[along])
  M   <- one/(one%e%one)
  z <- x-  M %e% x
  einstein.tensor( z,mark(z),by=names(M)) %e% M
  
}

A <- to.tensor(rnorm(4000),c(i=2,j=2,s=1000))

dim(A)["s"]
Var(A,"s")
delta.tensor(c(i=2,j=2))

Var(A,"s")-delta.tensor(c(i=2,j=2))



