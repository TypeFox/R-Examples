
library(giRaph)

G1 <- new("incidenceList", E=list(u(1,2), d(1,3),
                             u(3), d(2,5),
                             d(2,5),d(3,c(1,4),5),u(2,4,5),
                             d(c(3,4),c(2,1)),r(1,5)),
                             V=5:10)

g <- new("anyGraph",incidenceList=G1)
gs<- as(g,"simpleGraph")

n <- 10

X.empty     <- new("adjacencyMatrix",matrix(0,n,n))
X.full      <- new("adjacencyMatrix",matrix(1,n,n))

set.seed(90)
X <- matrix(0,n,n)
for (i in 2:n) X[i,sample(1:n,rbinom(1,n,.4))] <- 1
diag(X) <- 0

diag(X.empty)<-0
diag(X.full)<-0

## X <- as(G1,"adjacencyMatrix")
X <- new("adjacencyMatrix", X)
A <- as(X,"adjacencyList") 
I <- as(X,"incidenceMatrix")
G <- as(X,"incidenceList")

s1 <- new("simpleGraph",adjacencyMatrix=X)
s2 <- new("simpleGraph",adjacencyList=A)
s3 <- new("simpleGraph",incidenceMatrix=I)
s4 <- new("simpleGraph",incidenceList=G)

### the following has been commented out and should be fixed

### change representation
# s5 <- s4
# adjacencyMatrix(s5) <- adjacencyMatrix(s5)

### work with 2 representations
#s3 <- "adjacencyMatrix<-"(s3,force=FALSE,adjacencyMatrix(s3))  
#s3 <- "adjacencyMatrix<-"(s3,force=FALSE,as(G1,"adjacencyMatrix")) 

### coercion
#g4 <- as(s4,"anyGraph")


#is(s4) # gives all superclasses
#is(s4,"generalGraph")

#areTheSame(s1,s3) 
#areTheSame(s1,s4)
#areTheSame(s3,s4)

#efficiency <- function(X) {
#  ## X: adjacencyMatrix
#  A <- as(X,"adjacencyList")
#  I <- as(X,"incidenceMatrix")
#  G <- as(X,"incidenceList")

#  res <- matrix(NA,4,4)
#  nm <- c("G","I","A","X")
#  dimnames(res) <- list(nm,nm)

#  res[1,1] <- object.size(G)
#  res[2,2] <- object.size(I)
#  res[3,3] <- object.size(A)
#  res[4,4] <- object.size(X)

#  cat("I->G:")
#  res[2,1] <- system.time( as(I,"incidenceList")  ,gcFirst=TRUE)[3]
#  cat(res[2,1],"\nI->A:")
#  res[2,3] <- system.time( as(I,"adjacencyList")  ,gcFirst=TRUE)[3]
#  cat(res[2,3],"\nI->X:")
#  res[2,4] <- system.time( as(I,"adjacencyMatrix")  ,gcFirst=TRUE)[3]

#  cat(res[2,4],"\nX->G:")
#  res[4,1] <- system.time( as(X,"incidenceList")  ,gcFirst=TRUE)[3]
#  cat(res[4,1],"\nX->I:")
#  res[4,2] <- system.time( as(X,"incidenceMatrix")  ,gcFirst=TRUE)[3]#
#  cat(res[4,2],"\nX->A:")
#  gc()
#  res[4,3] <- system.time( as(X,"adjacencyList")  ,gcFirst=TRUE)[3]

#  cat(res[4,3],"\nG->I:")
#  res[1,2] <- system.time( as(G,"incidenceMatrix")  ,gcFirst=TRUE)[3]#
#  cat(res[1,2],"\nG->A:")
#  gc()
#  res[1,3] <- system.time( as(G,"adjacencyList")  ,gcFirst=TRUE)[3]
#  cat(res[1,3],"\nG->X:")
#  res[1,4] <- system.time( as(G,"adjacencyMatrix")  ,gcFirst=TRUE)[3]

#  cat(res[1,4],"\nA->G:")
#  res[3,1] <- system.time( as(A,"incidenceList")  ,gcFirst=TRUE)[3]
#  cat(res[3,1],"\nA->I:")
#  res[3,2] <- system.time( as(A,"incidenceMatrix")  ,gcFirst=TRUE)[3]
#  cat(res[3,2],"\nA->X:")
#  res[3,4] <- system.time( as(A,"adjacencyMatrix")  ,gcFirst=TRUE)[3]
#  cat(res[3,4],"\n")


#  cat( ncol(I), "Vertices, ", nrow(I)," Edges\n")
#  res  
#}

#efficiency(X.full)
#efficiency(X)
#efficiency(X.empty)

#S <- new("simpleGraph",adjacencyMatrix=X)
#M <- as(S,"multiGraph")
#H <- as(S,"generalGraph")
#G <- as(S,"anyGraph")

##system.time( as(M,"generalGraph")  ,gcFirst=TRUE)
##system.time( as(M,"anyGraph")  ,gcFirst=TRUE)
##system.time( as(H,"anyGraph")  ,gcFirst=TRUE)
