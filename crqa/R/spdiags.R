### both .R and .f code written by Prof. John C. Nash  <nashjc@uottawa.ca>

.packageName <- 'crqa'

spdiags <- function(A) { # A is a matrix

    m <- dim(A)[[1]]
    n <- dim(A)[[2]]
    k <- min(m, n) # length of diagonals
    if (m > n) { # tall matrix
        Adata <- as.vector(t(A)) # convert to vector BY ROWS
    } else { # fat or square matrix (m <= n)
        Adata <- as.vector(A) # convert to vector BY COLUMNS 
   }
#   print(Adata)
    la<-length(Adata)
    na<-m*n+2*k*(k-1)
    Adata <- c(Adata, rep(0,(na-la)))
    nb<-(m+n)*k
    nd<-m+n-1
    jb<-0
    Bdata<-rep(0,nb)
    d<-rep(0,nd)
    tv<-rep(0,k)
    tres<-.Fortran("jspd",m=as.integer(m), n=as.integer(n), k=as.integer(k),
                   Adata=as.double(Adata), jb=as.integer(jb), 
                   Bdata=as.double(Bdata), d=as.integer(d), tv=as.double(tv),
                   na=as.integer(na), nb=as.integer(nb), nd=as.integer(nd) )
                                        #   print(str(tres))
    jb<-tres$jb
    d<-tres$d[1:jb]
    Bdata<-tres$Bdata[1:(jb*k)]
#   print(Bdata)
    if (m > n) d <- -d # reset index
    B <- matrix(Bdata, nrow=k, byrow=FALSE) # convert to matrix form
    result<-list(B=B, d=d)
}
