library(partitions)

# Function f() tests the various next...() functions by creating a
# matrix by successively cbind()ing the next() partition; then
# comparing this matrix with the one produced by the vectorized
# function.  Note that the test matrices are produced by repeated
# calls to cbind(): not very efficient, but OTOH all three functions
# are tested [first...() , next...(), and islast...()].

f <- function(n1,n2,n3,n4,n5,n6,n7){
  a <- firstpart(n1)
  out <- cbind(a)
  while(!islastpart(a)){
    a <- nextpart(a)
    out <- cbind(out,a)
  }
  stopifnot(all(out == parts(n1)))
  
  a <- firstdiffpart(n1)
  out <- cbind(a)
  while(!islastdiffpart(a)){
    a <- nextdiffpart(a)
    out <- cbind(out,a)
  }
  stopifnot(all(out == diffparts(n1)))
  
  a <- firstrestrictedpart(n2,n3,FALSE)
  out <- cbind(a)
  while(!islastrestrictedpart(a)){
    a <- nextrestrictedpart(a)
    out <- cbind(out,a)
  }
  stopifnot(all(out == restrictedparts(n2,n3,FALSE)))
  
  a <- firstrestrictedpart(n2,n3,TRUE)
  out <- cbind(a)
  while(!islastrestrictedpart(a)){
    a <- nextrestrictedpart(a)
    out <- cbind(out,a)
  }
  stopifnot(all(out == restrictedparts(n2,n3,TRUE)))
  
  
  a <- firstblockpart(f=1:n4 , n=n5 , include.fewer=TRUE)
  out <- cbind(a)
  while(!islastblockpart(a, f=1:n4, n=n5, TRUE)){
    a <- nextblockpart(a, f=1:n4 , n=n5, TRUE)
    out <- cbind(out,a)
  }
  stopifnot(all(out == blockparts(1:n4, n5, TRUE)))



  a <- firstblockpart(f=1:n4 , n=n5 , include.fewer=FALSE)
  out <- cbind(a)
  while(!islastblockpart(a, f=1:n4 , FALSE)){
    a <- nextblockpart(a,f=1:n4, n=n5,FALSE)
    out <- cbind(out,a)
  }
  stopifnot(all(out == blockparts(1:n4, n5, FALSE)))


  
  a <- firstblockpart(f=1:n4 , n=NULL)
  out <- cbind(a)
  while(!islastblockpart(a, f=1:n4 , n=NULL)){
    a <- nextblockpart(a, f=1:n4 , n=NULL)
    out <- cbind(out,a)
  }
  stopifnot(all(out == blockparts(f=1:n4,n=NULL)))



  a <- firstcomposition(n6 , m=NULL , include.zero=TRUE)
  out <- cbind(a)
  while(!islastcomposition(a , FALSE , include.zero=TRUE)){
    a <- nextcomposition(a , FALSE , include.zero=TRUE)
    a <- c(a,rep(0,n6-length(a)))
    out <- cbind(out,a)
  }

  stopifnot(all(out == compositions(n6 , m=NULL, include.zero=TRUE)))


  a <- firstcomposition(n6 , m=NULL , include.zero=FALSE)
  out <- cbind(a)
  while(!islastcomposition(a , FALSE , include.zero=FALSE)){
    a <- nextcomposition(a , FALSE , include.zero=FALSE)
    a <- c(a,rep(0,n6-length(a)))
    out <- cbind(out,a)
  }
  stopifnot(all(out == compositions(n6 , m=NULL, include.zero=FALSE)))


    
  a <- firstcomposition(n6 , m=n7 , include.zero=FALSE)
  out <- cbind(a)
  while(!islastcomposition(a , TRUE, include.zero=FALSE)){
    a <- nextcomposition(a , TRUE , include.zero=FALSE)
    out <- cbind(out,a)
  }
  stopifnot(all(out == compositions(n6 , m=n7, include.zero=FALSE)))


  a <- firstcomposition(n6 , m=n7 , include.zero=TRUE)
  out <- cbind(a)
  while(!islastcomposition(a , TRUE , include.zero=TRUE)){
    a <- nextcomposition(a , TRUE , include.zero=TRUE)
    out <- cbind(out,a)
  }
  stopifnot(all(out == compositions(n6 , m=n7, include.zero=TRUE)))

}


f(7, 11, 3, 4, 4, 5, 4)
f(8, 10, 4, 4, 3, 5, 5)
