library(partitions)

# a little function for equality:
minmax <- function(x){all.equal(min(x), max(x))}

# size of tests:
N <- 20

#first check that the parts sum correctly:

f <- function(n){
  if (length(n) > 1) {
    return(sapply(n, match.fun(sys.call()[[1]])))
  }

  m <- ceiling(n/2)
  
  minmax(apply(           parts(n  ) ,2,sum))          &
  minmax(apply(       diffparts(n  ) ,2,sum))          &
  minmax(apply( restrictedparts(n,m) ,2,sum))          &
  minmax(apply(conjugate(          parts(n  )),2,sum)) &
  minmax(apply(conjugate(      diffparts(n  )),2,sum)) &
  minmax(apply(conjugate(restrictedparts(n,m)),2,sum)) 
}

stopifnot(all(f(3:N)))          
  

# now check that parts() has rep(1,n) as the last partition:
g <- function(n){
  if (length(n) > 1) {
    return(sapply(n, match.fun(sys.call()[[1]])))
  }
  jj <- parts(n)
  minmax(c(1,jj[,ncol(jj)]))
}

stopifnot(all(g(1:N)))

# now some spot checks:

# Andrews, page232:
stopifnot(all.equal(R(5,12), 13))

# a couple of values taken from A&S, table 21.5, page 836:
stopifnot(all.equal(P(100), 190569292))
stopifnot(all.equal(Q(100), 444793))


# some checks of conjugation.  Because conjugation is a bijection, the
# conjugate of an (unrestricted) enumeration of partitions should be
# an enumeration of unrestricted partitions, possibly in a different
# order.  The rather longwinded compare() below checks this
# explicitly; see the online helppage for order() for an explanation
# of the use of do.call(); here it is used to provide an unambiguous
# ordering for the partitions.  Note also the need to coerce to a
# matrix so that we can set the dimnames to NULL to allow identical()
# to work as desired.


compare <- function(n){
   if (length(n) > 1) {
        return(sapply(n, match.fun(sys.call()[[1]])))
    }
  jj <- parts(n)

  a <- as.data.frame(as.matrix(t(jj)))
  b <- as.data.frame(as.matrix(t(conjugate(jj))))

  a <- as.matrix(a[do.call(order,a),])
  b <- as.matrix(b[do.call(order,b),])
  
  dimnames(a) <- NULL
  dimnames(b) <- NULL

  return(identical(a,b))
}

stopifnot(all(compare(1:15)))





# some tests of durfee() and conjugate(), from p28 of Andrews:
a <- c(7,7,5,4,4,2,1)
stopifnot(identical(conjugate(a),as.integer(c(7,6,5,5,3,2,2))))


# now verify that the conjugate of an unequal partition is of the form
# c(rep(p,a0),rep(p-1,a1), ... , rep(p-i,ai), ... , rep(1,ap)),
# that is, something like  4 3 2 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 ---
# this being  "conjugate(diffparts(20))[,30]".

 f <- function(x){ #note second and third tests are the same
   (min(x)  %in%  0:1)            &
   all(diff(x) %in% -1:0)         &   
   all(diff(x[cumsum(rle(x)$lengths)]) == -1)  
 }
   
 g <- function(n){all(apply(conjugate(diffparts(n)),2,f))}
                  
 stopifnot(c(g(10),g(11),g(20)))


# Now verify that S() is independent of the order of y:
jj <- S(rep(1:4,each=2),5)
stopifnot(jj == 474)
stopifnot(jj == S(rep(1:4,each=2),5))



# Now some tests on compositions():
jj <- compositions(7)
stopifnot(all(apply(jj,2,sum)==7))

# test that the bug has been corrected:
stopifnot(all(apply(compositions(15,4,TRUE ),2,sum)==15))
stopifnot(all(apply(compositions(15,4,FALSE),2,sum)==15))




# some tests of setparts:
f <- function(jj){all(apply(setparts(jj),2,table)==jj)}
stopifnot(f(c(3,2)))
stopifnot(f(c(3,1)))
stopifnot(f(c(3,2,1)))


# some tests of the comptobin(); not run because it needs the elliptic package:
f <- function(n){
  jj <- as.integer(mobius(seq_len(n))==0)
  return(
         identical(jj , comptobin(bintocomp(jj,TRUE ))) &
         identical(jj , comptobin(bintocomp(jj,FALSE)))
         )
}
if(FALSE){
  f(1000)
}
