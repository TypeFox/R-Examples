library(magic)
n <- 10


# first check minmax (NB: includes checking all zeros):
stopifnot(all(sapply(0:10,function(i){minmax(rep(i,9))})))


# test magic() for magicness, standardness, and normality for magic(3)..magic(n):
stopifnot(is.magic   (magic(3:n)))
stopifnot(is.standard(magic(3:n)))
stopifnot(is.normal  (magic(3:n)))

# now test some of the specific algorithms:
stopifnot(is.magic(strachey(1:n)))
stopifnot(is.magic(lozenge (1:n)))
stopifnot(is.magic(magic.4n(1:n)))

stopifnot(sapply(1:n,function(i){is.square.palindromic(circulant(i))}))


# now test for magic.2np1() giving a generalized circulant:
stopifnot(sapply(1:n,function(i){is.circulant(magic.2np1(i)%%(2*i+1),c(2,-1))}))


# Now test that is.diagonally.correct() in fact extracts the correct elements,
# using a function that returns just the first element:
test.corners <- function(i){
  a <- magic(i)
  identical(a[c(1,i),c(1,i)],is.diagonally.correct(a,func=function(x){x[1]},g=TRUE)$diag.sums)
}

stopifnot(all(sapply(3:n,test.corners)))


# Now check that, in a 3x3x3 magic cube, the second element of each diagonal is the same:
f <- function(x){x[2]}
stopifnot(is.diagonally.correct(magiccube.2np1(1),func=f,boolean=FALSE,give=FALSE))


# Now check that the first eigenvalue of a magic square is indeed
# equal to its magic constant.

# First, define a wrapper to ensure that eigen() returns an integer:

eigen.wrap <- function(M){as.integer(round(Re(eigen(M,FALSE,TRUE)$values)))}

f <- function(i){minmax(c(eigen.wrap(magic(i))[1],magic.constant(i)))}
stopifnot(sapply(3:n,f))

# Now check that the sum of eigenvalues 2,...,n of a magic square is zero:
f <- function(i){minmax(c(1,1+sum(eigen.wrap(magic(i))[-1])))}
stopifnot(sapply(3:n,f))


# Check hudson() for 6n+1 and 6n-1:
stopifnot(sapply(c(6*(1:n)+1,6*(1:n)-1),function(i){is.magic(hudson(i))}))

# Check magichypercube.4n() for a range of dimensions and orders:
stopifnot(apply(expand.grid(m=1:2,d=2:4),1,function(x){is.magichypercube(magichypercube.4n(x[1],x[2]))}))

## Check magiccube.2np1():
stopifnot(sapply(1:n,function(i){is.magichypercube(magiccube.2np1( i))}))

## Sundry tests for transf;

## is transf(a,0) == a?
stopifnot(sapply(3:n , function(i){a <- magic(i);identical(a,transf(a,0))}))


## NB: following two tests *removed* following redefinition of
## "equal", that is eq(), or %eq%.  The _old_ definition was to put
## the square into Frenicle standard form, then compare.  The _new_
## definition considers the square directly.

## is transf(a,X) equal (ie eq()) to "a" for different X?
#stopifnot(sapply(3:n , function(i){a <- magic(i);eq(a,transf(a,i%%8  ))}))
#stopifnot(sapply(3:n , function(i){a <- magic(i);eq(a,transf(a,i%%8+1))}))


data(magiccubes)
stopifnot(unlist(lapply(magiccubes,is.magichypercube)))

data(Ollerenshaw)
stopifnot(is.mostperfect(Ollerenshaw))

data(cube2)
stopifnot(is.magichypercube(cube2))

data(hendricks)
stopifnot(is.perfect(hendricks))

data(perfectcube5)
stopifnot(is.perfect(perfectcube5))

data(perfectcube6)
stopifnot(is.perfect(perfectcube6))

data(Frankenstein)
stopifnot(is.perfect(Frankenstein))

# Comment out the line below because it takes too long
#stopifnot(apply(magic.8(),3,is.magic))

## Now check magic.product() works:
f <- function(x){is.magic(magic.product(x[1],x[2]))}
stopifnot(apply(expand.grid(3:5,3:5),1,f))


## Now check some identities for adiag():
a <- matrix(1:6,2,3)
a2 <- matrix(1,2,2)
a3 <- matrix(1,3,3)

x <- 0
dim(x) <- rep(1,7)

stopifnot(identical(dim(adiag(x,x,x)),rep(3:3,7)))
stopifnot(identical(adiag(a,t(a)),t(adiag(t(a),a))))
stopifnot(identical(adiag(1,1,1,1,1),diag(5)))
stopifnot(identical(adiag(a2,a2),kronecker(diag(2),a2)))
stopifnot(identical(adiag(a3,a3,a3),kronecker(diag(3),a3)))
stopifnot(identical(adiag(matrix(1,0,5),matrix(1,5,0),pad=1:5), kronecker(t(rep(1,5)),1:5)))

# Now some more tests.
# First, set the dimension.  Feel free to change this!
n <- 6 

# Check that pad value is correctly used:
a <- array(43,rep(2,n))
stopifnot(minmax(adiag(a,43,pad=43)))

# Check that adiag() plays nicely with subsums():
a <- array(1,rep(1,n))
stopifnot(minmax(subsums(adiag(a,a),2)))

# And another test:
a <- array(1,rep(1,n))
jj1 <- subsums(adiag(a,0,a),2,wrap=F)
x <- array(1,rep(2,n))
jj2 <- adiag(x,x)
stopifnot(identical(jj1,jj2))



# Now test adiag() for associativity:

jj1 <- array(seq_len(2^n),rep(2,n))
jj2 <- array(seq_len(3^n),rep(3,n))
jj3 <- array(seq_len(4^n),rep(4,n))

f <- function(x,y,z){stopifnot(identical(adiag(adiag(x,y),z),adiag(x,adiag(y,z))))} 
f(jj1,jj2,jj3)
f(jj2,jj3,jj1)
f(jj1,jj1,jj1)
f(jj3,jj3,jj3)



# Now some tests for is.circulant():
a <- array(0,rep(2,10))
a[1] <- a[1024] <- 1
stopifnot(is.circulant(a))

# "break" a by changing just one (randomly chosen) element:
a[1,1,1,1,2,1,2,1,1,1] <- 1
stopifnot(!is.circulant(a))

# Now test arev() with some tests:
a <- array(1:32,rep(2,5))
stopifnot(identical(as.vector(arev(a)),rev(a)))


jj <- as.vector(magic(19))[seq_len(360)]
stopifnot(identical(arev(array(jj,3:6)) , array(rev(jj),3:6)))


b <-  c(TRUE,FALSE,TRUE,FALSE,TRUE)
stopifnot(identical(a,arev(arev(a,b),b)))


stopifnot(identical(a[,2:1,,,],arev(a,2)))
stopifnot(identical(arev(a,c(2:4)),a[,2:1,2:1,2:1,]))


# now some tests of arot():
stopifnot(identical(arot(arot(a)),arot(a,2)))
stopifnot(identical(arot(arot(arot(a))),arot(a,3)))
b <- c(2,4)
stopifnot(identical(arot(arot(arot(a,p=b),p=b),p=b),arot(a,p=b,3)))
stopifnot(identical(arot(a,2),arev(a,1:2)))



#now some tests of shift:
stopifnot(identical(c(as.integer(10),1:9),shift(1:10)))
stopifnot(identical(shift(1:10,-2),c(3:10,1:2)))
stopifnot(identical(magic(4),ashift(ashift(ashift(ashift(magic(4)))))))
stopifnot(identical(ashift(ashift(ashift(magiccube.2np1(1)))),magiccube.2np1(1)))
a <- array(1:24,2:4)
stopifnot(identical(a,ashift(a,dim(a))))


stopifnot(is.magichypercube(ashift(magichypercube.4n(1))))
stopifnot(is.semimagichypercube(ashift(magichypercube.4n(1),1:3)))
stopifnot(is.semimagichypercube(ashift(magichypercube.4n(1,d=5),c(1,2,3,2,1))))


zero.extent <- array(0,c(3,0,3,2,3))
stopifnot(is.standard(zero.extent))


# Now test subsums.  With wrap=TRUE, a matrix with identical
# entries should have identical subsums:
a <- array(1,c(2,2,2,2,3,2,2))
stopifnot(minmax(subsums(a,2)))

# Now sundry tests of apltake(), apldrop(), apad(), arev() etc:
f1 <-
  function(a){
    zero <- as.integer(0)
    identical(ashift(a,dim(a)),a) &
    identical(apltake(a,dim(a)),a) &
    identical(apltake(a,-dim(a)),a) &
    identical(apldrop(a,dim(a)*0),a) &
    identical(apldrop(adiag(a,zero),-1+dim(a)*0),a) &
    identical(apldrop(adiag(a,a),dim(a)),a) &
    identical(apltake(adiag(a,a),dim(a)),a) &
    identical(apltake(adiag(a,a),-dim(a)),a) &
    identical(apldrop(apad(a,dim(a)),-dim(a)),a) &
    identical(apldrop(apad(a,dim(a),method="mirror"),dim(a)),arev(a)) &
    identical(apldrop(apad(a,dim(a)),dim(a)),array(do.call("[",c(list(a),as.list(dim(a)))),dim(a))) &
    identical(a,apldrop(ashift(adiag(a,a),dim(a)),dim(a)))    
  }


stopifnot(
          f1(magichypercube.4n(m=1,d=4)) &
          f1(array(1:24,2:4))            &
          f1(array(1:64,rep(2,7)))       &
          f1(matrix(1:30,5,6))
          )

             
    

# Some tests of do.index():
f1 <- function(x){as.integer(sum(x))}
f2 <- function(a){
  stopifnot(identical(do.index(a,f1),arow(a,1)+arow(a,2)+arow(a,3)+arow(a,4)))
}

f2(array(0L,c(2,3,4,5)))
f2(array(0L,c(3,5,4,2)))



# Some tests of the incidence functionality:

n <- 7
f <- function(a){
  is.latin(a)                                           &
  is.latin(unincidence(aperm(incidence(a),c(3,1,2))))   &
  is.latin(unincidence(aperm(incidence(a),c(3,2,1))))   &
  is.latin(unincidence(aperm(incidence(a),c(1,3,2))))
}

stopifnot(sapply(sapply(2:n , latin) , f))





#Some tests of the antimagic functions:

f <- function(x){is.sam(sam(x[1],x[2]))}
jj <- as.matrix(which(lower.tri(matrix(0,n,n)),arr.ind=TRUE))

stopifnot(all(apply(jj,1,f)))



#Some tests of fnsd():

a <- array(1:24,c(1,1,1,1,2,1,3,1,1,4))

stopifnot(length(fnsd(a,0))==0)
stopifnot(fnsd(a)==5)
stopifnot(all(fnsd(a,2) == c(5,7)))


#Some tests of recurse():



f <- function(p,i){
 stopifnot(all(recurse(start=recurse(p,i),p,-i) == seq_along(p)))
}


f(magic(10),0)
f(magic(11),1)
f(magic(12),2)
f(magic(13),3)
f(magic(14),4)




#Some multiplicative magic square tests; see Javier Cilleruelo and
#Florian Luca 2010, "On multiplicative magic squares", The Electronic
#Journal of Combinatorics", vol 17. #N8

f <- function(n,m){
  matrix(c(
           (n+2)*(m+0), (n+3)*(m+3), (n+1)*(m+2), (n+0)*(m+1),
           (n+1)*(m+1), (n+0)*(m+2), (n+2)*(m+3), (n+3)*(m+0),
           (n+0)*(m+3), (n+1)*(m+0), (n+3)*(m+1), (n+2)*(m+2),
           (n+3)*(m+2), (n+2)*(m+1), (n+0)*(m+0), (n+1)*(m+3)
           ),nrow=4,ncol=4,byrow=TRUE)
}
           
stopifnot(all(apply(as.matrix(expand.grid(seq_len(n),seq_len(n))),1,function(x){is.magic(f(x[1],x[2]),func=prod)})))



stopifnot(all(sapply(panmagic.6np1(seq_len(n)),is.panmagic)))
stopifnot(all(sapply(panmagic.6np1(seq_len(n)),is.panmagic)))
