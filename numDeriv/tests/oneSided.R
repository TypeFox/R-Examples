# test one-sided derivatives

 library(numDeriv)

fuzz <- 1e-8

##### scalar argument, scalar result (case 1)#####
 f <- function(x) if(x<=0) sin(x) else  NA 
################################################## 

## grad

err <- 1.0 - grad(f, x=0,    method="simple", side=-1) 
if( fuzz < err )    stop("grad case 1 method simple one-sided test 1 failed.")

if( ! is.na(grad(f, x=0, method="simple", side=1)))  stop("grad case 1 method simple one-sided test 2 failed.")  

err <- 1.0 - grad(f, x=0,    method="Richardson", side=-1)
if( fuzz < err ) stop("grad case 1 method Richardson one-sided test 1 failed.")

# print(grad(sin, x=-0.5, method="Richardson")         , digits=16)  # 0.8775825618862814
# print(grad(sin, x=-0.5, method="Richardson", side=-1), digits=16)  # 0.8775807270501326

err <- 0.8775807270501326 - grad(sin, x=-0.5, method="Richardson", side=-1)
if( fuzz < err ) stop("grad case 1 method Richardson one-sided test 2 failed.")


## jacobian

err <- 1.0 - jacobian(f, x=0,    method="simple", side= -1)
if( fuzz < err ) stop("jacobian case 1 method simple one-sided test failed.")

err <- 1.0 - jacobian(f, x=0,    method="Richardson", side= -1)
if( fuzz < err ) stop("jacobian case 1 method Richardson one-sided test 1 failed.")

if( ! is.na(jacobian(f, x=0, method="Richardson", side= 1))) stop("jacobian case 1 method Richardson one-sided test 2 failed.")



##### vector argument, vector result (case 3)#####
 f <- function(x) if(x[1]<=0) sin(x) else  c(NA, sin(x[-1]))
################################################## 

## grad

err <- 1.0 -  grad(f, x=c(0,0), method="simple", side=c(-1, -1))  #  1 1
if( fuzz < max(err) )  stop("grad case 3 method simple one-sided test 1 failed.")

err <- 1.0 -  grad(f, x=c(0,0), method="simple", side=c(-1,  1))  #  1 1
if( fuzz < max(err) )  stop("grad case 3 method simple one-sided test 2 failed.")

err <- 1.0 -  grad(f, x=c(0,0), method="simple", side=c(-1, NA))  #  1 1
if( fuzz < max(err) )  stop("grad case 3 method simple one-sided test 3 failed.")

err <- 1.0 -  grad(f, x=c(0,0), method="simple", side=c( 1,  1))  #  NA 1
if( fuzz < err[2] )    stop("grad case 3 method simple one-sided test 4 failed.")
if(!is.na( err[1]) )   stop("grad case 3 method simple one-sided test 4b failed.")


err <- 1.0 -  grad(f, x=c(0,0), method="Richardson", side=c(-1, -1)) #   1 1 
if( fuzz < max(err) )  stop("grad case 3 method Richardson one-sided test 1 failed.")

err <- 1.0 -  grad(f, x=c(0,0), method="Richardson", side=c(-1,  1)) #  1 1
if( fuzz < max(err) )  stop("grad case 3 method Richardson one-sided test 2 failed.")

err <- 1.0 -  grad(f, x=c(0,0), method="Richardson", side=c(-1, NA)) #  1 1
if( fuzz < max(err) )  stop("grad case 3 method Richardson one-sided test 3 failed.")

## jacobian

err <- 1.0 - jacobian(f, x=0,    method="simple", side= -1)
if( fuzz < err ) stop("jacobian case 3 method simple one-sided test failed.")

err <- 1.0 - jacobian(f, x=0,    method="Richardson", side= -1)
if( fuzz < err ) stop("jacobian case 3 method Richardson one-sided test 1 failed.")

if( ! is.na(jacobian(f, x=0, method="Richardson", side= 1))) stop("jacobian case 3 method Richardson one-sided test 2 failed.")



##### vector argument, scalar result (case 2)#####
 f <- function(x) if(x[1]<=0) sum(sin(x)) else  NA
################################################## 
 
## grad

err <- 1.0 - grad(f, x=c(0,0), method="simple", side=c(-1, -1))  #  1 1
if( fuzz < max(err) )  stop("grad case 2 method simple one-sided test 1 failed.")

err <- 1.0 - grad(f, x=c(0,0), method="simple", side=c(-1,  1))  #  1 1
if( fuzz < max(err) )  stop("grad case 2 method simple one-sided test 2 failed.")

err <- 1.0 - grad(f, x=c(0,0), method="simple", side=c(-1, NA))  #  1 1
if( fuzz < max(err) )  stop("grad case 2 method simple one-sided test 3 failed.")

err <- 1.0 - grad(f, x=c(0,0), method="simple", side=c( 1,  1))  #  NA 1
if( fuzz < err[2] )  stop("grad case 2 method simple one-sided test 4 failed.")
if(!is.na( err[1]) ) stop("grad case 2 method simple one-sided test 4b failed.")


err <- 1.0 - grad(f, x=c(0,0), method="Richardson", side=c(-1, -1)) #  1 1 
if( fuzz < max(err) )  stop("grad case 2 method Richardson one-sided test 1 failed.")

err <- 1.0 - grad(f, x=c(0,0), method="Richardson", side=c(-1,  1)) #  1 1
if( fuzz < max(err) )  stop("grad case 2 method Richardson one-sided test 2 failed.")

err <- 1.0 - grad(f, x=c(0,0), method="Richardson", side=c(-1, NA)) #  1 1
if( fuzz < max(err) )  stop("grad case 2 method Richardson one-sided test 3 failed.")

## jacobian

err <- 1.0 - jacobian(f, x=0,    method="simple", side= -1)
if( fuzz < err ) stop("jacobian case 2 method simple one-sided test failed.")

err <- 1.0 - jacobian(f, x=0,    method="Richardson", side= -1)
if( fuzz < err ) stop("jacobian case 2 method Richardson one-sided test 1 failed.")

if( ! is.na(jacobian(f, x=0, method="Richardson", side= 1))) stop("jacobian case 2 method Richardson one-sided test 2 failed.")



