require(hypergeo)
#  First some test values for z, generated randomly.  These values are
#  selected so as to use each of the subfunctions of hypergeom():

z <-
  c(
    0.28 - 0.21i,
   -0.79 - 0.40i,
    0.56 + 0.05i,
    2.13 + 0.68i, 
   -0.43 - 1.47i,
    1.23 + 0.48i
    )

# First test:  some randomly chosen values for A,B,C; my package:
my.ans <- hypergeo(A=1.21,B=1.443,C=1.88,z=z)

# Compare with Maple's:
maple.ans <-
  c(
    1.26242625279896547050 - 0.33515733250598455101i  ,
    0.55218026726586346626 - 0.11826012518395685586i  ,
    2.08656513655255552950 + 0.21074089910422408972i  ,
   -0.71334498434004598167 + 0.59647479085505445853i  ,
    0.36911468617705947291 - 0.35906488952504313903i  ,
   -0.44806924103752606401 + 1.91140611055324833040i
    )



stopifnot(all(Mod(my.ans-maple.ans) < 1e-9))


# Now verify eqn 15.1.15, p556:
"eqn15.1.15a" <- function(a,z){
  hypergeo(A=a, B=1-a, C=3/2, z=sin(z)^2,tol=1e-9)
}

"eqn15.1.15b" <- function(a,z){sin((2*a-1)*z)/((2*a-1)*sin(z))}

Mod(eqn15.1.15a(0.2,z)-eqn15.1.15b(0.2,z)) # Observe the mismatch for
                                           # element 4.  BUT
                                           # hypergeo() agrees with
                                           # Maple, which gives
                                           # 1.038571815 - 0.1598166166i



"eqn15.2.10" <- function(A,B,C,z){    # page 558
  (C-A)*hypergeo(A-1,B,C,z) + (2*A-C-A*z+B*z)*hypergeo(A,B,C,z) + A*(z-1)*hypergeo(A+1,B,C,z)
}

stopifnot(all(Mod(eqn15.2.10(0.1,0.44,0.611,z))<1e-6))


## Now some special elementary cases, on page 556.  Here, "rhs" means
## "right hand side" and "lhs" means "left hand side".

"f15.1.3" <- function(z){
  lhs <- hypergeo(1,1,2,z)
  rhs <- -log(1-z)/z
  return(rhs-lhs)
}

"f15.1.4" <- function(z){
  lhs <- hypergeo(1/2,1,3/2,z^2)
  rhs <- 0.5*log((1+z)/(1-z))/z
  return(rhs-lhs)
}

"f15.1.5" <- function(z){
  lhs <- hypergeo(1/2,1,3/2,-z^2)
  rhs <- atan(z)/z
  return(rhs-lhs)
}

"f15.1.6a" <- function(z){
  lhs <- hypergeo(1/2,1/2,3/2,z^2)
  rhs <-  sqrt(1-z^2)*hypergeo(1,1,3/2,z^2)
  return(rhs-lhs)
}

"f15.1.6b" <- function(z){
  lhs <- hypergeo(1/2,1/2,3/2,z^2)
  rhs <- asin(z)/z
  return(rhs-lhs)
}

"f15.1.7a" <- function(z){
  lhs <- hypergeo(1/2,1/2,3/2,-z^2)
  rhs <- sqrt(1+z^2)*hypergeo(1,1,3/2,-z^2)
  return(rhs-lhs)
}

"f15.1.7b" <- function(z){
  lhs <- hypergeo(1/2,1/2,3/2,-z^2)
  rhs <- log(z+sqrt(1+z^2))/z
  return(rhs-lhs)
}

f.all <- function(z){
  cbind(
    f15.1.3 (z),
    f15.1.4 (z),
    f15.1.5 (z),
    f15.1.6a(z),
    f15.1.6b(z),
    f15.1.7a(z),
    f15.1.7b(z)
    )
}
    
stopifnot(max(Mod(f.all(z))) < 1e-10)



# Below, jjR means value obtained from R via the package, and jjM means the value given by maple.

# Following test fails sometimes.  It passes on my mac, fails on the linuxbox:
jjR <- genhypergeo_contfrac_single(U=0.2 , L=c(9.9,2.7,8.7) , z=1+10i)
jjM <- 1.0007289707983569879 + 0.86250714217251837317e-2i
stopifnot(Mod(jjR-jjM)<1e-10)

# Test hypergeo_cover1():
jjR <- hypergeo(pi,pi/2,3*pi/2-4, z=0.1+0.2i)  # ie negative m;  ie f15.3.12()
jjM <- 0.53745229690249593045 + 1.8917456473240515664i
stopifnot(Mod(jjR-jjM)<1e-10)

jjR <- hypergeo(pi,pi/2,3*pi/2-4, z=10.1+0.2i)  # another negative m
jjM <- 0.31486642443024026933e-2 + 0.10505111398350790590e-2i
stopifnot(Mod(jjR-jjM)<1e-10)

jjR <- hypergeo(pi,pi/2,3*pi/2, z=0.1+0.2i)  # m=0 (ie 15.3.10)
jjM <- 1.0654685003741342889 +0.24452141417139649656i
stopifnot(Mod(jjR-jjM)<1e-10)

jjR <- hypergeo(pi,pi/2,3*pi/2+4, z=10.1+0.2i)  # This is positive m (15.3.11)
jjM <- -0.29639970263878733845 - 0.34765230143995441172i
stopifnot(Mod(jjR-jjM)<1e-10)

jjR <- hypergeo(pi,pi/2,3*pi/2+4, z=10.1+0.2i)  # m>0
jjM <- -0.29639970263878733845 -0.34765230143995441172i
stopifnot(Mod(jjR-jjM)<1e-10)

jjM <- -0.90818366414720846181e-2 - 0.10746858256201734833i
jjR <- hypergeo(pi,pi/2,3*pi/2, z=10.1+0.2i)  # m=0
stopifnot(Mod(jjR-jjM)<1e-10)


# Test hypergeo_cover2():
jjM <- -0.15888831928748121465e-5 + 0.40339599711492215912e-4i
jjR <- hypergeo(pi,pi+2, 1.1 , 1+10i)
stopifnot(Mod(jjR-jjM)<1e-10)


# Test hypergeo_cover3()
jjM <- -0.24397135980533720308e-1 + 0.28819643319432922231i
jjR <- hypergeo(pi, 1.4, pi+4,1+6i)
stopifnot(Mod(jjR-jjM)<1e-10)

jjM <- -0.10592465301475818414e-1 - 0.15993048891187879153e-1i
jjR <- hypergeo(pi , pi+1 , pi + 3 , 1+6i)
stopifnot(Mod(jjR-jjM)<1e-10)

# Test for hypergeo_taylor():
jjR <- hypergeo(pi,-4,2.2,1+5i)
jjM <- 1670.8287595795885335 - 204.81995157365381258i
stopifnot(Mod(jjR-jjM)<1e-10)


# quick test for a bug reported by Igor Kojanov

options(warn=2)
ignore <- hypergeo(1,2,3,0)
ignore <- hypergeo(1,1.64,2.64,-0.1111)


# another test for a bug reported by John Ormerod, following my
# ill-considered change from gamma(x)*gamma(y)/gamma(z) to
# exp(lgamma(x)+lgamma(y)-lgamma(z)).  Now using the much superior
# .f3() and .f4() notation.

## MMA> N[Hypergeometric2F1[525/100,1,65/10,501/1000],30]
## Out[6]= 1.70239432012007391092082702795


jjR <- hypergeo(5.25,1,6.5,0.501)
jj_Mathematica <- 1.70239432012007391092082702795
stopifnot(Mod(jjR-jj_Mathematica) < 1e-10)






 ## another test for a typo I corrected:


 A <- 1
 B <- 3
 m <- 2
jj_R1 <- hypergeo(A,B,A+B+m,0.9+0.01i)
jj_R2 <- f15.3.11(A,B,m,0.9+0.01i)

mma_string <- 'N[Hypergeometric2F1[1,3,1+3+2,9/10+I/100],30]'
jj_mathematica <- 2.04925816767572859287575551415 + 0.03163091033158608113987207436i
stopifnot(abs(jj_R1-jj_mathematica) < 1e-10)
stopifnot(abs(jj_R2-jj_mathematica) < 1e-10)

