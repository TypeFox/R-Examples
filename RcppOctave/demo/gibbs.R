
## Gibbs sampler in R and Octave
##
## Compare with RcppGibbs/ example in Rcpp package

## Simple Gibbs Sampler Example
## Adapted from Darren Wilkinson's post at
## http://darrenjw.wordpress.com/2010/04/28/mcmc-programming-in-r-python-java-and-c/
##
## Sanjog Misra and Dirk Eddelbuettel, June-July 2011

suppressMessages(library(compiler))
suppressMessages(library(rbenchmark))
suppressMessages(library(RcppOctave))

## This is Darren Wilkinsons R code (with the corrected variance)
## But we are returning only his columns 2 and 3 as the 1:N sequence
## is never used below
Rgibbs <- function(N,thin) {
    mat <- matrix(0,ncol=2,nrow=N)
    x <- 0
    y <- 0
    for (i in 1:N) {
        for (j in 1:thin) {
            x <- rgamma(1,3,y*y+4)
            y <- rnorm(1,1/(x+1),1/sqrt(2*(x+1)))
        }
        mat[i,] <- c(x,y)
    }
    mat
}

## We can also try the R compiler on this R function
RCgibbs <- cmpfun(Rgibbs)

Mgibbs <- OctaveFunction('
  function mat = Mgibbs(N, thin)
    mat = zeros(N, 2);
    x = 0;
    y = 0;
    for i = 1:N
      for j = 1:thin
        x = randg(3) / (y*y+4);
        y = randn(1)*1/sqrt(2*(x+1)) + 1/(x+1);
      end
      mat(i,:) = [ x, y ];
    end
  end
')

## Octave docs:
##     `gamma (a, b)' for `a > -1', `b > 0'
##               r = b * randg (a)

## CORRECTION:  set.seed(42); .O$gam(2,3)
##         ==   set.seed(42); rgamma(1,2,1,3)
##         ==   set.seed(42); rgamma(1,2,1)*3
##         ==   set.seed(42); o_rgamma(1,1,2,3)
# set.seed(42); .O$gam(47,88); set.seed(42); rgamma(1,47,1,88);  o_rgamma(1,1,47,88)
# set.seed(42); .O$gam(47,88); set.seed(42); rgamma(1,47,1,88); set.seed(42); o_rgamma(1,1,47,88)

set.seed(42)
matR <- Rgibbs(1000,10)
set.seed(42)
matO <- Mgibbs(1000,10)
stopifnot(all.equal(matR, matO))

##print(summary(matR))
##print(summary(matO))

## also use rbenchmark package
N <- 1000
thn <- 100
res <- benchmark(Rgibbs(N, thn),
                 RCgibbs(N, thn),
                 Mgibbs(N, thn),
                 columns=c("test", "replications", "elapsed",
                           "relative", "user.self", "sys.self"),
                 order="relative",
                 replications=10)
print(res)

