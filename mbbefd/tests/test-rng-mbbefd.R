library(mbbefd)

testfunc <- function(x)
  c(summary(x), sd=sd(x), tl=etl(x))
extensive <- TRUE
extensive <- FALSE


# test invalid param
n <- 5
a <- 0 
b <- -1/2
mbbefd:::rmbbefdCpp(n, a, b)
mbbefd:::rmbbefdR(n, a, b)
g <- 1/2 
b <- 3

mbbefd:::rMBBEFDCpp(n, g, b)
mbbefd:::rMBBEFDR(n, g, b)


#test of MBBEFD(a,b) distribution

n <- 10
a <- 0 
b <- 1/2

mbbefd:::rmbbefdCpp(n, a, b)
mbbefd:::rmbbefdR(n, a, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rmbbefdCpp(n, a, b)))
  print(testfunc(mbbefd:::rmbbefdR(n, a, b)))
}

a <- 1/2 
b <- 1
n <- 10

mbbefd:::rmbbefdCpp(n, a, b)
mbbefd:::rmbbefdR(n, a, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rmbbefdCpp(n, a, b)))
  print(testfunc(mbbefd:::rmbbefdR(n, a, b)))
}

a <- -1/2 
b <- 3
n <- 10

mbbefd:::rmbbefdCpp(n, a, b)
mbbefd:::rmbbefdR(n, a, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rmbbefdCpp(n, a, b)))
  print(testfunc(mbbefd:::rmbbefdR(n, a, b)))
}

a <- Inf
b <- 1/3
n <- 10

mbbefd:::rmbbefdCpp(n, a, b)
mbbefd:::rmbbefdR(n, a, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rmbbefdCpp(n, a, b)))
  print(testfunc(mbbefd:::rmbbefdR(n, a, b)))
}




#test of MBBEFD(g,b) distribution

n <- 10
g <- 1 
b <- 1/2

mbbefd:::rMBBEFDCpp(n, g, b)
mbbefd:::rMBBEFDR(n, g, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rMBBEFDCpp(n, g, b)))
  print(testfunc(mbbefd:::rMBBEFDR(n, g, b)))
}

n <- 10
g <- 2 
b <- 0

mbbefd:::rMBBEFDCpp(n, g, b)
mbbefd:::rMBBEFDR(n, g, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rMBBEFDCpp(n, g, b)))
  print(testfunc(mbbefd:::rMBBEFDR(n, g, b)))
}


n <- 10
g <- 2 
b <- 1/2

mbbefd:::rMBBEFDCpp(n, g, b)
mbbefd:::rMBBEFDR(n, g, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rMBBEFDCpp(n, g, b)))
  print(testfunc(mbbefd:::rMBBEFDR(n, g, b)))
}

n <- 10
g <- 2
b <- 1

mbbefd:::rMBBEFDCpp(n, g, b)
mbbefd:::rMBBEFDR(n, g, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rMBBEFDCpp(n, g, b)))
  print(testfunc(mbbefd:::rMBBEFDR(n, g, b)))
}

n <- 10
g <- 2 
b <- 3

mbbefd:::rMBBEFDCpp(n, g, b)
mbbefd:::rMBBEFDR(n, g, b)
if(extensive)
{
  n <- 1e6
  print(testfunc(mbbefd:::rMBBEFDCpp(n, g, b)))
  print(testfunc(mbbefd:::rMBBEFDR(n, g, b)))
}