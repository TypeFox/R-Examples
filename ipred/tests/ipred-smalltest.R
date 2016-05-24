
library(ipred)

# check if SdiffKM computes
# 
# int_start^stop  (exp(-h*t) - c)^2 dt 
# 
# in the correct way 

# low-level interface needed
myfoo <- function(times, prob, h, window=0.0001) {
 .Call("SdiffKM", as.double(c(0, times)), 
        as.double(c(prob[1], prob)), as.double(c(h,
        window)), PACKAGE = "ipred")
}

# to compare with
mexp <- function(start, stop, haz, c=0) {
  foo <- function(t)
    exp(-2*haz*t)/(-2*haz) - 2*c*exp(-haz*t)/(-haz) + c^2*t
  foo(stop) - foo(start)
}
 

times <- seq(from=0.01, to=8, by=0.01)

for (myc in c(0,0.5,0.9)) {
  for (h in c(1,2,3)) {
    prob <- rep(myc, length(times))   
    a <- round(mexp(0, max(times), h, c=myc),7)
    b <- round(myfoo(times, prob, h), 7)
    stopifnot(all.equal(a,b))
  }
}
