suppressMessages(library("circular"))
# ?median.circular


  med <- median.circular(circular(c(0, pi/2, pi, 3/2*pi)))

## expect:
## median is equal to NA
  stopifnot(
    is.na(med)
  )

## expect:
## all the points are minimizers of the function
  stopifnot(
    all.equal(attr(med, "medians"), c(0, pi/2, pi, 3/2*pi), tol = 2e-7)
  )

