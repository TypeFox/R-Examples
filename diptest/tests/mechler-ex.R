library(diptest)
## These are from
## the 217-readme.doc file that explains the bug fixed by
## Ferenc Mechler (fmechler@med.cornell.edu). [5/Sep/2002]
##
ex1 <- c(0.0198, 0.0198, 0.1961, 0.2898, 0.3184, 0.3687,
         0.4336, 0.4987, 0.5661, 0.6530, 0.7476, 0.8555)

ex2 <- c(0.0198, 0.1961, 0.2898, 0.3184, 0.3687, 0.4336,
         0.4987, 0.5661, 0.6530, 0.7476, 0.8555, 0.9912)

## Multiply them by 10000 here:

(D1 <- dip(10000*ex1, full=TRUE, debug=2))
str(D1, digits = 10, vec.len = 12)

(D2 <- dip(10000*ex2, full=TRUE, debug=2))
str(D2, digits = 10, vec.len = 12)
