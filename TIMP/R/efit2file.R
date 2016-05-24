efit2file <- 
function (filename, skip = 2, numcol, nrows = vector()) 
{
    xx <- as.matrix(scan(filename, skip = skip, na.strings=c("NA","NAN", 
       "NaNQ")))
    if (length(nrows) == 0) 
        nrows <- length(xx)/numcol
    dim(xx) <- c(numcol, sum(nrows))
    xx <- t(xx)
    cnt <- 1
    for (i in 1:length(nrows)) {
        write.table(xx[cnt:(cnt + nrows[i] - 1), ], row.names = FALSE, 
            col.names = FALSE, file = paste("plainmat", i, filename, 
                sep = "_"))
        cnt <- cnt + nrows[i]
        cat("Wrote the file", paste("plainmat", i, filename, 
            sep = "_"), "\n")
    }
}

# either copy this file into R or use
# source("efit2file.R")
# in R, when the file is in your working directory
#
# call the function with for example:
# 
# efit2file("SADSrc685nm5S22nm1010.efit", numcol=8, nrows=c(256,256,256))
#
# where the file "SADSrc685nm5S22nm1010.efit" contains 3 sets of spectra,
# with each set having 256 rows and 8 columns.   
#
# if nrows is not given, then just 1 matrix will be written to file,
# that includes everything in the efit file, for example with: 
# efit2file("SADSrc685nm5S22nm1010.efit", numcol=8)
