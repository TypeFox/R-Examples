# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

hz2bark <- function(f){
    
    if(!is.numeric(f) || f < 0)
      stop("frequencies have to be non-negative")

    # Inverse of Hynek's formula (see bark2hz)
    6 * asinh(f/600)
}

