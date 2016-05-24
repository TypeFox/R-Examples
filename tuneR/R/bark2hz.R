# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

bark2hz <- function(z){
    
    if(!is.numeric(z) || z < 0)
      stop("frequencies have to be non-negative")

    # Hynek's formula (taken from rasta/audspec.c)
    600 * sinh(z/6)
}

