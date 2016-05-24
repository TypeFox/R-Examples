# The ISOpureR package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION: ISOpure.util.matlab_log.R #############################################################
# 
# Input variables: 
#    x: a scalar or matrix/vector which may contain negative or complex number 
# 
# Output variables:
#    - log(x) if all entries of x >0
#    - for complex or negative input, x, where x = a + bi, the log function returns:
#      log(z) = log(abs(z)) + 1i*atan2(b,a)
#      where atan(b,a) is on the half-closed interval, (-pi, pi], as for the Matlab 
#      log function.
#
# Note: this function was made to match the Matlab log function in its output on
# negative entries

ISOpure.util.matlab_log <- function(x){

	indices_negative <- NULL;
	if (!is.complex(x)) {	
		indices_negative <- which(x<0);
	}
	
	# if there are any negative values in the vector/matrix, then take complex log
	if (length(indices_negative>0)) {
		# the positive, NA, NaN, or (+)Inf values, you can just take regular log of
		x[-indices_negative] <- log(x[-indices_negative]);
		x[indices_negative] <- as.complex(x[indices_negative]); 
		temp_imaginary <- atan2(Im(x[indices_negative]), Re(x[indices_negative]));
		x[indices_negative] <- complex(real = log(abs(x[indices_negative])) , imaginary = temp_imaginary);
	} else {
		# if there are no negative values in the vector/matrix, then regular log is fine
		x <- log(x);
	}

	return(x);

}
