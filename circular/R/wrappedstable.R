#############################################################
#                                                           #
#   rwrappedstable function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 15, 2007                                  #
#   Copyright (C) 2007 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-3                                           #
#############################################################

rwrappedstable <- function(n,  scale=1, index, skewness, control.circular=list()) {
    dc <- control.circular
    result <- rstable(n=n, scale=scale, index=index, skewness=skewness) %% (2 * pi)
    result <- conversion.circular(circular(result), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    return(result)
}

