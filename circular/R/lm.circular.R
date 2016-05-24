
#############################################################
#                                                           #
#   lm.circular function                                    #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 27, 2005                                   #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

lm.circular <- function(..., type=c("c-c", "c-l")) {
    type <- match.arg(type)
    if (type=="c-c") {
        lm.circular.cc(...)
    } else {
        lm.circular.cl(...)
    }
}
