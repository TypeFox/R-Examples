################################################################################
# Random KNN fitted Functions                                                  #
# File:   plot.R                                                               #
# Author: Shengqiao Li                                                         #
# Date:   June 24, 2008 (initial)                                              #
################################################################################

fitted.rknn<-
function (object, ...) 
{
    if (!inherits(object, "rknn")) 
        stop(deparse(substitute(object)), " is not a rknn object")
    object$pred
}
################################################################################
