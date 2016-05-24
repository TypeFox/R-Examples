#
# Functions for evaluation of objective function
# 
# zsum method definition
#   Evaluation of objective function
setGeneric("zsumlp", function(o, x=0, y=0, p=2) standardGeneric("zsumlp"))

# Definition of zsum for class loca.p with lp norm
setMethod("zsumlp", "loca.p", function(o, x=0, y=0, p=2)
   {
     if (p>=1) return(sum(o@w*(abs(o@x-x)^p+abs(o@y-y)^p)^(1/p)))
     else stop(paste(p, gettext("is not a valid value for p, use 1 <= p", domain = "R-orloca")))
   }
)

# zsumlpgra method definition
#   Evaluation of the gradient
setGeneric("zsumlpgra", function(o, x=0, y=0, p=2, partial=F) standardGeneric("zsumlpgra"))

# Definition os zsumlpgra for class loca.p with lp
setMethod("zsumlpgra", "loca.p", function(o, x=0, y=0, p=2, partial=F)
   {
     if (p>=1) {
       n<- o@w*(abs(x-o@x)^p+abs(y-o@y)^p)^(1/p-1)
       c(sum(sign(x-o@x)*abs(x-o@x)^(p-1)*n, na.rm=partial), sum(sign(y-o@y)*abs(y-o@y)^(p-1)*n, na.rm=partial))
     }
     else stop(paste(p, gettext("is not a valid value for p, use 1 <= p", domain = "R-orloca")))
   }
)
