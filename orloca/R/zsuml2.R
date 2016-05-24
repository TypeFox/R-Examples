#
# Functions for evaluation of objective function
# 
# zsuml2 method definition
#   Evaluation of objective function
setGeneric("zsuml2", function(o, x=0, y=0) standardGeneric("zsuml2"))

# Definition of zsum for class loca.p
setMethod("zsuml2", "loca.p", function(o, x=0, y=0)
   {
   sum(o@w*sqrt((o@x-x)^2+(o@y-y)^2))
   }
)

# zsuml2gra method definition
#   Evaluation of the gradient
setGeneric("zsuml2gra", function(o, x=0, y=0, partial=F) standardGeneric("zsuml2gra"))

# Definition os zsumgra for class loca.p
setMethod("zsuml2gra", "loca.p", function(o, x=0, y=0, partial=F)
   {
   n<- o@w/sqrt((x-o@x)^2+(y-o@y)^2)
   c(sum((x-o@x)*n, na.rm=partial), sum((y-o@y)*n, na.rm=partial))
   }
)
