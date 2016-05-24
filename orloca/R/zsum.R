#
# Functions for evaluation of objective function
# 
# zsum method definition
#   Evaluation of objective function
setGeneric("zsum", function(o, x=0, y=0, lp=numeric(0)) standardGeneric("zsum"))

# Definition of zsum for class loca.p
setMethod("zsum", "loca.p", function(o, x=0, y=0, lp=numeric(0))
   {
     if (length(lp) == 0) return(zsuml2(o=o, x=x, y=y))
     else if (lp >= 1) return(zsumlp(o=o, x=x, y=y, p=lp))
     else stop(paste(lp, gettext("is not a valid value for lp, use 1 <= lp", domain = "R-orloca")))
   }
)

# zsumgra method definition
#   Evaluation of the gradient
setGeneric("zsumgra", function(o, x=0, y=0, lp=numeric(0), partial=F) standardGeneric("zsumgra"))

# Definition os zsumgra for class loca.p
setMethod("zsumgra", "loca.p", function(o, x=0, y=0, lp=numeric(0), partial=F)
   {
     if (length(lp) == 0) return(zsuml2gra(o=o, x=x, y=y, partial=partial))
     else if (lp >= 1) return(zsumlpgra(o=o, x=x, y=y, p=lp, partial=partial))
     else stop(paste(lp, gettext("is not a valid value for lp, use 1 <= lp", domain = "R-orloca")))
   }
)
