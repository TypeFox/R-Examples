f.args.from.info <- function(info){
##
## EXTRACTS VALID haplin ARGUMENTS FROM AN info OBJECT
##
## NOTE: ASSUMES ALL PARAMETERS ARE ON LEVEL 2 IN info, EXCEPT filename
##
#
## EXTRACT STANDARD HAPLIN ARGUMENTS WITH DEFAULTS
.formals <- formals(haplin)

.formals$data <- NULL
.formals$pedIndex <- NULL


#
## FLATTEN info OBJECT, ONE LEVEL
.tmp <- unlist(info, recursive = F)
#
## EXTRACT NAMES FROM ONE LEVEL DOWN
.names <- sapply(info, names)
.names <- c("filename", unlist(.names[-1]))
#
## CHECK NAME LENGTH
if(length(.names) != length(.tmp)) stop()
#
## USE NAMES FROM ONE LEVEL DOWN
names(.tmp) <- .names
#
## PICK ONLY THOSE THAT ARE ARGUMENTS TO HAPLIN
.tmp <- .tmp[names(.formals)]
#
## FOR SAFETY'S SAKE
.test <- f.check.pars(.tmp, .formals)
#
##
return(.tmp)
}
