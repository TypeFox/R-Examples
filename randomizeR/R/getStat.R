# --------------------------------------------
# Generic functions for an issue
# --------------------------------------------

# Get statistics for an issue
# 
# Generic function which calculates statistics dependent on a \code{randSeq}, 
# an \code{issue}, and an \code{endp} object (optional). 
#
# @aliases getStat
# 
# @param randSeq object of the class randSeq.
# @param issue object of the class issue.
# @param endp object of the class endpoint (optional).
#
# @return
# Returns a vector of the statistics for the issue. For every randomization 
# sequence one value is calculated.
# 
# @seealso
# \code{\link{assess}}
#
# @name getStat
NULL

# @rdname getStat
# 
setGeneric("getStat", function(randSeq, issue, endp) standardGeneric("getStat"))