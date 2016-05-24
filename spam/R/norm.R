# This is file ../spam/R/norm.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


########################################################################

#norm <- function(x, type = "sup", ...){
#  typ <- charmatch(tolower(type), c("sup",'l1',"frobenius","hs"))
#  if (is.na(typ))          stop("undefined norm '",type,"'.",call.=FALSE)

#  switch(typ,
#         max(abs(x)),
#         sum(abs(x)),
#         sqrt(sum(x^2)),sqrt(sum(x^2))
#         )
#}
norm.spam <- function(x, type = "m", ...){
  typ <- substr(tolower(type),1,1)

  if (typ %in% c("o", "1")) {
    return( max( colSums(abs(x))))
  }
  if (typ %in% c("i")) {
    return( max( rowSums(abs(x))))
  }
  if (typ %in% c("f", "h")) {
    return( sqrt(sum(x@entries^2)))
  }
  if (typ %in% c("m","s")) {
    return( max(abs(x@entries)) )
  }
  
  stop("undefined norm '",type,"'.",call.=FALSE)

}

setMethod("norm",signature(x="spam",type="character"), 
          function(x, type, ...) norm.spam(x, type))

setMethod("norm",signature(x="spam",type="missing"),
          function(x, type, ...) norm.spam(x, type="O"))

setMethod("norm", signature(x = "numeric", type = "character"),
          function(x, type, ...) base::norm(as.matrix(x), type))

setMethod("norm", signature(x = "numeric", type = "missing"),
          function(x, type, ...) base::norm(as.matrix(x),  type="O"))

setMethod("norm", signature(x = "matrix", type = "character"),
	  function(x, type, ...) base::norm(x, type))
setMethod("norm", signature(x = "matrix", type = "missing"),
	  function(x, type, ...) base::norm(x, type='o'))


