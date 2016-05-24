model.matrix.varComp <-
function(object, what=c('fixed','random', 'varcov', 'X', 'K', 'Z'), ...)
{
# S3 method for getting design matrix of object of class varComp. 
# what: 'fixed' or 'X': Fixed effect design matrix. 
#       'random' or 'Z': random effect design matrix. NOTE: This is an equivalent design matrix, not necessarily the one supplied to varComp. 
#       'K': Returns the kernal matrix. 

  what=match.arg(what)
  what=switch(what, random='Z', fixed='X', varcov='K', what)
  if(what=='X'){
	if('X'%in%names(object)) object$X else stop("X matrix is not recored")
  }else if (what=='K' || what=='Z'){
	if(!('K'%in%names(object))) stop("K matrices are not recored")
	if(what=='Z') lapply(object$K, cholRoot) else object$K
  } 
}
