################################################################################
## FOR INTERNAL USE                                                           ##
## FORMULAUNIMVA: a method to create a list of m univariate formulas          ##
## from a formula with multivariate response of dimension m                   ##
################################################################################

formulaUnimva <- function( formula, var.subset, split.x=FALSE,
intercept=0, allow.noresp=FALSE ){

if(intercept==0){
  formula <- update(formula,  ~ . + 0)
} else {
  formula <- update(formula,  ~ . + 1)
}

if(split.x){
    foo.new <- extend.x.formula(formula, return.interaction=FALSE,
       extend.term=TRUE)
} else foo.new <- formula

terms.foo <- terms(foo.new)

term.labels <- attr(terms.foo,"term.labels")

formulasUni.mvad <- list()

resp  <- attr(terms.foo,"response")

if(resp==1) {
  respname  <- foo.new[[2]]
  datay     <- eval(respname)
  respname  <- deparse(respname, width.cutoff =500)
} else {

  if(!allow.noresp) stop("no formulaUnimva's can be build: formula has no response")
  formulasUni.mvad[[1]] <- foo.new
  return(formulasUni.mvad)

}

labelsyp  <- labels(datay)[[2]]

p         <- NCOL(datay)
if(p==0) stop("no formulaUnimva's can be build: formula has no response")

if(missing(var.subset)) {
	var.subset <- 1:p
} else {
  # Change logical var.subset to numerical var.subset, if necessary. Note that
  # NA values are logical as well, but should be excluded here.
  if(is.logical(var.subset) & any(!is.na(var.subset)))
    var.subset <- which(var.subset[!is.na(var.subset)])
  if((max (var.subset)>p) | (p<1)) stop("'var.subset' is not valid")
}
   
# replace response if necessary
if(p > 1) {
    j=1
    for (i in var.subset){
	tmp <- paste("(",respname,")[,",i,"]", sep="")			#Edited by Stephen
	#fooi <- reformulate(term.labels,response=tmp)			#Edited by Stephen
	fooi <- as.formula(paste(tmp,'~',c(terms.foo[[3]])))      

	if(intercept==0) { 
		formulasUni.mvad[[j]] <- update(fooi,  ~ . + 0)
      } else {
		formulasUni.mvad[[j]] <- update(fooi,  ~ . + 1)
	}

      j=1+j
    }
} else formulasUni.mvad[[1]] <- foo.new

names(formulasUni.mvad)<- labelsyp[var.subset]

class(formulasUni.mvad) <- c("formulaUnimva", "list")

return( formulasUni.mvad)

}


is.formulaUnimva <- function(x) {
  inherits(x, "formulaUnimva")
}


setMethod("show", "formulaUnimva",
function(object) {
obj<-object@g
show(obj)
})


