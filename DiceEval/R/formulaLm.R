formulaLm <- function(X,Y){

	X <- data.frame(X)
	p <- dim(X)[2]
	Y <- data.frame(Y)  

  F <- paste(names(Y),"~",names(X)[1],sep="")
  if (p>1){
    for (i in 2:p){
	F <- paste(F,"+",names(X)[i],sep="")
    }
  }
  return(as.formula(F))
}