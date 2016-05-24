RelComb <-
function(Combination , Delta , Crossing=matrix(0,nrow=0,ncol=0) , ParentPop=rep(0,0) , ShowIdentifiable=TRUE){

if (class(Delta) != "list"){
stop("Delta has to be a list")
}

if (length(Delta) != 9 && length(Delta) != 15){
stop("There is a problem with the number of relatedness coefficients (length of Delta)")
}

if (class(Combination)=="character"){
	if (length(Delta)==9){
		Combination <- switch(Combination,
			"simple relatedness" = c(1,0,1/2,0,1/2,0,1/2,1/4,0),
			"double relatedness" = c(1,0,0,0,0,0,1,0,0),
			"first inbreeding" = c(1,1,1,1,0,0,0,0,0),
			"second inbreeding" = c(1,1,0,0,1,1,0,0,0),
			"double inbreeding" = c(1,1,0,0,0,0,0,0,0))
	}else{
		Combination <- switch(Combination,
			"simple relatedness" = c(1,0,1/2,1/2,0,1/2,1/2,0,1/2,1/2,1/4,1/4,1/4,1/4,0),
			"double relatedness" = c(1,0,0,0,0,0,0,0,1,1,0,0,0,0,0),
			"first inbreeding" = c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0),
			"second inbreeding" = c(1,1,0,0,0,1,1,1,0,0,0,0,0,0,0),
			"double inbreeding" = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0))
	}
}


if (length(Combination) != length(Delta)){
stop("The length of Combination differs from the number of possible relatedness coefficients (length of Delta)")
}

NbIBD <- length(Delta)
NbHyb <- dim(Delta[[1]])[1]

if (nrow(Crossing)==0){
	Crossing <- matrix(1:(2*NbHyb),ncol=2,byrow=T)
}else{
	if (nrow(Crossing) != nrow(Delta[[1]])){
	stop("There is a problem between the number of individuals in Delta and in Crossing")
	}
	
	if (length(ParentPop)==0){
	ParentPop <- rep(1,max(Crossing))
	}else{
		if (length(ParentPop) != max(Crossing)){
		stop("There is a problem between the number of parents in Crossing and in ParentPop")
		}
	}
}



if (NbIBD==9){    
Ker <- cbind(c(0,1,0,-1,0,-1,-1,2,0))
} else {
Ker <- cbind(c(0,-1,0,0,1,0,0,1,0,1,0,-1,-1,0,0),c(0,0,0,0,0,0,0,0,-1,1,1,-1,-1,1,0))
}

Sum <- Reduce('+' , lapply(1:NbIBD , function(x) Combination[x]*Delta[[x]]))

if (prod(colSums(Ker*Combination)==0)==0){
print('The combination is not identifiable')
if (ShowIdentifiable==TRUE){
Couple <- combn(1:NbHyb , 2 , simplify=F)
Vec <- sapply(Couple , function(x) .ConditionIdentifiabilityCouple(Crossing[x,],ParentPop[c(Crossing[x,])]))
mat <- matrix(NA,NbHyb,NbHyb)
MatIdentifiability <- .MatTriSup(mat,Vec)
diag(MatIdentifiability) <- 1
Res <- Sum * MatIdentifiability
}else{
Res <- Sum
}
}else{
print('The combination is identifiable')
Res <- Sum
}

return(Res)
}
