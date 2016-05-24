

costLS <- function(sw, o, target=seq_along(sw),...){
   -sum(sw[o]*target)
}


# costED <- function(sw, o,node,...){ 
	# #node should be se[1,i]:se[4,i]
  # costPL(sw,o[node])
# }


costED <- function(sw, o,node,se,...){ 
   sw[o[se[2,node]],o[se[3,node]]]
}



costPL <- function(sw, o,...){
  .Call(cpl,sw,o)
}

costLPL <- function(sw, o,target=(nrow(sw)-1):1,...){ 	
  .Call(clpl,sw,o,as.numeric(target))
}


costBAR <- function(sw, o,target=max(2,floor(nrow(sw)/5)),...){
  .Call(cbar,sw,o,as.integer(target))
}

costARc <- function(sw, o,target=nrow(sw)-1,...){
   if (is.matrix(target))
   .Call(carct,sw,o,as.numeric(target))
  else .Call(cbar,sw,o,as.integer(target))
}





AR_target <- function(n) {
 #generates target matrix for ARc cost function.
  mat <- matrix(0L, nrow=n, ncol=n)
  targ <- n-abs(col(mat)-row(mat))
  targ[targ<0] <- 0L
  diag(targ) <- 0L
  mode(targ) <- "integer"
  targ
}






defaultcostArg <- function(costfn,sw){
	if (identical(costfn,costBAR)) as.integer(max(2,floor(nrow(sw)/5)))
	else if (identical(costfn,costARc)) nrow(sw)-1
	else if (identical(costfn,costLPL)) as.numeric((nrow(sw)-1):1)
	else if (identical(costfn,costLS)) seq_along(sw)
	else NULL
 }


# #cost functions from criterion function in seriation package

# cost_criterion <- function(d,o,method="Neumann_stress"){
	# a <- criterion(d[o,o], method=method)[[1]]
	# if (method == "Gradient_raw" || method == "Gradient_weighted")
	 # a <- -a
	# a
# }


