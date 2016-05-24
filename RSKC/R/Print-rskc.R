print.summary.rskc <-
function(x,...){
     if (is.null(x$L1)) {sparse<-FALSE; x$L1<-NULL} else{ sparse<-TRUE}
    sizeV<-NULL;C<-x$labels;n<-x$N;ncl<-x$ncl
	for ( i in 1 : ncl){
		sizeV[i]<-sum(C==i)
		}
    # Input
    if (x$miss){missin<-"YES"}else{missin<-"NO"}
    cat("\nInput:\n\nMissing values?",missin,
        "\n#obs=",x$N," #feature=",x$p,
 	    "\nL1=",x$L1," nstart=",x$nstart,
 	    "\nscaling=",x$scaling," correlation=",x$correlation
        ,"\nalpha=",x$alpha) 

    cat("\n\nResult:")
    # Result output
    if (sparse){    
 	cat("\n\nwbss:",x$WBSS[length(x$WBSS)])}
    else{
        x$weights<-rep(1,x$N)
        cat("\n\nwwss:",x$WWSS[length(x$WWSS)])
        }
    cat("\ncases trimmed in the squared weighted Euclidean dist:",x$oW,
 	    "\ncases trimmed in the squared Euclidean dist:", x$oE,
 	    "\n#non-zero weights:",sum(!(x$weights==0)),
 	    "\n",ncl,"clusters of sizes",paste(sizeV, collapse = ", "),"\n"
         )
        if (x$N<100) cat("Cluster labels:", x$labels,"\n")
 }
 
 	    		 	        		
 		

summary.rskc<-function(object,...){
	class(object)<-"summary.rskc"
	return(object)
	}

	
print.rskc<-function(x,...){
    sizeV<-NULL;C<-x$labels;n<-x$N;ncl<-x$ncl
    uniC <- unique(C)
    for ( i in 1 : ncl){
		sizeV[i]<-sum(C==uniC[i])
		}
    if (is.null(x$L1)){sparse<-FALSE; x$L1="NULL"}else{ sparse<-TRUE}
    if (is.character(x$oW)) x$oW<-NULL
    if (is.character(x$oE)) x$oE<-NULL
    # Input
    cat("\nInput:",
        "\n#obs=",x$N," #feature=",x$p,
 	"\nL1=",x$L1," alpha=",x$alpha) 

    cat("\n\nResult:")
    # Result output
    if (sparse){    
 	cat("\nwbss:",x$WBSS[length(x$WBSS)])}
    else{
        x$weights<-rep(1,x$N)
        cat("\n\nwwss:",x$WWSS[length(x$WWSS)])
        }
    trimmed<-union(x$oW,x$oE)
    if (length(trimmed)<20) cat("\ntrimmed cases:",trimmed)
      cat("\n#non-zero weights:",sum(!(x$weights==0)),
 	 "\n",ncl,"clusters of sizes",paste(sizeV, collapse = ", "),"\n"
         )
        if (x$N<100) cat("Cluster labels:", x$labels,"\n")
	}
