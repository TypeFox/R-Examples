peaks.in.clusters <- function(result, solution, new.order=1:solution){

    n.errors <- dim(result$measures.synthetic.peaks)[1]
    n.level <- dim(result$measures.synthetic.peaks)[2]
    cluster.rer <- change.order.clusters(result$cluster.assignment[[solution]], new.order)
    
    error.cluster.pos <- array(dim=c(n.errors,n.level))
    dimnames(error.cluster.pos) <- list(result$error.names, paste("level",1:n.level))
    cat("level &", paste(1:6, collapse="&"), "\\\\\n")
    cat("\\hline\n")
    for(error in 1:n.errors){
       if(result$use.som){
	       error.cluster.pos[error,] <-  cluster.rer$cluster[1+result$som.pos.errors[error,,1]+result$som.pos.errors[error,,2]*result$som.dim[1]]
	       cat(result$error.names[error], "&", paste( LETTERS[error.cluster.pos[error,]]  , collapse="&"), "\\\\\n")
       } else {
	       for(level in 1:n.level){
		    error.cluster.pos[error,level] <- 
		      which.min(apply(cluster.rer$centers, 1, FUN=eD, 
				      y= result$measures.uniform.synthetic.peaks[error,level,]))
	       }

	       cat(result$error.names[error], "&", paste(LETTERS[error.cluster.pos[error,]] , collapse="&"), "\\\\\n")
       }
    }

    
    cat("\n\nCluster & Error & Level\\\\\n")
    for(cluster in 1:solution){
       cluster.printed=FALSE
       for(error in 1:n.errors){
          level <- which(error.cluster.pos[error,]==cluster)
          if(length(level)>0){
              if(!cluster.printed){
                 cat(LETTERS[cluster])
                 cluster.printed=TRUE
              }
              cat(sep="", " & ", result$error.names[error],  " & ",
              paste(level, collaps=" "),  "\\\\\n")
          }
       }
    }

   return(error.cluster.pos)
}
