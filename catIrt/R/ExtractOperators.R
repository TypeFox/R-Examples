#####
# 1 # (EXTRACTING)
#####

# Making sure that selecting elements does not change the class:
"[.brm" <- function(x, ...){
  r <- NextMethod("[")
  class(r) <- class(x)
  return(r)
} # END [.brm FUNCTION

"[.grm" <- function(x, ...){
  r <- NextMethod("[")
  class(r) <- class(x)
  return(r)
} # END [.grm FUNCTION


#####
# 2 # (REPLACING)
#####

# Making sure that assignment elements does not change the class.
"[<-.brm" <- function(x, ..., value){
  r <- NextMethod("[<-")
  class(r) <- class(x)
  return(r)
} # END [<-.brm FUNCTION

"[<-.grm" <- function(x, ..., value){
  r <- NextMethod("[<-")
  class(r) <- class(x)
  return(r)
} # END [<-.grm FUNCTION	  


#####
# 3 # (SUBSETTING)
#####

sel.prm <- function(p, u, N, J, K){
	
	sel <- function(pi, ui, J) pi[ cbind(ui, 1:J) ]
	
	if(N == 1){
	  lik <- sel(p, drop(u), J)
	} else if( is.null( dim(u) ) & (J > 1) ){
	  p   <- split.data.frame(p, rep(1:N, each = K))
	  lik <- do.call(rbind, lapply(p, sel, u = u, J = J))
	} else{
	  p   <- split.data.frame(p, rep(1:N, each = K))
	  lik <- do.call(rbind, lapply(seq_along(p), function(i) sel(p[[i]], u[i, ], J)))
	} # END ifelse STATEMENTS
	
	return(lik)
	
} # END sel.prm FUNCTION