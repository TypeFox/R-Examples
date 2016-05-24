#' Adjust records in a \code{data.frame}
#'
#' A convenient wrapper around \code{\link{adjust}} that loops over all records in a
#' \code{data.frame}
#'
#' @param E a \code{\link[editrules]{editmatrix}}
#' @param dat a \code{data.frame}
#' @param adjust a \code{nrow(dat) x ncol(dat)} boolean array, indicating which fields 
#'        must be adjusted.
#' @param w a vector of length \code{ncol(dat)} or array of size \code{adjust} with adjustment weights.
#' @param verbose  print progress to console
#' @param ... extra options, passed through to \code{\link{adjust}}
#'
#'
#' @section Details:
#' This function is not written to be especiallty speedy or memory-efficient, but to offer a
#' convenient interface to adjusting a \code{data.frame} of records.
#'
#'
#'
#' @seealso \code{\link{adjust}} 
#'
#' @return An object of class \code{adjustedRecords}
#' @export
adjustRecords <- function(E, dat, adjust=array(TRUE,dim=dim(dat)), w=rep(1,ncol(dat)), verbose=FALSE, ... ){
	if (is.vector(w)){ 
		stopifnot(length(w) == ncol(dat))
		w <- t(array(w,dim=dim(dat)[2:1]))
		colnames(w) <- names(dat)
	}
		
   stopifnot(
      all(dim(adjust) == dim(dat)),
      all(getVars(E) %in% names(dat)),
	   all_finite(w),
		is.logical(adjust),
		sum(is.na(adjust))==0,
      all(w>0),
		all(dim(w) == dim(dat))
   )
	
   nm <- names(dat)
   if ( is.null(colnames(adjust)) ) colnames(adjust) <- nm
   if ( is.null(colnames(w)) ) colnames(w) <- nm 

   B <- blocks(E)
   status = NULL
   for ( i in 1:length(B) ){
      if (verbose ) cat(sprintf("adjusting block %4d of %4d\n",i, length(B)))
      e <- B[[i]]
		
      vars <- nm[nm %in% getVars(e)]
      adj <- adjustBlock(e, dat[vars], adjust[,vars,drop=FALSE], w[,vars,drop=FALSE], verbose=verbose, ...) 
      dat[vars] <- adj$adjusted
      status <- status %++% adj$status 
   }

   structure(list(adjusted=dat, status=status), class="adjustedRecords")

}



adjustBlock <- function(E, dat, adjust, w, verbose, ...){

	out <- t(dat)
	n <- nrow(dat)
	acc <- numeric(n)
   obj <- numeric(n)
	tpl <- getDuration(proc.time())
	dur <- array(0,dim=c(n,length(tpl)))
	colnames(dur) <- names(tpl)
	nit <- numeric(n)
	status <- new_status(rep(NA,n))

	for ( i in 1:nrow(dat) ){
		if( verbose ) cat( sprintf("\r          record %4d / %d ",i,n))
		r <- do.call(c,as.list(dat[i,]))
		J <- adjust[i,]
		if (!any(J)) next
	
		e <- reduce(substValue(E,names(r)[!J],r[!J]))
		y <- adjust(e, r[J], w = w[i,J],...)
		out[J,i]    <- y$x
		acc[i]      <- y$accuracy
      obj[i]      <- y$objective
		dur[i,]     <- getDuration(y$duration)
		status[i]   <- y$status
		nit[i]      <- y$niter
		status[i]   <- y$status
	}
   if(verbose) cat("\n")

	list(
		adjusted = as.data.frame(t(out)),
		status = data.frame(
			accuracy = acc,
         objective = obj,
		   niter    = nit,
			status   = status,
			dur
		)
	)
}






