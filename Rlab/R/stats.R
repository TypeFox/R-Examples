"stats" <-

function(x, by)



{



	# coerce into matrix or list format



	if(!missing(by)) {



		x <- cat.to.list(c(x), by)



		# if(!no.print) {



			# cat("Data from ", as.character(substitute(x)),



			#	"is grouped based on the categories in ",



			#	as.character(substitute(by)), fill = TRUE)



		# }



	}



	if(!is.list(x) & !is.matrix(x))



		x <- matrix(x, ncol = 1.)



	if(is.list(x)) {



		ncol <- length(names(x))



		out <- matrix(NA, ncol = ncol, nrow = length(describe(



			)))



		dimnames(out) <- list(describe(), names(x))



		for(j in (1.:ncol)) {



			if(is.numeric(x[[j]])) {



				out[, j] <- describe(x[[j]])



			}



		}



		return(out)



	}



	if(is.matrix(x)) {



		nc <- ncol(x)



		out <- matrix(NA, ncol = nc, nrow = length(describe(



			)))



		dimnames(out) <- list(describe(), dimnames(x)[[2.]])



		for(j in (1.:nc)) {



			out[, j] <- describe(x[, j])



		}



		return(out)



	}



}

