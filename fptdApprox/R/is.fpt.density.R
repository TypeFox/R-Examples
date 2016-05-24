is.fpt.density <-
function (obj) 
{
    	if (inherits(obj, "fpt.density") & is.list(obj) & (length(obj) == 3L)){
		if (all(sapply(obj[1:2], is.numeric)) & all(is.element(names(obj), c("x", "y", "y.x0"))) & all(is.element(names(attributes(obj)), 
          	c("names", "Call", "Steps", "cumIntegral", "skips", "CPUTime", "summary.fptl", "class")))){
			if (is.vector(obj$x) & is.vector(obj$y) & is.call(attr(obj, "Call")) & is.summary.fptl(attr(obj, "summary.fptl")) & is.matrix(attr(obj, "Steps")) & 
			is.numeric(attr(obj, "Steps")) & is.numeric(attr(obj, "cumIntegral")) & is.list(attr(obj, "skips")) & is.matrix(attr(obj, "CPUTime")) & 
			is.numeric(attr(obj, "CPUTime"))){
				if ((length(obj$x) == length(obj$y)) & all(diff(attr(obj, "cumIntegral")) >= 0) & (length(attr(obj, "cumIntegral")) == length(attr(obj, "skips"))) & 
				(length(attr(obj, "cumIntegral")) == nrow(attr(obj, "CPUTime"))) & (nrow(attr(obj, "Steps")) >= nrow(attr(obj, "CPUTime"))) & 
                  	(ncol(attr(obj, "Steps")) == 3) & (ncol(attr(obj, "CPUTime")) == 2) & all(sapply(attr(obj, "skips"), is.integer))) {
					Args <- as.list(attr(obj, "Call"))
					f <- as.character(Args[[1]])
					if (is.element(f, c("Approx.fpt.density", "Approx.cfpt.density"))){
						if (f == "Approx.fpt.density"){
							if (!is.null(obj$y.x0)){
								if (is.numeric(obj$y.x0) & is.matrix(obj$y.x0)){
									if ((length(obj$x) != nrow(obj$y.x0)) | (ncol(obj$y.x0) != (length(attr(obj, "summary.fptl")))) 
										| any(obj$y.x0 < 0L)) return(FALSE)
								}
								else return(FALSE)
							}
							if (is.numeric(Args$id)) m <- 1L
							else{
								m <- Args$m
								if (is.null(m)) m <- 100
							}
						}
						else{						
							if (!is.null(obj$y.x0)) return(FALSE)
                  				m <- 1
						}
						skip <- attr(obj, "skips")
						if (any(unlist(skip) > m)) return(FALSE)
						index <- 1:length(skip)
						indexjumps <- which(sapply(skip, identical, 1:m))
                    			h <- attr(obj, "Steps")[index, 3]
                    			if (length(indexjumps) > 0L) h[indexjumps] <- attr(obj, "Steps")[indexjumps, 2] - attr(obj, "Steps")[indexjumps, 1]
			  			lowerend <- attr(obj, "Steps")[index, 1]					
                    			x1 <- lowerend[1]
						x2 <- x1 + h[1]
						lowerend[1] <- x2
						lowerend <- lowerend + h
						if (lowerend[1] > attr(obj, "Steps")[1, 2]){
							lowerend <- lowerend[-1]
							index <- index[-1]
							h <- h[-1]
						}					
						z <- c(x1, x2, unlist(mapply(seq, lowerend, attr(obj, "Steps")[index, 2], by = h)))
						return(all(obj$x == z) & all(obj$y >= 0L))
					}
				}
                	}
		}
	}
     return(FALSE)

}
