"na.gam.replace" <-
function(frame)
{
	vars <- names(frame)
	if(!is.null(resp <- attr(attr(frame, "terms"), "response"))) {
		vars <- vars[ - resp]
		x <- frame[[resp]]
		pos <- is.na(x)
		if(any(pos)) {
			frame <- frame[!pos,  , drop = FALSE]
			warning(paste(sum(pos), 
				"observations omitted due to missing values in the response"
				))
		}
	}
	for(j in vars) {
		x <- frame[[j]]
		pos <- is.na(x)
		if(any(pos)) {
			if(length(levels(x))) {
				xx <- as.character(x)
				xx[pos] <- "NA"
				x <- factor(xx, exclude = NULL)
			}
			else if(is.matrix(x)) {
				ats <- attributes(x)
				w <- !pos
				x[pos] <- 0
				n <- nrow(x)
				TT <- array(1, c(1, n))
				xbar <- (TT %*% x)/(TT %*% w)
				xbar <- t(TT) %*% xbar
				x[pos] <- xbar[pos]
				attributes(x) <- ats
			}
			else {
				ats <- attributes(x)
				x[pos] <- mean(x[!pos])
				attributes(x) <- ats
			}
			frame[[j]] <- x
		}
	}
	frame
}
