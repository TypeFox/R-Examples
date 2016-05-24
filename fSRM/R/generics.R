#' @title Rerun a fSRM model with new parameters
#' @description
#' Rerun a fSRM model with new parameters
#'
#' @method update fSRM
#' @export

#' @param object A fSRM object.
#' @param evaluate Set to TRUE.
#' @param ... Other parameters (currently not used)
update.fSRM <- function(object, evaluate=TRUE, ...) {
	call <- object$call
    if(is.null(call)) stop("need an fRSM object with call slot")

    extras <- match.call(expand.dots = FALSE)$...

    if(length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for(a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if(any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate) {
        eval(call, parent.frame())
    }
    else call
}



#' @title Predict new cases based on a fitted fSRM model
#' @description
#' Predict new cases based on a fitted fSRM model
#'
#' @method predict fSRM
#' @export

#' @param object A fSRM object.
#' @param newdata A data frame with exactly the same structure as the data frame on which the original fSRM object is based on.
#' @param ... Other parameters (currently not used)
predict.fSRM <- function(object, newdata, ...) {
	
	# TODO: This snippet is also in fSRM function --> refactor
	# Restructure data format from long to wide
	fam0 <- list()
	for (v in c(object$var.id, object$add.variable)) {
		fam0[[v]] <- dcast(newdata[, c(object$var.id, object$actor.id, object$partner.id, object$group.id, object$add.variable)], formula(paste(object$group.id, "~", object$actor.id, "+", object$partner.id)), value.var=v)
		colnames(fam0[[v]])[-1] <- paste(colnames(fam0[[v]])[-1], v, sep="_")
	}

	fam <- merge.rec(fam0, by=object$group.id)
	
	#print(str(object$data))
	#print(str(fam))
	predict(object$fit, newdata=fam)
}


#' @title Plot an fSRM-object, two types
#' @description
#' This function provides two types of plots:
#' 1) Plot the relative variances of an fSRM-object (default)
#' 2) Plot the mean decomposition for each dyad (set \code{means=TRUE})
#'
#' @method plot fSRM
#' @export
#' @importFrom scales percent
#' @importFrom gridExtra grid.arrange
#' @import ggplot2

#' @param x A fSRM object.
#' @param ... Other parameters (currently not used)
#' @param means If FALSE, the relative variances are plotted. If TRUE, the mean structure is plotted.
#' @param bw Black/white plotting?
#' @param onlyStable In case of variance plots: Should only the partitioning of the \emph{stable} variance (without error) be plotted?
#' @examples
#' \dontrun{
#' data(two.indicators)

#' # 4 persons, 1 indicator
#' f4.1 <- fSRM(dep1 ~ actor.id*partner.id | family.id, two.indicators, means=TRUE)
#' f4.1
#' plot(f4.1)
#' plot(f4.1, bw=TRUE)
#' plot(f4.1, means=TRUE)
#' 
#' # 4 persons, 2 indicators
#' f4.2 <- fSRM(dep1/dep2 ~ actor.id*partner.id | family.id, two.indicators, means=TRUE)
#' f4.2
#' plot(f4.2)
#' plot(f4.2, bw=TRUE)
#' plot(f4.2, means=TRUE)
#' }
plot.fSRM <- function(x, ..., means=FALSE, bw=FALSE, onlyStable=FALSE) {
	
	# plot relative percentages
	if (means == FALSE) {
		if (is.null(x$group)) {
			p1 <- plot_relvar(x, bw=bw, onlyStable=onlyStable, ...)
			return(p1)
		} else {
			p1a <- plot_relvar(x, bw=bw, onlyStable=onlyStable, group=x$groupnames[1], ...) + ggtitle(paste("Group", x$groupnames[1]))
			p1b <- plot_relvar(x, bw=bw, onlyStable=onlyStable, group=x$groupnames[2], ...) + ggtitle(paste("Group", x$groupnames[2]))
			grid.arrange(p1a, p1b)
		}
		
	}
	
	# plot mean structure
	if (means == TRUE) {
		if (is.null(x$group)) {
			
			p1 <- plot_meanstruc(x, group="")
			return(p1)
			
		} else {
			p1a <- plot_meanstruc(x, group=x$groupnames[1]) + ggtitle(paste("Group", x$groupnames[1]))
			p1b <- plot_meanstruc(x, group=x$groupnames[2]) + ggtitle(paste("Group", x$groupnames[2]))
			
			# Equate y range
			R1 <- ggplot_build(p1a)
			R2 <- ggplot_build(p1b)
			y_a <- R1$panel$ranges[[1]]$y.range
			y_b <- R2$panel$ranges[[1]]$y.range
			p1a <- p1a + coord_cartesian(ylim=c(min(y_a, y_b), max(y_a, y_b)))
			p1b <- p1b + coord_cartesian(ylim=c(min(y_a, y_b), max(y_a, y_b)))
			grid.arrange(p1a, p1b)
		}
	}
}
