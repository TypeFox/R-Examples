#' Prepare string for regular expressions 
#' (backslash for all non-letter and non-digit characters)
#' 
#' @export
#' @param text A text string (smooth term label) that needs to be converted 
#' to a regular expression. 
#' @return A regular expression string.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' # Model for illustrating coefficients:
#' m0 <- bam(Y ~ s(Time) + s(Subject, bs='re') 
#' + s(Time, Subject, bs='re'), data=simdat)
#' 
#' # get all coefficients:
#' coef(m0)
#' # to get only the Subject intercepts:
#' coef(m0)[grepl(convertNonAlphanumeric("s(Subject)"), names(coef(m0)))]
#' # to get only the Subject slopes:
#' coef(m0)[grepl(convertNonAlphanumeric("s(Time,Subject)"), names(coef(m0)))]
#'
#' @family Utility functions
convertNonAlphanumeric <- function(text){
	return( gsub("([^a-zA-Z0-9])", "\\\\\\1", as.character(text)) )
}





#' Compare the formulas of two models and return the difference(s). 
#' 
#' @export
#' @import mgcv
#' @param model1 A fitted regression model (using lm, glm, gam, or bam).
#' @param model2 A fitted regression model (using lm, glm, gam, or bam).
#' @return A list with model terms that are not shared by both models.
#' @author Jacolien van Rij
#' @examples 
#' data(simdat)
#'
#' # Fit simple GAM model:
#' gam1 <- bam(Y ~ s(Time), data=simdat)
#' gam2 <- bam(Y ~ Group+s(Time), data=simdat)
#' diff_terms(gam1, gam2)
#'
#' @family Utility functions
diff_terms <- function(model1, model2){
	fa <- c()
	fb <- c()
	get_formula <- function(x){
		if("gam" %in% class(x)){
			return( attr(terms(x$formula), "term.labels") )
		}else{
			stop(sprintf("This function does not work for model of class %s.", class(x)))
		}
	}
	fa <- get_formula(model1)
	fb <- get_formula(model2)
	d1 <- fa[!fa %in% fb]
	d2 <- fb[!fb %in% fa]
	out <- list()
	out[[deparse(substitute(model1))]] <- d1
	out[[deparse(substitute(model2))]] <- d2
	return(out)
}





#' Return the regions in which the smooth is significantly different from zero.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param mean A vector with smooth predictions.
#' @param se A vector with the standard error on the smooth predictions.
#' @param xVals Optional vector with x values for the smooth. 
#' When \code{xVals} is provided, the regions are returned in terms of x-
#' values, otherwise as indices.
#' @param f A number to multiply the \code{se} with, to convert the \code{se} 
#' into confidence intervals. Use 1.96 for 95\% CI and 2.58 for 99\%CI.
#' @param as.vector Logical: whether or not to return the data points as 
#' vector, or not. Default is FALSE, and a list with start and end points will
#'  be returned.
#' @return The function returns a list with start points of each region 
#' (\code{start}) and end points of each region (\code{end}). The logical 
#' \code{xVals} indicates whether the returned values are on the x-scale 
#' (TRUE) or indices (FALSE).
#' @examples
#' data(simdat)
#' 
#' # Use aggregate to calculate mean and standard deviation per timestamp:
#' avg <- aggregate(simdat$Y, by=list(Time=simdat$Time),
#'     function(x){c(mean=mean(x), sd=sd(x))})
#' head(avg)
#' # Note that column x has two values in each row:
#' head(avg$x)
#' head(avg$x[,1])
#' 
#' # Plot line and standard deviation:
#' emptyPlot(range(avg$Time), c(-20,20), h0=0)
#' plot_error(avg$Time, avg$x[,'mean'], avg$x[,'sd'], 
#'    shade=TRUE, lty=3, lwd=3)
#'
#' # Show difference with 0:
#' x <- find_difference(avg$x[,'mean'], avg$x[,'sd'], xVals=avg$Time)
#' # Add arrows:
#' abline(v=c(x$start, x$end), lty=3, col='red')
#' arrows(x0=x$start, x1=x$end, y0=-5, y1=-5, code=3, length=.1, col='red')
#' @author Jacolien van Rij
#'
#' @family Utility functions
find_difference <- function(mean, se, xVals = NULL,f=1,
    as.vector=FALSE) {
    if (length(mean) != length(se)) {
        stop("The vectors mean and se are not equal in length.")
    } else {
        ub <- mean + f*se
        lb <- mean - f*se
        
        n <- which(!(ub >= 0 & lb <= 0))
        if(as.vector){
            if (length(n) == 0) {
                return(rep(FALSE, length(mean)))
            } else {
                out <- rep(FALSE, length(mean))
                out[n] <- TRUE
                return(out)
            }
        }else{
            if (length(n) == 0) {
                return(NULL)
            } else {
                n_prev <- c(NA, n[1:(length(n) - 1)])
                n_next <- c(n[2:length(n)], NA)
                if (!is.null(xVals) & (length(xVals) == length(mean))) {
                    return(list(start = xVals[n[which(is.na(n - n_prev) | (n - n_prev) > 1)]], end = xVals[n[which(is.na(n_next - 
                      n) | (n_next - n) > 1)]], xVals = TRUE))
                } else {
                    return(list(start = n[which(is.na(n - n_prev) | (n - n_prev) > 1)], end = n[which(is.na(n_next - n) | (n_next - 
                      n) > 1)], xVals = FALSE))
                }
                
            }
        }
    }
}
 





#' Return n neighbors around given indices.
#' 
#' @export
#' @import grDevices
#' @import graphics
#' @param el A numeric vector.
#' @param n Number indicating how many points around the elements of \code{el} 
#' need to be selected.
#' @param max The maximum value of the returned elements.
#' @return A vector with the elements of x surrounded by n points.
#' @examples
#' vectorIndices <- 1:1000
#' indOutliers <- c(2,10, 473, 359, 717, 519)
#' fn3 <- find_n_neighbors(indOutliers, n=3, max=max(vectorIndices))
#' fn20 <- find_n_neighbors(indOutliers, n=20, max=max(vectorIndices))
#'
#' # check fn3:
#' print(fn3)
#'
#' # Plot:
#' emptyPlot(c(-10,1000), c(-1,1), h0=0, v0=indOutliers)
#' points(fn3, rep(.5, length(fn3)), pch='*')
#' points(fn20, rep(-.5, length(fn20)), pch='*')
#' @author Jacolien van Rij
#' @family Utility functions
find_n_neighbors <- function(el, n, max) {
    if (length(el) > 0) {
        new_el <- sort(unique(unlist(lapply(el, FUN = function(x) {
            return(sort(unique(c(x, (x - n):x, x:(x + n)))))
        }))))
        new_el <- new_el[new_el >= 1 & new_el <= max]
        return(new_el)
    } else {
        return(NULL)
    }
}





#' Return the value (or the element with the value) closest to zero. 
#' 
#' @export
#' @import stats
#' @param x A numeric vector.
#' @param element Logical: whether or not to return the value (FALSE, default) 
#' or the index (TRUE).
#' @return The value or index of the element closest to zero (absolute 
#' minimum).
#' @examples
#' (test <- seq(-25,25, by=3))
#' min(test[test>0])
#' max(test[test<0])
#' min(abs(test))
#' findAbsMin(test)
#' @author Jacolien van Rij
#' @family Utility functions
findAbsMin <- function(x, element = FALSE) {
    abs_x <- abs(x)
    el <- max(which(abs_x == min(abs_x)))
    if (element) {
        return(el)
    } else {
        return(x[el])
    }
}





#' Return the number of decimal places.
#' 
#' @export
#' @import stats
#' @param x A numeric vector.
#' @return Number of decimals
#' 
#' @author Based on http://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r
#' @family Utility functions
#' 
getDec <- function(x){
	dec <- 0
	dec <- sapply(x, function(a){
		if ((a %% 1) != 0) {
	        nchar(unlist(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE))[[2]])
	    } else {
	        return(0)
	    }
	})
	return(dec)
}





#' Function for rounding and/or segmenting a range.
#' 
#' @export
#' @import stats
#' @param x A numeric vector.
#' @param dec Number of decimal points for rounding using function 
#' \code{\link[base]{round}}. Applied after argument 
#' \code{step}. If NULL (default), no rounding is applied.
#' @param step Round the 
#' @param n.seg Numeric value, number of values in the equally spaced sequence. 
#' Default is 2 (min, max).
#' @return vector, range of equally spaced sequence.
#' @examples
#' zlim <- c(-2.5, 3.01)
#' # does not change anything:
#' getRange(zlim)
#' # create a range of 5 numbers: 
#' # (basically just using seq )
#' getRange(zlim, n.seg=5)
#' # rounds the numbers:
#' getRange(zlim, dec=0)
#' getRange(zlim, n.seg=5, dec=0)
#' # extreme values are multiplications of 5
#' # that contains zlim values:
#' getRange(zlim, step=5)
#' getRange(zlim, step=5, n.seg=5)
#' # similar, but not the same:
#' getRange(zlim, n.seg=5, dec=0)
#' getRange(zlim, n.seg=5, step=1)
#' # combining:
#' getRange(zlim, n.seg=5, step=1, dec=0)
#' 
#' @author Jacolien van Rij
#' @family Utility functions
#' 
getRange <- function(x, dec=NULL, step=NULL, n.seg=2){
	vals <- seq(min(x), max(x), length=n.seg)
    if (!is.null(step)){
        vals <- seq(floor(min(x)/step)*step, ceiling(max(x)/step)*step, length=n.seg)
    }
    if (!is.null(dec)){
        vals <- round(vals, dec)
    }
    return(vals)
}





#' Sort split by grouping predictor.
#' 
#' @export
#' @description Function uses \code{\link[base]{sort.list}} to return indices
#' of of a vector, sorted per group.
#' @param x A vector to be sorted.
#' @param group A names list that specify the different groups to split the 
#' data.
#' @param decreasing Logical: whether or not the sort order should be 
#' decreasing.
#' @return Indices indicating the order of vector x per group.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' range(simdat$Y)
#' ind <- group_sort(simdat$Y, 
#'     group=list(Group=simdat$Group, Trial=simdat$Trial))
#' head(simdat[ind,])
#' @seealso \code{\link[base]{sort.list}}
#' @family Utility functions
group_sort <- function(x, group=NULL, decreasing=FALSE){
    if(is.null(group)){
        return(sort.list(as.numeric(x), decreasing=decreasing))
    }else{
        el <- which(is.na(x))
        tmp <- data.frame(x=x, i=1:length(x))
        x.split <- split(tmp, f=group, drop=TRUE)
        x.split <- as.vector(unlist(lapply(x.split, 
            function(x){
                x[sort.list(as.numeric(x$x), decreasing=decreasing),'i']
            })))
        if(length(el) > 0){
            x.split <- c(x.split, el[!el %in% x.split])
        }
        return(x.split)
    }
}





#' Combine list values as string.
#'
#' @export
#' @param x A vector with the names or numbers of list elements to be combined.
#' @param inputlist A (named) list with information, e.g., graphical parameter settings.
#' @return String
#' @family Utility functions
list2str <- function(x, inputlist) {
    out <- c()
    for(i in x){
        name.i <- NULL
        val.i  <- NULL
        if(is.numeric(i)){
            if(i > 0 & i <= length(inputlist)){
                name.i <- sprintf("el%.0f", i)
                val.i  <- inputlist[[i]]
            }
        }else if(i %in% names(inputlist)){
            name.i <- i
            val.i  <- inputlist[[i]]
        }
        if(! is.null(name.i)){
            if(inherits(val.i, c("numeric", "logical"))){
                out <- c(out, sprintf("%s=c(%s)", name.i, paste(val.i, collapse=",")))
            }else if(inherits(val.i, c("character", "factor"))){
                out <- c(out, sprintf("%s=c(%s)", name.i, paste(sprintf("'%s'",val.i), collapse=",")))
            }else{
                warning(sprintf("Class %s is not supported, element %s is ignored.",
                    class(name.i)[1], name.i))
            } 
        }
    }
    return(paste(out, collapse=", "))
}





#' Return indices of data that were not fitted by the model.
#' 
#' @export
#' @param model A fitted regression model (using lm, glm, gam, or bam).
#' @return The indices of the data that were not fitted by the model.
#' @author Jacolien van Rij
#' @examples 
#' data(simdat)
#'
#' # Add missing values:
#' set.seed(123)
#' simdat[sample(nrow(simdat), size=20),]$Y <- NA
#' # Fit simple linear model:
#' lm1 <- lm(Y ~ Time, data=simdat)
#' na.el <- missing_est(lm1)
#' length(na.el)
#'
#' @family Utility functions
missing_est <- function(model){
	if("lm" %in% class(model)){
		el <- unique(model$na.action)
		if(length(el) > 0){
			return(sort(el))
		}else{
			return(el)
		}
	}else{
		stop("This method is currently only implemented for lm, glm, gam, and bam models.")
	}
}





#' Move a vector n elements forward or backward.
#' 
#' @export
#' @param x A vector.
#' @param n Number indicating how many steps the vector should shift forward 
#' (N > 0) or backward (n < 0).
#' @param na_value The value to replace the empty cells with (e.g., the first 
#' or last points). Defaults to NA.
#' @return A vector with the same length of \code{x}, all moved \code{n} steps.
#' @examples
#' data(simdat)
#' (test <- simdat[1:20,])
#' test$Y.prev <- move_n_point(test$Y)
#' test$change <- test$Y - test$Y.prev
#' test$Y.post5 <- move_n_point(test$Y, n=-5)
#'
#' emptyPlot(nrow(test), range(test$Y))
#' lines(test$Y)
#' lines(test$Y.prev, col='red')
#' lines(test$Y.post5, col='blue')
#'
#' @author Jacolien van Rij
#' @family Utility functions
move_n_point <- function(x, n = 1, na_value = NA) {
    x <- as.vector(x)
    
    if (length(x) > abs(n)) {
        if (n > 0) {
            return(c(rep(na_value, n), x[1:(length(x) - n)]))
        } else {
            return(c(x[(abs(n) + 1):(length(x))], rep(na_value, abs(n))))
        }
    } else if (length(x) == abs(n)) {
        return(NA)
    } else {
        return(NULL)
    }
} 





#' Calculate standard error of the mean.
#' 
#' @export
#' @param x A vector.
#' @return Standard Error of the mean.
#' @family Utility functions
se <- function(x){
    s <- sd(x)
    if (is.na(s)){
        warning("Problem in calculating SD.")
        return(NA)
    }
    else{ return(sd(x)/sqrt(length(x))) } 
}





#' Print a descriptive summary of a data frame.
#' 
#' @export
#' @description The function prints a summary of the data. 
#' Similar to the function \code{\link[utils]{str}}, but easier readable.
#' @param data A data frame.
#' @param print Logical: whether or not to print the summary.
#' @param n Number: maximum number of values being mentioned in the summary.
#' If NULL all values are being mentioned. Defaults to 10.
#' @return Optionally returns a named list with info.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' summary_data(simdat)
#' @family Utility functions
summary_data <- function(data, print=TRUE, n=10){    
    
    labelColumns <- function(x, data){
        out <- NULL
        cn <- ifelse(is.numeric(x), colnames(data)[x], x)
        cl <- class(data[,cn])
        if(inherits(data[,cn],'factor')){
            vals <- sort(unique(as.character(data[,x])))
            n.cur <- length(vals)+1
            if(!is.null(n)){
                n.cur <- n
            }
            if(length(vals)>n.cur){
                out <- sprintf("factor with %d values; set to the value(s): %s, ...", 
                    length(vals),
                    paste( vals[1:n.cur], collapse=", ") )  
            }else{
                out <- sprintf("factor; set to the value(s): %s.", 
                    paste( vals, collapse=", ") )
            }
        }else if(inherits(data[,cn],'numeric')){
            if(length(unique(data[,x])) > 2){
                out <- sprintf("numeric predictor; with %d values ranging from %f to %f.", 
                    length(unique(data[,x])),
                    min(data[,x], na.rm=TRUE), max(data[,x], na.rm=TRUE)) 
            }else{
                out <- sprintf("numeric predictor; set to the value(s): %s.", paste(unique(data[,x]), collapse=", ") )
            }
        }else if(inherits(data[,cn], "matrix")){
            if(length(unique(data[,x])) > 2){
                out <- sprintf("a matrix predictor; with %d values ranging from %f to %f.", 
                    length(unique(data[,x])),
                    min(data[,x], na.rm=TRUE), max(data[,x], na.rm=TRUE))
            }else{
                out <- sprintf("matrix predictor; set to the value(s): %s.", paste(unique(data[,x]), collapse=", ") ) 
            }
        }else{
            vals <- sort(unique(data[,x]))
            n.cur <- length(vals)+1
            if(!is.null(n)){
                n.cur <- n
            }
            if(length(vals)>n.cur){
                out <- sprintf("%s vector with %d values; set to the value(s): %s, ...", 
                    class(data[,cn])[1],
                    length(vals),
                    paste( vals[1:n.cur], collapse=", ") )
            }else{
                out <- sprintf("%s vector; set to the value(s): %s.", 
                    class(data[,cn])[1],
                    paste( vals, collapse=", ") )
            }
        }
        return(out)
    }
    mysummary <- sapply(colnames(data), function(x){labelColumns(x, data)})
    if(print){
        print_summary(mysummary)
    }  
    invisible(mysummary)
}
#' Print a named list of strings, output from \code{\link{summary_data}}.
#' 
#' @export
#' @param sumlist Named list, output of \code{\link{summary_data}}.
#' @param title Optional, text string that will be printed as title.
#' @author Jacolien van Rij
#' @family Utility functions
print_summary <- function(sumlist, title=NULL){
    if(is.null(title)){
        cat("Summary:\n")
    }else{
        cat(title, "\n")
    }
    for(x in names(sumlist)){
        cat("\t*", x, ":", sumlist[[x]], "\n")
    }
}





#' Label timestamps as timebins of a given binsize.
#' 
#' @export
#' @description Function for calculating timebins.
#' @param x Numerical vector with timestamp information.
#' @param binsize Size of the timebin, measured in the same units (often ms) 
#' as \code{x}.
#' @param pos Numerical value that determines the label of the binsize 
#' as proportion of the binsize. A value of 0 will provide the minimum 
#' timestamp within the bin as label, a value of 1 will provide the maximum 
#' value within the bin as label. Defaults to 0.5, the center of the bin.
#' @return Anumerical vector of the same size as \code{x} with timebin 
#' information.
#' @author Jacolien van Rij
#' @examples
#' data(simdat)
#' # grouping Time values in bins:
#' simdat$Timebin <- timeBins(simdat$Time, 200)
#' head(simdat)
#' 
#' # different labels:
#' simdat$Timebin2 <- timeBins(simdat$Time, 200, pos=0)
#' head(simdat)
#' @family Utility functions
timeBins <- function(x, binsize, pos=.5){
  return( ( floor( x / binsize )+ pos ) * binsize ) 
}





