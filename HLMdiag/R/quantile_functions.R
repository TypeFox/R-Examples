#' Constructing a normal Q-Q plot
#'
#' This function will construct a normal Q-Q plot within the
#' \code{ggplot2} framework. It combines the functionality of \code{qqnorm} and
#' \code{qqline}.
#' 
#' 
#' @param x a numeric vector
#' @param line the method used to fit a reference line. If no reference line is
#' desired, leave the value as \code{NULL}. \code{line = "rlm"} will use robust
#' regression to fit a reference line. \code{line = "quantile"} will fit a line
#' through the first and third quartiles. These options are the same as those
#' given for the \code{qqPlot} function in the \code{car} package.
#' @param ... other arguments to be passed to \code{qplot()}
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @keywords hplot
#' @seealso \code{\link{qqnorm}}, \code{\link{qqline}}
ggplot_qqnorm <- function(x, line = NULL, ...){
  p <- ppoints(x)
  theory <- qnorm(p = p)
  yp <- sort(x)
  ret <- qplot(x = theory, y = yp, xlab = "Theoretical Quantiles", 
               ylab = "Sample Quantiles", ...)
  if(!is.null(line)){
    if(line == "quantile"){
      line.info <- qqlineInfo(x)
      ret <- ret + geom_abline(intercept = as.numeric(line.info[1]), 
                               slope = as.numeric(line.info[2]))
    }
    if(line == "rlm"){
      line.info <- coef(rlm(yp ~ theory))
      ret <- ret + geom_abline(intercept = as.numeric(line.info[1]), 
                               slope = as.numeric(line.info[2]))
    }
  }
  return(ret)
}



# Adding a line to a normal quantile-quantile plot
# This function will calculate the equation of the line that can be
# used to help read the normal quantile plot of interest.
qqlineInfo <- function(x){
	yp <- quantile(x, c(0.25, 0.75))
	theory <- qnorm(p = c(0.25, 0.75))
	slope <- diff(yp)/diff(theory)
	intercept <- yp[1L] - slope * theory[1L]
	return(c(intercept = as.numeric(intercept), slope = as.numeric(slope)))
}

#' Overlaying normal Q-Q plots
#'
#' This function will overlay multiple normal Q-Q plots on the same plot. This
#' will be particulary useful when comparing the distribution between groups.
#' In this situation, significantly different slopes would indicate the normal 
#' distributions for the groups do not share a common standard deviation.
#'
#' @param x a numeric vector from which quantiles will be calculated
#' @param group a vector indicating group membership for each value in \code{x}.
#' @param line the method used to fit reference lines. If no reference lines are desired,
#' leave the value as \code{NULL}. \code{line = "rlm"} will use robust regression to fit
#' reference lines. \code{line = "quantile"} will fit lines through the first and third quartiles.
#' @param ... other arguments to be passed to \code{ggplot}
#' @param alpha_point alpha value specified for the points
#' @param alpha_line alpha value specified for the lines
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @references 
#' Hilden-Minton, J. A. (1995) 
#' Mulilevel Diagnostics for Mixed and Hierarchical Linear Models,
#' Ph.D. thesis, University of California Los Angeles. 
#' @export
#' @keywords hplot
group_qqnorm <- function(x, group, line = NULL, alpha_point = 1, alpha_line = 1, ...){
	
  p <- theory <- yp <- intercept <- slope <- NULL # Make codetools happy
  
  # Finding the slope and intercept for each group
	qq.list <- split(x, group)
	if(is.null(line) == FALSE){
		if(line == "quantile"){
			qq.coefs <- lapply(qq.list, qqlineInfo)
		}
		if(line == "rlm"){
			qq.coefs <- lapply(qq.list, function(j){
				p <- ppoints(j)
				theory <- qnorm(p = p)
				yp <- sort(j)
				coef(rlm(yp ~ theory))
				})
		}
	qq.coefs <- do.call('rbind', qq.coefs)
	colnames(qq.coefs) <- c("intercept", "slope")
	qq.coefs <- data.frame(qq.coefs, group = row.names(qq.coefs))
	}
	
	# Defining the quantiles of interest for each group
	group.quant <- data.frame(x = x, group = group)
	group.quant <- ddply(group.quant, .(group), transform, p = ppoints(x), yp = sort(x))
	group.quant <- ddply(group.quant, .(group), transform, theory = qnorm(p = p))
	
	# Plotting
	qq <- ggplot(data = group.quant, mapping = aes(x = theory, y = yp), ...) + 
    geom_point(alpha = alpha_point)
	if(is.null(line) == FALSE){
		qq <- qq + geom_abline(aes(intercept = intercept, slope = slope), alpha = alpha_line, data = qq.coefs)
	}
	qq <- qq + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")
	
	return(qq)
}
