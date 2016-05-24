
#' Generic pointsOnBoxplot function. Calls pointsOnBoxplot.default or pointsOnBoxplot.formula
#' @param x a vector of numeric values to be passed on
#' @param ... further arguments for the methods, such as a vector of categories 'y' for the default method 
#' @author Jason Waddell
#' @seealso \code{\link{pointsOnBoxplot.default}} \code{\link{pointsOnBoxplot.formula}}
#' @export
#' @examples 
#' # Examples run in the formula and default methods
#' x2 <- runif(50, 0, 10); 
#' table(customRound(x2, roundTo = 0.5))
#' boxplot(x2)
#' pointsOnBoxplot(x2, pch = 19, roundTo = 0.5)
#' 
#' # Set up input data
#' x <- c(sample(1:5, size = 25, replace = TRUE), rpois(25, lambda = 4))
#' colVec <- c(rep("olivedrab", 10), rep("red", 5), rep("goldenrod", 15), 
#'     rep("red", 15), rep("olivedrab", 5))
#' y <- rep(c("Awesome Rats", "Stupid Rats"), each = 25)
#' y2 <- rep(c("Open", "Analytics"), 25)
#' 
#' x2 <- c(1, 2, 2, 3, 3, 1, 1, 1, 4, 5)
#' y3 <- c(rep("A", 5), rep("B", 5))
#' levels(y3) <- c("A", "B", "C")
#' 
#' boxplot(x ~ y, horizontal = TRUE)
#' pointsOnBoxplot(x ~ y, horizontal = TRUE)
#' 
#' boxplot(x ~ y)
#' pointsOnBoxplot(x = x, y = y, col = colVec, pch = 19, cex = 2)
#' 
#' boxplot(x ~ y + y2)
#' pointsOnBoxplot(x ~ y + y2, col = colVec, pch = 19, cex = 2)
pointsOnBoxplot <- function(x, ...){
	UseMethod("pointsOnBoxplot")
}





#' Draw Points on Top of a Boxplot using Appropriate Shifting 
#' @param x vector of numeric values that were used to create boxplots
#' @param y vector of values representing a categorical variable
#' @param totalSpread total spread of point plotting range within a boxplot. Defaults to 0.3 so that points plot between 0.85 and 1.15
#' @param roundTo optional rounding interval. For example, if given roundTo = 0.25, all numeric x values will be rounded to the nearest quarter
#' @param ... further parameters to be passed to the points function
#' @param horizontal logical indicating if the boxplots should be horizontal; default FALSE means vertical boxes.
#' @return points are drawn to the current device
#' @author Jason Waddell
#' @method pointsOnBoxplot default
#' @export
#' @examples 
#' # Examples run in the formula and default methods
#' x2 <- runif(50, 0, 10); 
#' table(customRound(x2, roundTo = 0.5))
#' boxplot(x2)
#' pointsOnBoxplot(x2, pch = 19, roundTo = 0.5)
#' 
#' # Set up input data
#' x <- c(sample(1:5, size = 25, replace = TRUE), rpois(25, lambda = 4))
#' colVec <- c(rep("olivedrab", 10), rep("red", 5), rep("goldenrod", 15), 
#'     rep("red", 15), rep("olivedrab", 5))
#' y <- rep(c("Awesome Rats", "Stupid Rats"), each = 25)
#' y2 <- rep(c("Open", "Analytics"), 25)
#' 
#' x2 <- c(1, 2, 2, 3, 3, 1, 1, 1, 4, 5)
#' y3 <- c(rep("A", 5), rep("B", 5))
#' levels(y3) <- c("A", "B", "C")
#' 
#' boxplot(x ~ y, horizontal = TRUE)
#' pointsOnBoxplot(x ~ y, horizontal = TRUE)
#' 
#' boxplot(x ~ y)
#' pointsOnBoxplot(x = x, y = y, col = colVec, pch = 19, cex = 2)
#' 
#' boxplot(x ~ y + y2)
#' pointsOnBoxplot(x ~ y + y2, col = colVec, pch = 19, cex = 2)
pointsOnBoxplot.default <- function(x = NULL, y = NULL, totalSpread = 0.3, 
		roundTo = NULL, horizontal = FALSE,	...){
	
	
	if(!is.null(roundTo))
		x <- customRound(x, roundTo)
	
	if(is.null(y)){
		
		vec <- x; 
		maxRep <- max(table(vec))	
		xLeft = if(maxRep == 2) 1 - totalSpread/4  else 1 - totalSpread/2 
		xRight = if(maxRep == 2) 1 + totalSpread/4 else 1 + totalSpread/2
		center <- (xLeft+xRight)/2 
		span <- abs(xRight - xLeft)
		space <- span/(maxRep-1)	
		
		xLocations <- numeric(length = length(vec))
		
		for(k in seq_along(table(vec))){
			numRep <- table(vec)[k]
			tempLoc <- findLocations(n = numRep, space = space, center = center)
			
			tempIndex <- which(vec == sort(unique(vec))[k])
			
			for(j in seq_along(tempIndex))
				xLocations[tempIndex[j]] <- sort(tempLoc)[j]	
			
		}
		globalX <- xLocations
		
	} else 
	
	{ # !is.null(y)
		y <- as.factor(y)
		nCat <- nlevels(y)
		
		globalX <- numeric(length = length(x))
		
		for(i in 1:nCat){
			
			index <- which(y == levels(y)[i])
			vec <- x[index]; catSub <- y[index]
			
			# todo add if statement here
			
			maxRep <- max(table(vec))	
			xLeft = if(maxRep == 2) 1 - totalSpread/4 + (i-1) else 1 - totalSpread/2 + (i-1)
			xRight = if(maxRep == 2) 1 + totalSpread/4 + (i-1) else 1 + totalSpread/2 + (i-1)
			center <- (xLeft+xRight)/2 
			span <- abs(xRight - xLeft)
			space <- span/(maxRep-1)	
			
			xLocations <- numeric(length = length(vec))
			
			for(k in seq_along(table(vec))){
				numRep <- table(vec)[k]
				tempLoc <- findLocations(n = numRep, space = space, center = center)
				
				tempIndex <- which(vec == sort(unique(vec))[k])
				
				for(j in seq_along(tempIndex))
					xLocations[tempIndex[j]] <- sort(tempLoc)[j]	
				
			}
			
			globalX[index] <- xLocations
		}
	}
	
	if(horizontal){
		points(x = x, y = globalX,  ...)
	} else {
		points(x = globalX, y = x,  ...)
	}
}






#' Draw Points on Top of a Boxplot using Appropriate Shifting 
#' @param formula a formula of the form a ~ b (+ c, etc.), where a is a numeric vector and all other variables are categorical
#' @param data an optional input parameter of a data.frame containing the variables used in the formula
#' @param ... further arguments to be passed to pointsOnBoxplot.default
#' @param na.action parameter specifying how to handle missingness
#' @author Jason Waddell
#' @method pointsOnBoxplot formula
#' @export
pointsOnBoxplot.formula <- function(formula, data = NULL, ..., na.action = NULL){
	if (missing(formula) || (length(formula) != 3L))
		stop("'formula' missing or incorrect")
	m <- match.call(expand.dots = FALSE)
	if (is.matrix(eval(m$data, parent.frame())))
		m$data <- as.data.frame(data)
	m$... <- NULL
	m$na.action <- na.action
	m[[1L]] <- as.name("model.frame")
	mf <- eval(m, parent.frame())
	response <- attr(attr(mf, "terms"), "response")
	pointsOnBoxplot(mf[[response]], interaction(mf[-response]), ...)
	
	#return(c(mf, response, interaction(mf[-response]) ))
}






