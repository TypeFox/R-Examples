confusion <- function (x, ...)
	UseMethod("confusion")

## TODO: implement weights
.confusion <- function (classes, labels, useNA, prior, ...)
{
	## useNA can be "no", "always" or "ifany", but with the later value
	## one takes the risk to get non square matrix if there are NAs in only
	## on vector of classes => change to "no" or "always", depending if there
	## are missing data or not
	if (useNA == "ifany")
		if (any(is.na(classes))) useNA <- "always" else useNA <- "no"
	res <- table(classes, dnn = labels, useNA = useNA)
	total <- sum(res)
	truePos <- sum(diag(res))
	row.freqs <- rowSums(res)	
	
	## Additional data as attributes
	attr(res, "row.freqs") <- row.freqs
	attr(res, "col.freqs") <- colSums(res)
	attr(res, "levels") <- levels(classes[1, ]) # These are *initial* levels!
	## Final levels may differ if there are empty levels, or NAs!
	attr(res, "prior") <- row.freqs # Initial prior are row.freqs
	attr(res, "stats") <- c(total = total, truepos = truePos,
		error = 1 - (truePos / total))
	
	## This is a confusion object, inheriting from table
	class(res) <- c("confusion", "table")
	
	## Do we rescale the confusion matrix?
	if (!missing(prior)) prior(res) <- prior
	
	res
}
	
confusion.default <- function (x, y = NULL, vars = c("Actual", "Predicted"),
labels = vars, merge.by = "Id", useNA = "ifany", prior, ...)
{	
	## If the object is already a 'confusion' object, return it
	if (inherits(x, "confusion")) {
		if (!missing(y))
			warning("you cannot provide 'y' when 'x' is a 'confusion' object")
		## Possibly rescale it
		if (!missing(prior)) prior(x) <- prior		
		return(x)
	}
	
	## Idem if there is a 'confusion' attribute and no y
	conf <- attr(x, "confusion")
	if (!is.null(conf) && missing(y)) {
		## Possibly reweight it
		if (!missing(prior)) prior(conf) <- prior
		return(conf)
	}	

	## Reworks and check arguments
	vars <- as.character(vars)
	if (length(vars) != 2)
		stop("You must provide exactly 2 strings for 'vars'")
	merge.by <- as.character(merge.by)
	
	## There are three possibilities:
	## 1) A single data frame => use vars
	if (missing(y)) {
		## Special case of a data frame or list of two factors: keep as it is
		if (is.list(x) && length(x) == 2 && is.null(vars)) {
			clCompa <- as.data.frame(x)
			if (missing(labels)) labels <- names(clCompa)
		} else {
			x <- as.data.frame(x)
			## Check that vars exist
			if (is.null(names(x)) || !all(vars %in% names(x)))
				stop("'vars' are not among column names of 'x'")
			## Check that levels of two vars do match
			lev1 <- levels(x[[vars[1]]])
			lev2 <- levels(x[[vars[2]]])
			if (!all(lev1 == lev2)) {
				## If difference is only in the order of both levels, reorder #2
				if (!all(sort(lev1) == sort(lev2))) {
					stop("levels of the two variables in 'x' do not match")
				} else x[[vars[2]]] <- factor(as.character(x[[vars[2]]]),
					levels = lev1)
			}
			clCompa <- data.frame(class1 = x[[vars[1]]], class2 = x[[vars[2]]])
		}
	} else { # y is provided
		## 2) Two vectors of factors (must have same length/same levels)
		if (is.factor(x) && is.factor(y)) {
			## Check length match
			if (length(x) != length(x))
				stop("lengths of 'x' and 'y' are not the same")
			## Check levels match
			lev1 <- levels(x)
			lev2 <- levels(y)
			if (!all(lev1  == lev2)) {
				## If difference is only in the order of both levels, reorder #2
				if (!all(sort(lev1)  == sort(lev1))) {
					stop("'x' and 'y' levels do not match")
				} else y <- factor(as.character(y), levels = lev1)
			}
			clCompa <- data.frame(class1 = y, class2 = x)
		} else {
			## 3) Two data frames => merge first, then use vars
			## Check vars exist
			if (is.null(names(x)) || !(vars[1] %in% names(x)))
				stop("first item of 'vars' is not among names of 'x'")
			if (is.null(names(y)) || !(vars[2] %in% names(y)))
				stop("second item of 'vars' is not among names of 'y'")
			## Check that levels of two vars do match
			lev1 <- levels(x[[vars[1]]])
			lev2 <- levels(y[[vars[2]]])
			if (!all(lev1  == lev2)) {
				## If difference is only in the order of both levels, reorder #2
				if (!all(sort(lev1)  == sort(lev2))) {
					stop("levels of the variables in 'x' and 'y' do not match")
				} else x[[vars[2]]] <- factor(as.character(x[[vars[2]]]),
					levels = lev1)
			}
			## Merge data according to merge.by
			clCompa <- merge(y[, c(vars[2], merge.by)],
				x[, c(vars[1], merge.by)], by = merge.by)
			nc <- ncol(clCompa)
			clCompa <- clCompa[, c(nc - 1, nc)]
			## Are there common objects left?
			if (!nrow(clCompa)) stop("no common objects between 'x' and 'y'")
		}
	}
	
	## Construct the confusion object
	if (missing(prior)) {
		.confusion(classes = clCompa, labels = labels, useNA = useNA, ...)
	} else {
		.confusion(classes = clCompa, labels = labels, useNA = useNA,
			prior = prior, ...)
	}
}

confusion.mlearning <- function (x, y = response(x),
labels = c("Actual", "Predicted"), useNA = "ifany", prior, ...) {
	## Check labels
	labels <- as.character(labels)
	if (length(labels) != 2)
		stop("You must provide exactly 2 character strings for 'labels'")
	
	## Extract class2 by using predict on the mlearning object
	class2 <- predict(x, ...)
	
	## Check that both variables are of same length and same levels
	if (length(y) != length(class2))
		stop("lengths of 'x' and 'y' are not the same")
	lev1 <- levels(y)
	lev2 <- levels(class2)
	if (!all(lev1  == lev2)) {
		## If difference is only in the order of both levels, reorder #2
		if (!all(sort(lev1)  == sort(lev2))) {
			stop("levels of 'x' and 'y' do not match")
		} else class2 <- factor(as.character(class2), levels = lev1)
	}
	
	## Construct the confusion object
	if (missing(prior)) {
		.confusion(data.frame(class1 = y, class2 = class2),
			labels = labels, useNA = useNA, ...)
	} else {
		.confusion(data.frame(class1 = y, class2 = class2),
			labels = labels, useNA = useNA, prior = prior, ...)
	}
}

prior <- function (object, ...)
	UseMethod("prior")

prior.confusion <- function (object, ...)
	attr(object, "prior")

`prior<-`<- function (object, ..., value)
	UseMethod("prior<-")

`prior<-.confusion`<- function (object, ..., value)
{
	rsums <- rowSums(object)
	if (!length(value)) { # value is NULL or of zero length
		## Reset prior to original frequencies
		value <- attr(object, "row.freqs")
		res <- round(object / rsums * value)
	
	} else if (is.numeric(value)) { # value is numeric
		
		if (length(value) == 1) { # value is a single number
			if (is.na(value) || !is.finite(value) || value <= 0)
				stop("value must be a finite positive number")
			res <- object / rsums * as.numeric(value)
		
		} else { # value is a vector of numerics
			## It must be either of the same length as nrow(object) or of
			## levels(objects)
			l <- length(value)
			n <- names(value)
			l2 <- levels(object)
			
			if (l == nrow(object)) {
				## If the vector is named, check names and possibly reorder it
				if (length(n))
					if (all(n %in% rownames(object))) {
						value <- value[rownames(object)]
					} else stop("Names of the values do not match levels in the confusion matrix")
			
			} else if (l == length(l2)) {
				## Assume names as levels(object), if they are not provides
				if (!length(n)) names(value) <- n <- l2
				
				## If the vector is named, check names match levels
				if (length(n))
					if (all(n %in% l2)) {
						## Extract levels used in the confusion matrix
						value <- value[rownames(object)]
					} else stop("Names of the values do not match levels in the confusion matrix")

			} else stop("length of 'value' do not match the number of levels in the confusion matrix")	
			
			res <- object / rsums * as.numeric(value)
		}
		
	} else stop("value must be a numeric vector, a single number or NULL")
	
	attr(res, "prior") <- value
	## Take care to rows with no items! => put back zeros!
	res[rsums == 0] <- 0
	res
}

print.confusion <- function (x, sums = TRUE, error.col = sums, digits = 0,
sort = "ward", ...)
{
	## General stats on the confusion matrix
	Stats <- attr(x, "stats")
	Error <- round(Stats["error"] * 100, 1)
	cat(Stats["total"], " items classified with ", Stats["truepos"],
		" true positives (error rate = ", Error, "%)\n",
		sep = "")
	row.freqs <- attr(x, "row.freqs")
	if (!all(attr(x, "prior") == row.freqs)) {
		cat("with initial row frequencies:\n")
		print(row.freqs)
		cat("Rescaled to:\n")
	}
	
	## Print the confusion matrix itself
	X <- x
	class(X) <- "table"

    n <- ncol(X)

	## Do we sort items?
	if (length(sort) && !is.na(sort) && sort != FALSE && sort != "") {
		## Grouping of items
		confuSim <- X + t(X)
		confuSim <- 1 - (confuSim / sum(confuSim) * 2)
		confuDist <- structure(confuSim[lower.tri(confuSim)], Size = n,
			Diag = FALSE, Upper = FALSE, method = "confusion", call = "",
			class = "dist")
		order <- hclust(confuDist, method = sort)$order
		X <- X[order, order]
	}
	
	## Change row and column names to a more compact representation
	nbrs <- formatC(1:ncol(X), digits = 1, flag = "0")
	colnames(X) <- nbrs
	rownames(X) <- paste(nbrs, rownames(X))
	
	## Add sums?
	if (isTRUE(as.logical(sums))) {
		## Calculate error (%)
		ErrorTot <- (1 - (sum(diag(x)) / sum(x))) * 100
		Errors <- as.integer(round(c((1 - diag(X) / apply(X, 1, sum)) * 100,
			ErrorTot), 0))
		## ... and add row and column sums
		X <- addmargins(X, FUN = list(`(sum)` = sum), quiet = TRUE)
	} else Errors <- as.integer(round((1 - diag(X) / apply(X, 1, sum)) * 100, 0))
	
	## Add class errors?
	if (isTRUE(as.logical(error.col))) {
		X <- as.table(cbind(X, `(FNR%)` = Errors))
		dn <- dimnames(X)
		names(dn) <- names(dimnames(x))
		dimnames(X) <- dn
	}
	print(round(X, digits))
	
	## Return the original object invisibly
	invisible(x)
}

plot.confusion <- function (x, y = NULL,
type = c("image", "barplot", "stars", "dendrogram"), stat1 = "Recall",
stat2 = "Precision", names, ...)
{
	if (is.null(y)) type <- match.arg(type)[1] else type <- "stars"
	if (missing(names)) names <- c(substitute(x), substitute(y))
	res <- switch(type,
		image = confusionImage(x, y, ...),
		barplot = confusionBarplot(x, y, ...),
		stars = confusionStars(x, y, stat1 = stat1, stat2 = stat2, names, ...),
		dendrogram = confusionDendrogram(x, y, ...),
		stop("'type' must be 'image', 'barplot', 'stars' or 'dendrogram'"))
	invisible(res)
}

## Representation of the confusion matrix
confusionImage <- function (x, y = NULL, labels = names(dimnames(x)),
sort = "ward", numbers = TRUE, digits = 0, mar = c(3.1, 10.1, 3.1, 3.1),
cex = 1, asp = 1, colfun, ncols = 41, col0 = FALSE, grid.col = "gray", ...)
{
	if (!inherits(x, "confusion"))
        stop("'x' must be a 'confusion' object")

	if (!is.null(y))
		stop("cannot use a second classifier 'y' for this plot")
	
	## Default labels in case none provided
	if (is.null(labels)) labels <- c("Actual", "Predicted")
	
	## Default color function
	## (greens for correct values, reds for errors, white for zero)
	if (missing(colfun)) colfun <- function (n, alpha = 1, s = 0.9, v = 0.9) {
		if ((n <- as.integer(n[1L])) <= 0) return(character(0L))
		## Initial (red) and final (green) colors with white in between
		cols <- c(hsv(h = 0, s = s, v = v, alpha = alpha),   # Red
				  hsv(h = 0, s = 0, v = v, alpha = alpha),   # White
				  hsv(h = 2/6, s = s, v = v, alpha = alpha)) # Green
		## Use a color ramp from red to white to green
		return(colorRampPalette(cols)(n))
	}
	
    n <- ncol(x)

	## Do we sort items?
	if (length(sort) && !is.na(sort) && sort != FALSE && sort != "") {
		## Grouping of items
		confuSim <- x + t(x)
		confuSim <- 1 - (confuSim / sum(confuSim) * 2)
		confuDist <- structure(confuSim[lower.tri(confuSim)], Size = n,
			Diag = FALSE, Upper = FALSE, method = "confusion", call = "",
			class = "dist")
		order <- hclust(confuDist, method = sort)$order
		x <- x[order, order]
	}
	
	## Recode row and column names for more compact display
	colnames(x) <- names2 <- formatC(1:n, digits = 1, flag = "0")
	rownames(x) <- names1 <- paste(rownames(x), names2)
	
	## Transform for better colorization
	## (use a transfo to get 0, 1, 2, 3, 4, 7, 10, 15, 25+)
	confuCol <- x
	confuCol <- log(confuCol + .5) * 2.33
	confuCol[confuCol < 0] <- if (isTRUE(as.logical(col0))) 0 else NA
	confuCol[confuCol > 10] <- 10
	
	## Negative values (in green) on the diagonal (correct IDs)
	diag(confuCol) <- -diag(confuCol)	
	
	## Make an image of this matrix
	opar <- par(no.readonly = TRUE)
	on.exit(par(opar))
	par(mar = mar, cex = cex)
	image(1:n, 1:n, -t(confuCol[nrow(confuCol):1, ]), zlim = c(-10, 10),
		asp = asp, bty = "n", col = colfun(ncols), xaxt = "n", yaxt = "n",
		xlab = "", ylab = "", ...)
	
	## Indicate the actual numbers
	if (isTRUE(as.logical(numbers))) {
		confuTxt <- as.character(round(x[n:1, ], digits = digits))
		confuTxt[confuTxt == "0"] <- ""
		text(rep(1:n, each = n), 1:n, labels = confuTxt)
	}
	
	## Add the grid
	if (length(grid.col)) {
		abline(h = 0:n + 0.5, col = grid.col)
		abline(v = 0:n + 0.5, col = grid.col)
	}
	
	## Add the axis labels
	axis(1, 1:n, labels = names2, tick =  FALSE, padj = 0)
	axis(2, 1:n, labels = names1[n:1], tick =  FALSE, las = 1, hadj = 1)
	axis(3, 1:n, labels = names2, tick =  FALSE)
	axis(4, 1:n, labels = names2[n:1], tick =  FALSE, las = 1, hadj = 0)
	
	## Add labels at top-left
	if (length(labels)) {
		if (length(labels) != 2) stop("You must provide two labels")
		mar[2] <- 1.1
		par (mar = mar, new = TRUE)
		plot(0, 0, type = "n", xaxt = "n", yaxt = "n", bty = "n")
		mtext(paste(labels, collapse = " // "), adj = 0, line = 1, cex = cex)
	}
	
	## Return the confusion matrix, as displayed, in text format
	invisible(x)
}

## Confusion barplot with recall and precision in green bars
## TODO: various bar rescaling possibilities!!!
confusionBarplot <- function (x, y = NULL,
col = c("PeachPuff2", "green3", "lemonChiffon2"), mar = c(1.1, 8.1, 4.1, 2.1),
cex = 1, cex.axis = cex, cex.legend = cex,
main = "F-score (precision versus recall)", numbers = TRUE, min.width = 17, ...)
{
    if (!inherits(x, "confusion"))
        stop("'x' must be a 'confusion' object")
		
	if (!is.null(y))
		stop("cannot use a second classifier 'y' for this plot")
	
	## F-score is 2 * recall * precision / (recall + precision), ... but also
	## F-score = TP / (TP + FP/2 + FN/2). We represent this in a barplot
	TP <- tp <- diag(x)
	FP <- fp <- colSums(x) - tp
	FN <- fn <- rowSums(x) - tp
	## In case we have missing data...
	fn[is.na(tp)] <- 50
    fp[is.na(tp)] <- 50
    tp[is.na(tp)] <- 0

	## We scale these values, so that the sum fp/2 + tp + fn/2 makes 100
	scale <- fp/2 + tp + fn/2
    res <- matrix(c(fp/2 / scale * 100, tp / scale * 100, fn/2 / scale * 100),
		ncol = 3)
    colnames(res) <- c("FPcontrib", "Fscore", "FNcontrib") # In %
    Labels <- names(attr(x, "col.freqs"))
    
	## The graph is ordered in decreasing F-score values
	pos <- order(res[, 2], decreasing = TRUE)
    res <- res[pos, ]
    FN <- FN[pos]
    FP <- FP[pos]
    TP <- TP[pos]
    Labels <- Labels[pos]
    l <- length(FN)
    
	## Plot the graph
	omar <- par("mar")
    on.exit(par(omar))
    par(mar = mar)
    ## The barplot
	barplot(t(res), horiz = TRUE, col = col, xaxt = "n", las = 1, space = 0,
		main = main, ...)
    ## The line that shows where symmetry is
	lines(c(50, 50), c(0, l), lwd = 1)
    
	## Do we add figures into the plot?
	if (isTRUE(as.logical(numbers))) {
		## F-score is written in the middle of the central bar
		xpos <- res[, 1] + res[, 2] / 2
		text(xpos, 1:l - 0.5, paste("(", round(res[, 2]), "%)", sep = ""),
			adj = c(0.5, 0.5), cex = cex)
		
		## Add the number of FP and FN to the left and right, respectively
		text(rep(1, l), 1:l - 0.5, round(FP), adj = c(0, 0.5), cex = cex)
		text(rep(99, l), 1:l - 0.5, round(FN), adj = c(1, 0.5), cex = cex)
	}

    ## Add a legend (if cex.legend is not NULL)
	if (length(cex.legend)) {
		legend(50, l * 1.05, legend = c("False Positives",
			"2*TP (F-score %)", "False Negatives"), cex = cex.legend, xjust = 0.5, yjust = 1,
			fill = col, bty = "n", horiz = TRUE)
	}
    
	## Add axes if cex.axis is not NULL
	if (length(cex.axis))
		axis(2, 1:l - 0.5, tick = FALSE, las = 1, cex.axis = cex.axis,
			labels = Labels)
    
    invisible(res)
}

## TODO: check the box around the legend
confusionStars <- function(x, y = NULL, stat1 = "Recall", stat2 = "Precision",
names, main, col = c("green2", "blue2", "green4", "blue4"), ...)
{
    ## Check objects
	if (!inherits(x, "confusion"))
        stop("'x' must be a 'confusion' object")
    if (!is.null(y) && !inherits(x, "confusion"))
        stop("'y' must be NULL or a 'confusion' object")
	
	## Check stats
	SupportedStats <- c("Recall", "Precision", "Specificity",
        "NPV", "FPR", "FNR", "FDR", "FOR")
    if (!stat1 %in% SupportedStats)
        stop("stats1 must be one of Recall, Precision, Specificity, NPV, FPR, FNR, FDR, FOR")
    if (!stat2 %in% SupportedStats)
        stop("stats2 must be one of Recall, Precision, Specificity, NPV, FPR, FNR, FDR, FOR")
	
	## Choose colors TODO: add a colors argument!
	Blue <- topo.colors(16)
    Green <- terrain.colors(16)
    Stat <- summary(x)
    if (!is.null(y)) { # Comparison of two confusion matrices
		Stat2 <- summary(y)
		Data <- data.frame(Stat2[, stat1], Stat[, stat1], Stat[, stat2],
			Stat2[, stat2])
		Data <- rbind(Data, rep(0, 4))
		colnames(Data) <- paste(rep(c(stat1, stat2), each = 2), c(2, 1, 1, 2))
		if (missing(main)) { # Calculate a suitable title
			if (missing(names)) {
				names <- c(substitute(x), substitute(y))
			} else if (length(names) != 2)
				stop("you must provide two nmaes for the two compared classifiers")
			names <- as.character(names)
			main <- paste("Groups comparison (1 =", names[1], ", 2 =",
			names[2], ")")
		}
		if (length(col) >= 4) {
			col <- col[c(3, 1, 2, 4)]
		} else stop("you must provide four colors for the two statistics and the two classifiers")
	} else { # Single confusion matrix
		Data <- data.frame(Stat[, stat1], Stat[, stat2])
		Data <- rbind(Data, rep(0, 2))
		colnames(Data) <- c(stat1, stat2)
		if (missing(main))
			main <- paste("Groups comparison")
		if (length(col) >= 2) {
			col <- col[1:2]
		} else stop("you must provide two colors for the two statistics")
	}
    rownames(Data) <- c(rownames(Stat), " ")
	## Note: last one is empty box for legend
	
	## Save graph parameters and restore on exit
	opar <- par(no.readonly = TRUE)
	on.exit(par(opar))
	
	## Calculate key location
	kl <- stars(Data, draw.segments = TRUE, scale = FALSE,
		len = 0.8, main = main, col.segments = col, plot = FALSE, ...)
	kcoords <- c(max(kl[, 1]), min(kl[, 2]))
	kspan <- apply(kl, 2, min) / 1.95
	
	## Draw the plot	
	res <- stars(Data, draw.segments = TRUE, scale = FALSE, key.loc = kcoords,
		len = 0.8, main = main, col.segments = col, ...)
	
	## Draw a rectangle around key to differentiate it from the rest
	rect(kcoords[1] - kspan[1], kcoords[2] - kspan[2], kcoords[1] + kspan[1],
		kcoords[2] + kspan[2])
	
	res
}

## Representation of the confusion matrix as a dendrogram
confusionDendrogram <- function (x, y = NULL, labels = rownames(x),
sort = "ward", main = "Groups clustering", ...)
{
    ## Check objects
	if (!inherits(x, "confusion"))
        stop("'x' must be a 'confusion' object")
	if (!is.null(y))
		stop("cannot use a second classifier 'y' for this plot")	
	
    ## Transform the confusion matrix into a symmetric matrix
    ConfuSim <- x + t(x)
    ConfuSim <- 1 - (ConfuSim / sum(ConfuSim) * 2)
	
	
	## Create the structure of a "dist" object
    ConfuDist <- structure(ConfuSim[lower.tri(ConfuSim)], Size = nrow(x),
        Diag = FALSE, Upper = FALSE, method = "confusion", call = "",
        class = "dist")
    
	## method :"ward", "single", "complete", "average", "mcquitty",
	## "median" or "centroid"
    HC <- hclust(ConfuDist, method = as.character(sort)[1])
    plot(HC, labels = labels, main = main, ...)
    
	invisible(HC)
}

## Table with stats per groupe precision, recall, etc
summary.confusion <- function(object, type = "all", sort.by = "Fscore",
decreasing = TRUE, ...)
{
    ## Check objects
	if (!inherits(object, "confusion"))
        stop("'object' must be a 'confusion' object")
	
	## General parameters
    ## Number of groups
    Ngp <- nrow(object)
    
    ## Total : TP + TN + FP + FN
    Tot <- sum(object)
    
    ## TP : True positive item : All items on diagonal
    TP <- diag(object)
    
    ## TP + TN : sum of diagonal = All correct identification
    TP_TN <- sum(TP)
    
    ## TP + FP : sum of columns : Automatic classification
    TP_FP <- colSums(object)
    
    ## TP + FN : sum of rows : Manual classification
    TP_FN <- rowSums(object)
    
    ## FP : False positive items
    FP <- TP_FP - TP    

    ## FN : False negative item
    FN <- TP_FN - TP

    ## TN : True Negative = Total - TP - FP - FN
    TN <- rep(Tot, Ngp) - TP - FP - FN

    ## The 8 basic ratios
    ## Recall = TP / (TP + FN) = 1 - FNR
    Recall <- TP / (TP_FN)

    ## Specificity = TN / (TN + FP) = 1 - FPR
    Specificity <- TN / (TN + FP)

    ## Precision = TP / (TP + FP) = 1 - FDR
    Precision <- TP / (TP_FP)
    
    ## NPV : Negative predicted value = TN / (TN + FN) = 1 - FOR
    NPV <- TN / (TN + FN)
    
    ## FPR : False positive rate = 1 - Specificity = FP / (FP + TN) 
    FPR <- FP / (FP + TN) #1 - Specificity
    
    ## FNR : False negative rate = 1 - Recall = FN / (TP + FN)
    FNR <- FN / (TP + FN) #1 - Recall

    ## FDR : False Discovery Rate = 1 - Precision = FP / (TP + FP)
    FDR <- FP / (TP_FP) #1 - Precision
    
    ## FOR : False omission rate = 1 - NPV = FN / (FN + TN)
    FOR <- FN / (FN + TN) #1 - NPV

    ## The 4 ratios of ratios
    ## LRPT = Likelihood Ratio for Positive Tests = Recall / FPR = Recall /
	## (1 - Specificity)
    LRPT <- Recall / (FPR)
    
    ## LRNT = Likelihood Ratio for Negative Tests = FNR / Specificity =
	## (1 - Recall) / Specificity
    LRNT <- FNR / (Specificity)
    
    ## LRPS : Likelihood Ratio for Positive Subjects = Precision / FOR =
	## Precision / (1 - NPV)
    LRPS <- Precision / (FOR)
    
    ## LRNS : Likelihood Ratio Negative Subjects = FDR / NPV = (1 - Precision) /
	## (1 - FOR)
    LRNS <- FDR / (NPV)
    
	## Additional statistics
    ## F-score = F-measure = F1 score = Harmonic mean of Precision and recall
    Fscore <- 2 * ((Precision * Recall) / (Precision + Recall))
	## F-score is also TP/(TP + (FP + FN) / 2). As such, if TP is null but
	## at least one of FP or FN is not null, F-score = 0. In this case,
	## as both Recall and Precision equal zero, we got NaN => do the correction!
	Fscore[is.nan(Fscore)] <- 0
    
    ## Balanced accuracy = (Sensitivity + Specificity) / 2
    BalAcc <- (Recall + Specificity) / 2

    ## MCC : Matthews correlation coefficient
    Sum1 <- TP + FP
    Sum2 <- TP + FN
    Sum3 <- TN + FP
    Sum4 <- TN + FN
    Denominator <- sqrt(Sum1 * Sum2 * Sum3 * Sum4)
    ZeroIdx <- Sum1 == 0 | Sum2 == 0 | Sum3 == 0 | Sum4 == 0
    if (any(ZeroIdx)) Denominator[ZeroIdx] <- 1
    MCC <- ((TP * TN) - (FP * FN)) / Denominator
    
    ## Chisq : Significance
    Chisq <- (((TP * TN) - (FP * FN))^2 * (TP + TN + FP + FN)) /
		((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    ## Automatic classification - Manual classification
    Auto_Manu <- TP_FP - TP_FN
    
    ## Bray-Curtis dissimilarity index
    Bray <- abs(Auto_Manu) / (sum(TP_FP) + sum(TP_FN))
	
	## General statistics
    ## Error = 1 - Accuracy = 1 - ((TP + TN) / (TP + TN + FP + FN))
    Error <- 1 - (TP_TN / Tot)
	## Micro and macro-averaged F-score
	meanRecall <- sum(Recall, na.rm = TRUE) / Ngp
	meanPrecision <- sum(Precision, na.rm = TRUE) / Ngp
	Fmicro <- 2 * meanRecall * meanPrecision / (meanRecall + meanPrecision)
	Fmacro <- sum(Fscore, na.rm = TRUE) / Ngp
    
    ## Take care to avoid missing data for data frame rownames!
	nms <- names(Fscore)
	nms[is.na(nms)] <- "<NA>"
	names(Fscore) <- nms
	
	## Create a data frame with all results
	res <- data.frame(
	   Fscore = Fscore, Recall = Recall, Precision = Precision,
	   Specificity = Specificity, NPV = NPV, FPR = FPR, FNR = FNR, FDR = FDR,
	   FOR = FOR, LRPT = LRPT, LRNT = LRNT, LRPS = LRPS, LRNS = LRNS, 
	   BalAcc = BalAcc, MCC = MCC, Chisq = Chisq, Bray = Bray, Auto = TP_FP,
	   Manu = TP_FN, A_M = Auto_Manu, TP = TP, FP = FP, FN = FN, TN = TN)

    lev <- rownames(object)
    lev[is.na(lev)] <- "<NA>"
	rownames(res) <- lev
	
	## Sort the table in function of one parameter... by default Fscore
	if (length(sort.by) && sort.by != FALSE) {
		if (sort.by %in% names(res)) {
			ord <- order(res[, sort.by], decreasing = decreasing)
			res <- res[ord, ]
			lev <- lev[ord]
		} else warning("wrong sort.by: ignored and no sort performed")
	}
	
	## What type of results should we return?
	if (length(type) && type != "all") {
		okType <- type[type %in% names(res)]
		if (!length(okType)) stop("Wrong type specified")
		if (length(okType) < length(type))
			warning("one or more wrong types are ignored")
		res <- res[, okType]
		## If the data are reduced to a numeric vector, reinject names
		## and return only this vector
		if (is.numeric(res)) {
			res <- as.numeric(res)
			names(res) <- lev
			attr(res, "stat.type") <- okType
			return(res)
		}
	}

    attr(res, "stats") <- attr(object, "stats")
    attr(res, "stats.weighted") <-
		c(error = Error, Fmicro = Fmicro, Fmacro = Fmacro)
    
    class(res) <- c("summary.confusion", class(res))
	res
}

print.summary.confusion <- function (x, ...)
{
	## General stats on the confusion matrix
	Stats <- attr(x, "stats")
	Error <- round(Stats["error"] * 100, 1)
	cat(Stats["total"], " items classified with ", Stats["truepos"],
		" true positives (error = ", Error, "%)\n",
		sep = "")
	cat("\nGlobal statistics on reweighted data:\n")
	Stats2 <- attr(x, "stats.weighted")
	cat("Error rate: ", round(Stats2["error"] * 100, digits = 1),
		"%, F(micro-average): ", round(Stats2["Fmicro"], digits = 3),
		", F(macro-average): ", round(Stats2["Fmacro"], digits = 3), "\n\n",
		sep = "")
	X <- x
	class(X) <- class(X)[-1]
	print(X)
	
	invisible(x)
}
