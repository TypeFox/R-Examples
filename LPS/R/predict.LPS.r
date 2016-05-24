## Linear Predictor Score prediction
## Author : Sylvain Mareschal <maressyl@gmail.com>
rain <- function(n) rainbow(n, v=0.8)
excl <- function(n) c("#FFCC00", "#333399", "#993333", "#66CC00", "#CC99FF", "#000000", "#FFFFFF", "#99CCFF")[1:n]
predict.LPS <- function(
		object,
		newdata,
		type = c("class", "probability", "score"),
		method = c("Wright", "Radmacher", "exact"),
		threshold = 0.9,
		na.rm = TRUE,
		subset = NULL,
		col.lines = "#FFFFFF",
		col.classes = c("#FFCC00", "#1144CC"),
		plot = FALSE,
		side = NULL,                           # heat.map()
		cex.col = NA,                          # heat.map()
		cex.row = NA,                          # heat.map()
		mai.left = NA,                         # heat.map()
		mai.bottom = NA,                       # heat.map()
		mai.right = 1,                         # heat.map()
		mai.top = 0.1,                         # heat.map()
		side.height = 1,                       # heat.map()
		side.col = NULL,                       # heat.map()
		col.heatmap = heat(),                  # heat.map()
		zlim = "0 centered",                   # heat.map()
		norm = c("rows", "columns", "none"),   # heat.map()
		norm.robust = FALSE,                   # heat.map()
		customLayout = FALSE,                  # heat.map()
		getLayout = FALSE,                     # heat.map()
		...                                    # ignored
	) {
	# Arguments
	type <- match.arg(type)
	method <- match.arg(method)
	norm <- match.arg(norm)
	
	# Plot layout
	if(isTRUE(plot) && !isTRUE(customLayout)) {
		# Layout heights
		heights <- 3
		if(!is.null(side))    heights <- c(lcm(ncol(side)*side.height), heights)
		if(type == "class") { heights <- c(lcm(2*side.height), heights)
		} else              { heights <- c(1, heights)
		}
		
		# Layout widths
		widths <- 1L
		
		# Layout matrix
		if(is.null(side)) { mat <- matrix(c(2,1), ncol=1)
		} else            { mat <- matrix(c(3,1,2), ncol=1)
		}
	} else {
		# No layout call
		mat <- as.integer(NA)
		heights <- as.integer(NA)
		widths <- as.integer(NA)
	}
	
	# Stop returning layout
	if(isTRUE(getLayout)) {
		return(list(mat=mat, heights=heights, widths=widths))
	}
	
	if(is.vector(newdata)) {
		# Use provided scores
		score <- newdata
		
		# Subsetting
		if(!is.null(subset)) score <- score[ subset ]
	} else {
		# Expression matrix
		expr <- as.matrix(newdata[, names(object$coeff) ])
		
		# Row subsetting
		if(!is.null(subset)) expr <- expr[ subset , , drop=FALSE ]
		
		# LPS computation
		score <- apply(t(expr) * object$coeff, 2, sum, na.rm=na.rm)
	}
	
	if(type == "score") {
		# Score only
		out <- score
	} else {
		# Probabilities
		if(method == "Radmacher") {
			# Radmacher (closest mean)
			P1 <- as.double(abs(score - object$means[1]) < abs(score - object$means[2]))
			P2 <- 1 - P1
		} else if(method == "Wright") {
			# Bayesian probabilities (original LPS)
			D1 <- dnorm(score, mean=object$means[1], sd=object$sd[1])
			D2 <- dnorm(score, mean=object$means[2], sd=object$sd[2])
			P1 <- D1 / (D1 + D2)
			P2 <- 1 - P1
		} else if(method == "exact") {
			# Exact probabilities
			if(object$means[1] < object$means[2]) { gLow <- 1L; gHigh <- 2L		
			} else                                { gLow <- 2L; gHigh <- 1L
			}
			D <- list(double(length(score)), double(length(score)))
			for(i in 1:length(score)) {
				D[[ gLow ]][i] <- sum(object$scores[[ gLow ]] >= score[i]) / length(object$scores[[ gLow ]])
				D[[ gHigh ]][i] <- sum(object$scores[[ gHigh ]] <= score[i]) / length(object$scores[[ gHigh ]])
			}
			P1 <- D[[1]] / (D[[1]] + D[[2]])
			P2 <- 1 - P1
		}
	}
	
	if(type == "probability") {
		# Probability matrix
		out <- cbind(P1, P2)
		colnames(out) <- object$classes
		rownames(out) <- names(score)
	} else if(type == "class") {
		# Class
		out <- rep(as.character(NA), length(P1))
		names(out) <- names(score)
		out[ P1 > threshold ] <- object$classes[1]
		out[ P2 > threshold ] <- object$classes[2]
		out[ P1 > threshold & P2 > threshold ] <- paste(object$classes, collapse=" & ")
	}
	
	# Plot
	if(isTRUE(plot)) {
		# Layout call
		if(!isTRUE(customLayout)) {
			layout(mat=mat, heights=heights)
			on.exit(layout(1))
		}
		
		# Order by raw t statistics
		expr <- expr[ order(score) , order(object$t) ]
		
		# Annotation and heatmap
		out <- c(
			list(out),
			heat.map(
				expr = expr,
				customLayout = TRUE,
				cex.col = cex.col,
				cex.row = cex.row,
				mai.left = mai.left,
				mai.bottom = mai.bottom,
				mai.right = mai.right,
				mai.top = 0.1,
				side = side,
				side.height = side.height,
				side.col = side.col,
				col.heatmap = col.heatmap,
				zlim = zlim,
				norm = norm,
				norm.robust = norm.robust
			)
		)
		
		# Print raw t statistics
		axis(side=4, at=(1:ncol(expr) - 1L) / (ncol(expr) - 1L), labels=sprintf("%+0.3f", object$t[ colnames(expr) ]), las=2, cex.axis=out$cex.row, tick=FALSE)
		
		# Sample groups
		if(type != "score") {
			Y1 <- which(head(P1[ order(score) ] > threshold, -1) != tail(P1[ order(score) ] > threshold, -1))
			Y2 <- which(head(P2[ order(score) ] > threshold, -1) != tail(P2[ order(score) ] > threshold, -1))
			abline(v=c(Y1-0.5, Y2-0.5)/(nrow(expr) - 1L), lwd=2, col=col.lines)
			box()
		}
		
		# Feature groups
		x <- tail(which(sort(object$coeff) <= 0), 1)-0.5
		abline(h=x/(ncol(expr) - 1L), lwd=2, col=col.lines)
		
		# Prediction plot
		par(mai=c(0, out$mai.left, mai.top, mai.right))
		if(type == "score") {
			# Score
			plot(x=1:nrow(expr), xlim=c(0.5, nrow(expr)+0.5), y=score[ order(score) ], type="o", pch=16, cex=0.5, xaxs="i", yaxs="i", xpd=NA, xaxt="n", yaxt="n", xlab="", ylab="Score")
		} else if(type == "probability") {
			# Probabilities
			plot(x=1:nrow(expr), xlim=c(0.5, nrow(expr)+0.5), ylim=0:1, y=P1[ order(score) ], type="o", pch=16, cex=0.5, xaxs="i", yaxs="i", xpd=NA, xaxt="n", xlab="", ylab="Probability", col=col.classes[1])
			par(new=TRUE)
			plot(x=1:nrow(expr), xlim=c(0.5, nrow(expr)+0.5), ylim=0:1, y=P2[ order(score) ], type="o", pch=16, cex=0.5, xaxs="i", yaxs="i", xpd=NA, xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col=col.classes[2])
			abline(h=threshold, lty="dotted")
			legend(x="right", inset=0.01, bg="#EEEEEE", lty="solid", pch=16, col=col.classes, legend=object$classes)
		} else if(type == "class") {
			# Classes
			plot(x=NA, y=NA, xlim=c(0.5, nrow(expr)+0.5), ylim=0:1, xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="Class")
			k <- factor(out[[1]][ order(score) ])
			levels(k) <- col.classes
			rect(xleft=(1:nrow(expr))-0.5, xright=(1:nrow(expr))+0.5, ybottom=0, ytop=1, col=as.character(k))
			text(y=0.5, labels=object$classes[ which.min(object$means) ], adj=0, col="white", font=2, x=1)
			text(y=0.5, labels=object$classes[ which.max(object$means) ], adj=1, col="white", font=2, x=nrow(expr))
		}
	}
	
	return(out)
}

