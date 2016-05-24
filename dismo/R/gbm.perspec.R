#
# gbm.perspec version 2.9 April 2007
# J Leathwick/J Elith
#
# takes a gbm boosted regression tree object produced by gbm.step and
# plots a perspective plot showing predicted values for two predictors
# as specified by number using x and y
# values for all other variables are set at their mean by default
# but values can be specified by giving a list consisting of the variable name
# and its desired value, e.g., c(name1 = 12.2, name2 = 57.6)


gbm.perspec <- function(gbm.object, 
    x = 1,                # the first variable to be plotted
    y = 2,                # the second variable to be plotted
    pred.means = NULL,    # allows specification of values for other variables
    x.label = NULL,       # allows manual specification of the x label
    x.range = NULL,       # manual range specification for the x variable
    y.label = NULL,       # and y la seminar committeebel
    z.label = "fitted value", #default z label
    y.range = NULL,       # and the y
    z.range = NULL,       # allows control of the vertical axis
    leg.coords = NULL,  #can specify coords (x, y) for legend
    ticktype = "detailed",# specifiy detailed types - otherwise "simple"
    theta = 55,           # rotation 
    phi=40,               # and elevation 
    smooth = "none",      # controls smoothing of the predicted surface
    mask = FALSE,         # controls masking using a sample intensity model
    perspective = TRUE,   # controls whether a contour or perspective plot is drawn
    ...)                  # allows the passing of additional arguments to plotting routine
                           # useful options include shade, ltheta, lphi for controlling illumination
                           # and cex for controlling text size - cex.axis and cex.lab have no effect
{

	if (! requireNamespace('gbm') ) { stop('you need to install the gbm package to use this function') }
	requireNamespace('splines')
#get the boosting model details

	gbm.call <- gbm.object$gbm.call
	gbm.x <- gbm.call$gbm.x
	n.preds <- length(gbm.x)
	gbm.y <- gbm.call$gbm.y
	pred.names <- gbm.call$predictor.names
	family = gbm.call$family

# and now set up range variables for the x and y preds

	have.factor <- FALSE

	x.name <- gbm.call$predictor.names[x]
	if (is.null(x.label)) {
		x.label <- gbm.call$predictor.names[x]
	}    
	y.name <- gbm.call$predictor.names[y]
	if (is.null(y.label)) {
		y.label <- gbm.call$predictor.names[y]
	}
	data <- gbm.call$dataframe[ , gbm.x, drop=FALSE]  
	n.trees <- gbm.call$best.trees

# if marginal variable is a vector then create intervals along the range

	if (is.vector(data[,x])) {
		if (is.null(x.range)) {
			x.var <- seq(min(data[,x],na.rm=T),max(data[,x],na.rm=T),length = 50)
		} else {
			x.var <- seq(x.range[1],x.range[2],length = 50) 
		}
	} else {
		x.var <- names(table(data[,x]))
		have.factor <- TRUE
	}
	if (is.vector(data[,y])) {
		if (is.null(y.range)) {
			y.var <- seq(min(data[,y],na.rm=T),max(data[,y],na.rm=T),length = 50)
		}  else {y.var <- seq(y.range[1],y.range[2],length = 50)}
	} else {
		y.var <- names(table(data[,y]))
		if (have.factor) {  #check that we don't already have a factor 
			stop("at least one marginal predictor must be a vector!")
		} else {have.factor <- TRUE}
	}
	pred.frame <- expand.grid(list(x.var,y.var))
	names(pred.frame) <- c(x.name,y.name)
	pred.rows <- nrow(pred.frame)
	
#make sure that the factor variable comes first

	if (have.factor) {
		if (is.factor(pred.frame[,2])) {  # swap them about
			pred.frame <- pred.frame[,c(2,1)]
			x.var <- y.var
		}
	}
	
	j <- 3
# cycle through the predictors
# if a non-target variable find the mean
	for (i in 1:n.preds) {
		if (i != x & i != y) {
			if (is.vector(data[,i])) {
				m <- match(pred.names[i],names(pred.means))
				if (is.na(m)) {
					pred.frame[,j] <- mean(data[,i],na.rm=T)
				} else pred.frame[,j] <- pred.means[m]
			}
			if (is.factor(data[,i])) {
				m <- match(pred.names[i],names(pred.means))
				temp.table <- table(data[,i])
				if (is.na(m)) {
					pred.frame[,j] <- rep(names(temp.table)[2],pred.rows)
				} else {
					pred.frame[,j] <- pred.means[m]
				}
				pred.frame[,j] <- factor(pred.frame[,j],levels=names(temp.table))
			}
			names(pred.frame)[j] <- pred.names[i]
			j <- j + 1
		}
	}

#
# form the prediction
#
#assign("pred.frame", pred.frame, pos=1)
	prediction <- gbm::predict.gbm(gbm.object,pred.frame,n.trees = n.trees, type="response")
#assign("prediction", prediction, pos=1, immediate =T)

# model smooth if required

	if (smooth == "model") {
		pred.glm <- glm(prediction ~ ns(pred.frame[,1], df = 8) * ns(pred.frame[,2], df = 8), data=pred.frame,family=poisson)
		prediction <- fitted(pred.glm)
	}

# report the maximum value and set up realistic ranges for z

	max.pred <- max(prediction)
	message("maximum value = ",round(max.pred,2),"\n")

	if (is.null(z.range)) {
		if (family == "bernoulli") {
			z.range <- c(0,1)
		} else if (family == "poisson") {
			z.range <- c(0,max.pred * 1.1)
		} else {
			z.min <- min(data[,y],na.rm=T)
			z.max <- max(data[,y],na.rm=T)
			z.delta <- z.max - z.min
			z.range <- c(z.min - (1.1 * z.delta), z.max + (1.1 * z.delta))
		}
	}

# now process assuming both x and y are vectors

	if (have.factor == FALSE) {

# form the matrix

		pred.matrix <- matrix(prediction,ncol=50,nrow=50)

# kernel smooth if required

		if (smooth == "average") {  #apply a 3 x 3 smoothing average
			pred.matrix.smooth <- pred.matrix
			for (i in 2:49) {
				for (j in 2:49) {
					pred.matrix.smooth[i,j] <- mean(pred.matrix[c((i-1):(i+1)),c((j-1):(j+1))])
				}
			}
			pred.matrix <- pred.matrix.smooth
		}

# mask out values inside hyper-rectangle but outside of sample space

		if (mask) {
			mask.trees <- gbm.object$gbm.call$best.trees
			point.prob <- gbm::predict.gbm(gbm.object[[1]],pred.frame, n.trees = mask.trees, type="response")
			point.prob <- matrix(point.prob,ncol=50,nrow=50)
			pred.matrix[point.prob < 0.5] <- 0.0
		}
#
# and finally plot the result
#
		if (!perspective) {
			image(x = x.var, y = y.var, z = pred.matrix, zlim = z.range)
		} else {
			persp(x=x.var, y=y.var, z=pred.matrix, zlim= z.range,      # input vars
			xlab = x.label, ylab = y.label, zlab = z.label,   # labels
			theta=theta, phi=phi, r = sqrt(10), d = 3,               # viewing pars
			ticktype = ticktype, mgp = c(4,1,0), ...) #
		}
	}
	if (have.factor) {
  # we need to plot values of y for each x
		factor.list <- names(table(pred.frame[,1]))
		n <- 1
#add this bit so z.range still works as expected:
		if (is.null(z.range)) {
			vert.limits <- c(0, max.pred * 1.1)
		} else {
			vert.limits <- z.range
		}

		plot(pred.frame[pred.frame[,1]==factor.list[1],2],
		prediction[pred.frame[,1]==factor.list[1]],
		type = 'l', 
		#ylim = c(0, max.pred * 1.1), 
		ylim = vert.limits,
		xlab = y.label,
		ylab = z.label, ...)
		for (i in 2:length(factor.list)) { 
			#factor.level in factor.list) {
			factor.level <- factor.list[i]
			lines(pred.frame[pred.frame[,1]==factor.level,2],
			prediction[pred.frame[,1]==factor.level], lty = i)
		}

# now draw a legend
		if(is.null(leg.coords)){
			x.max <- max(pred.frame[,2])
			x.min <- min(pred.frame[,2])
			x.range <- x.max - x.min
			x.pos <- c(x.min + (0.02 * x.range),x.min + (0.3 * x.range))

			y.max <- max(prediction)
			y.min <- min(prediction)
			y.range <- y.max - y.min
			y.pos <- c(y.min + (0.8 * y.range),y.min + (0.95 * y.range))
			legend(x = x.pos, y = y.pos, factor.list, lty = c(1:length(factor.list)), bty = "n")
		} else {
			legend(x = leg.coords[1], y = leg.coords[2], factor.list, lty = c(1:length(factor.list)), bty = "n")
		}
	}
}


