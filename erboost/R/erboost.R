# The code is a modified version of gbm library originally written by Greg Ridgeway. See
# 
# Ridgeway, G. (2007). Generalized boosted models: A guide to the gbm package. R pack-
# age vignette. http://cran.r-project.org/web/packages/gbm.

.onAttach <- function(libname, pkgname)
     packageStartupMessage("Loaded erboost ",
                           utils::packageVersion(pkgname, libname),"\n")


predict.erboost <- function(object,newdata,n.trees,
                        single.tree = FALSE,
                        ...)
{
   if(!is.null(object$Terms))
   {
      x <- model.frame(terms(reformulate(object$var.names)),
                       newdata,
                       na.action=na.pass)
   }
   else
   {
      x <- newdata
   }

   cRows <- nrow(x)
   cCols <- ncol(x)

   for(i in 1:cCols)
   {
      if(is.factor(x[,i]))
      {
         j <- match(levels(x[,i]), object$var.levels[[i]])
         if(any(is.na(j)))
         {
            stop(paste("New levels for variable ",
                        object$var.names[i],": ",
                        levels(x[,i])[is.na(j)],sep=""))
         }
         x[,i] <- as.numeric(x[,i])-1
      }
   }

   x <- as.vector(unlist(x))
   if(missing(n.trees) || any(n.trees > object$n.trees))
   {
      n.trees[n.trees>object$n.trees] <- object$n.trees
      warning("Number of trees not specified or exceeded number fit so far. Using ",paste(n.trees,collapse=" "),".")
   }
   i.ntree.order <- order(n.trees)
      
   predF <- .Call("erboost_pred",
                  X=as.double(x),
                  cRows=as.integer(cRows),
                  cCols=as.integer(cCols),
                  n.trees=as.integer(n.trees[i.ntree.order]),
                  initF=object$initF,
                  trees=object$trees,
                  c.split=object$c.split,
                  var.type=as.integer(object$var.type),
                  single.tree = as.integer(single.tree),
                  PACKAGE = "erboost")
   
   if(length(n.trees)>1) 
   {
      predF <- matrix(predF,ncol=length(n.trees),byrow=FALSE)
      colnames(predF) <- n.trees
      predF[,i.ntree.order] <- predF
   }
   
   if(!is.null(attr(object$Terms,"offset")))
   {
      warning("predict.erboost does not add the offset to the predicted values.")
   }

   return(predF)
}


plot.erboost <- function(x,
                     i.var=1,
                     n.trees=x$n.trees,
                     continuous.resolution=100,
                     return.grid = FALSE,
                     ...)
{
   if(all(is.character(i.var)))
   {
      i <- match(i.var,x$var.names)
      if(any(is.na(i)))
      {
         stop("Plot variables not used in erboost model fit: ",i.var[is.na(i)])
      } else
      {
         i.var <- i
      }
   }

   if((min(i.var)<1) || (max(i.var)>length(x$var.names)))
   {
      warning("i.var must be between 1 and ",length(x$var.names))
   }
   if(n.trees > x$n.trees)
   {
      warning(paste("n.trees exceeds the number of trees in the model, ",x$n.trees,
                    ". Plotting using ",x$n.trees," trees.",sep=""))
      n.trees <- x$n.trees
   }

   if(length(i.var) > 3)
   {
     warning("plot.erboost creates up to 3-way interaction plots.\nplot.erboost will only return the plotting data structure.")
     return.grid = TRUE
   }

   # generate grid to evaluate erboost model
   grid.levels <- vector("list",length(i.var))
   for(i in 1:length(i.var))
   {
      # continuous
      if(is.numeric(x$var.levels[[i.var[i]]]))
      {
         grid.levels[[i]] <- seq(min(x$var.levels[[i.var[i]]]),
                                 max(x$var.levels[[i.var[i]]]),
                                 length=continuous.resolution)
      }
      # categorical or ordered
      else
      {
         grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]],
                                               levels=x$var.levels[[i.var[i]]]))-1
      }
   }
   X <- expand.grid(grid.levels)
   names(X) <- paste("X",1:length(i.var),sep="")

   # evaluate at each data point
   X$y <- .Call("erboost_plot",
                X=as.double(data.matrix(X)),
                cRows=as.integer(nrow(X)),
                cCols=as.integer(ncol(X)),
                i.var=as.integer(i.var-1),
                n.trees=as.integer(n.trees),
                initF=as.double(x$initF),
                trees=x$trees,
                c.splits=x$c.splits,
                var.type=as.integer(x$var.type),
                PACKAGE = "erboost")

   # transform categorical variables back to factors
   f.factor <- rep(FALSE,length(i.var))
   for(i in 1:length(i.var))
   {
      if(!is.numeric(x$var.levels[[i.var[i]]]))
      {
         X[,i] <- factor(x$var.levels[[i.var[i]]][X[,i]+1],
                         levels=x$var.levels[[i.var[i]]])
         f.factor[i] <- TRUE
      }
   }

   if(return.grid)
   {
      names(X)[1:length(i.var)] <- x$var.names[i.var]
      return(X)
   }


   # create the plots
   if(length(i.var)==1)
   {
      if(!f.factor)
      {
         j <- order(X$X1)
         plot(X$X1,X$y,
              type="l",
              xlab=x$var.names[i.var],
              ylab=paste("f(",x$var.names[i.var],")",sep=""),...)
      }
      else
      {
         plot(X$X1,
              X$y,
              xlab=x$var.names[i.var],
              ylab=paste("f(",x$var.names[i.var],")",sep=""),
              ...)
      }
   }
   else if(length(i.var)==2)
   {
      if(!f.factor[1] && !f.factor[2])
      {
         levelplot(y~X1*X2,data=X,
                     xlab=x$var.names[i.var[1]],
                     ylab=x$var.names[i.var[2]],...)
      }
      else if(f.factor[1] && !f.factor[2])
      {
         xyplot(y~X2|X1,data=X,
                xlab=x$var.names[i.var[2]],
                ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                type="l",
                ...)
      }
      else if(!f.factor[1] && f.factor[2])
      {
         xyplot(y~X1|X2,data=X,
                xlab=x$var.names[i.var[1]],
                ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                type="l",
                ...)
      }
      else
      {
         stripplot(y~X1|X2,data=X,
                   xlab=x$var.names[i.var[2]],
                   ylab=paste("f(",x$var.names[i.var[1]],",",x$var.names[i.var[2]],")",sep=""),
                   ...)
      }
   }
   else if(length(i.var)==3)
   {
      i <- order(f.factor)
      X.new <- X[,i]
      X.new$y <- X$y
      names(X.new) <- names(X)

      # 0 factor, 3 continuous
      if(sum(f.factor)==0)
      {
         X.new$X3 <- equal.count(X.new$X3)
         levelplot(y~X1*X2|X3,data=X.new,
                     xlab=x$var.names[i.var[i[1]]],
                     ylab=x$var.names[i.var[i[2]]],...)
      }
      # 1 factor, 2 continuous
      else if(sum(f.factor)==1)
      {
         levelplot(y~X1*X2|X3,data=X.new,
                     xlab=x$var.names[i.var[i[1]]],
                     ylab=x$var.names[i.var[i[2]]],...)
      }
      # 2 factors, 1 continuous
      else if(sum(f.factor)==2)
      {
         xyplot(y~X1|X2*X3,data=X.new,
                type="l",
                xlab=x$var.names[i.var[i[1]]],
                ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                ...)
      }
      # 3 factors, 0 continuous
      else if(sum(f.factor)==3)
      {
         stripplot(y~X1|X2*X3,data=X.new,
                   xlab=x$var.names[i.var[i[1]]],
                   ylab=paste("f(",paste(x$var.names[i.var[1:3]],collapse=","),")",sep=""),
                   ...)
      }
   }
}


erboost.more <- function(object,
                     n.new.trees = 3000,
                     data = NULL,
                     weights = NULL,
                     offset = NULL,
                     verbose = NULL)
{
   if(is.null(object$Terms) && is.null(object$data))
   {
      stop("The erboost model was fit using erboost.fit (rather than erboost) and keep.data was set to FALSE. erboost.more cannot locate the dataset.")
   }
   else if(is.null(object$data) && is.null(data))
   {
      stop("keep.data was set to FALSE on original erboost call and argument 'data' is NULL")
   }
   else if(is.null(object$data))
   {
      m <- eval(object$m, parent.frame())
      Terms <- attr(m, "terms")
      a <- attributes(Terms)

      y <- as.vector(model.extract(m, "response"))
      offset <- model.extract(m,"offset")
      x <- model.frame(delete.response(Terms),
                       data,
                       na.action=na.pass)

      w <- weights
      if(length(w)==0) w <- rep(1, nrow(x))
      w <- w*length(w)/sum(w) # normalize to N

      if(is.null(offset) || (offset==0))
      {
         offset <- NA
      }
      Misc <- NA

      # create index upfront... subtract one for 0 based order
      x.order <- apply(x[1:object$nTrain,,drop=FALSE],2,order,na.last=FALSE)-1
      x <- data.matrix(x)
      cRows <- nrow(x)
      cCols <- ncol(x)
   }
   else
   {
      y           <- object$data$y
      x           <- object$data$x
      x.order     <- object$data$x.order
      offset      <- object$data$offset
      Misc        <- object$data$Misc
      w           <- object$data$w
      cRows <- length(y)
      cCols <- length(x)/cRows
   }

   if(is.null(verbose))
   {
      verbose <- object$verbose
   }
   x <- as.vector(x)

   erboost.obj <- .Call("erboost",
                    Y = as.double(y),
                    Offset = as.double(offset),
                    X = as.double(x),
                    X.order = as.integer(x.order),
                    weights = as.double(w),
                    Misc = as.double(Misc),
                    cRows = as.integer(cRows),
                    cCols = as.integer(cCols),
                    var.type = as.integer(object$var.type),
                    var.monotone = as.integer(object$var.monotone),
                    distribution = as.character(object$distribution$name),
                    n.trees = as.integer(n.new.trees),
                    interaction.depth = as.integer(object$interaction.depth),
                    n.minobsinnode = as.integer(object$n.minobsinnode),
                    shrinkage = as.double(object$shrinkage),
                    bag.fraction = as.double(object$bag.fraction),
                    nTrain = as.integer(object$nTrain),
                    fit.old = as.double(object$fit),
                    n.cat.splits.old = length(object$c.splits),
                    n.trees.old = as.integer(object$n.trees),
                    verbose = as.integer(verbose),
                    PACKAGE = "erboost")
   names(erboost.obj) <- c("initF","fit","train.error","valid.error",
                       "oobag.improve","trees","c.splits")

   erboost.obj$initF         <- object$initF
   erboost.obj$train.error   <- c(object$train.error, erboost.obj$train.error)
   erboost.obj$valid.error   <- c(object$valid.error, erboost.obj$valid.error)
   erboost.obj$oobag.improve <- c(object$oobag.improve, erboost.obj$oobag.improve)
   erboost.obj$trees         <- c(object$trees, erboost.obj$trees)
   erboost.obj$c.splits      <- c(object$c.splits, erboost.obj$c.splits)

   # cv.error not updated when using erboost.more
   erboost.obj$cv.error      <- object$cv.error

   erboost.obj$n.trees        <- length(erboost.obj$trees)
   erboost.obj$distribution   <- object$distribution
   erboost.obj$train.fraction <- object$train.fraction
   erboost.obj$shrinkage      <- object$shrinkage
   erboost.obj$bag.fraction   <- object$bag.fraction
   erboost.obj$var.type       <- object$var.type
   erboost.obj$var.monotone   <- object$var.monotone
   erboost.obj$var.names      <- object$var.names
   erboost.obj$interaction.depth <- object$interaction.depth
   erboost.obj$n.minobsinnode    <- object$n.minobsinnode
   erboost.obj$nTrain            <- object$nTrain
   erboost.obj$Terms             <- object$Terms
   erboost.obj$var.levels        <- object$var.levels
   erboost.obj$verbose           <- verbose

   if(!is.null(object$data))
   {
      erboost.obj$data <- object$data
   }
   else
   {
      erboost.obj$data <- NULL
      erboost.obj$m <- object$m
   }

   class(erboost.obj) <- "erboost"
   return(erboost.obj)
}


erboost.fit <- function(x,y,
                    offset = NULL,
                    misc = NULL,
                    distribution = list(name = "expectile", alpha = 0.5),
                    w = NULL,
                    var.monotone = NULL,
                    n.trees = 3000,
                    interaction.depth = 3,
                    n.minobsinnode = 10,
                    shrinkage = 0.001,
                    bag.fraction = 0.5,
                    train.fraction = 1.0,
                    keep.data = TRUE,
                    verbose = TRUE,
                    var.names = NULL,
                    response.name = NULL)
{
   cRows <- nrow(x)
   cCols <- ncol(x)
   if(is.null(var.names))
   {
      if(is.matrix(x))
      {
         var.names <- colnames(x)
      }
      else if(is.data.frame(x))
      {
         var.names <- names(x)
      }
      else
      {
         var.names <- paste("X",1:cCols,sep="")
      }
   }
   if(is.null(response.name))
   {
      response.name <- "y"
   }

   # check dataset size
   if(cRows*train.fraction*bag.fraction <= 2*n.minobsinnode+1)
   {
      stop("The dataset size is too small or subsampling rate is too large: cRows*train.fraction*bag.fraction <= n.minobsinnode")
   }

   if(nrow(x) != length(y))
   {
      stop("The number of rows in x does not equal the length of y.")
   }

   if(interaction.depth < 1)
   {
      stop("interaction.depth must be at least 1.")
   }

   if(length(w)==0) w <- rep(1, cRows)
   else if(any(w < 0)) stop("negative weights not allowed")
   w <- w*length(w)/sum(w) # normalize to N

   if(any(is.na(y))) stop("Missing values are not allowed in the response, ",
                          response.name)

   if(is.null(offset) || (offset==0))
   {
      offset <- NA
   }
   else if(length(offset) != length(y))
   {
      stop("The length of offset does not equal the length of y.")
   }
   Misc <- NA

   # setup variable types
   var.type <- rep(0,cCols)
   var.levels <- vector("list",cCols)
   for(i in 1:length(var.type))
   {
      if(all(is.na(x[,i])))
      {
         stop("variable ",i,": ",var.names[i]," has only missing values.")
      }
      if(is.ordered(x[,i]))
      {
         var.levels[[i]] <- levels(x[,i])
         x[,i] <- as.numeric(x[,i])-1
         var.type[i] <- 0
      }
      else if(is.factor(x[,i]))
      {
         if(length(levels(x[,i]))>1024)
            stop("erboost does not currently handle categorical variables with more than 1024 levels. Variable ",i,": ",var.names[i]," has ",length(levels(x[,i]))," levels.")
         var.levels[[i]] <- levels(x[,i])
         x[,i] <- as.numeric(x[,i])-1
         var.type[i] <- max(x[,i],na.rm=TRUE)+1
      }
      else if(is.numeric(x[,i]))
      {
         var.levels[[i]] <- quantile(x[,i],prob=(0:10)/10,na.rm=TRUE)
      }
      else
      {
         stop("variable ",i,": ",var.names[i]," is not of type numeric, ordered, or factor.")
      }

      # check for some variation in each variable
      if(length(unique(var.levels[[i]])) == 1)
      {
         warning("variable ",i,": ",var.names[i]," has no variation.")
      }
   }

   nTrain <- as.integer(train.fraction*cRows)

   if(is.character(distribution)) distribution <- list(name=distribution)
   if(!("name" %in% names(distribution)))
   {
      stop("The distribution is missing a 'name' component, for example list(name=\"expectile\")")
   }
   supported.distributions <-
      c("expectile")
   # check potential problems with the distributions
   if(!is.element(distribution$name,supported.distributions))
   {
      stop("Distribution ",distribution$name," is not supported")
   }
   if(!is.numeric(y))
   {
      stop("The response ",response.name," must be numeric. Factors must be converted to numeric")
   }
   if(distribution$name == "expectile")
   {
      if(length(unique(w)) > 1)
      {
         stop("This version of package for the expectile regression lacks a weighted expectile. For now the weights must be constant.")
      }
      if(is.null(distribution$alpha))
      {
         stop("For expectile regression, the distribution parameter must be a list with a parameter 'alpha' indicating the expectile, for example list(name=\"expectile\",alpha=0.95).")
      } else
      if((distribution$alpha<0) || (distribution$alpha>1))
      {
         stop("alpha must be between 0 and 1.")
      }
      Misc <- c(alpha=distribution$alpha)
   }

   # create index upfront... subtract one for 0 based order
   x.order <- apply(x[1:nTrain,,drop=FALSE],2,order,na.last=FALSE)-1

   x <- as.vector(data.matrix(x))
   predF <- rep(0,length(y))
   train.error <- rep(0,n.trees)
   valid.error <- rep(0,n.trees)
   oobag.improve <- rep(0,n.trees)

   if(is.null(var.monotone)) var.monotone <- rep(0,cCols)
   else if(length(var.monotone)!=cCols)
   {
      stop("Length of var.monotone != number of predictors")
   }
   else if(!all(is.element(var.monotone,-1:1)))
   {
      stop("var.monotone must be -1, 0, or 1")
   }
   fError <- FALSE
   erboost.obj <- .Call("erboost",
                    Y=as.double(y),
                    Offset=as.double(offset),
                    X=as.double(x),
                    X.order=as.integer(x.order),
                    weights=as.double(w),
                    Misc=as.double(Misc),
                    cRows=as.integer(cRows),
                    cCols=as.integer(cCols),
                    var.type=as.integer(var.type),
                    var.monotone=as.integer(var.monotone),
                    distribution=as.character(distribution$name),
                    n.trees=as.integer(n.trees),
                    interaction.depth=as.integer(interaction.depth),
                    n.minobsinnode=as.integer(n.minobsinnode),
                    shrinkage=as.double(shrinkage),
                    bag.fraction=as.double(bag.fraction),
                    nTrain=as.integer(nTrain),
                    fit.old=as.double(NA),
                    n.cat.splits.old=as.integer(0),
                    n.trees.old=as.integer(0),
                    verbose=as.integer(verbose),
                    PACKAGE = "erboost")
   names(erboost.obj) <- c("initF","fit","train.error","valid.error",
                       "oobag.improve","trees","c.splits")

   erboost.obj$bag.fraction <- bag.fraction
   erboost.obj$distribution <- distribution
   erboost.obj$interaction.depth <- interaction.depth
   erboost.obj$n.minobsinnode <- n.minobsinnode
   erboost.obj$n.trees <- length(erboost.obj$trees)
   erboost.obj$nTrain <- nTrain
   erboost.obj$response.name <- response.name
   erboost.obj$shrinkage <- shrinkage
   erboost.obj$train.fraction <- train.fraction
   erboost.obj$var.levels <- var.levels
   erboost.obj$var.monotone <- var.monotone
   erboost.obj$var.names <- var.names
   erboost.obj$var.type <- var.type
   erboost.obj$verbose <- verbose
   erboost.obj$Terms <- NULL

   if(keep.data)
   {
      erboost.obj$data <- list(y=y,x=x,x.order=x.order,offset=offset,Misc=Misc,w=w)
   }
   else
   {
      erboost.obj$data <- NULL
   }

   class(erboost.obj) <- "erboost"
   return(erboost.obj)
}


erboost <- function(formula = formula(data),
                distribution = list(name = "expectile", alpha = 0.5),
                data = list(),
                weights,
                var.monotone = NULL,
                n.trees = 3000,
                interaction.depth = 3,
                n.minobsinnode = 10,
                shrinkage = 0.001,
                bag.fraction = 0.5,
                train.fraction = 1.0,
                cv.folds=0,
                keep.data = TRUE,
                verbose = TRUE)
{
   mf <- match.call(expand.dots = FALSE)    
   m <- match(c("formula", "data", "weights", "offset"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf$na.action <- na.pass
   mf[[1]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   Terms <- attr(mf, "terms")

   y <- model.response(mf, "numeric")
   w <- model.weights(mf)
   offset <- model.offset(mf)
   
   var.names <- attributes(Terms)$term.labels
   x <- model.frame(terms(reformulate(var.names)),
                    data,
                    na.action=na.pass)

   # get the character name of the response variable
   response.name <- as.character(formula[[2]])
   if(is.character(distribution)) distribution <- list(name=distribution)

   cv.error <- NULL
   if(cv.folds>1)
   {
      i.train <- 1:floor(train.fraction*length(y))
      cv.group <- sample(rep(1:cv.folds, length=length(i.train)))
      cv.error <- rep(0, n.trees)
      for(i.cv in 1:cv.folds)
      {
         if(verbose) cat("CV:",i.cv,"\n")
         i <- order(cv.group==i.cv)
         erboost.obj <- erboost.fit(x[i.train,,drop=FALSE][i,,drop=FALSE], 
                            y[i.train][i],
                            offset = offset[i.train][i],
                            distribution = distribution,
                            w = ifelse(w==NULL,NULL,w[i.train][i]),
                            var.monotone = var.monotone,
                            n.trees = n.trees,
                            interaction.depth = interaction.depth,
                            n.minobsinnode = n.minobsinnode,
                            shrinkage = shrinkage,
                            bag.fraction = bag.fraction,
                            train.fraction = mean(cv.group!=i.cv),
                            keep.data = FALSE,
                            verbose = verbose,
                            var.names = var.names,
                            response.name = response.name)
         cv.error <- cv.error + erboost.obj$valid.error*sum(cv.group==i.cv)
      }
      cv.error <- cv.error/length(i.train)
   }

   erboost.obj <- erboost.fit(x,y,
                      offset = offset,
                      distribution = distribution,
                      w = w,
                      var.monotone = var.monotone,
                      n.trees = n.trees,
                      interaction.depth = interaction.depth,
                      n.minobsinnode = n.minobsinnode,
                      shrinkage = shrinkage,
                      bag.fraction = bag.fraction,
                      train.fraction = train.fraction,
                      keep.data = keep.data,
                      verbose = verbose,
                      var.names = var.names,
                      response.name = response.name)
   erboost.obj$Terms <- Terms
   erboost.obj$cv.error <- cv.error
   erboost.obj$cv.folds <- cv.folds

   return(erboost.obj)
}


erboost.perf <- function(object,
            plot.it=TRUE,
            oobag.curve=FALSE,
            overlay=TRUE,
            method)
{
   smoother <- NULL

   if((method == "OOB") || oobag.curve)
   {
      if(object$bag.fraction==1)
         stop("Cannot compute OOB estimate or the OOB curve when bag.fraction=1")
      if(all(!is.finite(object$oobag.improve)))
         stop("Cannot compute OOB estimate or the OOB curve. No finite OOB estimates of improvement")
      x <- 1:object$n.trees
      smoother <- loess(object$oobag.improve~x,
                        enp.target=min(max(4,length(x)/10),50))
      smoother$y <- smoother$fitted
      smoother$x <- x

      best.iter.oob <- x[which.min(-cumsum(smoother$y))]
      best.iter <- best.iter.oob
   }

   if(method == "OOB")
   {
      warning("OOB generally underestimates the optimal number of iterations although predictive performance is reasonably competitive. Using cv.folds>0 when calling erboost usually results in improved predictive performance.")
   }

   if(method == "test")
   {
      best.iter.test <- which.min(object$valid.error)
      best.iter <- best.iter.test
   }

   if(method == "cv")
   {
      if(is.null(object$cv.error))
         stop("In order to use method=\"cv\" erboost must be called with cv.folds>1.")
      if(length(object$cv.error) < object$n.trees)
         warning("cross-validation error is not computed for any additional iterations run using erboost.more().")
      best.iter.cv <- which.min(object$cv.error)
      best.iter <- best.iter.cv
   }

   if(!is.element(method,c("OOB","test","cv")))
      stop("method must be cv, test, or OOB")

   if(plot.it)
   {
      par(mar=c(5,4,4,4)+.1)
      ylab <- "Expectile loss"
      if(object$train.fraction==1)
      {
         ylim <- range(object$train.error)
      }
      else
      {
         ylim <- range(object$train.error,object$valid.error)
      }

      plot(object$train.error,
           ylim=ylim,
           type="l",
           xlab="Iteration",ylab=ylab)

      if(object$train.fraction!=1)
      {
         lines(object$valid.error,col="red")
      }
      if(method=="cv")
      {
         lines(object$cv.error,col="green")
      }

      if(oobag.curve)
      {
         if(overlay)
         {
            par(new=TRUE)
            plot(smoother$x,
                 cumsum(smoother$y),
                 col="blue",
                 type="l",
                 xlab="",ylab="",
                 axes=FALSE)
            axis(4,srt=0)
            at <- mean(range(smoother$y))
            mtext(paste("OOB improvement in",ylab),side=4,srt=270,line=2)
            abline(h=0,col="blue",lwd=2)
            if(!is.na(best.iter)) abline(v=best.iter,col="blue",lwd=2)
         }

         plot(object$oobag.improve,type="l",
              xlab="Iteration",
              ylab=paste("OOB change in",ylab))
         lines(smoother,col="red",lwd=2)
         abline(h=0,col="blue",lwd=1)

         abline(v=best.iter,col="blue",lwd=1)
      }
   }

   return(best.iter)
}


relative.influence <- function(object,
                               n.trees)
{
   get.rel.inf <- function(obj)
   {
      lapply(split(obj[[6]],obj[[1]]),sum) # 6 - Improvement, 1 - var name
   }

   temp <- unlist(lapply(object$trees[1:n.trees],get.rel.inf))
   rel.inf.compact <- unlist(lapply(split(temp,names(temp)),sum))
   rel.inf.compact <- rel.inf.compact[names(rel.inf.compact)!="-1"]

   # rel.inf.compact excludes those variable that never entered the model
   # insert 0's for the excluded variables
   rel.inf <- rep(0,length(object$var.names))
   i <- as.numeric(names(rel.inf.compact))+1
   rel.inf[i] <- rel.inf.compact

   return(rel.inf=rel.inf)
}

erboost.loss <- function(y,f,w,offset,dist,baseline)
{
   if(!is.na(offset))
   {
      f <- offset+f
   }
   switch(dist,
          gaussian = weighted.mean((y - f)^2,w) - baseline,
          bernoulli = -2*weighted.mean(y*f - log(1+exp(f)),w) - baseline,
          laplace = weighted.mean(abs(y-f),w) - baseline,
          adaboost = weighted.mean(exp(-(2*y-1)*f),w) - baseline,
          poisson = -2*weighted.mean(y*f-exp(f),w) - baseline,
          stop(paste("Distribution",dist,"is not yet supported for method=permutation.test.erboost")))
}

permutation.test.erboost <- function(object,
                                 n.trees)
{
   # get variables used in the model
   i.vars <- sort(unique(unlist(lapply(object$trees[1:n.trees],
                                       function(x){unique(x[[1]])}))))
   i.vars <- i.vars[i.vars!=-1] + 1
   rel.inf <- rep(0,length(object$var.names))

   if(!is.null(object$data))
   {
      y           <- object$data$y
      os          <- object$data$offset
      Misc        <- object$data$Misc
      w           <- object$data$w
      x <- matrix(object$data$x,ncol=length(object$var.names))
      object$Terms <- NULL # this makes predict.erboost take x as it is
   }
   else
   {
      stop("Model was fit with keep.data=FALSE. permutation.test.erboost has not been implemented for that case.")
   }

   # the index shuffler
   j <- sample(1:nrow(x))
   for(i in 1:length(i.vars))
   {
      x[ ,i.vars[i]]  <- x[j,i.vars[i]]

      new.pred <- predict.erboost(object,newdata=x,n.trees=n.trees)
      rel.inf[i.vars[i]] <- erboost.loss(y,new.pred,w,os,
                                     object$distribution$name,
                                     object$train.error[n.trees])

      x[j,i.vars[i]] <- x[ ,i.vars[i]]
   }

   return(rel.inf=rel.inf)
}


summary.erboost <- function(object,
                        cBars=length(object$var.names),
                        n.trees=object$n.trees,
                        plotit=TRUE,
                        order=TRUE,
                        method=relative.influence,
                        normalize=TRUE,
                        ...)
{
   if(n.trees < 1)
   {
      stop("n.trees must be greater than 0.")
   }
   if(n.trees > object$n.trees)
   {
      warning("Exceeded total number of erboost terms. Results use n.trees=",object$n.trees," terms.\n")
      n.trees <- object$n.trees
   }

   rel.inf <- method(object,n.trees)
   rel.inf[rel.inf<0] <- 0

   if(order)
   {
      i <- order(-rel.inf)
   }
   else
   {
      i <- 1:length(rel.inf)
   }
   if(cBars==0) cBars <- min(10,length(object$var.names))
   if(cBars>length(object$var.names)) cBars <- length(object$var.names)

   if(normalize) rel.inf <- 100*rel.inf/sum(rel.inf)

   if(plotit)
   {
      barplot(rel.inf[i[cBars:1]],
              horiz=TRUE,
              col=rainbow(cBars,start=3/6,end=4/6),
              names=object$var.names[i[cBars:1]],
              xlab="Relative influence",...)
   }
   return(data.frame(var=object$var.names[i],
                     rel.inf=rel.inf[i]))
}

