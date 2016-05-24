setMethod("rFLXM", signature(model="FLXM", components="list"),
          function(model, components, class, ...) {
            y <- NULL
            for (l in seq_along(components)) {
              yl <- as.matrix(rFLXM(model, components[[l]], ...))
              if (is.null(y))  y <- matrix(NA, nrow = length(class), ncol = ncol(yl))
              y[class == l,] <- yl[class==l,]
            }
            y 
          })

setMethod("rFLXM", signature(model = "FLXMRglm", components="FLXcomponent"),
          function(model, components, ...) {
            family <- model@family
            n <- nrow(model@x)
            if(family == "gaussian") {
              sigma <- components@parameters$sigma
              y <- rnorm(n, mean = components@predict(model@x, ...), sd = sigma)
            }
            else if (family == "binomial") {
              dotarg = list(...)
              if ("size" %in% names(dotarg)) 
                size <- dotarg$size
              else {
                if (nrow(model@y)!=n) stop("no y values - specify a size argument")
                size <- rowSums(model@y)
              }
              parms <- components@parameters
              y <- rbinom(n, prob = components@predict(model@x, ...), size=size)
              y <- cbind(y, size - y)
            }
            else if (family == "poisson") {
              y <- rpois(n, lambda = components@predict(model@x, ...))
              
            }
            else if (family == "Gamma") {
              shape <- components@parameters$shape
              y <- rgamma(n, shape = shape, scale = components@predict(model@x, ...)/shape)
            }
            else stop("family not supported")
            y
          })

setMethod("rFLXM", signature(model = "FLXMRglmfix", components="list"),
          function(model, components, class, ...) {
            k <- sum(model@nestedformula@k)
            n <- nrow(model@x)/k
            y <- matrix(NA, nrow = length(class), ncol = ncol(model@y))
            model.sub <- as(model, "FLXMRglm")
            for (l in seq_len(k)) {
              rok <- (l-1)*n + seq_len(n)
              model.sub@x <- model@x[rok, as.logical(model@design[l,]), drop=FALSE]
              model.sub@y <- model@y[rok,,drop=FALSE]
              yl <- as.matrix(rFLXM(model.sub, components[[l]], ...))
              y[class==l,] <- yl[class==l,]
            }
            y
          })

rmvbinom <- function(n, size, prob) sapply(prob, function(p) rbinom(n, size, p))
rmvbinary <- function(n, center) sapply(center, function(p) rbinom(n, 1, p))

setMethod("rFLXM", signature(model = "FLXMC", components = "FLXcomponent"),
          function(model, components, class, ...) {
            rmvnorm <- function(n, center, cov) mvtnorm::rmvnorm(n = n, mean = center, sigma = cov)
            dots <- list(...)
            FUN <- paste("r", model@dist, sep = "")
            args <- c(n = nrow(model@x), dots, components@parameters)
            return(do.call(FUN, args))
          })

