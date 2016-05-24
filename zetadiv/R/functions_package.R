.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is zetadiv 0.1")
  packageStartupMessage("Warning: package \"zetadiv\" was built under R.3.0.2")
}


#' Zeta diversity decline
#'
#' Computes zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), from 1 to 'order'.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param orders  Range of number of assemblages or sites for which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta-diversity is computed for each number of assemblages or sites.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @param sd  Boolean value (TRUE or FALSE) indicating if the standard deviation of each zeta-diversity value must be plotted.
#' @return \code{Zeta.decline} returns a list containing the following components:
#' @return \item{zeta.order}{The number of assemblages or sites for which the zeta-diversity was computed.}
#' @return \item{zeta.val}{The zeta-diversity values.}
#' @return \item{zeta.val.sd}{The zeta-diversity standard deviation values.}
#' @return \item{zeta.exp}{Object of class "\code{lm}", containing the output of the exponential regression.}
#' @return \item{zeta.pl}{Object of class "\code{lm}", containing the output of the power law regression.}
#' @return \item{aic}{AIC values for \code{zeta.exp} and \code{zeta.pl}.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.order}}
#' @examples
#' data(BCI.spec.coarse)
#' xy <- BCI.spec.coarse[1:2]
#' data.spec <- BCI.spec.coarse[3:308]
#' 
#' dev.new(width = 12, height = 4)
#' zeta <- Zeta.decline(data.spec, orders = 1:5, sam = 100)
#' zeta
#' 
#' ##########
#' 
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' 
#' dev.new(width = 12, height = 4)
#' zeta.species <- Zeta.decline(data.species, orders = 1:5, sam = 100)
#' zeta.species
#' 
#' @export
#'
Zeta.decline <- function(data.spec, orders = 1:10, sam = 1000, plot = TRUE, sd = TRUE){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.",  sep = ""))
  }
  if(max(orders)>dim(data.spec)[1]){
  	stop("Wrong value for \"orders\": the maximum value must be equal or lower than the number of sites.")
  }
    
  x <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()
  
  for(j in orders){
    if(choose(x, j)>sam){
      u <- rep(NA, sam)
      for(z in 1:sam){
        samp <- sample(1:x, j, replace = FALSE)
        u[z] <- sum(apply(data.spec[samp, ], 2, prod))
      }
    }else{
      u <- rep(NA, choose(x, j))
      samp <- combn(1:x, j)
      for(z in 1:dim(samp)[2]){
        u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
      }
    }
        
    zeta.val[j]<-mean(u)
    zeta.val.sd[j] <- sd(u)
  }
  
  
  ##create a single list for output
  zeta <- list()
  zeta$zeta.order <- orders
  zeta$zeta.val <- zeta.val
  zeta$zeta.val.sd <- zeta.val.sd
    
  
  ##regression - exponential
  zeta.val.log <- log10(zeta.val)
  zeta.val.log[which(is.infinite(zeta.val.log))] <- NA
  zeta.exp <- lm(zeta.val.log ~ c(orders), na.action = na.omit)
  zeta$zeta.exp <- zeta.exp
  
  
  ##regression - power law
  zeta.pl <- lm(zeta.val.log ~ log10(c(orders)), na.action = na.omit)
  zeta$zeta.pl <- zeta.pl
  
  zeta$aic <- AIC(zeta$zeta.exp, zeta$zeta.pl)
  
  
  ##Plot zeta and regressions
  if(plot == TRUE){
    par(mfrow = c(1, 3))
    if (sd == TRUE){
      plot(orders, zeta.val, xlab = "Zeta order", ylab = "Zeta-diversity", pch = 20, ylim = c(0, zeta.val[1] + zeta.val.sd[1]), main = "Zeta diversity decline")
      lines(orders, zeta.val)
      ##sd of zeta as error bars
      for(i in orders){
        suppressWarnings(arrows(i, zeta.val[i], i, zeta.val[i] + zeta.val.sd[i], angle = 90, length = 0.1))
        suppressWarnings(arrows(i, zeta.val[i], i, zeta.val[i] - zeta.val.sd[i], angle = 90, length = 0.1))
      }
    }else{
      plot(orders, zeta.val, xlab = "number of sites", ylab = "zeta-diversity", pch = 20, ylim = c(0, zeta.val[1]))
      lines(orders, zeta.val)
    }
    plot(orders, zeta.val, log = "y", pch = 20, xlab = "Zeta order", ylab = "Zeta-diversity", main = "Exponential regression")
    lines(orders, 10^predict.lm(zeta.exp, data.frame(orders)))
    plot(orders, zeta.val, log = "xy", pch = 20, xlab = "Zeta order", ylab = "Zeta-diversity", main = "Power law regression")
    lines(orders, 10^predict.lm(zeta.pl, data.frame(orders)))
  }
  
  return(zeta)

}




#' Zeta diversity for a specific number of assemblages or sites
#'
#' Computes zeta diversity, the number of species shared by multiple assemblages, for a specific order (number of assemblages or sites).
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta-diversity is computed.
#' @return \code{Zeta.order}  returns a list containing the following components:
#' @return \item{zeta.order}{The number of assemblages or sites for which the zeta-diversity was computed.}
#' @return \item{zeta.val}{The zeta-diversity values.}
#' @return \item{zeta.val.sd}{The zeta-diversity standard deviation values.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline}}
#' @examples
#' data(BCI.spec.coarse)
#' xy <- BCI.spec.coarse[1:2]
#' data.spec <- BCI.spec.coarse[3:308]
#' 
#' zeta <- Zeta.order(data.spec, order = 3, sam=100)
#' zeta
#' 
#' ##########
#' 
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' 
#' zeta.species <- Zeta.order(data.species, order = 3, sam=100)
#' zeta.species
#' 
#' @export
Zeta.order <- function(data.spec, order = 1, sam = 1000){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
  	stop("Wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
      
  x <- dim(data.spec)[1]
    
  if(choose(x, order)>sam){
    u <- rep(NA, sam)
    for(z in 1:sam){
      samp <- sample(1:x, order, replace = FALSE)
      u[z] <- sum(apply(data.spec[samp, ], 2, prod))
    }
  }else{
    u <- rep(NA, choose(x, order))
    samp <- combn(1:x, order)
    for(z in 1:dim(samp)[2]){
      u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
    }
  }
  
  
  zeta.val <- mean(u)
  zeta.val.sd <- sd(u)
  
  
  
  zeta.order <- list()
  zeta.order$val <- zeta.val
  zeta.order$sd <- zeta.val.sd
    
  return(zeta.order)

}




#' Sensitivity analysis for the sample size of zeta
#'
#' Computes zeta diversity for a given order (number of assemblages or sites) for a range of sample sizes, to assess the sensitivity to this parameter.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam.min  Minimum number of samples for which the zeta-diversity is computed.
#' @param sam.max  Maximum number of samples for which the zeta-diversity is computed.
#' @param sam.incr  Increment of the number of samples for which the zeta-diversity is computed. It also serves as the initial number.
#' @param reps  Number of replicates of zeta-diversity computations for each sample size
#' @param display  Boolean value (TRUE or FALSE) indicating if the current value of the sample size must be displayed. Acts as a counter.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted as a boxplot of the zeta-diversity distributions for each sample size
#' @param notch  Boolean value (TRUE or FALSE) indicating if the notches must be plotted in the boxplot.
#' @return \code{Zeta.sam.sensitivity} returns a matrix with \code{(sam.max-sam.min)/sam.incr} columns and \code{reps} rows.
#' @details Note that the execution of \code{Zeta.sam.sensitivity} can be quite lengthy, because of the number of replicates needed.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline}}, \code{\link{Zeta.order}}
#' @examples
#' \donttest{
#' #Note that the sensitivity analyses in the following two examples are quite long to run, 
#' #typically around 10 minutes for the first example and 1-2 minutes for the second.
#'
#' data(BCI.spec.coarse)
#' xy <- BCI.spec.coarse[1:2]
#' data.spec <- BCI.spec.coarse[3:308]
#' 
#' dev.new()
#' zeta.sens <- Zeta.sam.sensitivity(data.spec, order = 3, sam.min = 250, sam.max = 1000,
#'     sam.incr = 250, reps = 20, display = TRUE, plot = TRUE, notch = TRUE)
#' zeta.sens
#' 
#' ##########
#' 
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' 
#' dev.new()
#' zeta.sens.species <- Zeta.sam.sensitivity(data.species, order = 3, sam.min = 250, sam.max = 1000, 
#'     sam.incr = 250, reps = 20, plot = TRUE, notch = TRUE)
#' zeta.sens.species
#' }
#'
#' @export
Zeta.sam.sensitivity <- function(data.spec, order = 1, sam.min = 500, sam.max = 2000, sam.incr = 500, reps = 20, display = TRUE, plot = TRUE, notch = TRUE){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
  	stop("Wrong value for \"order\": it must be equal or lower than the number of sites.")
  	}
  
  x <- dim(data.spec)[1]
  zeta.val <- matrix(NA, reps, length(seq(sam.min, sam.max, sam.incr)))
  zeta.val.sd <- matrix(NA, reps, length(seq(sam.min, sam.max, sam.incr)))

  i.sam <- 0
  for(sam in seq(sam.min, sam.max, sam.incr)){
    if(display == TRUE){print(sam)}
    i.sam <- i.sam + 1
    for(i in 1:reps){
      u <- rep(NA, sam)
      for(j in 1:order){
        for(z in 1:sam){
          samp <- sample(1:x, j, replace = FALSE)
          u[z] <- sum(apply(data.spec[samp, ], 2, prod))
        }
        zeta.val[i, i.sam]<-mean(u)
        zeta.val.sd[i, i.sam] <- sd(u)
      }
    }
  }
  
  if (plot == TRUE){
    boxplot(zeta.val, notch = notch, names = seq(sam.min, sam.max, sam.incr), xlab = "number of samples", ylab = paste("zeta ", order, sep = ""), main = "Distributions of zeta diversities for different number of samples")
  }
  
  zeta.sens <- zeta.val
  
  return(zeta.sens)
  
}






#' Linear regression of zeta for a set of environmental variables and distance
#'
#' Computes the linear regression of zeta diversity for a given order (number of assemblages or sites) against a set of environmental variables and distances between sites.
#' @param xy  Site coordinates
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param data.env  Sites-by-variable data frame, with sites as rows and environmental variables as columns.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta-diversity is computed.
#' @param confint.level  Percentage for the confidence intervals of the coefficient from the linear regression.
#' @param standard  Boolean parameter indicating if the spatial distances and differences in environmental variables should be standardized between 0 and 1.
#' @param method  Indicates how to combine the pairwise differences and distances for more than 3 sites. Method can be "\code{mean}" or "\code{max}".
#' @param silent  Boolean parameter indicating if warnings must be printed
#' @return \code{Zeta.lm} returns a list containing the following components:
#' @return \item{model}{Object of class "\code{lm}", containing the output of the linear regression of zeta over the environmental variables.}
#' @return \item{confint}{The confidence intervals for the coefficients.}
#' @return \item{vif}{The variance inflation factors for all the variables.}
#' @details The environmental variables can be numeric or factorial. 
#' @details If \code{order = 1}, the variables are used as such in the regression, and factorial variables must be dummy for the output of the regression to be interpretable. 
#' @details For numeric variables, if \code{order>1} the pairwise difference between sites is computed and combined according to \code{method}. For factorial variables, the distance corresponds to the number of unique values over the number of assemblages of sites specified by \code{order}. 
#' @details If xy = NULL, \code{Zeta.lm} only uses environmental variables in the regression. Otherwise, it also computes and uses euclidian distance (average or maximum distance between multiple sites, depending on the parameters \code{method}) as an explanatory variable.
#' @details Zeta is regressed against the differences of values of the environmental variables divided by the maximum difference for each variable, to be rescaled between 0 and 1. If \code{xy != NULL}, distances between sites are also divided by the maximum distance.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline}}, \code{\link{Zeta.order}}, \code{\link{Zeta.gam}}
#' @examples
#' data(BCI.spec.coarse)
#' xy <- BCI.spec.coarse[1:2]
#' data.spec <- BCI.spec.coarse[3:308]
#' data(BCI.env.coarse)
#' data.env <- BCI.env.coarse[10:15]
#' 
#' zeta.lm <- Zeta.lm(data.spec, data.env, sam = 100, order = 3)
#' zeta.lm
#' dev.new()
#' plot(zeta.lm$model)
#' 
#' zeta.lm <- Zeta.lm(data.spec, data.env, xy, sam = 100, order = 3)
#' zeta.lm
#' dev.new()
#' plot(zeta.lm$model)
#' 
#' ##########
#' 
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' data(Marion.env)
#' data.env <- Marion.env[3:4]
#' 
#' zeta.lm.species <- Zeta.lm(data.species, data.env, sam = 100, order = 3)
#' zeta.lm.species
#' dev.new()
#' plot(zeta.lm.species$model)
#' 
#' zeta.lm.species <- Zeta.lm(data.species, data.env, xy, sam = 100, order = 2)
#' zeta.lm.species
#' dev.new()
#' plot(zeta.lm.species$model)
#' 
#' @export
Zeta.lm <- function(data.spec, data.env, xy = NULL, order = 1, sam = 1000, confint.level = 0.95, standard = TRUE, method = "mean", silent = FALSE){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(class(data.env) != "data.frame"){
  	stop(paste(deparse(substitute(data.env)), " is a ", class(data.env), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
  	stop("Wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  
  if (length(setdiff(sapply(data.env, class), c("factor", "numeric")))>0){
  	stop("Variables must be numeric or factor")
  }
  if(order == 1 & !is.null(xy)){
  	stop("Cannot include distance for order = 1")
  }
  if(silent == FALSE & order == 1 & length(intersect(sapply(data.env, class), "factor"))>0){
  	warning("factor variables should be dummy for order = 1")
  } 
  
  x <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()
  
  
  if(choose(x, order)>sam){
    u <- rep(NA, sam)
    data.var <- as.data.frame(matrix(NA, sam, dim(data.env)[2]))
    distance <- rep(NA, sam)
    
    for(z in 1:sam){
      samp <- sample(1:x, order, replace = FALSE)
      u[z] <- sum(apply(data.spec[samp, ], 2, prod))
      
      if(order == 1){
        data.var[z, ] <- data.env[samp, ]
      }else {
        
        
        fac <- which(sapply(data.env, class) == "factor")
        num <- which(sapply(data.env, class) == "numeric")
        
        if(order>2){
          if(length(num)>1){
            if(method == "mean"){
              data.var[z, num] <- apply(apply(data.env[samp, num], 2, dist), 2, mean)  ##computing directly the mean measure of euclidian distance
            }else if(method == "max"){
              data.var[z, num] <- apply(apply(data.env[samp, num], 2, dist), 2, max)
            }else{
              stop("Error: method for distance does not exist")
            }
          }else if(length(num)>0){
            data.var[z, num] <- mean(dist(data.env[samp, num]))
          }
        }else{
          if(length(num)>1){
            data.var[z, num] <- apply(data.env[samp, num], 2, dist)
          }else if(length(num)>0){
            data.var[z, num] <- dist(data.env[samp, num])
          }
        }
        if(length(fac)>1){
          data.var[z, fac] <- apply(data.env[samp, fac], 2, function(x){length(unique(x))}) - 1
        }else if (length(fac)>0){
          data.var[z, fac] <- length(unique(data.env[samp, fac])) - 1
        }
  
        if (!is.null(xy)){
          if(method == "mean"){
            distance[z] <- mean(dist(xy[samp, ]))
          }else if(method == "max"){
            distance[z] <- max(dist(xy[samp, ]))
          }else{
            stop("Error: method for distance does not exist")
          }
        }
      }
    }
    
  }else{
    u <- rep(NA, choose(x, order))
    data.var <- as.data.frame(matrix(NA, choose(x, order), dim(data.env)[2]))
    distance <- rep(NA, choose(x, order))
    samp <- combn(1:x, order)
    
    for(z in 1:dim(samp)[2]){
      u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
      
      if(order == 1){
        data.var[z, ] <- data.env[samp[, z], ]
      }else {
        
        
        fac <- which(sapply(data.env, class) == "factor")
        num <- which(sapply(data.env, class) == "numeric")
        
        if(order>2){
          if(length(num)>1){
            if(method == "mean"){
              data.var[z, num] <- apply(apply(data.env[samp[, z], num], 2, dist), 2, mean)  ##computing directly the mean measure of euclidian distance
            }else if(method == "max"){
              data.var[z, num] <- apply(apply(data.env[samp[, z], num], 2, dist), 2, max)
            }else{
              stop("Error: method for distance does not exist")
            }
          }else if(length(num)>0){
            data.var[z, num] <- mean(dist(data.env[samp[, z], num]))
          }
        }else{
          if(length(num)>1){
            data.var[z, num] <- apply(data.env[samp[, z], num], 2, dist)
          }else if(length(num)>0){
            data.var[z, num] <- dist(data.env[samp[, z], num])
          }
        }
        if(length(fac)>1){
          data.var[z, fac] <- apply(data.env[samp[, z], fac], 2, function(x){length(unique(x))}) - 1
        }else if (length(fac)>0){
          data.var[z, fac] <- length(unique(data.env[samp[, z], fac])) - 1
        }
  
        if (!is.null(xy)){
          if(method == "mean"){
            distance[z] <- mean(dist(xy[samp[, z], ]))
          }else if(method == "max"){
            distance[z] <- max(dist(xy[samp[, z], ]))
          }else{
            stop("Error: method for distance does not exist")
          }
        }
      }
    }
  }
      
  ##rescale the environmental variables or distances between 0 and 1
  if(standard == TRUE){
    if(order>1){
      data.var <- data.var / matrix(rep(apply(data.var, 2, max), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.env)[2], byrow = T)
    }else{
      num <- which(sapply(data.env, class) == "numeric")
      if(length(num)>1){
        data.var[, num] <- (data.var[, num] - matrix(rep(apply(data.var[, num], 2, min), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.var[, num])[2], , byrow = T)) / matrix(rep((apply(data.var[, num], 2, max) - apply(data.var[, num], 2, min)), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.var[, num])[2], byrow = T)
      }else{
        data.var[, num] <- (data.var[, num] - min(data.var[, num])) / (max(data.var[, num]) - min(data.var[, num]))
      }
    }
    if(!is.null(xy)){
    	distance <- distance / max(distance)
    }
  }
  names(data.var) <- names(data.env)
    
  if(is.null(xy)){data.tot <- data.var}else{data.tot <- cbind(data.var, distance)}
  
  zeta.val <- u
  
  zeta.lm.model <- glm(zeta.val ~ ., data = data.tot)
  zeta.lm.confint <- suppressMessages(confint(zeta.lm.model, level = confint.level))
  if(dim(data.env)[2]>1){zeta.lm.vif <- car::vif(zeta.lm.model)}else{zeta.lm.vif <- NA}
  
  zeta.lm <- list()
  zeta.lm$model <- zeta.lm.model
  zeta.lm$confint <- zeta.lm.confint
  zeta.lm$vif <- zeta.lm.vif
  
  return(zeta.lm)
  
  
  
}


#' Generalized additive model of zeta for a set of environmental variables and distance
#'
#' Computes a generalized additive model of zeta diversity for a given order (number of assemblages or sites) against a set of environmental variables and distances between sites.
#' @param xy  Site coordinates
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param data.env  Sites-by-variable data frame, with sites as rows and environmental variables as columns.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta-diversity is computed.
#' @param standard  Boolean parameter indicating if the spatial distances and differences in environmental variables should be standardized between 0 and 1.
#' @param method  Indicates how to combine the pairwise differences and distances for more than 3 sites. Method can be "\code{mean}" or "\code{max}".
#' @return \code{Zeta.gam} returns an object of class "\code{gam}", containing the output of the generalized additive model of zeta over the environmental variables.
#' @details If \code{order = 1}, the environmental variables must be numeric and are used as such in the regression. If \code{order>1}, the environmental variables can be numeric or factorial. 
#' @details For numeric variables, if \code{order>1} the pairwise difference between sites is computed and combined according to \code{method}. For factorial variables, the distance corresponds to the number of unique values over the number of assemblages of sites specified by \code{order}. 
#' @details If xy = NULL, \code{Zeta.gam} only uses environmental variables in the regression. Otherwise, it also computes and uses euclidian distance (average or maximum distance between multiple sites, depending on the parameters \code{method}) as an explanatory variable.
#' @details Zeta is regressed against the differences of values of the environmental variables divided by the maximum difference for each variable, to be rescaled between 0 and 1. If \code{xy != NULL}, distances between sites are also divided by the maximum distance.
#' @details The \code{gam} function uses knots to generate the splines, whose number depends on the degrees of freedom and returns an error if the number of unique values of a variable is too low with respect to the number of knots. This situation is likely to occur for factorial variables with few levels.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline}}, \code{\link{Zeta.order}}, \code{\link{Zeta.lm}}
#' @examples
#' data(BCI.spec.coarse)
#' xy <- BCI.spec.coarse[1:2]
#' data.spec <- BCI.spec.coarse[3:308]
#' data(BCI.env.coarse)
#' data.env <- BCI.env.coarse[10:15]
#' 
#' zeta.gam <- Zeta.gam(data.spec, data.env, sam = 100, order = 3)
#' summary(zeta.gam)
#' dev.new()
#' plot(zeta.gam)
#' 
#' zeta.gam <- Zeta.gam(data.spec, data.env, xy, sam = 100, order = 3)
#' summary(zeta.gam)
#' dev.new()
#' plot(zeta.gam)
#' 
#' ##########
#' 
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' data(Marion.env)
#' data.env <- Marion.env[3]
#' 
#' zeta.gam.species <- Zeta.gam(data.species, data.env, sam = 100, order = 3)
#' summary(zeta.gam.species)
#' dev.new()
#' plot(zeta.gam.species)
#' 
#' zeta.gam.species <- Zeta.gam(data.species, data.env, xy, sam = 100, order = 2)
#' summary(zeta.gam.species)
#' dev.new()
#' plot(zeta.gam.species)
#' 
#' @export
Zeta.gam <- function(data.spec, data.env, xy = NULL, order = 1, sam = 1000, standard = TRUE, method = "mean"){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(class(data.env) != "data.frame"){
  	stop(paste(deparse(substitute(data.env)), " is a ", class(data.env), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
  	stop("Wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  
  if (length(setdiff(sapply(data.env, class), c("factor", "numeric")))>0){
  	stop("Variables must be numeric or factor")
  }
  if(order == 1 & !is.null(xy)){
  	stop("Cannot include distance for order = 1")
  }
  if(order == 1 & length(setdiff(sapply(data.env, class), "numeric"))>0){
  	stop("Environmental variables must be numeric for order = 1")
  }
  if(length(intersect(sapply(data.env, class), "factor"))>0){
  	warning("Factor variables with few levels are likely to prevent the gam from running")
  }
    
  x <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()
  
  if(choose(x, order)>sam){
    u <- rep(NA, sam)
    data.var <- as.data.frame(matrix(NA, sam, dim(data.env)[2]))
    distance <- rep(NA, sam)
    
    for(z in 1:sam){
      samp <- sample(1:x, order, replace = FALSE)
      u[z] <- sum(apply(data.spec[samp, ], 2, prod))
      
      if(order == 1){
        data.var[z, ] <- data.env[samp, ]
      }else {
        
        
        fac <- which(sapply(data.env, class) == "factor")
        num <- which(sapply(data.env, class) == "numeric")
        
        if(order>2){
          if(length(num)>1){
            if(method == "mean"){
              data.var[z, num] <- apply(apply(data.env[samp, num], 2, dist), 2, mean)  ##computing directly the mean measure of euclidian distance
            }else if(method == "max"){
              data.var[z, num] <- apply(apply(data.env[samp, num], 2, dist), 2, max)
            }else{
              stop("Error: method for distance does not exist")
            }
          }else if(length(num)>0){
            data.var[z, num] <- mean(dist(data.env[samp, num]))
          }
        }else{
          if(length(num)>1){
            data.var[z, num] <- apply(data.env[samp, num], 2, dist)
          }else if(length(num)>0){
            data.var[z, num] <- dist(data.env[samp, num])
          }
        }
        if(length(fac)>1){
          data.var[z, fac] <- apply(data.env[samp, fac], 2, function(x){length(unique(x))}) - 1
        }else if (length(fac)>0){
          data.var[z, fac] <- length(unique(data.env[samp, fac])) - 1
        }
  
        if (!is.null(xy)){
          if(method == "mean"){
            distance[z] <- mean(dist(xy[samp, ]))
          }else if(method == "max"){
            distance[z] <- max(dist(xy[samp, ]))
          }else{
            stop("Error: method for distance does not exist")
          }
        }
      }
    }
    
  }else{
    u <- rep(NA, choose(x, order))
    data.var <- as.data.frame(matrix(NA, choose(x, order), dim(data.env)[2]))
    samp <- combn(1:x, order)
    distance <- rep(NA, choose(x, order))
    
    for(z in 1:dim(samp)[2]){
      u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
      
      if(order == 1){
        data.var[z, ] <- data.env[samp[, z], ]
      }else {
        
        
        fac <- which(sapply(data.env, class) == "factor")
        num <- which(sapply(data.env, class) == "numeric")
        
        if(order>2){
          if(length(num)>1){
            if(method == "mean"){
              data.var[z, num] <- apply(apply(data.env[samp[, z], num], 2, dist), 2, mean)  ##computing directly the mean measure of euclidian distance
            }else if(method == "max"){
              data.var[z, num] <- apply(apply(data.env[samp[, z], num], 2, dist), 2, max)
            }else{
              stop("Error: method for distance does not exist")
            }
          }else if(length(num)>0){
            data.var[z, num] <- mean(dist(data.env[samp[, z], num]))
          }
        }else{
          if(length(num)>1){
            data.var[z, num] <- apply(data.env[samp[, z], num], 2, dist)
          }else if(length(num)>0){
            data.var[z, num] <- dist(data.env[samp[, z], num])
          }
        }
        if(length(fac)>1){
          data.var[z, fac] <- apply(data.env[samp[, z], fac], 2, function(x){length(unique(x))}) - 1
        }else if (length(fac)>0){
          data.var[z, fac] <- length(unique(data.env[samp[, z], fac])) - 1
        }
  
        if (!is.null(xy)){
          if(method == "mean"){
            distance[z] <- mean(dist(xy[samp[, z], ]))
          }else if(method == "max"){
            distance[z] <- max(dist(xy[samp[, z], ]))
          }else{
            stop("Error: method for distance does not exist")
          }
        }
      }
    }
  }

  
  
  ##rescale the environmental variables or distances between 0 and 1
  if(standard == TRUE){
    if(order>1){
      data.var <- data.var / matrix(rep(apply(data.var, 2, max), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.env)[2], byrow = T)
    }else{
      num <- which(sapply(data.env, class) == "numeric")
      if(length(num)>1){
        data.var[, num] <- (data.var[, num] - matrix(rep(apply(data.var[, num], 2, min), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.var[, num])[2], , byrow = T)) / matrix(rep((apply(data.var[, num], 2, max) - apply(data.var[, num], 2, min)), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.var[, num])[2], byrow = T)
      }else{
        data.var[, num] <- (data.var[, num] - min(data.var[, num])) / (max(data.var[, num]) - min(data.var[, num]))
      }
    }
    if(!is.null(xy)){
    	distance <- distance / max(distance)
    }
  }
  names(data.var) <- names(data.env)
  
  zeta.val <- u
  
  if(is.null(xy)){data.tot <- data.var}else{data.tot <- cbind(data.var, distance)}
  
  ##create formula to be used in gam
  xnam <- names(data.tot)
  fm <- as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = ") + s("), ")"))
  zeta.gam <- mgcv::gam(fm, data = data.tot)
  
  return(zeta.gam)
  
  
  
}




#' Zeta distance decay for a specific number of assemblages or sites
#'
#' Computes the distance decay of zeta diversity for a specific order (number of assemblages or sites), using linear regression.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta-diversity is computed.
#' @param confint.level  Percentage for the confidence intervals of the coefficient from the linear regression
#' @param method  Indicates which distance to consider for more than 3 sites. Method can be "\code{mean}" or "\code{max}"
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @return \code{Zeta.ddecay} returns a list containing the following components:
#' @return \item{coefs}{An object of class "\code{lm}" corresponding to the linear regression over distance for the number of assemblages or sites specified in 'order'.}
#' @return \item{confint}{The confidence intervals for the coefficients from the linear regression.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline}}, \code{\link{Zeta.order}}, \code{\link{Zeta.ddecays}}
#' @examples
#' data(BCI.spec.coarse)
#' xy <- BCI.spec.coarse[1:2]
#' data.spec <- BCI.spec.coarse[3:308]
#' 
#' dev.new()
#' zeta.ddecay <- Zeta.ddecay(xy, data.spec, sam = 100, order = 3, confint.level = 0.95)
#' 
#' ##########
#' 
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' 
#' dev.new()
#' zeta.ddecay.species <- Zeta.ddecay(xy, data.species, sam = 100, order = 3, confint.level = 0.95)
#' 
#' @export
Zeta.ddecay <- function(xy, data.spec, order = 2, sam = 1000, confint.level = 0.95, method = "mean", plot = TRUE){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
  	stop("Wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  
  x <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()
    
  if(choose(x, order)>sam){
    u <- rep(NA, sam)
    distance <- rep(NA, sam)
    
    for(z in 1:sam){
      samp <- sample(1:x, order, replace = FALSE)
      u[z] <- sum(apply(data.spec[samp, ], 2, prod))
      
      if(order == 1){
        stop("Error: distance decay cannot be computed for zeta 1")
      }else {
        if(method == "mean"){
          distance[z] <- mean(dist(xy[samp, ]))
        }else if(method == "max"){
          distance[z] <- max(dist(xy[samp, ]))
        }else{
          stop("Error: method does not exist")
        }
      }
    }
    
  }else{
    u <- rep(NA, choose(x, order))
    distance <- rep(NA, choose(x, order))
    samp <- combn(1:x, order)
    
    for(z in 1:dim(samp)[2]){
      u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
      
      if(order == 1){
        stop("Error: distance decay cannot be computed for zeta 1")
      }else {
        if(method == "mean"){
          distance[z] <- mean(dist(xy[samp[, z], ]))
        }else if(method == "max"){
          distance[z] <- max(dist(xy[samp[, z], ]))
        }else{
          stop("Error: method does not exist")
        }
      }
    }
  }
    
  

      
  zeta.val <- u
  
  zeta.ddcay.lm <- lm(zeta.val ~ distance)
  
  zeta.ddcay <- list()
  zeta.ddcay$lm <- zeta.ddcay.lm
  zeta.ddcay$confint <- suppressMessages(confint(zeta.ddcay.lm, level = confint.level))
  
  if(plot == TRUE){
    plot(distance, zeta.val, xlab = "Distance", ylab = paste("Zeta ", order, sep = ""), pch = 16)
    prd <- predict.lm(zeta.ddcay.lm, data.frame(distance = sort(distance)), interval = c("confidence"), 
level = 0.95, type = "response")
    lines(sort(distance), prd[, 1], col = "red", lwd = 2)
    lines(sort(distance), prd[, 2], col = "red", lty = 2, lwd = 2)
    lines(sort(distance), prd[, 3], col = "red", lty = 2, lwd = 2)
  }
  
  return(zeta.ddcay)

  
}




#' Zeta distance decay for a range of numbers of assemblages or sites
#'
#' Computes the distance decay of zeta diversity for a range of orders (number of assemblages or sites), using linear regressions.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param orders  Range of number of assemblages or sites at which zeta diversity is computed. All the orders must be striclty greater than 1.
#' @param sam  Number of samples for which the zeta-diversity is computed.
#' @param confint.level  Percentage for the confidence intervals of the coefficient from the linear regression.
#' @param method  Indicates which distance to consider for more than 3 sites. Method can be "\code{mean}" or "\code{max}"
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @return \code{Zeta.ddecays} returns a list containing the following components:
#' @return \item{coefs}{A vector of the coefficients from the linear regressions over distance for the numbers of sites specified by 'order'.}
#' @return \item{confint}{The confidence intervals for the coefficients from the linear regressions.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline}}, \code{\link{Zeta.order}}, \code{\link{Zeta.ddecay}}
#' @examples
#' data(BCI.spec.coarse)
#' xy <- BCI.spec.coarse[1:2]
#' data.spec <- BCI.spec.coarse[3:308]
#' 
#' dev.new()
#' zeta.ddecays <- Zeta.ddecays(xy, data.spec, sam = 100, orders = 2:5, 
#'     plot = TRUE, confint.level = 0.95)
#' 
#' ##########
#' 
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' 
#' dev.new()
#' zeta.ddecays.species <- Zeta.ddecays(xy, data.species, sam = 100, orders = 2:5, 
#'    plot = TRUE, confint.level = 0.95)
#' 
#' @export
Zeta.ddecays <- function(xy, data.spec, orders = 2:10, sam = 1000, confint.level = 0.95, method = "mean", plot = TRUE){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(max(orders)>dim(data.spec)[1]){
  	stop("Wrong value for \"orders\": the maximum value must be equal or lower than the number of sites.")
  }
  
  if(length(which(orders<= 1))>0){stop("orders must be striclty greater than 1")}
  
  zeta.ddecays.coefs <- rep(NA, length(orders))
  zeta.ddecays.confint <- matrix(NA, length(orders), 2)
  ii <- 0
  for (i in orders){
    ii <- ii + 1
    temp <- Zeta.ddecay(xy, data.spec, sam = sam, order = i, confint.level = confint.level, plot = FALSE)
    zeta.ddecays.coefs[ii] <- coef(temp$lm)[2]
    zeta.ddecays.confint[ii, ] <- temp$confint[2, ]
    
  }
  
  zeta.ddcays <- list()
  zeta.ddcays$coefs <- zeta.ddecays.coefs
  zeta.ddcays$confint <- zeta.ddecays.confint
  
  if (plot == TRUE){
    plot(orders, zeta.ddcays$coefs, pch = 16, ylim = c(min(0, range(zeta.ddcays$confint)[1]), max(0, range(zeta.ddcays$confint)[2])), xlab = "Number of sites", ylab = "Slope", main = "Distance decay of zeta-diversity")
    lines(orders, zeta.ddcays$coefs)
    suppressWarnings(arrows(x0 = c(orders), y0 = zeta.ddcays$coefs, x1 = c(orders), y1 = zeta.ddcays$confint[, 1], angle = 90, length = 0.2))
    suppressWarnings(arrows(x0 = c(orders), y0 = zeta.ddcays$coefs, x1 = c(orders), y1 = zeta.ddcays$confint[, 2], angle = 90, length = 0.2))
    lines(0:(max(orders) + 2), rep(0, max(orders) + 3), lty = 2)
  }
  
  return(zeta.ddcays)
  
}





#' Zeta diversity scaling with sample grain using hierarchically increases in grain size
#'
#' Computes zeta diversity scaling with sample grain for a specific order (number of assemblages or sites), increasing grain by hierarchically nesting of regularly spaced sites.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param n  Vector of aggregation indices: regularly spaced sites are grouped as n[i] x n[i] sites.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta-diversity is computed.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @details The sites (plots or quadrates) are incrementally aggregated as nearest neighbouring groups of 4, 9, etc. sites, using a nested approach, starting from the lowest x and y, to increase the grain. The sites can be spatially contiguous or discontiguous, as long as they are regularly spaced. This function is not suitable for irregularly spaced sites.
#' @return \code{Zeta.scale.regular} returns a vector \code{zeta.scale.reg} containing the zeta-diversity values for each grain.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
#' @seealso \code{\link{Zeta.decline}}, \code{\link{Zeta.order}}, \code{\link{Zeta.scale.min.dist}}
#' @examples
#' data(BCI.spec.fine)
#' xy <- BCI.spec.fine[1:2]
#' data.spec <- BCI.spec.fine[3:308]
#' 
#' dev.new()
#' ##sam = 25 is used here for fast execution, but a higher value is advised
#' zeta.scale.reg <- Zeta.scale.regular(xy, data.spec, n = 1:3, order = 3, sam = 25)
#' @export
Zeta.scale.regular <- function(xy, data.spec, n, order = 1, sam = 1000, plot = TRUE){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
  	stop("Wrong value for \"order\": it must be equal or lower than the number of sites.")
  	
  }
  
  n <- sort(n)
  n2 <- n
  
  zeta.scale.reg <- rep(NA, length(n))
  
  names <- c(names(xy), names(data.spec))
  
  if(n[1] == 1){
      zeta.scale.reg[1] <- Zeta.order(data.spec, order = order, sam = sam)$val
      n2 <- n[2:length(n)]
  }
  
  ##sort data according to the plots coordinates
  xy <- xy[order(xy$x, xy$y), ]
  data.spec <- data.spec[order(xy$x, xy$y), ]
  
  ##compute the scale dependence for the specified grains
  for (nn in 1:length(n2)){
    
    Ux <- sort(unique(xy$x))
    Uy <- sort(unique(xy$y))
    
    max.x <- Ux[length(Ux) - length(Ux)%%(n2[nn]^2)]
    max.y <- Uy[length(Uy) - length(Uy)%%(n2[nn]^2)]
    
    
    if(length(which(xy$x>max.x | xy$y>max.y))>0){
      xy2 <- xy[-which(xy$x>max.x | xy$y>max.y), ]
      data.spec2 <- data.spec[-which(xy$x>max.x | xy$y>max.y), ]
    }else{
      xy2 <- xy
      data.spec2 <- data.spec
    }
    
    
    data2 <- data.frame(setNames(replicate(length(names), numeric(0), simplify = F), names(data)))

    
    
    for (i in seq(1, length(unique(xy2$x)), n2[nn])){
      for (j in seq(1, length(unique(xy2$y)), n2[nn])){
        temp.xy <- xy2[which(xy2$x %in% Ux[i:(i + n2[nn] - 1)] & xy2$y %in% Uy[j:(j + n2[nn] - 1)]), ]
        temp.spec <- data.spec2[which(xy2$x %in% Ux[i:(i + n2[nn] - 1)] & xy2$y %in% Uy[j:(j + n2[nn] - 1)]), ]
        
        
        ##compute the mean coordinates, the unions of species, and the mean of environmental variables
        temp.xy <- apply(temp.xy, 2, mean)
        temp.spec <- (apply(temp.spec, 2, sum)>0) * 1
        
        ## add the vector to the new coarser dataset
        data2 <- rbind(data2, as.vector(c(temp.xy, temp.spec)))
    
      }
    }
    

    
    names(data2) <- names
    
    ##compute the zeta diversity of the new grain
    if(n[1] == 1){
        zeta.scale.reg[nn + 1] <- Zeta.order(data2[3:(2 + dim(data.spec)[2])], order = order, sam = sam)$val
    }else{
        zeta.scale.reg[nn] <- Zeta.order(data2[3:(2 + dim(data.spec)[2])], order = order, sam = sam)$val
    }
    
  }
  
  if(plot == TRUE){
    plot(n^2, zeta.scale.reg, xlab = "Grain", ylab = paste("Zeta ", order, sep = ""), main = "Zeta-Scale Relationship", pch = 16)
    lines(n^2, zeta.scale.reg)
  }
  
  return(zeta.scale.reg)
  
}






#' Zeta diversity scaling with sample grain dependency based on minimum distance
#'
#' Computes zeta diversity scaling with sample grain for a specific order (number of assemblages or sites), increasing grain by sequentially adding sites based on the minimum distance between them.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param n  Vector of aggregation indices: regularly spaced sites are grouped as n[i] x n[i] sites.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param reorder  Number of times the sites are rearranged and grouped together for the computation of zeta.
#' @param sam  Number of samples for which the zeta-diversity is computed.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @param sd  Boolean value (TRUE or FALSE) indicating if the standard deviation must be plotted for each grain.
#' @details The nearest neighbouring sites (plots, quadrates, or areas of varying shapes) are grouped as spatial clusters of 2, 3, 4, etc. sites, based on the minimum distance between them. Since the procedure is based on the relative distance between sites, the site order can have an impact on the output. The procedure is therefore performed 'reorder' times, for which sites are randomly reordered each time, and the mean zeta is computed. This function is suitable for both regularly and irregularly spaced sites, contiguous or non contiguous. For regularly spaced sites, the use of \code{\link{Zeta.scale.regular}} is recommended.  
#' @return \code{zeta.scale.min.dist} returns a matrix containing the zeta-diversity values over the '\code{reorder}' computations, for each grain.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
#' @seealso \code{\link{Zeta.decline}}, \code{\link{Zeta.order}}, \code{\link{Zeta.scale.regular}}
#' @examples
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' 
#' dev.new()
#' zeta.scale.irreg.species <- Zeta.scale.min.dist(xy, data.species, n = 1:3, order = 3,
#'     reorder = 3, sam = 50)
#' 
#' @export
Zeta.scale.min.dist <- function(xy, data.spec, n, order = 1, reorder = 100, sam = 1000, plot = TRUE, sd = TRUE){
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
  	stop("Wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  
  n <- sort(n)
  n2 <- n
  
  zeta.scale.irreg <- matrix(NA, reorder, length(n))
  
  names <- c(names(xy), names(data.spec))
  
  if (n[1] == 1){
    zeta.scale.irreg[, 1] <- rep(Zeta.order(data.spec, order = order, sam = sam)$val, reorder)
    n2 <- n[2:length(n)]
  }
  
  ##compute the scale dependence for the specified grains
  for (nn in 1:length(n2)){
    
    zeta.scale.temp <- rep(NA, reorder)
    
    ##repeat 'reorder times'
    for(reord in 1:reorder){
      
      ##randomize the plot orders
      xy <- xy[sample(nrow(xy)), ]
      data.spec <- data.spec[sample(nrow(xy)), ]
      
      #pairwise distance
      D <- as.matrix(dist(xy))
      D[which(D == 0)] <- NA
      
      data2 <- data.frame(setNames(replicate(length(names), numeric(0), simplify = F), names(data)))
      
      
      i.plot <- 1
      for(i in 1:floor(dim(xy)[1] / n2[nn])){
        
        while(length(which(!is.na(D[, i.plot]))) == 0){
          i.plot <- i.plot + 1
        }
        
        ##select the plots in the group
        temp.xy <- xy[c(i.plot, order(D[, i.plot])[1:(n2[nn] - 1)]), ]
        temp.spec <- data.spec[c(i.plot, order(D[, i.plot])[1:(n2[nn] - 1)]), ]
        
        ##compute the mean coordinates, the unions of species, and the mean of environmental variables
        temp.xy <- apply(temp.xy, 2, mean)
        temp.spec <- (apply(temp.spec, 2, sum)>0) * 1
        
        ## add the vector to the new coarser dataset
        data2 <- rbind(data2, as.vector(c(temp.xy, temp.spec)))
        
        ##set the column corresponding to the plots to NA in the distance matrix to avoid using them in the following
        DD <- D
        D[, order(DD[, i.plot])[1:(n2[nn] - 1)]] <- NA
        D[order(DD[, i.plot])[1:(n2[nn] - 1)], ] <- NA
        D[, i.plot] <- NA
        D[i.plot, ] <- NA
        
      }
      
      names(data2) <- names
      
      ##compute the zeta diversity of the new grain
      if (n[1] == 1){
	      zeta.scale.irreg[reord, nn + 1] <- Zeta.order(data2[3:(2 + dim(data.spec)[2])], order = order, sam = sam)$val
	  }else{
	      zeta.scale.irreg[reord, nn] <- Zeta.order(data2[3:(2 + dim(data.spec)[2])], order = order, sam = sam)$val
	  }
      
    }
  }
  
  
  if(plot == TRUE){
    plot(n, apply(zeta.scale.irreg, 2, mean), ylim = c(apply(zeta.scale.irreg, 2, mean)[1] - apply(zeta.scale.irreg, 2, sd)[1], apply(zeta.scale.irreg, 2, mean)[length(n)] + apply(zeta.scale.irreg, 2, sd)[length(n)]), xlab = "Grain", ylab = paste("Zeta ", order, sep = ""), main = "Zeta-Scale Relationship", pch = 16)
    lines(n, apply(zeta.scale.irreg, 2, mean))
    if(sd == TRUE){
      for(i in n){
        suppressWarnings(arrows(n[i], apply(zeta.scale.irreg, 2, mean)[i], n[i], apply(zeta.scale.irreg, 2, mean)[i] + apply(zeta.scale.irreg, 2, sd)[i], angle = 90, length = 0.1))
        suppressWarnings(arrows(n[i], apply(zeta.scale.irreg, 2, mean)[i], n[i], apply(zeta.scale.irreg, 2, mean)[i] - apply(zeta.scale.irreg, 2, sd)[i], angle = 90, length = 0.1))
      }
    }
  }
  
  return(zeta.scale.irreg)
  
}





#' Variation partitioning for zeta diversity
#'
#' Variation partitioning of zeta diversity for a specific order (number of assemblages or sites) over distance and environmental variables.
#' @param xy  Site coordinates
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param data.env  Sites-by-variable data frame, with sites as rows and environmental variables as columns.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta-diversity is computed.
#' @param method Indicates which distance to consider for more than 3 sites. Method can be "\code{mean}" or "\code{max}"
#' @param standard  Boolean parameter indicating if the spatial distances and differences in environmental variables should be standardized between 0 and 1.
#' @return \code{Zeta.varpart} returns a list containing the following components:
#' @return \item{varpart}{Object of class "\code{varpart}", containing the output of the variation partitioning of zeta over the distance and environmental variables.}
#' @return \item{signif.abc}{Object of class "\code{anova}", showing the significance of the partitioning of zeta over all variables.}
#' @return \item{signif.ab}{Object of class "\code{anova}", showing the significance of the partitioning of zeta over distance.}
#' @return \item{signif.bc}{Object of class "\code{anova}", showing the significance of the partitioning of zeta over the environmental variables.}
#' @return \item{signif.a}{Object of class "\code{anova}", showing the significance of the partitioning of zeta over distance, removing the part relative to the environmental variables.}
#' @return \item{signif.c}{Object of class "\code{anova}", showing the significance of the partitioning of zeta over the environmental variables, removing the part relative to distance.}
#' @details This function calls function \code{varpart} from package \{vegan\} and apply it to zeta for a specific order (number of assemblages or sites). The variation partitioning is based on adjusted R squared to account for the differences in numbers of explanatory variables.
#' @details The environmental variables can be numeric or factorial, and \code{order} must be greater than 1. 
#' @details For numeric variables, the pairwise difference between sites is computed and combined according to \code{method}. For factorial variables, the distance corresponds to the number of unique values over the number of assemblages of sites specified by \code{order}. 
#' @details If xy = NULL, \code{Zeta.lm} only uses environmental variables in the regression. Otherwise, it also computes and uses euclidian distance (average or maximum distance between multiple sites, depending on the parameters \code{method}) as an explanatory variable.
#' @details Zeta is regressed against the differences of values of the environmental variables divided by the maximum difference for each variable, to be rescaled between 0 and 1. If \code{xy != NULL}, distances between sites are also divided by the maximum distance.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Borcard, D., Legendre, P. & Drapeau, P. (1992). Partialling out the spatial component of ecological variation. \emph{Ecology} 73, 1045-1055.
#' @references Legendre, P. &  Legendre, L.F. (2012). \emph{Numerical ecology}, 3rd English edition. Elsevier Science BV, Amsterdam.
#' @seealso \code{\link{Zeta.decline}}, \code{\link{Zeta.order}}, \code{\link{Zeta.lm}}, \code{\link{Zeta.gam}}
#' @examples
#' data(BCI.spec.coarse)
#' xy <- BCI.spec.coarse[1:2]
#' data.spec <- BCI.spec.coarse[3:308]
#' data(BCI.env.coarse)
#' data.env <- BCI.env.coarse[10:15]
#' 
#' zeta.varpart <- Zeta.varpart(xy, data.spec, data.env, order = 3, sam = 100)
#' zeta.varpart
#' dev.new()
#' plot(zeta.varpart$varpart)
#' dev.new()
#' pie.neg(zeta.varpart$varpart$part$indfract[, 3], density = c(4, 0, 8, -1), 
#'     angle = c(90, 0, 0, 0), labels = c("distance", "b", "environment", "d"), radius = 0.9)
#' 
#' ##########
#' 
#' data(Marion.species)
#' xy <- Marion.species[1:2]
#' data.species <- Marion.species[3:33]
#' data(Marion.env)
#' data.env <- Marion.env[3:4]
#' 
#' zeta.varpart.species <- Zeta.varpart(xy, data.species, data.env, order = 3, sam = 100)
#' zeta.varpart.species
#' dev.new()
#' plot(zeta.varpart.species$varpart)
#' dev.new()
#' pie.neg(zeta.varpart.species$varpart$part$indfract[, 3], density = c(4, 0, 8, -1), 
#'     angle = c(90, 0, 0, 0), labels = c("distance", "b", "environment", "d"), radius = 0.9)
#' 
#' @export
Zeta.varpart <- function(xy, data.spec, data.env, order = 2, sam = 1000, method = "mean", standard = TRUE){
  
  if(order<2){
  	stop("Cannot compute variation partitioning for order<2")
  }
  if(order>dim(data.spec)[1]){
  	stop("Wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  
  if(class(data.spec) != "data.frame"){
  	stop(paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(class(data.env) != "data.frame"){
  	stop(paste(deparse(substitute(data.env)), " is a ", class(data.env), ". It must be a data frame.", sep = ""))
  }
  
  if (length(setdiff(sapply(data.env, class), c("factor", "numeric")))>0){
  	stop("Variables must be numeric or factor")
  }
  
  
  x <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()
  
  
  if(choose(x, order)>sam){
    u <- rep(NA, sam)
    data.var <- as.data.frame(matrix(NA, sam, dim(data.env)[2]))
    distance <- rep(NA, sam)
    
    for(z in 1:sam){
      samp <- sample(1:x, order, replace = FALSE)
      u[z] <- sum(apply(data.spec[samp, ], 2, prod))
      
      if(order == 1){
        data.var[z, ] <- data.env[samp, ]
      }else {
        
        
        fac <- which(sapply(data.env, class) == "factor")
        num <- which(sapply(data.env, class) == "numeric")
        
        if(order>2){
          if(length(num)>1){
            if(method == "mean"){
              data.var[z, num] <- apply(apply(data.env[samp, num], 2, dist), 2, mean)  ##computing directly the mean measure of euclidian distance
            }else if(method == "max"){
              data.var[z, num] <- apply(apply(data.env[samp, num], 2, dist), 2, max)
            }else{
              stop("Error: method for distance does not exist")
            }
          }else if(length(num)>0){
            data.var[z, num] <- mean(dist(data.env[samp, num]))
          }
        }else{
          if(length(num)>1){
            data.var[z, num] <- apply(data.env[samp, num], 2, dist)
          }else if(length(num)>0){
            data.var[z, num] <- dist(data.env[samp, num])
          }
        }
        if(length(fac)>1){
          data.var[z, fac] <- apply(data.env[samp, fac], 2, function(x){length(unique(x))}) - 1
        }else if (length(fac)>0){
          data.var[z, fac] <- length(unique(data.env[samp, fac])) - 1
        }
  
        if (!is.null(xy)){
          if(method == "mean"){
            distance[z] <- mean(dist(xy[samp, ]))
          }else if(method == "max"){
            distance[z] <- max(dist(xy[samp, ]))
          }else{
            stop("Error: method for distance does not exist")
          }
        }
      }
    }
    
  }else{
    u <- rep(NA, choose(x, order))
    data.var <- as.data.frame(matrix(NA, choose(x, order), dim(data.env)[2]))
    distance <- rep(NA, choose(x, order))
    samp <- combn(1:x, order)
    
    for(z in 1:dim(samp)[2]){
      u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
      
      if(order == 1){
        data.var[z, ] <- data.env[samp[, z], ]
      }else {
        
        
        fac <- which(sapply(data.env, class) == "factor")
        num <- which(sapply(data.env, class) == "numeric")
        
        if(order>2){
          if(length(num)>1){
            if(method == "mean"){
              data.var[z, num] <- apply(apply(data.env[samp[, z], num], 2, dist), 2, mean)  ##computing directly the mean measure of euclidian distance
            }else if(method == "max"){
              data.var[z, num] <- apply(apply(data.env[samp[, z], num], 2, dist), 2, max)
            }else{
              stop("Error: method for distance does not exist")
            }
          }else if(length(num)>0){
            data.var[z, num] <- mean(dist(data.env[samp[, z], num]))
          }
        }else{
          if(length(num)>1){
            data.var[z, num] <- apply(data.env[samp[, z], num], 2, dist)
          }else if(length(num)>0){
            data.var[z, num] <- dist(data.env[samp[, z], num])
          }
        }
        if(length(fac)>1){
          data.var[z, fac] <- apply(data.env[samp[, z], fac], 2, function(x){length(unique(x))}) - 1
        }else if (length(fac)>0){
          data.var[z, fac] <- length(unique(data.env[samp[, z], fac])) - 1
        }
  
        if (!is.null(xy)){
          if(method == "mean"){
            distance[z] <- mean(dist(xy[samp[, z], ]))
          }else if(method == "max"){
            distance[z] <- max(dist(xy[samp[, z], ]))
          }else{
            stop("Error: method for distance does not exist")
          }
        }
      }
    }
  }
      
  ##rescale the environmental variables or distances between 0 and 1
  if(standard == TRUE){
    if(order>1){
      data.var <- data.var / matrix(rep(apply(data.var, 2, max), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.env)[2], byrow = T)
      
    }else{
      num <- which(sapply(data.env, class) == "numeric")
      if(length(num)>1){
        data.var[, num] <- (data.var[, num] - matrix(rep(apply(data.var[, num], 2, min), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.var[, num])[2], , byrow = T)) / matrix(rep((apply(data.var[, num], 2, max) - apply(data.var[, num], 2, min)), min(choose(x, order), sam)), min(choose(x, order), sam), dim(data.var[, num])[2], byrow = T)
      }else{
        data.var[, num] <- (data.var[, num] - min(data.var[, num])) / (max(data.var[, num]) - min(data.var[, num]))
      }
    }
    if(!is.null(xy)){
    	distance <- distance / max(distance)
    }
  }
  names(data.var) <- names(data.env)
  
  if(!is.null(xy)){
  	distance <- distance / max(distance)
  }
  zeta.val <- u
  
  data.tot <- cbind(data.var, distance)
  distance <- data.frame(distance)
  
  
  zeta.varpart <- list()
  
  ##variation partitioning
  zeta.varpart$varpart <- vegan::varpart(zeta.val, distance, data.var)
  
  ##significance of [distance + environment]
  zeta.varpart$signif.abc <- anova(vegan::rda(zeta.val ~ ., data = data.tot))
  
  ##significance of [distance -- a + b]
  zeta.varpart$signif.ab <- anova(vegan::rda(zeta.val ~ ., data = distance))
  
  ##significance of [environment -- b + c]
  zeta.varpart$signif.bc <- anova(vegan::rda(zeta.val ~ ., data = data.var))
  
  ##significance of [distance only -- a]
  zeta.varpart$signif.a <- anova(vegan::rda(X = zeta.val, Y = distance, Z = data.var))
  
  ##significance of [environment only -- c]
  zeta.varpart$signif.c <- anova(vegan::rda(X = zeta.val, Y = data.var, Z = distance))
  
  
    
  return(zeta.varpart)
  
  
  
}



#' Pie Chart, considering negative values as zeros
#'
#' Plots a pie chart, considering negative values as zeros, for the purpose of illustrating variation partitioning.
#' @param x  a vector of non-negative numerical quantities. The values in x are displayed as the areas of pie slices.
#' @param labels one or more expressions or character strings giving names for the slices. Other objects are coerced by as.graphicsAnnot. For empty or NA (after coercion to character) labels, no label nor pointing line is drawn.
#' @param edges  the circular outline of the pie is approximated by a polygon with this many edges.
#' @param radius  the pie is drawn centered in a square box whose sides range from -1 to 1. If the character strings labeling the slices are long it may be necessary to use a smaller radius.
#' @param clockwise  logical indicating if slices are drawn clockwise or counter clockwise (i.e., mathematically positive direction), the latter is default.
#' @param init.angle  number specifying the starting angle (in degrees) for the slices. Defaults to 0 (i.e., '3 o'clock') unless clockwise is true where init.angle defaults to 90 (degrees), (i.e., '12 o'clock').
#' @param density  the density of shading lines, in lines per inch. The default value of NULL means that no shading lines are drawn. Non-positive values of density also inhibit the drawing of shading lines.
#' @param angle  the slope of shading lines, given as an angle in degrees (counter-clockwise).
#' @param col  a vector of colors to be used in filling or shading the slices. If missing a set of 6 pastel colours is used, unless density is specified when par("fg") is used.
#' @param border,lty  (possibly vectors) arguments passed to polygon which draws each slice.
#' @param main  an overall title for the plot.
#' @param warning Boolean value. Set to FALSE to avoid displaying a warning if some values are negative and set to 0.
#' @param ...  graphical parameters can be given as arguments to pie. They will affect the main title and labels only.
#' @details This function is identical to the function \code{\link{pie}} in \{graphics\}, except that it considers all negative values as zeros, to allow for plotting variation partitioning outputs. The original \code{\link{pie}} function returns an error when negative values are present. However, variation partitioning can return negative values, which can then be treated as zeros (Legendre & Legendre, 2008). This function allows direct use of the results from \code{\link{Zeta.varpart}} without editing the data.
#' @seealso \code{\link{pie}}, \code{\link{Zeta.varpart}}
#' @references  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988). \emph{The new S language}. Wadsworth & Brooks/Cole.
#' @references  Cleveland, W. S. (1985). \emph{The elements of graphing data}. Wadsworth: Monterey, CA, USA.
#' @references Legendre, P. &  Legendre, L.F. (2012). \emph{Numerical ecology}, 3rd English edition. Elsevier Science BV, Amsterdam.
#' @examples
#' pie.neg(rep(1, 24), col = rainbow(24), radius = 0.9)
#' @export
pie.neg <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
    init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
    col = NULL, border = NULL, lty = NULL, main = NULL, warning = TRUE, ...) 
{
    if (!is.numeric(x) || any(is.na(x))) 
      stop("'x' values must be numeric.")
    if (sum(x < 0)>0){
      x[which(x<0)] <- 0
      warning("Negative values set to 0.")
  }
    
    
    if (is.null(labels)) 
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x) / sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L] / pin[2L]) * xlim
    else ylim <- (pin[2L] / pin[1L]) * ylim
    dev.hold()
    on.exit(dev.flush())
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "lightblue", "mistyrose", "lightcyan", 
                "lavender", "cornsilk")
        else par("fg")
    if (!is.null(col)) 
        col <- rep_len(col, nx)
    if (!is.null(border)) 
        border <- rep_len(border, nx)
    if (!is.null(lty)) 
        lty <- rep_len(lty, nx)
    angle <- rep(angle, nx)
    if (!is.null(density)) 
        density <- rep_len(density, nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi / 180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
        P <- t2xy(mean(x[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
            text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
                adj = ifelse(P$x < 0, 1, 0), ...)
        }
    }
    title(main = main, ...)
    invisible(NULL)
}






####################
##DATA DESCRIPTION##
####################


#' Barro Colorado Island 50 Hectare Quadrat Environmental Dataset at Coarse Scale
#' 
#' Measurements of soil composition and elevation in 50, 100m x 100m contiguous quadrats.
#' 
#' 
#' The data set contains the following variables:
#' 
#' \itemize{
#' \item{x}: x-position in meters
#' \item{y}: y-position in meters
#' \item{Al}: Aluminium concentration
#' \item{B}: Boron concentration
#' \item{Ca}: Calcium concentration
#' \item{Cu}: Copper concentration
#' \item{Fe}: Iron concentration
#' \item{K}: Potassium concentration
#' \item{Mg}: Magnesium concentration
#' \item{Mn}: Manganese concentration
#' \item{P}: Phosphorous concentration
#' \item{Zn}: Zinc concentration
#' \item{N}: Nitrogen concentration
#' \item{pH}: Interpolated pH content in the surface soil
#' \item{elevation}: mean elevation
#' }
#'
#' Location: Panama -- 9 degrees 9' 4.5" N, 79 degrees 51' 19.08" W
#' 
#' Data owner: Center for Tropical Forest Science of the Smithsonian Tropical Research Institute
#' @name BCI.env.coarse
#' @usage data(BCI.env.coarse)
#' @docType data
#' @references Hubbell, S.P., Condit, R., and Foster, R.B. (2005). Barro Colorado Forest Census Plot Data. \url{http://ctfs.arnarb.harvard.edu/webatlas/datasets/bci/}
#' @references Condit, R. (1998). \emph{Tropical Forest Census Plots}. Springer-Verlag and R. G. Landes Company, Berlin, Germany, and Georgetown, Texas. 
#' @references Hubbell, S.P., R.B. Foster, S.T. O'Brien, K.E. Harms, R. Condit, B. Wechsler, S.J. Wright, and S. Loo de Lao. (1999). Light gap disturbances, recruitment limitation, and tree diversity in a neotropical forest. \emph{Science} 283, 554-557.
#' @references \url{http://ctfs.arnarb.harvard.edu/webatlas/datasets/bci/soilmaps/BCIsoil.html}
#' @keywords data
#' @format A data frame with 50 rows (quadrats) and 15 columns (environmental variables).
"BCI.env.coarse"


#' Barro Colorado Island 50 Hectare Quadrat Environmental Dataset at Fine Scale
#' 
#' Measurements of soil composition and elevation in 1250, 20m x 20m contiguous quadrats.
#' 
#' 
#' The data set contains the following variables:
#' 
#' \itemize{
#' \item{x}: x-position in meters
#' \item{y}: y-position in meters
#' \item{Al}: Aluminium concentration
#' \item{B}: Boron concentration
#' \item{Ca}: Calcium concentration
#' \item{Cu}: Copper concentration
#' \item{Fe}: Iron concentration
#' \item{K}: Potassium concentration
#' \item{Mg}: Magnesium concentration
#' \item{Mn}: Manganese concentration
#' \item{P}: Phosphorous concentration
#' \item{Zn}: Zinc concentration
#' \item{N}: Nitrogen concentration
#' \item{pH}: Interpolated pH content in the surface soil
#' \item{elevation}: mean elevation
#' }
#'
#' Location: Panama -- 9 degrees 9' 4.5" N, 79 degrees 51' 19.08" W
#' 
#' Data owner: Center for Tropical Forest Science of the Smithsonian Tropical Research Institute
#' @name BCI.env.fine
#' @usage data(BCI.env.fine)
#' @docType data
#' @references Hubbell, S.P., Condit, R., and Foster, R.B. (2005). Barro Colorado Forest Census Plot Data. \url{http://ctfs.arnarb.harvard.edu/webatlas/datasets/bci/}
#' @references Condit, R. (1998). \emph{Tropical Forest Census Plots}. Springer-Verlag and R. G. Landes Company, Berlin, Germany, and Georgetown, Texas. 
#' @references Hubbell, S.P., R.B. Foster, S.T. O'Brien, K.E. Harms, R. Condit, B. Wechsler, S.J. Wright, and S. Loo de Lao. (1999). Light gap disturbances, recruitment limitation, and tree diversity in a neotropical forest. \emph{Science} 283, 554-557.
#' @references \url{http://ctfs.arnarb.harvard.edu/webatlas/datasets/bci/soilmaps/BCIsoil.html}
#' @keywords data
#' @format A data frame with 1250 rows (quadrats) and 15 columns (environmental variables)
"BCI.env.fine"


#' Barro Colorado Island 50 Hectare Quadrat Species Presence-Absence Dataset at Coarse Scale
#' 
#' Inventory of species presence-absence in 50, 100m x 100m contiguous quadrats.
#' 
#'
#' Location: Panama -- 9 degrees 9' 4.5" N, 79 degrees 51' 19.08" W
#' 
#' Data owner: Center for Tropical Forest Science of the Smithsonian Tropical Research Institute
#' @name BCI.spec.coarse
#' @usage data(BCI.spec.coarse)
#' @docType data
#' @references Hubbell, S.P., Condit, R., and Foster, R.B. (2005). Barro Colorado Forest Census Plot Data. \url{http://ctfs.arnarb.harvard.edu/webatlas/datasets/bci/}
#' @references Condit, R. (1998). \emph{Tropical Forest Census Plots}. Springer-Verlag and R. G. Landes Company, Berlin, Germany, and Georgetown, Texas. 
#' @references Hubbell, S.P., R.B. Foster, S.T. O'Brien, K.E. Harms, R. Condit, B. Wechsler, S.J. Wright, and S. Loo de Lao. (1999). Light gap disturbances, recruitment limitation, and tree diversity in a neotropical forest. \emph{Science} 283, 554-557.
#' @keywords data
#' @format A data frame with 50 rows (quadrats) and 308 columns (species).
"BCI.spec.coarse"


#' Barro Colorado Island 50 Hectare Quadrat Species Presence-Absence Dataset at Fine Scale
#' 
#' Inventory of species presence-absence in 1250, 20m x 20m contiguous quadrats.
#' 
#'
#' Location: Panama -- 9 degrees 9' 4.5" N, 79 degrees 51' 19.08" W
#' 
#' Data owner: Center for Tropical Forest Science of the Smithsonian Tropical Research Institute
#' @name BCI.spec.fine
#' @usage data(BCI.spec.fine)
#' @docType data
#' @references Hubbell, S.P., Condit, R., and Foster, R.B. (2005). Barro Colorado Forest Census Plot Data. \url{http://ctfs.arnarb.harvard.edu/webatlas/datasets/bci/}
#' @references Condit, R. (1998). \emph{Tropical Forest Census Plots}. Springer-Verlag and R. G. Landes Company, Berlin, Germany, and Georgetown, Texas. 
#' @references Hubbell, S.P., R.B. Foster, S.T. O'Brien, K.E. Harms, R. Condit, B. Wechsler, S.J. Wright, and S. Loo de Lao. (1999). Light gap disturbances, recruitment limitation, and tree diversity in a neotropical forest. \emph{Science} 283, 554-557.
#' @keywords data
#' @format A data frame with 1250 rows (quadrats) and 308 columns (species).
"BCI.spec.fine"



#' Marion Island Species Presence-Absence Dataset
#' 
#' Inventory of springtails and mite species presence-absence in 12 plots (4 transects and 3 altitudes) on Marion Island.
#' 
#' The data set contains the following variables:
#' 
#' \itemize{
#' \item{x}: x-position in meters in UTM37G projection
#' \item{y}: y-position in meters in UTM37G projection
#' \item{columns 3-24}: mite species presence absence
#' \item{columns 25-33}: springtail species presence absence
#' }
#'
#' Location: Marion Island -- 46 degrees 53' 34.2" S, 37 degrees 45' 02.3" E
#' 
#' Data owner: Melodie A. McGeoch
#' @name Marion.species
#' @usage data(Marion.species)
#' @docType data
#' @references Nyakatya, M.J. & McGeoch, M.A. (2008). Temperature variation across Marion Island associated with a keystone plant species (\emph{Azorella selago} Hook. (Apiaceae)). Polar Biology, 31, 139-151.
#' @references McGeoch, M.A., Le Roux, P.C., Hugo, E.A. & Nyakatya, M.J. (2008). Spatial variation in the terrestrial biotic system. The Prince Edward Islands: Land-Sea Interactions in a Changing World (ed. by S.L. Chown and P.W. Froneman), pp. 245-276. African SunMedia, Stellenbosch.
#' @keywords data
#' @format A data frame with 12 rows (plots) and 33 columns (species).
"Marion.species"

#' Marion Island Environmental Dataset
#' 
#' Geographic coordinates, altitude and island side (East, West) at 12 plots (4 transects and 3 altitudes) on Marion Island.
#' 
#' The data set contains the following variables:
#' 
#' \itemize{
#' \item{x}: x-position in meters in UTM37G projection
#' \item{y}: y-position in meters in UTM37G projection
#' \item{Altitude}: mean elevation
#' \item{Side}: cardinal (East or West) side of the island
#' }
#'
#' Location: Marion Island -- 46 degrees 53' 34.2" S, 37 degrees 45' 02.3" E
#' 
#' Data owner: Melodie A. McGeoch
#' @name Marion.env
#' @usage data(Marion.env)
#' @docType data
#' @references Nyakatya, M.J. & McGeoch, M.A. (2008). Temperature variation across Marion Island associated with a keystone plant species (\emph{Azorella selago} Hook. (Apiaceae)). Polar Biology, 31, 139-151.
#' @references McGeoch, M.A., Le Roux, P.C., Hugo, E.A. & Nyakatya, M.J. (2008). Spatial variation in the terrestrial biotic system. The Prince Edward Islands: Land-Sea Interactions in a Changing World (ed. by S.L. Chown and P.W. Froneman), pp. 245-276. African SunMedia, Stellenbosch.
#' @keywords data
#' @format A data frame with 12 rows (plots) and 4 columns (variables).
"Marion.env"





