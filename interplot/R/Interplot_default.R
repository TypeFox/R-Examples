#' Plot Conditional Coefficients in (Generalized) Linear Models with Interaction Terms
#' 
#' \code{interplot.default} is a method to calculate conditional coefficient estimates from the results of (generalized) linear regression models with interaction terms. 
#' 
#' @param m A model object including an interaction term, or, alternately, a data frame recording conditional coefficients.
#' @param var1 The name (as a string) of the variable of interest in the interaction term; its conditional coefficient estimates will be plotted.
#' @param var2 The name (as a string) of the other variable in the interaction term.
#' @param plot A logical value indicating whether the output is a plot or a dataframe including the conditional coefficient estimates of var1, their upper and lower bounds, and the corresponding values of var2.
#' @param hist A logical value indicating if there is a histogram of `var2` added at the bottom of the conditional effect plot.
#' @param var2_dt A numerical value indicating the frequency distibution of `var2`. It is only used when `hist == TRUE`. When the object is a model, the default is the distribution of `var2` of the model. 
#' @param point A logical value determining the format of plot. By default, the function produces a line plot when var2 takes on ten or more distinct values and a point (dot-and-whisker) plot otherwise; option TRUE forces a point plot.
#' @param sims Number of independent simulation draws used to calculate upper and lower bounds of coefficient estimates: lower values run faster; higher values produce smoother curves.
#' @param xmin A numerical value indicating the minimum value shown of x shown in the graph. Rarely used.
#' @param xmax A numerical value indicating the maximum value shown of x shown in the graph. Rarely used.
#' @param ercolor A character value indicating the outline color of the whisker or ribbon.
#' @param esize A numerical value indicating the size of the whisker or ribbon.
#' @param ralpha A numerical value indicating the transparency of the ribbon.
#' @param rfill A character value indicating the filling color of the ribbon.
#' @param ... Other ggplot aesthetics arguments for points in the dot-whisker plot or lines in the line-ribbon plots. Not currently used.
#' 
#' @details \code{interplot.default} is a S3 method from the \code{interplot}. It works on two classes of objects:
#' \itemize{
#'   \item Ordinary linear models (object class: \code{lm});
#'   \item Generalized linear models (object class: \code{glm}).
#'   }
#' 
#' Because the output function is based on \code{\link[ggplot2]{ggplot}}, any additional arguments and layers supported by \code{ggplot2} can be added with the \code{+}. 
#' 
#' @return The function returns a \code{ggplot} object.
#' 
#' @importFrom arm sim
#' @importFrom stats quantile
#' @import  ggplot2
#' 
#' 
#' 
#' @export

# S3 method for class 'lm' and 'glm'
interplot.default <- function(m, var1, var2, plot = TRUE, hist = FALSE, var2_dt = NA, point = FALSE, 
    sims = 5000, xmin = NA, xmax = NA, ercolor = NA, esize = 0.5, ralpha = 0.5, rfill = "grey70", 
    ...) {
    set.seed(324)
    
    m.class <- class(m)
    m.sims <- arm::sim(m, sims)
    
    
    ### For factor base terms###
    factor_v1 <- factor_v2 <- FALSE
    
    if (is.factor(eval(parse(text = paste0("m$model$", var1)))) & is.factor(eval(parse(text = paste0("m$model$", 
        var2))))) 
        stop("The function does not support interactions between two factors.")
    
    
    if (is.factor(eval(parse(text = paste0("m$model$", var1))))) {
        var1_bk <- var1
        var1 <- paste0(var1, eval(parse(text = paste0("m$xlevel$", var1))))
        factor_v1 <- TRUE
        ifelse(var1 == var2, var12 <- paste0("I(", var1, "^2)"), var12 <- paste0(var2, ":", var1)[-1])
        
        # the first category is censored to avoid multicolinarity
        for (i in seq(var12)) {
            if (!var12[i] %in% names(m$coef)) 
                var12[i] <- paste0(var1, ":", var2)[-1][i]
            if (!var12[i] %in% names(m$coef)) 
                stop(paste("Model does not include the interaction of", var1, "and", var2, "."))
        }
        
    } else if (is.factor(eval(parse(text = paste0("m$model$", var2))))) {
        var2_bk <- var2
        var2 <- paste0(var2, eval(parse(text = paste0("m$xlevel$", var2))))
        factor_v2 <- TRUE
        ifelse(var1 == var2, var12 <- paste0("I(", var1, "^2)"), var12 <- paste0(var2, ":", var1)[-1])
        
        # the first category is censored to avoid multicolinarity
        for (i in seq(var12)) {
            if (!var12[i] %in% names(m$coef)) 
                var12[i] <- paste0(var1, ":", var2)[-1][i]
            if (!var12[i] %in% names(m$coef)) 
                stop(paste("Model does not include the interaction of", var1, "and", var2, "."))
        }
        
    } else {
        ifelse(var1 == var2, var12 <- paste0("I(", var1, "^2)"), var12 <- paste0(var2, ":", var1))
        
        # the first category is censored to avoid multicolinarity
        for (i in seq(var12)) {
            if (!var12[i] %in% names(m$coef)) 
                var12[i] <- paste0(var1, ":", var2)[i]
            if (!var12[i] %in% names(m$coef)) 
                stop(paste("Model does not include the interaction of", var1, "and", var2, "."))
        }
    }
    
    ################### 
    
    
    if (factor_v2) {
        xmin <- 0
        xmax <- 1
        steps <- 2
    } else {
        if (is.na(xmin)) 
            xmin <- min(m$model[var2], na.rm = T)
        if (is.na(xmax)) 
            xmax <- max(m$model[var2], na.rm = T)
        
        steps <- eval(parse(text = paste0("length(unique(na.omit(m$model$", var2, ")))")))
        if (steps > 100) 
            steps <- 100  # avoid redundant calculation
    }
    
    coef <- data.frame(fake = seq(xmin, xmax, length.out = steps), coef1 = NA, ub = NA, lb = NA)
    coef_df <- data.frame(fake = numeric(0), coef1 = numeric(0), ub = numeric(0), lb = numeric(0), 
        model = character(0))
    
    if (factor_v1) {
        for (j in 1:(length(eval(parse(text = paste0("m$xlevel$", var1_bk)))) - 1)) {
            # only n - 1 interactions; one category is avoided against multicolinarity
            
            for (i in 1:steps) {
                coef$coef1[i] <- mean(m.sims@coef[, match(var1[j + 1], names(m$coef))] + coef$fake[i] * 
                  m.sims@coef[, match(var12[j], names(m$coef))])
                coef$ub[i] <- quantile(m.sims@coef[, match(var1[j + 1], names(m$coef))] + coef$fake[i] * 
                  m.sims@coef[, match(var12[j], names(m$coef))], 0.975)
                coef$lb[i] <- quantile(m.sims@coef[, match(var1[j + 1], names(m$coef))] + coef$fake[i] * 
                  m.sims@coef[, match(var12[j], names(m$coef))], 0.025)
            }
            
            if (plot == TRUE) {
                coef$value <- var1[j + 1]
                coef_df <- rbind(coef_df, coef)
                if (hist == TRUE) {
                  if (is.na(var2_dt)) {
                    var2_dt <- eval(parse(text = paste0("m$model$", var2)))
                  } else {
                    var2_dt <- var2_dt
                  }
                }
            } else {
                names(coef) <- c(var2, "coef", "ub", "lb")
                return(coef)
            }
        }
        coef_df$value <- as.factor(coef_df$value)
        interplot.plot(m = coef_df, hist = hist, var2_dt = var2_dt, point = point, ercolor = ercolor, esize = esize, ralpha = ralpha, rfill = rfill, ...) + facet_grid(. ~ value)
        
    } else if (factor_v2) {
        for (j in 1:(length(eval(parse(text = paste0("m$xlevel$", var2_bk)))) - 1)) {
            # only n - 1 interactions; one category is avoided against multicolinarity
            
            for (i in 1:steps) {
                coef$coef1[i] <- mean(m.sims@coef[, match(var1, names(m$coef))] + coef$fake[i] * 
                  m.sims@coef[, match(var12[j], names(m$coef))])
                coef$ub[i] <- quantile(m.sims@coef[, match(var1, names(m$coef))] + coef$fake[i] * 
                  m.sims@coef[, match(var12[j], names(m$coef))], 0.975)
                coef$lb[i] <- quantile(m.sims@coef[, match(var1, names(m$coef))] + coef$fake[i] * 
                  m.sims@coef[, match(var12[j], names(m$coef))], 0.025)
            }
            
            if (plot == TRUE) {
                coef$value <- var2[j + 1]
                coef_df <- rbind(coef_df, coef)
                if (hist == TRUE) {
                  if (is.na(var2_dt)) {
                    var2_dt <- eval(parse(text = paste0("m$model$", var2)))
                  } else {
                    var2_dt <- var2_dt
                  }
                }
            } else {
                names(coef) <- c(var2, "coef", "ub", "lb")
                return(coef)
            }
        }
        coef_df$value <- as.factor(coef_df$value)
        interplot.plot(m = coef_df, hist = hist, var2_dt = var2_dt, point = point, ercolor = ercolor, 
            esize = esize, ralpha = ralpha, rfill = rfill, ...) + facet_grid(. ~ value)
        
        
    } else {
        ## Correct marginal effect for quadratic terms
        multiplier <- if (var1 == var2) 
            2 else 1
        
        for (i in 1:steps) {
            coef$coef1[i] <- mean(m.sims@coef[, match(var1, names(m$coef))] + multiplier * coef$fake[i] * 
                m.sims@coef[, match(var12, names(m$coef))])
            coef$ub[i] <- quantile(m.sims@coef[, match(var1, names(m$coef))] + multiplier * coef$fake[i] * 
                m.sims@coef[, match(var12, names(m$coef))], 0.975)
            coef$lb[i] <- quantile(m.sims@coef[, match(var1, names(m$coef))] + multiplier * coef$fake[i] * 
                m.sims@coef[, match(var12, names(m$coef))], 0.025)
        }
        
        if (plot == TRUE) {
            if (hist == TRUE) {
                if (is.na(var2_dt)) {
                  var2_dt <- eval(parse(text = paste0("m$model$", var2)))
                } else {
                  var2_dt <- var2_dt
                }
            }
            interplot.plot(m = coef, hist = hist, var2_dt = var2_dt, point = point, ercolor = ercolor, 
                esize = esize, ralpha = ralpha, rfill = rfill, ...)
        } else {
            names(coef) <- c(var2, "coef", "ub", "lb")
            return(coef)
        }
        
    }
    
} 
