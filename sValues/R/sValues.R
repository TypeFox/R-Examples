##### User functions ####

utils::globalVariables(c("variable", "coefficients"))

# New functions -----------------------------------------------------------
# to be exported soon

# Z/Chi s-values (to be added)
##' @importFrom stats setNames coef
zChi <- function(x){
  if(class(x) != "sValues") stop("Object must be of class 'sValues'")
  z   <- coef(x)[, "t_ols_all"]
  chi <- sqrt(Reduce("*", summary(x$all)$fstatistic[c("value", "numdf")]))
  setNames(z/chi, row.names(coef(x)))
}

# g-prior svalue (to be added)
g_s_values <- function(x, R2_bounds){
  if(class(x) != "sValues") stop("Object must be of class 'sValues'")
  zchi <- zChi(x)
  
  r1  <- min(R2_bounds)
  r2  <- max(R2_bounds)
  n   <- length(x$all$residuals)
  k   <- nrow(x$s_values)
  R2  <- summary(x$all)$adj.r.squared
  
  zchi*((2*(r1*r2)/((k*(1- R2))/n)) + r1 + r2)/(r2 - r1)
}




# Main function - sValues -------------------------------------------------



##' S-values: conventional model ambiguity measures
##' 
##' @description 
##' The function \code{sValues} performs the extreme bound analysis proposed by Leamer (2014) and 
##' discussed in Leamer (2015). 
##' For further details see the package vignette. 
##'
##' @param ... arguments passed to other methods. The first argument should be a \code{formula} followed by a \code{data.frame};
##' alternatively, as a shortcut, you can omit the \code{formula} and provide only a \code{matrix} or 
##' a \code{data.frame}: in that case, the function will automatically consider the first column as the dependent variable 
##' and the rest as the independent variables.
##' @return 
##' \code{sValues} returns an object a list of class "sValues" containing the main results of the analysis:
##' 
##' \itemize{ 
##' \item \code{info}: a \code{list} with the general information about the paramaters used in the analysis, such as the 
##' formula, the data, the bounds and favorite variables.
##' 
##' \item \code{simple}: a \code{list} with the results of the simple linear regressions for each variable.
##' 
##' \item \code{all}: the results of the linear regression with all variables.
##' 
##' \item \code{bayes}: a \code{list} with the results of the bayesian regression for each combination of the R2 bounds. 
##' Each bayesian regression includes the coefficient estimates, the variance-covariance matrix and the t-values.
##' 
##' \item \code{ext_bounds}: a \code{list} with the extreme bounds estimates for each combination of the R2 bounds.
##' 
##' \item \code{s_values}: a \code{data.frame} with the s_values for each combination of the R2 bounds.
##' }
##' 
##' @examples 
##' # growth regressions example
##' ## All variables, No favorites
##' data(economic_growth)
##' eg_sv <- sValues(GR6096 ~ ., data = economic_growth)
##' eg_sv # prints results
##' plot(eg_sv, R2_bounds = c(0.5, 1))
##' plot(eg_sv, type = "beta_plot", variable = "P60", error_bar = TRUE)
##' coefs_eg <- coef(eg_sv) # extract coefficients
##' coefs_eg
##' 
##' ##  only 14 variables
##' eg_sv_14 <-  sValues(GR6096 ~GDPCH60L + OTHFRAC + ABSLATIT + 
##'                        LT100CR + BRIT + GOVNOM1 + WARTIME + 
##'                        SCOUT + P60 + PRIEXP70 + OIL + 
##'                        H60 + POP1560 + POP6560, data = economic_growth)
##' eg_sv_14
##' coefs_eg_14 <- coef(eg_sv_14)
##' 
##' ## With 14 favorites among all variables
##' favorites <- c("GDPCH60L", "OTHFRAC", "ABSLATIT", "LT100CR", 
##'               "BRIT", "GOVNOM1", "WARTIME", "SCOUT", 
##'               "P60", "PRIEXP70", "OIL", "H60", 
##'               "POP1560", "POP6560")
##' eg_sv_fav <- sValues(GR6096 ~ ., data = economic_growth, R2_bounds = c(0.5, 1),
##'                     favorites = favorites, R2_favorites = c(0.4, 0.8))
##' eg_sv_fav
##' plot(eg_sv_fav, R2_bounds = c(0.5, 1))
##' plot(eg_sv_fav, type = "beta_plot", variable = "P60", error_bar = TRUE)
##' coefs_eg_fav <- coef(eg_sv_fav)
##' coefs_eg_fav
##' 
##' @seealso 
##' \code{\link{coef.sValues}} to extract coefficients or statistics; 
##' 
##' \code{\link{print.sValues}} for printing;
##' 
##' \code{\link{summary.sValues}} for summaries;
##' 
##' \code{\link{plot.sValues}} for plots.
##' 
##' @references
##'  
##' Leamer, E. (2014). S-values: Conventional context-minimal measures of the sturdiness of regression coefficients. Working Paper
##' 
##' Leamer, E. (2015). S-values and bayesian weighted all-subsets regressions. European Economic Review.
##' 
##' @export
### sValues ###
sValues <- function(..., R2_bounds = c(0.1, 0.5, 1), 
                    favorites = NULL, R2_favorites = NULL, scale = TRUE){
  UseMethod("sValues")
}

##' @param formula an object of the class \code{\link{formula}}: a symbolic description of the model to be fitted. 
##' @param data needed only when you pass a formula as first parameter. An object of the class \code{\link{data.frame}} containing the variables used in the analysis. 
##' @param R2_bounds a numeric vector with two or more R2 bounds to be considered in the analysis. The default values are
##'  \code{c(0.1, 0.5, 1)}, proposed by Leamer (2014).
##' @param favorites \emph{optional} - a character vector that specifies the "favorite" varibles to be used in the analysis.
##' These variables will have different lower and upper R2 bounds as defined in the \code{R_favorites} argument.
##' @param R2_favorites \emph{optional} - a numeric vector with two or more R2 bounds for the "favorite" variables.
##' @param scale should the variables be scaled/standardized to zero mean and unit variance? 
##' The default is \code{TRUE}. If your data is already scaled/standardized you should set this to \code{FALSE}.
##' @export
##' @name sValues
##' @importFrom stats lm
sValues.formula <- function(formula, data, R2_bounds = c(0.1, 0.5, 1), 
                            favorites = NULL, R2_favorites = NULL, scale = TRUE, ...){
  
  stopifnot(class(formula)=="formula",
            class(data) == "data.frame",
            length(R2_bounds) >= 2)
  
  if(!is.null(R2_favorites) && length(R2_bounds) != length(R2_favorites))
    stop("R2 bounds and R2 favorites must have the same length")
  
  #-- save call information --#
  info                   <- list(formula     = formula, 
                                 data        = deparse(substitute(data)),
                                 R2_bounds   = R2_bounds,
                                 favorites   = favorites,
                                 R2_favorites = R2_favorites)
  
  #-- organize R2 bounds --#
  
  #-- regular bounds
  R2_bounds              <- sort(R2_bounds)
  r2_combs               <- r2_combs(R2_bounds)
  
  #-- Favorites
  if(!is.null(favorites)){
    R2_favorites         <- sort(R2_favorites)
    r2_combs_fav         <- r2_combs(R2_favorites)
  }else{
    R2_favorites         <- vector("list", length(R2_bounds))
    r2_combs_fav         <- vector("list", length(r2_combs))
  }
  
  
  #-- model matrix and scale data --#
  new_df                 <- model_frame_scale(formula, data, scale = scale)
  
  #-- ols simple --#
  singles_formula        <- lapply(names(new_df)[-1],  function(i) formula(new_df[c(deparse(formula[[2]]),i)]))
  names(singles_formula) <- names(new_df)[-1]
  simple                 <- lapply(singles_formula, function(x) lm(x, new_df)) 
  
  #-- ols with all variables --#
  all                    <- lm(formula, new_df)
  bayes                  <- Map(function(x, y) .bayes(all, x, favorites,y), R2_bounds, R2_favorites)
  names(bayes)           <- paste0("R2_", R2_bounds) 
  
  #-- extreme bounds --#
  ext_bounds             <- Map(function(r2_bounds, r2_fav) lapply(1:(ncol(new_df)-1), 
                                                                   function(x) .bounds(all, x, r2_bounds, favorites, r2_fav)), 
                                r2_combs, r2_combs_fav)
  ext_bounds             <- lapply(ext_bounds, function(x) do.call("rbind", x))
  ext_bounds             <- lapply(ext_bounds, function(x){row.names(x) <- names(new_df)[-1];x})
  
  #-- s-values --#
  s_values               <- lapply(ext_bounds, function(x) apply(x, 1, .s_value))
  s_values               <- do.call("cbind", s_values)
  s_values               <- data.frame(s_values)
  row.names(s_values)    <- names(new_df)[-1]
  names(s_values)        <- paste0("s_", names(s_values))
  
  #-- return --#
  ret <- list(info = info,
              simple = simple, 
              all = all, 
              bayes = bayes, 
              ext_bounds = ext_bounds, 
              s_values = s_values)
  
  class(ret) <- "sValues"
  
  return(ret)
}



##' @param m an object of class \code{\link{matrix}} with the dependent variable in the first column 
##' followed by the covariates. 
##' The matrix must have column names.
##' @export
##' @name sValues
sValues.matrix <- function(m, R2_bounds = c(0.1, 0.5, 1), 
                           favorites = NULL, R2_favorites = NULL, scale = TRUE, ...){
  if(!length(colnames(m)>0)) stop("Matrix must have column names!")
  df <- as.data.frame(m)
  res <- sValues(df, R2_bounds = R2_bounds, favorites = favorites, 
                 R2_favorites = R2_favorites, scale = scale)
  res$info$data <- deparse(substitute(m))
  res
}

##' @param df an object of class \code{\link{data.frame}} with the dependent variable in the first column 
##' followed by the covariates.
##' @export
##' @name sValues
##' @importFrom stats as.formula
sValues.data.frame <- function(df, R2_bounds = c(0.1, 0.5, 1), 
                               favorites = NULL, R2_favorites = NULL, scale = TRUE, ...){
  formula <- as.formula(paste(names(df)[1], "~ ."))
  res <- sValues(formula = formula, data = df, R2_bounds = R2_bounds, 
                 favorites = favorites, R2_favorites = R2_favorites, scale = scale)
  res$info$data <- deparse(substitute(df))
  res
}




# Coef --------------------------------------------------------------------


##' Extract sValues Model Coefficients/Statistics 
##' 
##' @param object an object of class \code{\link{sValues}}.
##' @param type which coefficient/statistic to extract? Current options are "betas", "t_values", 
##' "s_values", "extreme_bounds" and "default". See details.
##' @param ... further arguments passed to or from other methods.
##' @details 
##' 
##' For the \code{coef} function, the default is to extract the beta coefficients, t-values and s-values. You can can get 
##' each one of those individually by setting \code{type} to either "betas", "t_values" or "s_values". 
##' You can also get the extreme bounds of the estimates by setting \code{type} to "extreme_bounds". 
##' Finally, you can set \code{type = "all"} to get everything. 
##' 
##' For each option of \code{coef}, there is an alternative helper function with the same name. 
##' That is, \code{coef(x, "betas")} is equivalent to \code{betas(x)}, or \code{coef(x, "extreme_bounds")} is equivalent
##' to \code{extreme_bounds(x)}. 
##' 
##' @return 
##' The function returns a \code{data.frame} with the estimates for each variable.
##' 
##' @examples
##' data(economic_growth)
##' eg_sv <- sValues(GR6096 ~ ., data = economic_growth)
##' eg_betas <- coef(eg_sv, "betas")
##' eg_t_values <- coef(eg_sv, "t_values")
##' eg_s_values <- coef(eg_sv, "s_values")
##' eg_ext_bounds <- coef(eg_sv, "extreme_bounds")
##' 
##' # get sturdy estimates for R2 bounds 0.5 - 1
##' eg_s_values[abs(eg_s_values[3]) > 1, 3, drop = FALSE]
##' 
##' 
##' @seealso 
##' \code{\link{summary.sValues}}.
##' 
##' @export
##' @name coef.sValues
coef.sValues <- function(object, type = "default", ...){
  switch(type, 
         betas = betas(object),
         t_values = t_values(object),
         s_values = s_values(object),
         extreme_bounds = extreme_bounds(object),
         default = cbind(betas(object), t_values(object), s_values(object))[order(abs(s_values(object)[,1]), decreasing = TRUE),],
         all = cbind(betas(object), t_values(object), s_values(object), extreme_bounds(object))[order(abs(s_values(object)[,1]), decreasing = TRUE),],
         stop("Invalid argument 'type'. See ?coef.Svalues for valid options."))
}


##' @export
##' @name coef.sValues
betas <- function(object){
  if(class(object) != "sValues") stop("Object must be of class 'sValues'")
  #simple
  s <- do.call("rbind", lapply(object$simple, function(object) object$coefficients[-1]))
  # bayes
  b <- do.call("cbind", lapply(object$bayes, function(object) object$coefficients))
  # all
  a <- coef(object$all)[-1]
  # combining
  data.frame(ols_simple = c(s), b, ols_all = a)
}


##' @export
##' @name coef.sValues
##' @importFrom stats coef
t_values <- function(object){
  if(class(object) != "sValues") stop("Object must be of class 'sValues'")
  #simple
  st <- do.call("rbind", lapply(object$simple, function(object) coef(summary(object))[-1,"t value"]))
  #bayes
  bt <- do.call("cbind", lapply(object$bayes, function(object) object$t))
  #all
  at <- coef(summary(object$all))[-1, "t value"]
  # combining
  data.frame(t_ols_simple = c(st), bt, t_ols_all = at) 
}


##' @export
##' @name coef.sValues
s_values <- function(object){
  if(class(object) != "sValues") stop("Object must be of class 'sValues'")
  object$s_values
}


##' @export
##' @name coef.sValues
extreme_bounds <- function(object){
  if(class(object) != "sValues") stop("Object must be of class 'sValues'")
  do.call("cbind", object$ext_bounds)
}


# Print -------------------------------------------------------------------



##' Succinct display of S-values results.
##'  
##' Succinct display of S-values results.
##' 
##' @param x an object of class \code{\link{sValues}}.
##' @param ... further arguments passed to or from other methods.
##' @param print.length how many variables to show in the screen?
##' This is used for pretty printing. The default is 6.
##' 
##'  @return 
##'  \code{NULL}
##'  
##' @examples 
##' data(economic_growth)
##' eg_sv <- sValues(GR6096 ~ ., data = economic_growth)
##' eg_sv
##' str(eg_sv) 
##'  
##' @export
print.sValues <- function(x, ..., print.length = 6){
  with(x$info,{
    cat("Data: ", data,",", "    Formula: ", deparse(formula), sep = "")
    cat("\nR2 bounds:", paste(R2_bounds, collapse=" - "))
    
    if(!is.null(favorites)){
      cat("\nFavorites:", text_break(favorites, print.length = print.length))
      cat("\nR2 favorites:", paste(R2_favorites, collapse = " - "))
    }
    
    cat("\n")
    cat("\nabs(S-value) > 1:\n")
    indices <- lapply(x$s_values, function(x) which(abs(x)>1))
#     indices2 <- which( abs(coef(x, "t_values")$t_ols_all)>2)
#     indices <- lapply(indices, function(x) intersect(x, indices2))
    e <- lapply(indices, function(i) row.names(x$s_values)[i])
    
    for(i in seq_along(e)){
      cat(" ", 
          paste0(gsub("s_(R2)_([0-9]\\.?[0-9]?)_([0-9]\\.?[0-9]?)", 
                      "\\1 (\\2, \\3)", names(e[i])), ":"), 
          if(length(e[[i]])> 0){
            text_break(e[[i]], print.length = print.length)
            }else{
              "None"
              }, "\n")
    }
    
    cat("\nabs(t-value) > 2:\n")
    indices<- lapply(coef(x, "t_values")[-1], function(x) which(abs(x)>2))
    e <- lapply(indices, function(i) row.names(x$s_values)[i])
    for(i in 1:(length(e)-1)){
      cat(" ", 
          paste0(gsub("t_b(ayes|ols)_([0-9].?[0-9]?)?|all", 
                      "B\\1ian (R2 = \\2)", names(e[i])), ":"), 
          if(length(e[[i]])> 0){
            text_break(e[[i]], print.length = print.length)
            }else{
              "None"
              }, "\n")
    }
    cat("  Unconstrained OLS:", text_break(e[[length(e)]],print.length = print.length), "\n")
  })
}


##' str sValues
##' 
##' \code{str} method for \code{sValues}. 
##' 
##' @param object an object of class \code{\link{sValues}}.
##' @param max.level maximal level of nesting which is 
##' applied for displaying nested structures. Default is 1.
##' @param ... further arguments passed to or from other methods.
##' @export
##' @importFrom utils str
str.sValues <- function(object, max.level = 1, ...){
  str(unclass(object), max.level = max.level, ...)
}


# Summary -----------------------------------------------------------------

##' summary sValues
##' 
##' For now, this function is equivalent to \code{\link{print.sValues}}.
##' 
##' 
##' @param object an object of class \code{\link{sValues}}.
##' @param ... further arguments passed to or from other methods. 
##' @export
summary.sValues <- function(object, ...){
  print(object)
}




# Plots -------------------------------------------------------------------


##' Plot method for S-values
##' 
##' @description 
##' 
##' Plot methods for objects of the class \code{\link{sValues}}.
##' 
##' @param x an object of class \code{\link{sValues}}.
##' @param type the type of the plot. Current options are \emph{t_s_plot} which returns
##' a scatterplot of s-values vs t-values for all coefficients and \emph{beta_plot} which returns
##' a plot of the different estimates for the coefficients. 
##' @param ... aditional arguments to be passed to the plot functions. See details.
##' 
##' @details 
##' 
##' Addional arguments:
##' 
##' \code{t_s_plot} 
##' \itemize{
##'     \item \code{R2_bounds}: a numeric vector of length two specifying which R2 bounds range to plot.
##' }
##' 
##' \code{beta_plot}
##' \itemize{
##'     \item \code{variables}: a character vector specifying which variables to plot. Default is "all".
##'     \item \code{error_bar}: should the error bars be plotted? Default is \code{FALSE}.
##'     \item \code{ext_bounds_shades}: should shades representing the extreme bounds be plotted? Default is \code{FALSE}.
##' }
##' 
##' @return 
##' It returns a \code{\link{ggplot}} object with the requested plot.
##' 
##' @examples
##' # growth regressions example
##' data(economic_growth)
##' eg_sv <- sValues(GR6096 ~ ., data = economic_growth)
##' plot(eg_sv, R2_bounds = c(0.5, 1))
##' plot(eg_sv, R2_bounds = c(0.1, 1))
##' plot(eg_sv, type = "beta_plot", variable = "OPENDEC1", error_bar = FALSE)
##' plot(eg_sv, type = "beta_plot", variable = "OPENDEC1", error_bar = TRUE)
##' 
##' @export
plot.sValues <- function(x, type = "t_s_plot", ...){
  switch(type,
         t_s_plot = t_s_plot(x, ...),
         beta_plot = beta_plot(x, ...),
         stop("Invalid argument 'type'. See ?plot.Svalues for valid options."))
}


##' @import ggplot2 reshape2
##' @importFrom stats coef
t_s_plot <- function(x, R2_bounds = NULL){
  data <- coef(x)
  if(!is.null(R2_bounds)){
    rangeR2 <- paste0("s_R2_",paste(R2_bounds, collapse = "_"))
    if(!any(names(data) %in% rangeR2)) stop("R2 bounds not in data")
  }else{
    rangeR2 <- names(data)[length(data)] 
    R2_bounds <- gsub(".*([0-9].[0-9])_([0-9])", "\\1 - \\2", rangeR2)
    warning(paste("R2 bounds", gsub(".*([0-9].[0-9])_([0-9])", "c(\\1, \\2)", rangeR2),
                  "chosen as default. If you want to change this, set the 'R2_bounds' parameter."))
  }
  
  rangeY <- range(data[,rangeR2])
  rangeY <- c(min(rangeY, -3), max(rangeY, 3))
  rangeX <- range(data[,"t_ols_all"])
  rangeX <- c(min(rangeX, -3), max(rangeX,3))
  ggplot(data, aes_string("t_ols_all",rangeR2)) + 
    annotate("rect", xmin=-2,xmax=2,ymin=-Inf,ymax=Inf, fill="gray", alpha=0.5) +
    annotate("rect", ymin=-1,ymax=1,xmin=-Inf,xmax=Inf,fill="gray", alpha=0.5) + 
    geom_point(shape=5, size=3) + geom_smooth(method="lm", se=F) + theme_bw() +
    ylab(paste("s-values\n R2 range", paste("(", paste(R2_bounds, collapse = " - "), ")"))) +
    xlab("t-values\n Unconstrained OLS regression") +
    ggtitle("S-values x t-values \nUncertain and Fragile Estimates Shaded") +
    ylim(rangeY*1.2)+
    xlim(rangeX*1.2)
}  


##' @import ggplot2 reshape2
##' @importFrom stats coef
beta_plot <- function(x, variables = "all", error_bar = FALSE, ext_bounds_shades = FALSE){
  
  betas <- .melt.coef(x, "betas")
  
  t_values <- .melt.coef(x, "t_values")
  
  ext_bounds <-coef(x, "extreme_bounds")
  ext_bounds$coefficients <- row.names(ext_bounds)
  ext_bounds<- do.call("rbind", replicate(ncol(coef(x, "betas")), ext_bounds, simplify = FALSE))
  
  melted <- cbind(betas, t_values, ext_bounds)
  
  if(variables[1] != "all"){ 
    if(!any(melted$coefficients %in% variables)) stop("Invalid variable names.")
    melted <- melted[melted$coefficients %in% variables,]
  }
  
  p <- ggplot(melted, aes(variable, betas, group=coefficients)) + 
    geom_line(aes(variable, betas, group=coefficients), col = "transparent")
  
  if(ext_bounds_shades){
    names <- grep("^R2", names(melted), value = TRUE)
    for(i in 1:(length(names)/2)){
      p <- p + geom_rect(aes_string(ymin=names[2*i - 1],ymax=names[2*i],xmin=-Inf,xmax=Inf), fill="gray", alpha=0.1)
    }
  }
  
  p <- p + geom_line(aes(variable, betas, group=coefficients)) +
    facet_wrap(~coefficients, scales="free_y") + 
    ylab("Coefficient Values") +
    xlab("Prior Expected R2") +
    ggtitle("Estimated Beta-Coefficients") + theme_bw() 
  
  if(error_bar)
    p <- p + geom_errorbar(aes(ymin=betas-1.96*(betas/t_values), 
                               ymax=betas+1.96*(betas/t_values)),
                           colour="black", width=.1, linetype = 2) +
      geom_hline(yintercept=0, col="blue", linetype="dotdash", size=1)
  p  
}

