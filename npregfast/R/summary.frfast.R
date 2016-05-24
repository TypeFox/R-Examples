#' Summarizing fits of \code{frfast} class
#' @aliases print.frfast
#' @description Takes a fitted \code{frfast} object produced by \code{frfast()}
#' and produces various useful summaries from it.
#' @param object a fitted \code{frfast} object as produced by \code{frfast()}.
#' @param \ldots additional arguments affecting the predictions produced.
#' @details \code{print.frfast} tries to be smart about \code{summary.frfast}.
#' @return \code{summary.frfast} computes and returns a list of summary 
#' information for a fitted \code{frfast} object.
#' \item{model}{type of model: nonparametric or allometric.}
#' \item{smooth}{type of smoother: kernel or splines.}
#' \item{h}{the kernel bandwidth smoothing parameter.}
#' \item{dp}{degree of the polynomial.}
#' \item{nboot}{number of bootstrap repeats.}
#' \item{kbin}{number of binning nodes over which the function is to be estimated.}
#' \item{n}{sample size.}
#' \item{fmod}{factor's levels.}
#' \item{coef}{if \code{model = "allo"}, coefficients of the model.}
#' 
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' 
#' @references 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#' 
#' 
#' @examples
#' library(npregfast)
#' data(barnacle)
#' 
#' # Nonparametric regression without interactions
#' fit <- frfast(DW ~ RC, data = barnacle, nboot = 100) 
#' fit
#' summary(fit)
#' 
#' # Nonparametric regression with interactions
#' fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 100)
#' fit2
#' summary(fit2)
#' 
#' # Allometric model
#' fit3 <- frfast(DW ~ RC, data = barnacle, model = "allo", nboot = 100)
#' fit3
#' summary(fit3)
#' 
#' @export



summary.frfast <- function(object = model, ...) {
  model <- object
  if (missing(model)) 
    stop("Argument 'model' is missing with no default")
  
  if (model$nmodel == 1) {
    
    cat("\nCall:\n")
    print(model$call)
    cat("", "\n")
    if (model$nmodel == 1) {
      m <- "Nonparametric"
    } else {
      m <- "Alometric"
    }
    cat("*********************************************", "\n")
    cat(m, "Model", "\n")
    cat("*********************************************", "\n")
   
    
    if (model$smooth == "kernel"){
      cat("\nType of nonparametric smoother: ",format(model$smooth), "\n")
     if (model$kernel == 1) 
      cat("Kernel: Epanechnikov \n")
    if (model$kernel == 2) 
      cat("Kernel: Triangular \n")
    if (model$kernel == 3) 
      cat("Kernel: Gaussian \n")
    
    if (model$nf != 1) {
      cat("Bandwidth:", format(c(model$h0, model$h), digits = 2), "\n")
    } else {
      cat("Bandwidth:", format(c(model$h0), digits = 2), "\n")
    }
    cat("Polynomial degree:", model$dp, "\n")
    cat("Number of bootstrap repeats:", model$nboot, "\n")
    cat("Number of binning nodes", model$kbin, "\n")
    cat("", "\n")
    cat("", "\n")
    cat("The number of data is: ", model$n, "\n")
    }else{
      cat("\nType of nonparametric smoother: ",format(model$smooth), "\n")
      cat("Number of bootstrap repeats:", model$nboot, "\n")
      cat("", "\n")
      cat("The number of data is: ", model$n, "\n")
    }
    # cat('The factor's levels are: ',etiquetas<-model$etiquetas, '\n')
    
    nf <- length(model$label)
    if (nf != 1) {
      cat("The factor's levels are: ", etiquetas <- model$label, "\n")
      for (factor in 1:nf) {
        if(model$smooth == "splines") nn <- length(model$xdata[model$fmod == etiquetas[factor]])
        if(model$smooth == "kernel") nn <- length(model$xdata[model$fmod == factor])
        cat("The number of data for the level", etiquetas[factor], "is:", nn, "\n")
      }
      cat("", "\n")
      cat("Summaries for the response variable (for each level): ")
      for (factor in c(1:nf)) {
        cat("", "\n")
        cat("Level", etiquetas[factor], ":", "\n")
        if(model$smooth == "kernel") print(summary(model$ydata[model$fmod == factor]))
        if(model$smooth == "splines") print(summary(model$ydata[model$fmod == etiquetas[factor]]))
      }
      
    } else {
      cat("", "\n")
      cat("Summaries for the response variable: \n")
      print(summary(model$ydata))
    }
    
    
  }
  
  
  if (model$nmodel == 2) {
    etiquetas <- model$label
    nf <- length(etiquetas)
    
    
    cat("\nCall:\n")
    print(model$call)
    cat("", "\n")
    if (model$nmodel == 1) {
      m <- "Nonparametric"
    } else {
      m <- "Allometric"
    }
    cat("*********************************************", "\n")
    cat(m, "Model", "\n")
    cat("*********************************************", "\n")
    
    
    
    
    cat("", "\n")
    if (nf == 1) {
      cat("Coefficients:", "\n")
    } else {
      cat("Coefficients (for each level):", "\n")
    }
    cat("", "\n")
    for (factor in c(1:nf)) {
      if (nf > 1) 
        cat("Level", etiquetas[factor], ":", "\n")
      aux <- c(model$a[factor], model$al[factor], model$au[factor], model$b[factor], 
               model$bl[factor], model$bu[factor])
      aux <- round(aux, 6)
      tabla <- matrix(aux, ncol = 3, nrow = 2, byrow = T)
      colnames(tabla) <- c("", "2.5 %", "97.5 %")
      rownames(tabla) <- c("a", "b")
      print(tabla)
      cat("", "\n")
    }
    cat("Adjusted R-squared: ", model$r2, "\n")
    cat("", "\n")
    cat("*********************************************", "\n")
    
    
    cat("", "\n")
    cat("", "\n")
    cat("Polynomial degree:", model$dp, "\n")
    cat("Number of bootstrap repeats:", model$nboot, "\n")
    cat("Number of binning nodes", model$kbin, "\n")
    cat("", "\n")
    cat("", "\n")
    cat("The number of data is: ", model$n, "\n")
    cat("The factor's levels are: ", etiquetas <- model$label, "\n")
    
    
    if (nf != 1) {
      for (factor in c(1:nf)) {
        cat("The number of data for the level", etiquetas[factor], "is:", 
            length(model$xdata[model$fmod == factor]), "\n")
      }
    }
    
    cat("", "\n")
    if (nf > 1) {
      cat("Summaries for the variable y (for each level): ")
    } else {
      "Summaries for the variable y:"
    }
    
    for (factor in c(1:nf)) {
      cat("", "\n")
      if (nf > 1) 
        cat("Level", etiquetas[factor], ":", "\n")
      
      print(summary(model$ydata[model$fmod == factor]))
    }
    
  }
} 
