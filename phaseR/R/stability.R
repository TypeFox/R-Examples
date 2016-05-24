stability <- function(deriv, y.star = NULL, parameters = NULL,
                      system = "two.dim", h = 1e-07, summary = TRUE){
  if (!(system %in% c("one.dim", "two.dim"))){
    stop("system must either be set to one.dim or two.dim")
  }
  if (is.null(y.star)){
    y.star <- locator(n = 1)
    if (system == "two.dim"){
      y.star <- c(y.star[1], y.star[2])
    }
    if (system == "one.dim"){
      y.star <- y.star[2]
    }
  }
  if (h <= 0){
    stop("h is less than or equal to zero")
  }
  if (!is.logical(summary)){
    stop("summary must either be set to TRUE or FALSE")
  }
  if (system == "one.dim"){
    discriminant <- as.numeric(deriv(0, y.star + h, parameters = parameters))/h
    if (discriminant > 0){
      classification <- "Unstable"
    }
    if (discriminant < 0){
      classification <- "Stable"
    }
    if (discriminant == 0){
      classification <- "Indeterminate"
    }
  }
  if (system == "two.dim"){
    df <- deriv(t = 0, y = y.star, parameters = parameters)
    dfx <- deriv(t = 0, y = c(y.star[1] + h, y.star[2]), 
                 parameters = parameters)
    dfy <- deriv(t = 0, y = c(y.star[1], y.star[2] + h), 
                 parameters = parameters)
    A <- (as.numeric(dfx[[1]][1]) - as.numeric(df[[1]][1]))/h
    B <- (as.numeric(dfy[[1]][1]) - as.numeric(df[[1]][1]))/h
    C <- (as.numeric(dfx[[1]][2]) - as.numeric(df[[1]][2]))/h
    D <- (as.numeric(dfy[[1]][2]) - as.numeric(df[[1]][2]))/h
    Delta <- A*D - B*C
    tr <- A + D
    discriminant <- tr^2 - 4*Delta
    jacobian <- matrix(c(A, B, C, D), 2, 2, byrow = TRUE)
    eigen <- eigen(jacobian)
    eigenvalues <- eigen$values
    eigenvectors <- eigen$vectors
    if (Delta < 0){
      classification <- "Saddle"
    }
    if (Delta == 0){
      classification <- "Indeterminate"
    }
    if (Delta > 0){
      if (discriminant > 0){
        if (tr < 0){
          classification <- "Stable node"
        }
        if (tr > 0){
          classification <- "Unstable node"
        }
      }
      if (discriminant < 0){
        if (tr < 0) {
          classification <- "Stable focus"
        }
        if (tr > 0){
          classification <- "Unstable focus"
        }
        if (tr == 0){
          classification <- "Centre"
        }
      }
    }
  }
  if (summary == TRUE){
    if (system == "one.dim"){
      cat("\nDiscriminant:", discriminant, "  Classification:", 
          classification)
    }
    if (system == "two.dim"){
      cat("\nT:", tr, "  Delta:", Delta, "  Discriminant:", 
        discriminant, "  Classification:", classification)
    }
  }
  output                <- list()
  output$classification <- classification
  if (system == "two.dim"){
    output$Delta        <- Delta
  }
  output$deriv          <- deriv
  output$discriminant   <- discriminant
  if (system == "two.dim"){
    output$eigenvalues  <- eigenvalues
    output$eigenvectors <- eigenvectors
  }
  output$h              <- h
  if (system == "two.dim"){
    output$jacobian     <- jacobian
  }
  output$parameters     <- parameters
  output$system         <- system
  if (system == "two.dim"){
    output$tr           <- tr
  }
  output$y.star         <- y.star
  return(output)
}