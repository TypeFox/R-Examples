#########################################
#                                       #
#  White's Test for Heteroskedasticity  #
#        Based on Doornik (1996)        #
#                                       #
#  Implemented by Sebastian Andersson   #
#  sebastian.andersson@statistik.uu.se  #
#                                       #
#     Currently supports VARs with:     #
#   - intercept                         #
#   - trend                             #
#   - intercept and trend               #
#                                       #
#     Support for:                      #
#   - neither intercept nor trend       #
#   - exogenous variables               #
#     is yet to be implemented.         #
#                                       #
#    Last Updated February 25, 2013     #
#                                       #
#########################################

# First we create a new class called "whitetest"
setClass("whitetest", representation("list"))

# Specify the appearance of the output
setMethod("show", "whitetest", function(object) {
  text1 <- "White's Test for Heteroskedasticity:"
  cat(paste("\n", text1, "\n", sep = ""))
  row <- paste(rep("=", nchar(text1)), collapse = "")
  cat(row, "\n")
  cat("\n")
  cat(" No Cross Terms\n")
  cat("\n")
  cat(" H0: Homoskedasticity\n")
  cat(" H1: Heteroskedasticity\n")
  cat("\n")
  cat(" Test Statistic:\n")
  cat("", sprintf("%.4f", object$statistic), "\n")
  cat("\n")
  cat(" Degrees of Freedom:\n")
  cat("", object$degrees, "\n")
  cat("\n")
  cat(" P-value:\n")
  cat("", sprintf("%.4f", object$p.value), "\n")
})

# Create the function for the test
whites.htest <- function(var.model) {
  
  if (class(var.model) != "varest") {
    stop("The VAR model must be of class varest.")
  }
  
  if ("exogen" %in% names(var.model$call) && class(var.model$call$exogen) != "NULL") {
    stop("The script does not yet support VAR models with exogenous variables.")
  }
  
  if ("season" %in% names(var.model$call) && class(var.model$call$season) != "NULL") {
    stop("The script does not yet support VAR models with seasonal variables.")
  }
  
  if (class(var.model$restrictions) != "NULL") {
    stop("The script does not yet support VAR models with restrictions.")
  }
  
  # Create a matrix with the residuals from the VAR model
  residual.matrix <- residuals(var.model) 
  
  # Create an empty matrix to store cross products of residuals
  res.products <- matrix(numeric(0), nrow = nrow(residual.matrix), ncol = (sum(1:ncol(residual.matrix))-ncol(residual.matrix)))
  counter <- 1 # A counter used in the loop

  # Save the lag order and dimension
  lag.order <- var.model$p
  var.dimension <- var.model$K

  # Two empty vectors
  namevector <- c()
  temp.namevector <- c()

  # A temporary matrix for cross products when i=j
  temp.matrix <- matrix(numeric(0), nrow = nrow(residual.matrix), ncol = ncol(residual.matrix))

  # Two loops which create matrices with residual cross products, one for i=j and one for i!=j
  for (i in 1:ncol(residual.matrix)) {
    for (j in i:ncol(residual.matrix)) {
      if (i==j) {
        temp.matrix[,i] <- residual.matrix[,i]*residual.matrix[,j]
        temp.namevector <- c(temp.namevector, paste(i, "*", j, sep=""))
      } else {
      res.products[,counter] <- residual.matrix[,i]*residual.matrix[,j]
      namevector <- c(namevector, paste(i, "*", j, sep=""))
      counter <- counter + 1
      }
    }
  }

  # Merge the two matrices. The columns in the new matrix will be i=j followed by i!=j
  res.products <- cbind(temp.matrix, res.products)
  colnames(res.products) <- colnames(res.products, do.NULL = FALSE)
  colnames(res.products) <- c(temp.namevector, namevector)

  # Subtract the mean from the residual cross product series
  for (i in 1:ncol(res.products)) {
    res.products[,i] <- res.products[,i]-mean(res.products[,i])
  }

  # The matrix rcov is:
  # - transpose of the residual products matrix times the residual products matrix
  # - then divided by the difference between the number of residuals per equation and the number of equations in the VAR
  # - lastly, the inverse is taken of this matrix
  rcov <- solve(t(res.products) %*% res.products/(nrow(residual.matrix)-ncol(residual.matrix)))

  # lagged.matrix is a matrix containing all the lagged series.
  if (var.model$type=="const") {
    lagged.matrix <- var.model$datamat[-ncol(var.model$datamat)]
  } else if (var.model$type=="both") {
    lagged.matrix <- var.model$datamat[-(ncol(var.model$datamat)-1)]
  } else {
    lagged.matrix <- var.model$datamat
  }
  lagged.matrix <- lagged.matrix[-1:-var.dimension]

  # num.lagged is the number of series in lagged.matrix
  num.lagged <- ncol(lagged.matrix)
  
  # Add the squares of the series in lagged.matrix to lagged.matrix
  for (i in 1:num.lagged) {
    lagged.matrix <- cbind(lagged.matrix, lagged.matrix[,i]^2)
    colnames(lagged.matrix)[num.lagged+i] <- paste("sq(", colnames(lagged.matrix)[i], ")", sep="")
  }
  
  # Update num.lagged to include the squares
  num.lagged <- ncol(lagged.matrix)

  #Create ucov and fill it with the residuals from the regression models when using lagged.matrix as the explanatory matrix and the squared VAR model residuals as response variables.
  ucov <- matrix(numeric(0), nrow = nrow(res.products), ncol = ncol(res.products))
  for (i in 1:ncol(res.products)) {
    ucov[,i] <- lm(res.products[,i] ~ ., data=lagged.matrix)$residuals
  }

  # Redefine ucov as its transpose times itself divided by the difference between n in the VAR residual series and the number of residual series
  ucov <- (t(ucov) %*% ucov)/(nrow(residual.matrix)-ncol(residual.matrix))
  
  # The matrix rucov is rcov times ucov
  rucov <- rcov %*% ucov

  # Compute the trace of rucov
  trace.cov <- sum(diag(rucov))

  # Compute the statistic
  stat <- unname(nrow(residual.matrix)*(ncol(res.products)-trace.cov))
  names(stat) <- "Test Statistic"
  
  # Compute the degrees of freedom
  dfree <- num.lagged*var.dimension*(var.dimension+1)/2
  names(dfree) <- "Degrees of Freedom"
  
  # Compute the p-value
  pval <- unname(1-pchisq(stat, dfree))
  names(pval) <- "p-value"
  
  call <- match.call()  
  names(call) <- "Call"
  
  result <- list(statistic = stat, degrees = dfree, p.value = pval, rcov = rcov, ucov = ucov, res.products = res.products, lagged.variables = lagged.matrix, call = call)
  new("whitetest", result)
  
}