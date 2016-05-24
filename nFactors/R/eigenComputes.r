eigenComputes <-
function(x, cor=TRUE, model="components", ...) {
 dataType <- eigenFrom(x)

 if (model == "components") {
  res <- switch(dataType,
   eigenvalues = as.vector(x),
   correlation = {if (cor == FALSE) eigen(x)$values           else  eigen(cov2cor(x))$values},
   data        = {if (cor == TRUE)  eigen(cor(x, ...))$values else  eigen(cov(x, ...))$values}
   )
  }
  
 if (model == "factors") {
  res <- switch(dataType,
   eigenvalues = as.vector(x),
   correlation = {if (cor == FALSE) eigen(corFA(x, method="ginv"))$values else   eigen(cov2cor(corFA(x, method="ginv")))$values},
   data        = {if (cor == TRUE)  eigen(corFA(cor(x, ...),method="ginv"))$values else  eigen(corFA(cov(x, ...),method="ginv"))$values}
   )
  }
 return(res)
 }
