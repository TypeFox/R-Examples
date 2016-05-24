#' @export
print.FuzzyPairwiseComparisonMatrix  <- function(x, ...){
  textMatrix = matrix(data = "",  nrow = nrow(x@fnMin), ncol = ncol(x@fnMin))

  for(i in 1:nrow(x@fnMin)){
    for(j in 1:ncol(x@fnMin)){
      textMatrix[i,j] = paste("(", round(x@fnMin[i,j], digits = 4), ";",
                              round(x@fnModal[i,j], digits = 4), ";",
                              round(x@fnMax[i,j], digits = 4), ")", sep = "")
    }
  }

  print(textMatrix)
}

#' @export
print.PairwiseComparisonMatrix  <- function(x, ...){
  textMatrix = matrix(data = "",  nrow = nrow(x@valuesChar), ncol = ncol(x@valuesChar))

  for(i in 1:nrow(x@valuesChar)){
    for(j in 1:ncol(x@valuesChar)){
      textMatrix[i,j] = x@valuesChar[i,j]
    }
  }

  print(textMatrix)
}

#' @export
print.Weights <- function(x, ...){
  print(round(x@weights, digits = 4))
}

#' @export
print.FuzzyWeights <- function(x, ...){

  textMatrix = matrix(data = "",  nrow = length(x@fnMin), ncol = 1)
  rowNames = c()
  for(i in 1:length(x@fnMin)){
    textMatrix[i,1] = paste("(", round(x@fnMin[i], digits = 4), ";",
                            round(x@fnModal[i], digits = 4), ";",
                            round(x@fnMax[i], digits = 4), ")", sep = "")

    rowNames = append(rowNames, paste("w",i, sep = ""))
  }
  rownames(textMatrix) <- rowNames
  print(textMatrix)
}

#' @export
print.FuzzyData <- function(x, ...){

  textMatrix = matrix(data = "",  nrow = nrow(x@fnMin), ncol = ncol(x@fnMin))

  for(i in 1:nrow(x@fnMin)){
    for(j in 1:ncol(x@fnMin)){
      textMatrix[i,j] = paste("(", round(x@fnMin[i,j], digits = 4), ";",
                                   round(x@fnModal[i,j], digits = 4), ";",
                                   round(x@fnMax[i,j], digits = 4), ")", sep = "")

    }
  }

  print(textMatrix)
}
