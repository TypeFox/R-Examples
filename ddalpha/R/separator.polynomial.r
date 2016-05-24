
.ddalpha.learn.polynomial <- function(ddalpha){
  # Separating (calculating extensions and normals)
  counter <- 1
  # Determining multi-class behaviour
  if (ddalpha$methodAggregation == "majority"){
    for (i in 1:(ddalpha$numPatterns - 1)){
      for (j in (i + 1):ddalpha$numPatterns){
        # Creating a classifier
        polynomial <- .polynomial_learn_C(ddalpha$maxDegree, 
                                   rbind(ddalpha$patterns[[i]]$depths, 
                                         ddalpha$patterns[[j]]$depths), 
                                   ddalpha$patterns[[i]]$cardinality, 
                                   ddalpha$patterns[[j]]$cardinality, 
                                   ddalpha$numChunks, ddalpha$seed)
 # DEBUG       
 if (F){
        print(polynomial$coefficients)
        print(GetEmpiricalRisk (polynomial$coefficients,
                          rbind(ddalpha$patterns[[i]]$depths, ddalpha$patterns[[j]]$depths), 
                          ddalpha$patterns[[i]]$cardinality, 
                          ddalpha$patterns[[j]]$cardinality))
      }  
        # Adding the classifier to the list of classifiers
        ddalpha$classifiers[[counter]] <- 
          list(
            index          = counter,
            index0         = i,
            index1         = j,
            polynomial     = polynomial$coefficients,
            degree         = polynomial$degree,
            axis           = polynomial$axis)
        
        counter <- counter + 1
      }
    }
    ddalpha$numClassifiers <- counter - 1
  }
  if (ddalpha$methodAggregation == "sequent"){
    for (i in 1:ddalpha$numPatterns){
      anotherClass <- NULL
      for (j in 1:ddalpha$numPatterns){
        if (j != i){
          anotherClass <- rbind(anotherClass, ddalpha$patterns[[j]]$depths)
        }
      }
      polynomial <- .polynomial_learn_C(ddalpha$maxDegree, rbind(ddalpha$patterns[[i]]$depths, anotherClass), 
                                        ddalpha$patterns[[i]]$cardinality, nrow(anotherClass), ddalpha$numChunks, ddalpha$seed)

      # Adding the classifier to the list of classifiers
      ddalpha$classifiers[[i]] <- 
        list(index          = counter, 
             index0         = i, 
             index1         = -1, 
             polynomial     = polynomial$coefficients,
             degree         = polynomial$degree,
             axis           = polynomial$axis)
    }
    ddalpha$numClassifiers <- ddalpha$numPatterns
  }
  
  return (ddalpha)
}

################################################################################
# Functions for intermediate calculations are presented below
################################################################################

nlm_optimize_r <- function(r_minCandidate, r_points, numClass1, numClass2){
  return ( nlm(CGetEmpiricalRiskSmoothed, r_minCandidate, r_points, numClass1, numClass2)$estimate )
}

.polynomial_learn_C <- function(maxDegree, data, numClass1, numClass2, numChunks, seed){
  points <- as.vector(t(data))
  numPoints <- numClass1 + numClass2
  dimension <- ncol(data)
  cardinalities <- c(numClass1, numClass2)
  upToPower <- maxDegree
  minFeatures <- 2
  maxExtDimension <- (factorial(dimension + maxDegree) / (factorial(dimension)*factorial(maxDegree))) - 1;
  

  res <- .C("PolynomialLearnCV", 
          as.double(points), 
          as.integer(numPoints), 
          as.integer(dimension), 
          as.integer(cardinalities),  
          as.integer(upToPower), 
          as.integer(numChunks), 
          as.integer(seed),
          degree = integer(1),
          axis = integer(1),
          polynomial=double(upToPower))
  
  degree <- res$degree
  axis <- res$axis
  polynomial <- res$polynomial[1:degree]
  
  return(list(coefficients = polynomial, axis = axis, degree = degree))
}

GetNumsErrors <- function(polynomial, depths, numClass1, numClass2){
  # Calculates the number of classification error for two classes on the
  # basis of given depths
  # 
  # Args:
  #   polynomial: Polynomial as a vector of coefficients starting with the 
  #               first degree (a0 = 0 always)
  #   depths:     nx2 matrix of depths, where each column contains the depths
  #               against the corresponding class
  #   numClass1:  Number of points belonging to the first class
  #   numClass2:  Number of points belonging to the second class
  # Returns:
  #   Vector containing number of errors of the points from the firts and 
  #   the second class
  degree <- length(polynomial)
  numErrors1 <- 0
  if(numClass1 != 0){
    for(i in 1:numClass1){      
      val <- depths[i,1]
      res <- 0
      for(j in 1:degree){res <- res + polynomial[j]*val^j}
      if(depths[i,2] > res){
        numErrors1 <- numErrors1 + 1
      }
    }
  }
  numErrors2 <- 0
  if(numClass2 != 0){
    for(i in (numClass1 + 1):(numClass1 + numClass2)){
      val <- depths[i,1]
      res <- 0
      for(j in 1:degree){res <- res + polynomial[j]*val^j}
      if(depths[i,2] < res){
        numErrors2 <- numErrors2 + 1
      }
    }
  }
  return(c(numErrors1, numErrors2))
}

GetEmpiricalRiskSmoothed <- function(polynomial, depths, numClass1, numClass2){
  res = (colSums(sapply(depths[,1], '^', (1:length(polynomial)))*polynomial) - depths[,2])*c(rep(-1, numClass1), rep(1, numClass2))
  risk = sum(1/(1 + exp(-100*(res))))
  return (risk/(numClass1 + numClass2))
}

GetEmpiricalRisk <- function(polynomial, depths, numClass1, numClass2){
  # Calculates the empirical risk for two classes on the basis of given depths
  # 
  # Args:
  #   polynomial: Polynomial as a vector of coefficients starting with the 
  #               first degree (a0 = 0 always)
  #   depths:     nx2 matrix of depths, where each column contains the depths
  #               against the corresponding class
  #   numClass1:  Number of points belonging to the first class
  #   numClass2:  Number of points belonging to the second class
  # Returns:
  #   Empirical risk
  risk1 <- 0
  degree <- length(polynomial)
  for(i in 1:numClass1){
    val <- depths[i,1]
    res <- 0
    for(j in 1:degree){res <- res + polynomial[j]*val^j}
    if(depths[i,2] > res){
      risk1 <- risk1 + 1
    }
  }
  risk2 <- 0
  for(i in (numClass1 + 1):(numClass1 + numClass2)){
    val <- depths[i,1]
    res <- 0
    for(j in 1:degree){res <- res + polynomial[j]*val^j}
    if(depths[i,2] < res){
      risk2 <- risk2 + 1
    }
  }
  risk <- (risk1 + risk2)/(numClass1 + numClass2)
  return(risk)
}
