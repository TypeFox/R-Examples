# Similarity Measures Version 1.4, 2015-02-06
# Created by Kevin Toohey, 2014, License: GPL-3

# This file contains functions to run and assist four different similarity 
# measures. The similarity measures included are: longest common 
# subsequence (LCSS), Frechet distance, edit distance and dynamic time 
# warping (DTW). Each of these similarity measures can be calculated from 
# two n-dimensional trajectories, both in matrix form.

# Please see the function and package help files for more information.

LCSSRatio <- function(traj1, traj2, pointSpacing=-1, pointDistance=20,
                      errorMarg=2, returnTrans=FALSE) {
  # Function to calculate the ratio of the longest common sub sequence to
  # the shortest trajectory.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   pointSpacing: An integer value of the maximum index difference between 
  #                 trajectory1 and trajectory2 allowed in the calculation. 
  #                 A negative value sets the point spacing to unlimited.
  #   pointDistance: A floating point number representing the maximum 
  #                  distance in each dimension allowed for points to be 
  #                  considered equivalent.
  #   errorMarg: A floating point error margin used to scale the accuracy 
  #              and speed of the calculation.
  #   returnTrans: A boolean value to allow the best translation found to 
  #                be returned as well as the LCSS value if set to true.
  #
  # Returns:
  #   A floating point value is returned. This represents the maximum LCSS 
  #   ratio obtained using the variables provided. If returnTrans is set to 
  #   TRUE, then the LCSS ratio and the translations are returned as a 
  #   vector. The first value of this vector is the LCSS ratio and the 
  #   translations follow directly afterwards. If a problem occurs, then a 
  #   string containing information about the problem is returned.

  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  # Calculating the number of points in each trajectory and the dimensions.
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # If a trajectory has no points then the ratio is 0.
  if (length1 == 0 | length2 == 0) {
    warning("At least one trajectory contains 0 points.")
    return(0.0)
  }
  # Calculating the ratio based on the shortest trajectory.
  length <- min(length1, length2)
  return(LCSS(traj1, traj2, pointSpacing, pointDistance,
              errorMarg, returnTrans) * 1.0 / length)
}

LCSS <- function(traj1, traj2, pointSpacing=-1, pointDistance=20,
                 errorMarg=2, returnTrans=FALSE) {
  # Function to calculate the longest common subsequence for two
  # given trajectories.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   pointSpacing: An integer value of the maximum index difference between 
  #                 trajectory1 and trajectory2 allowed in the calculation. 
  #                 A negative value sets the point spacing to unlimited.
  #   pointDistance: A floating point number representing the maximum 
  #                  distance in each dimension allowed for points to be 
  #                  considered equivalent.
  #   errorMarg: A floating point error margin used to scale the accuracy 
  #              and speed of the calculation.
  #   returnTrans: A boolean value to allow the best translation found to 
  #                be returned as well as the LCSS value if set to true.
  #
  # Returns:
  #   An integer value is returned. This represents the maximum LCSS 
  #   value obtained using the variables provided. If returnTrans is set 
  #   to TRUE, then the LCSS value and the translations are returned as a 
  #   vector. The first value of this vector is the LCSS value and the 
  #   translations follow directly afterwards. If a problem occurs, then a 
  #   string containing information about the problem is returned.

  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  # Calculating the number of points in each trajectory and the dimensions.
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # If a trajectory has no points then there are 0 similar points.
  if (length1 == 0 | length2 == 0) {
    warning("At least one trajectory contains 0 points.")
    return(0)
  }
  # If the dimension is 0 then the points are considered equal.
  if (dimensions == 0) {
    warning("The dimension is 0.")
    return(min(length1, length2))
  }
  # Setting the default point spacing if required.
  if (pointSpacing < 0) {
    pointSpacing <- length1
    if (length1 < length2) {
      pointSpacing <- length2
    }
  }
  # Calculating the subsets of translations.
  translations <- list()
  for (d in 1:dimensions) {
    translations[[d]] <- TranslationSubset(traj1[, d], traj2[, d],
                                         pointSpacing, pointDistance)
  }
  # Storing the most optimal translations and similarity found so far.
  similarity <- LCSSCalc(traj1, traj2, pointSpacing, pointDistance)
  optimalTrans <- rep(0.0, dimensions)
  similarity <- c(similarity, optimalTrans)
  # Calculating how many translation possibilities are skipped for
  # every one that is checked using the error margin given.
  spacing <- length(translations[[1]]) / (4 * pointSpacing / errorMarg)
  if (spacing < 1) {
    spacing <- 1
  } else if (spacing > (length(translations[[1]]) / 2.0)) {
    spacing <- length(translations[[1]]) / 2.0
  }
  spacing <- as.integer(spacing)
  # Running the LCSS algorithm on each of the translations to be checked.
  similarity <- SimLoop(traj1, traj2, pointSpacing, pointDistance, spacing,
                      similarity, translations, dimensions)
  # Returning the similarity and translations if requested.
  if (returnTrans == TRUE) {
    return(similarity)
  }
  # Returning the best similarity found.
  return(similarity[1])
}

SimLoop <- function(traj1, traj2, pointSpacing, pointDistance, spacing,
                    similarity, translations, dimensions, dimLeft=dimensions,
                    currentTrans=rep(0.0, dimensions)) {
  # Function to loop over and test the trajectories using the different
  # translations in each dimension (this should not be called directly).
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   pointSpacing: An integer value of the maximum index difference between 
  #                 trajectory1 and trajectory2 allowed in the calculation. 
  #                 A negative value sets the point spacing to unlimited.
  #   pointDistance: A floating point number representing the maximum 
  #                  distance in each dimension allowed for points to be 
  #                  considered equivalent.
  #   spacing: The integer spacing between each translation that will 
  #            be tested.
  #   similarity: A vector containing the current best similarity and 
  #               translations calculated.
  #   translations: A list of vectors containing the translations in 
  #                 each dimension.
  #   dimensions: An integer representing the number of dimensions being 
  #               used for the calculation.
  #   dimLeft: An integer number of dimensions which have not been 
  #            looped over yet.
  #   currentTrans: A vector containing the current translation being tested.
  #
  # Returns:
  #   Returns the current best LCSS value and the translations that created 
  #   this as a vector.

  # Testing each translation in this dimension.
  thisDim <- 1 + dimensions - dimLeft
  prevTrans <- NULL
  for (i in seq(spacing, length(translations[[thisDim]]), by=spacing)) {
    # The newest translation.
    currentTrans[thisDim] <- translations[[thisDim]][round(i)]
    # Skipping translations which have already been checked.
    if (!(isTRUE(all.equal(currentTrans[thisDim], prevTrans,
                           tolerance=.Machine$double.eps * 1000)))) {
      if (dimLeft > 1) {
        similarity <- SimLoop(traj1, traj2, pointSpacing, pointDistance,
                              spacing, similarity, translations,
                              dimensions, (dimLeft - 1), currentTrans)
      } else {
        # Running the LCSS algorithm on each of the translations to
        # be checked.
        newValue <- LCSSCalc(traj1, traj2, pointSpacing, pointDistance,
                             currentTrans)
        # Keeping the new similarity and translations if they are better than
        # the previous best.
        if (newValue > similarity[1]) {
          similarity[1] <- newValue
          for (d in 1:dimensions) {
            similarity[(d + 1)] <- currentTrans[d]
          }
        }
      }
      prevTrans <- currentTrans[thisDim]
    }
  }
  # Returning the vector containing the current best similarity
  # with translations.
  return(similarity)
}

LCSSCalc <- function(traj1, traj2, pointSpacing=-1, pointDistance=20,
                     trans=rep(0.0, (dim(traj1)[2]))) {
  # Function to calculate the LCSS of two trajectories using a set translation.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   pointSpacing: An integer value of the maximum index difference between 
  #                 trajectory1 and trajectory2 allowed in the calculation. 
  #                 A negative value sets the point spacing to unlimited.
  #   pointDistance: A floating point number representing the maximum 
  #                  distance in each dimension allowed for points to be 
  #                  considered equivalent.
  #   trans: A vector containing translations in each dimension to be applied 
  #          to trajectory2 in this calculation.
  #
  # Returns:
  #   An integer value is returned. This represents the maximum LCSS value 
  #   obtained using the variables provided. If a problem occurs, then a 
  #   string containing information about the problem is returned.
  
  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  # Calculating the number of points in each trajectory and the dimensions.
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # If a trajectory has no points then there are 0 similar points.
  if (length1 == 0 | length2 == 0) {
    warning("At least one trajectory contains 0 points.")
    return(0)
  }
  # If the dimension is 0 then the points are considered equal.
  if (dimensions == 0) {
    warning("The dimension is 0.")
    return(min(length1, length2))
  }
  # Setting the default point spacing if required.
  if (pointSpacing < 0) {
    pointSpacing <- length1
    if (length1 < length2) {
      pointSpacing <- length2
    }
  }
  # Rounding the point spacing if necessary.
  pointSpacing <- round(pointSpacing)

  distMatrix <- matrix(0, nrow=length1, ncol=length2)
  similarity <- 0
  for (row in 1:length1) {
    # Calculating the relevant columns for each row.
    minCol <- 1
    maxCol <- length2
    if (row > pointSpacing + 1) {
      minCol <- row - pointSpacing
    }
    if (row < length2 - pointSpacing) {
      maxCol <- row + pointSpacing
    }
    if (minCol <= maxCol) {
      for (col in minCol:maxCol) {
        newValue <- 0
        finalValue <- 0
        # Calculating the new LCSS value for the current two points.
        # Checking the diagonal.
        if (row != 1 & col != 1) {
          newValue <- distMatrix[row - 1, col - 1]
          finalValue <- newValue
        }
        # Checking below.
        if (row != 1) {
          below <- distMatrix[row - 1, col]
          if (below > finalValue) {
            finalValue <- below
          }
        }
        # Checking to the left.
        if (col != 1) {
          before <- distMatrix[row, col - 1]
          if (before > finalValue) {
            finalValue <- before
          }
        }
        # Checking if the current points can increment the LCSS.
        if (finalValue < newValue + 1) {
          checkPoint <- DistanceCheck(traj1[row, ], (traj2[col, ] + trans),
                                      pointDistance, dimensions)
          if (checkPoint) {
            newValue <- newValue + 1
            finalValue <- newValue
          }
        }
        # Updating the distance matrix.
        distMatrix[row, col] <- finalValue
        # Updating the similarity if a new maximum has been found.
        if (finalValue > similarity) {
          similarity <- finalValue
        }
      }
    }
  }
  # Returning the largest similarity.
  return(similarity)
}

LCSSRatioCalc <- function(traj1, traj2, pointSpacing=-1, pointDistance=20,
                     trans=rep(0.0, (dim(traj1)[2]))) {
  # A function to calculate the LCSS ratio between two trajectories using
  # a set translation.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   pointSpacing: An integer value of the maximum index difference between 
  #                 trajectory1 and trajectory2 allowed in the calculation. 
  #                 A negative value sets the point spacing to unlimited.
  #   pointDistance: A floating point number representing the maximum 
  #                  distance in each dimension allowed for points to be 
  #                  considered equivalent.
  #   trans: A vector containing translations in each dimension to be applied 
  #          to trajectory2 in this calculation.
  #
  # Returns:
  #   A floating point value is returned. This represents the maximum LCSS 
  #   ratio obtained using the variables provided. If a problem occurs, then 
  #   a string containing information about the problem is returned.
  
  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  # Calculating the number of points in each trajectory and the dimensions.
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # If a trajectory has no points then the ratio is 0.
  if (length1 == 0 | length2 == 0) {
    warning("At least one trajectory contains 0 points.")
    return(0.0)
  }
  # Calculating the ratio based on the shortest trajectory.
  length <- min(length1, length2)
  return(LCSSCalc(traj1, traj2, pointSpacing, pointDistance,
                  trans) * 1.0 / length)
}

LCSSTranslation <- function(traj1, traj2, pointSpacing=-1, pointDistance=20,
                            errorMarg=2) {
  # A function for returning the best translation calculated using the
  # LCSS method.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   pointSpacing: An integer value of the maximum index difference between 
  #                 trajectory1 and trajectory2 allowed in the calculation. 
  #                 A negative value sets the point spacing to unlimited.
  #   pointDistance: A floating point number representing the maximum 
  #                  distance in each dimension allowed for points to be 
  #                  considered equivalent.
  #   errorMarg: A floating point error margin used to scale the accuracy 
  #              and speed of the calculation.
  #
  # Returns:
  #   A vector of length n is returned containing the translation in each 
  #   dimension. If there is a problem this returns a string of information 
  #   instead.
  
  # Running the LCSS method.
  values <- LCSS(traj1, traj2, pointSpacing, pointDistance,
                 errorMarg, TRUE)
  # Returning the translations as a vector.
  return(values[2:(length(values))])
}

TranslationSubset <- function(traj1, traj2, pointSpacing, pointDistance) {
  # A function for calculating the subsets of translations to be tested
  # using the LCSS methods (this should not be called directly).
  #
  # Args:
  #   traj1: A vector containing one dimension of trajectory1.
  #   traj2: A vector containing one dimension of trajectory2.
  #   pointSpacing: An integer value of the maximum index difference between 
  #                 trajectory1 and trajectory2 allowed in the calculation.
  #   pointDistance: A floating point number representing the maximum 
  #                  distance in each dimension allowed for points to be 
  #                  considered equivalent.
  #
  # Returns:
  #   A vector of floating point numbers is returned containing the 
  #   translations calculated. This vector is sorted in ascending order.
  
  # Calculating the lengths of the trajectories.
  length1 <- length(traj1)
  length2 <- length(traj2)
  translations <- vector()
  for (row in 1:length1) {
    # Calculating the relevant columns for each row.
    minCol <- 1
    maxCol <- length2
    if (row > pointSpacing + 1) {
      minCol <- row - pointSpacing
    }
    if (row < length2 - pointSpacing) {
      maxCol <- row + pointSpacing
    }
    if (minCol <= maxCol) {
      for (col in minCol:maxCol) {
        # Adding the new translations calculated from the distance boundaries.
        translations <- c(translations, (traj1[row] - traj2[col] + 
                                           pointDistance))
        translations <- c(translations, (traj1[row] - traj2[col] - 
                                           pointDistance))
      }
    }
  }
  # Returning the translations as a sorted vector.
  return(sort(translations))
}

Frechet <- function(traj1, traj2, testLeash=-1.0) {
  # A function to calculate the Frechet distance between two trajectories.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   testLeash: A numeric leash value, which if positive, checks whether 
  #              the leash can be used between the two trajectories. If this 
  #              value is negative, then it is not used and the standard 
  #              calculation is performed.
  #
  # Returns:
  #   A floating point value representing the Frechet distance is returned. 
  #   If a test leash is given, then a boolean value is returned as true if 
  #   the leash was successful and false if not. If a problem occurs, then a 
  #   string containing information about the problem is returned.
  
  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  # Calculating lengths and dimensions.
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # If a length is 0 then a trajectory has no points.
  if (length1 == 0 | length2 == 0) {
    warning("At least one trajectory contains 0 points.")
    return("At least one trajectory contains 0 points.")
  }
  # If the dimension is 0 then the points are considered equal.
  if (dimensions == 0) {
    warning("The dimension is 0.")
    return(0.0)
  }
  # Calculation if either trajectory has only one point.
  if (length1 == 1 | length2 == 1) {
    leash <- SinglePointCalc(traj1, traj2)
    if (testLeash >= 0) {
      if (testLeash >= leash) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      return(leash)
    }
  }
  dist1 <- rep(0, length1 - 1)
  dist2 <- rep(0, length2 - 1)
  
  # Computing the distances between each point and the next for
  # both trajectories.
  for (point in 1:(length1 - 1)) {
    dist <- DistanceSq(traj1[point + 1, ], traj1[point, ], dimensions)
    dist1[point] <- sqrt(dist)
  }
  for (point in 1:(length2 - 1)) {
    dist <- DistanceSq(traj2[point + 1, ], traj2[point, ], dimensions)
    dist2[point] <- sqrt(dist)
  }
  

  # Checking the distances at the end points for the minimum leash
  # possibility.
  minLeashSq <- DistanceSq(traj1[1, ], traj2[1, ], dimensions)
  endDistSq <- DistanceSq(traj1[length1, ], traj2[length2, ], dimensions)
  if (minLeashSq < endDistSq) {
    minLeashSq <- endDistSq
  }
  leashList <- sqrt(minLeashSq)
  
  # Computing the squares of all the distances between the two trajectories.
  distSq12 <- matrix(0, nrow=length1, ncol=length2)
  for (point1 in 1:length1) {
    for (point2 in 1:length2) {
      dist <- DistanceSq(traj1[point1, ], traj2[point2, ], dimensions)
      distSq12[point1, point2] <- dist
      # Adding in these leash possibilties because they are critical points.
      if (dist > minLeashSq) {
        leashList <- c(leashList, sqrt(dist))
      }
    }
  }
  # If test leash is being used.
  if (testLeash >= 0) {
    return(FrechetCheck(traj1, traj2, testLeash, dist1, dist2, distSq12))
  }
  # Adding critical point leash possibilities to the leash list.
  for (point1 in 1:(length1 - 1)) {
    # Creating a unit vector in the direction of the next point from point1.
    unitV1 <- rep(0, times=dimensions)
    if (dist1[point1] != 0) {
      unitV1 <- 1.0 * (traj1[point1 + 1, ] - traj1[point1, ]) / 
        (dist1[point1])
    }
    for (point2 in 1:length2) {
      # Creating a vector from point1 to point 2.
      vect12 <- traj2[point2, ] - traj1[point1, ]
      # Dot product finds how far from point1 the closest point on the
      # line is.
      pointDist <- Dot(unitV1, vect12, dimensions)
      # Squaring for easy calculations
      pointDistSq <- pointDist * pointDist
      # The square of the distance between the line segment and the point.
      shortDist <- distSq12[point1, point2] - pointDistSq
      leashSq <- 0.0
      if (pointDist < 0) {
        # This case is accounted for earlier.
      } else if (pointDist > dist1[point1]) {
        leashSq <- distSq12[(point1 + 1), point2]
      } else {
        leashSq <- shortDist
      }
      # Adding the leash possibility to the list.
      if (leashSq > minLeashSq) {
        leashList <- c(leashList, sqrt(leashSq))
      }
      
    }
  }
  
  for (point2 in 1:(length2 - 1)) {
    # Creating a unit vector in the direction of the next point from point2.
    unitV1 <- rep(0, times=dimensions)
    if (dist2[point2] != 0) {
      unitV1 <- 1.0 * (traj2[point2 + 1, ] - traj2[point2, ]) / 
        (dist2[point2])
    }
    for (point1 in 1:length1) {
      # Creating a vector from point2 to point 1.
      vect12 <- traj1[point1, ] - traj2[point2, ]
      # Dot product finds how far from point2 the closest point on the
      # line is.
      pointDist <- Dot(unitV1, vect12, dimensions)
      # Squaring for easy calculations
      pointDistSq <- pointDist * pointDist
      # The square of the distance between the line segment and the point.
      shortDist <- distSq12[point1, point2] - pointDistSq
      leashSq <- 0.0
      if (pointDist < 0) {
        # This case is accounted for earlier.
      } else if (pointDist > dist2[point2]) {
        leashSq <- distSq12[point1, (point2 + 1)]
      } else {
        leashSq <- shortDist
      }
      # Adding the leash possibility to the list.
      if (leashSq > minLeashSq) {
        leashList <- c(leashList, sqrt(leashSq))
      }
      
    }
  }
  
  # Calculating the critical points where new passages may open.
  if (length1 > 3) {
    for (point2 in 1:(length2 - 1)) {
      # Creating a unit vector in the direction of the next point from point1.
      unitV2 <- rep(0, times=dimensions)
      if (dist2[point2] != 0) {
        unitV2 <- 1.0 * (traj2[point2 + 1, ] - traj2[point2, ]) / 
          (dist2[point2])
      }
      for (point1 in 2:(length1 - 2)) {
        # Creating a vector from point 2 to point 1.
        vect21 <- traj1[point1, ] - traj2[point2, ]
        # Dot product finds how far from point2 the closest point on the
        # line is.
        pointDist <- Dot(unitV2, vect21, dimensions)
        if (pointDist > 0) {
          # Squaring for easy calculations
          pointDistSq <- pointDist * pointDist
          # The square of the distance between the line segment and the point.
          shortDist <- distSq12[point1, point2] - pointDistSq
          # The second point where the passage opens up.
          for (newPoint in (point1 + 1):(length1 - 1)) {
            # Creating a vector from point 2 to newPoint.
            vect2new <- traj1[newPoint, ] - traj2[point2, ]
            # Dot product finds how far from point2 the closest point on the
            # line is.
            newPointDist <- Dot(unitV2, vect2new, dimensions)
            if (newPointDist < pointDist) {
              newPointDistSq <- newPointDist * newPointDist
              newShortDist <- distSq12[newPoint, point2] - newPointDistSq
              # The distance between the two closest points on the line.
              pointDiff <- pointDist - newPointDist
              # Finding the point where the passage opens.
              equivPoint <- (pointDiff * pointDiff + shortDist -
                               newShortDist) / (2.0 * pointDiff)
              if (equivPoint > 0 & equivPoint < dist2[point2]) {
                # Adding this leash to the list of possible leashes.
                leashSq <- newShortDist + equivPoint * equivPoint
                if (leashSq > minLeashSq) {
                  leashList <- c(leashList, sqrt(leashSq))
                }
              }
            }
          }
        }
      }
    }
  }
  if (length2 > 3) {
    for (point1 in 1:(length1 - 1)) {
      # Creating a unit vector in the direction of the next point from point2.
      unitV1 <- rep(0, times=dimensions)
      if (dist1[point1] != 0) {
        unitV1 <- 1.0 * (traj1[point1 + 1, ] - traj1[point1, ]) / 
          (dist1[point1])
      }
      for (point2 in 2:(length2 - 2)) {
        # Creating a vector from point 1 to point 2.
        vect12 <- traj2[point2, ] - traj1[point1, ]
        # Dot product finds how far from point1 the closest point on the
        # line is.
        pointDist <- Dot(unitV1, vect12, dimensions)
        if (pointDist > 0) {
          # Squaring for easy calculations
          pointDistSq <- pointDist * pointDist
          # The square of the distance between the line segment and the point.
          shortDist <- distSq12[point1, point2] - pointDistSq
          # The second point where the passage opens up.
          for (newPoint in (point2 + 1):(length2 - 1)) {
            # Creating a vector from point 1 to newPoint.
            vect1new <- traj2[newPoint, ] - traj1[point1, ]
            # Dot product finds how far from point1 the closest point on the
            # line is.
            newPointDist <- Dot(unitV1, vect1new, dimensions)
            if (newPointDist < pointDist) {
              newPointDistSq <- newPointDist * newPointDist
              newShortDist <- distSq12[point1, newPoint] - newPointDistSq
              # The distance between the two closest points on the line.
              pointDiff <- pointDist - newPointDist
              # Finding the point where the passage opens.
              equivPoint <- (pointDiff * pointDiff + shortDist -
                               newShortDist) / (2.0 * pointDiff)
              if (equivPoint > 0 & equivPoint < dist1[point1]) {
                # Adding this leash to the list of possible leashes.
                leashSq <- newShortDist + equivPoint * equivPoint
                if (leashSq > minLeashSq) {
                  leashList <- c(leashList, sqrt(leashSq))
                }
              }
            }
          }
        }
      }
    }
  }
  # Sorting the leash list and removing multiples.
  leashList <- sort(leashList)
  uniqLeash <- leashList[1]
  lastLeash <- uniqLeash
  for (item in leashList) {
    if (lastLeash != item) {
      lastLeash <- item
      uniqLeash <- c(uniqLeash, item)
    }
  }
  # Setting up binary search for the list.
  startSearch <- 1
  endSearch <- length(uniqLeash)
  # Making sure that the final leash is large enough.
  if (FrechetCheck(traj1, traj2, uniqLeash[endSearch], dist1, dist2,
                   distSq12)) {
    # Performing binary search on the list of leashes.
    while (startSearch < endSearch) {
      current <- (as.integer((endSearch - startSearch) / 2) + startSearch)
      if (FrechetCheck(traj1, traj2, uniqLeash[current], dist1, dist2,
                       distSq12)) {
        endSearch <- current
      } else {
        startSearch <- current + 1
      }
    }
    # Returning the shortest leash for the trajectories.
    return(uniqLeash[endSearch])
  } else {
    warning("The Frechet distance was unable to be found")
    return(-1)
  }
}

DistanceSq <- function(point1, point2, dimensions=length(point1)) {
  # A function to calculate the square of the distance between two points
  # with the same dimensions.
  #
  # Args:
  #   point1: An n dimensional vector representing point1.
  #   point2: An n dimensional vector representing point2.
  #   dimensions: An integer representing the number of dimensions being 
  #               used for the distance square calculation. This defaults 
  #               to the length of the first vector.
  #
  # Returns:
  #   A floating point value is returned, representing the square of the
  #   distance between the two points.
  
  dist <- 0.0
  # Adding on the square from each dimension.
  for (d in 1:dimensions) {
    dist <- dist + (point1[d] - point2[d])^2
  }
  return(dist)
}

DistanceCheck <- function(point1, point2, dist, dimensions=length(point1)) {
  # A function to check whether two points lie within some distance
  # in every dimension.
  #
  # Args:
  #   point1: An n dimensional vector representing point1.
  #   point2: An n dimensional vector representing point2.
  #   dist: A floating point number representing the maximum distance in 
  #         each dimension allowed for points to be considered equivalent.
  #   dimensions: An integer representing the number of dimensions being 
  #               checked. This defaults to the length of the first vector.
  #
  # Returns:
  #   A boolean value is returned. The value is true if the points are 
  #   within the distance in every dimension and false if not.
  
  # Initialising to TRUE which is returned if the check is successful.
  check <- TRUE
  for (d in 1:dimensions) {
    if (check) {
      newDist <- abs(point1[d] - point2[d])
      # If a the points are not within the distance in a dimension
      # then FALSE is returned.
      if (!(isTRUE(all.equal(newDist, dist,
                             tolerance=.Machine$double.eps * 1000)))
          & newDist > dist) {
        check <- FALSE
      }
    }
  }
  # Returning the boolean value.
  return(check)
}

SinglePointCalc <- function(traj1, traj2) {
  # A function to calculate the Frechet distance when one trajectory
  # contains only a single point (this should not be called directly).
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #
  # Returns:
  #   A floating point value representing the Frechet distance is returned.
  
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  leashSq <- -1.0
  # Calculating each leash possibility.
  if (length1 == 1) {
    for (point2 in 1:length2) {
      newLeashSq <- DistanceSq(traj1[1, ], traj2[point2, ], dimensions)
      # Keeping the new leash if it is longer than the previous.
      if (newLeashSq > leashSq) {
        leashSq <- newLeashSq
      }
    }
  } else if (length2 == 1) {
    for (point1 in 1:length1) {
      newLeashSq <- DistanceSq(traj1[point1, ], traj2[1, ], dimensions)
      # Keeping the new leash if it is longer than the previous.
      if (newLeashSq > leashSq) {
        leashSq <- newLeashSq
      }
    }
  }
  # Returning the leash.
  if (leashSq >= 0) {
    return(sqrt(leashSq))
  } else {
    warning("Error in single point trajectory calculation.")
    return("Error in single point trajectory calculation.")
  }
}

FrechetCheck <- function(traj1, traj2, leash, dist1, dist2, distSq12) {
  # A function to check whether a leash will work between two trajectories 
  # using the Frechet method (this should not be called directly).
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   leash: A numeric leash value to be checked by the function.
  #   dist1: A vector containing the distance between each successive two 
  #          points in trajectory1.
  #   dist2: A vector containing the distance between each successive two 
  #          points in trajectory2.
  #   distSq12: A matrix containing the distance between each pair of two 
  #             points where 1 point lies in trajectory1 and the other 
  #             in trajectory2.
  #
  # Returns:
  #   A boolean value is returned. A value of true is returned if the leash 
  #   is successful and false if not.
  
  # Setting up the dimensions, lengths and required arrays.
  leashSq <- leash * leash
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  left <- array(-1, dim=c(length1, (length2 - 1), 2))
  bottom <- array(-1, dim=c((length1 - 1), length2, 2))
  newLeft <- array(-1, dim=c(length1, (length2 - 1), 2))
  newBot <- array(-1, dim=c((length1 - 1), length2, 2)) 
  # Checking the leash works at the endpoints.
  if (leashSq < distSq12[1, 1] | leashSq < distSq12[length1, length2]) {
    return(FALSE)
  }
  # Calculating the freespace of the first trajectory with respect
  # to the second.
  for (point1 in 1:(length1 - 1)) {
    # Creating a unit vector in the direction of the next point from point1.
    unitV1 <- rep(0, times=dimensions)
    if (dist1[point1] != 0) {
      unitV1 <- 1.0 * (traj1[point1 + 1, ] - traj1[point1, ]) / 
        (dist1[point1])
    }
    for (point2 in 1:length2) {
      # Creating a vector from point1 to point 2.
      vect12 <- traj2[point2, ] - traj1[point1, ]
      # Dot product finds how far from point1 the closest point on the
      # line is.
      pointDist <- Dot(unitV1, vect12, dimensions)
      # Squaring for easy calculations
      pointDistSq <- pointDist * pointDist
      # The square of the distance between the line segment and the point.
      shortDist <- distSq12[point1, point2] - pointDistSq
      # If some part of the current line can be used by the leash.
      if (shortDist <= leashSq) {
        # Calculating the envelope along the line.
        envSize <- sqrt(leashSq - shortDist)
        envLow <- pointDist - envSize
        envHigh <- pointDist + envSize
        # If the whole line is within the envelope.
        if (envHigh >= dist1[point1] & envLow <= 0) {
          bottom[point1, point2, 1] <- 0.0
          bottom[point1, point2, 2] <- 1.0
        } else if (envHigh >= 0 & envLow <= 0) {
          # If the start of the line is within the envelope.
          bottom[point1, point2, 1] <- 0.0
          bottom[point1, point2, 2] <- 1.0 * envHigh / dist1[point1]
        } else if (envHigh >= dist1[point1] & envLow <= dist1[point1]) {
          # If the end of the line is within the envelope.
          bottom[point1, point2, 1] <- 1.0 * envLow / dist1[point1]
          bottom[point1, point2, 2] <- 1.0
        } else if (envHigh >= 0 & envLow <= dist1[point1]) {
          # If the envelope is completely within the line.
          bottom[point1, point2, 1] <- 1.0 * envLow / dist1[point1]
          bottom[point1, point2, 2] <- 1.0 * envHigh / dist1[point1]
        }
      }
    }
  }
  # Calculating the freespace of the second trajectory with respect
  # to the first.
  for (point2 in 1:(length2 - 1)) {
    # Creating a unit vector in the direction of the next point from point2.
    unitV1 <- rep(0, times=dimensions)
    if (dist2[point2] != 0) {
      unitV1 <- 1.0 * (traj2[point2 + 1, ] - traj2[point2, ]) / 
        (dist2[point2])
    }
    for (point1 in 1:length1) {
      # Creating a vector from point2 to point 1.
      vect12 <- traj1[point1, ] - traj2[point2, ]
      # Dot product finds how far from point2 the closest point on the
      # line is.
      pointDist <- Dot(unitV1, vect12, dimensions)
      # Squaring for easy calculations
      pointDistSq <- pointDist * pointDist
      # The square of the distance between the line segment and the point.
      shortDist <- distSq12[point1, point2] - pointDistSq
      # If some part of the current line can be used by the leash.
      if (shortDist <= leashSq) {
        # Calculating the envelope along the line.
        envSize <- sqrt(leashSq - shortDist)
        envLow <- pointDist - envSize
        envHigh <- pointDist + envSize
        # If the whole line is within the envelope.
        if (envHigh >= dist2[point2] & envLow <= 0) {
          left[point1, point2, 1] <- 0.0
          left[point1, point2, 2] <- 1.0
        } else if (envHigh >= 0 & envLow <= 0) {
          # If the start of the line is within the envelope.
          left[point1, point2, 1] <- 0.0
          left[point1, point2, 2] <- 1.0 * envHigh / dist2[point2]
        } else if (envHigh >= dist2[point2] & envLow <= dist2[point2]) {
          # If the end of the line is within the envelope.
          left[point1, point2, 1] <- 1.0 * envLow / dist2[point2]
          left[point1, point2, 2] <- 1.0
        } else if (envHigh >= 0 & envLow <= dist2[point2]) {
          # If the envelope is completely within the line.
          left[point1, point2, 1] <- 1.0 * envLow / dist2[point2]
          left[point1, point2, 2] <- 1.0 * envHigh / dist2[point2]
        }
      }
    }
  }
  # Setting up the new arrays to find the monotone freespace.
  newLeft[1, 1, 1] <- left[1, 1, 1]
  newLeft[1, 1, 2] <- left[1, 1, 2]
  newBot[1, 1, 1] <- bottom[1, 1, 1]
  newBot[1, 1, 2] <- bottom[1, 1, 2]
  # Setting the first line of the new left array.
  if (length2 > 2) {
    for (point2 in 2:(length2 - 1)) {
      if (newLeft[1, (point2 - 1), 2] == 1) {
        newLeft[1, point2, 1] <- left[1, point2, 1]
        newLeft[1, point2, 2] <- left[1, point2, 2]
      }
    }
  }
  # Setting the first line of the new bottom array.
  if (length1 > 2) {
    for (point1 in 2:(length1 - 1)) {
      if (newBot[(point1 - 1), 1, 2] == 1) {
        newBot[point1, 1, 1] <- bottom[point1, 1, 1]
        newBot[point1, 1, 2] <- bottom[point1, 1, 2]
      }
    }
  }
  # Calculating the monotone freespace
  for (point1 in 1:length1) {
    for (point2 in 1:length2) {
      if (point1 != length1 & point2 != 1) {
        # If the area is allowable from the freespace.
        if (bottom[point1, point2, 1] > -0.1) {
          # If the area can be entered from the left.
          if (newLeft[point1, (point2 - 1), 1] > -0.1) {
            # Setting the monotone freespace for these points.
            newBot[point1, point2, 1] <- bottom[point1, point2, 1]
            newBot[point1, point2, 2] <- bottom[point1, point2, 2]
          } else if (newBot[point1, (point2 - 1), 1] > -0.1) {
            # Otherwise using the new bottom array.
            # Setting the monotone freespace according to what is reachable.
            if (newBot[point1, (point2 - 1), 1] <= 
                  bottom[point1, point2, 1]) {
              newBot[point1, point2, 1] <- bottom[point1, point2, 1]
              newBot[point1, point2, 2] <- bottom[point1, point2, 2]
            } else if (newBot[point1, (point2 - 1), 1] <= 
                         bottom[point1, point2, 2]) {
              newBot[point1, point2, 1] <- newBot[point1, (point2 - 1), 1]
              newBot[point1, point2, 2] <- bottom[point1, point2, 2]
            }
          }
        }
      }
      if (point2 != length2 & point1 != 1) {
        # If the area is allowable from the freespace.
        if (left[point1, point2, 1] > -0.1) {
          # If the area can be entered from below.
          if (newBot[(point1 - 1), point2, 1] > -0.1) {
            # Setting the monotone freespace for these points.
            newLeft[point1, point2, 1] <- left[point1, point2, 1]
            newLeft[point1, point2, 2] <- left[point1, point2, 2]
          } else if (newLeft[(point1 - 1), point2, 1] > -0.1) {
            # Otherwise using the new left array.
            # Setting the monotone freespace according to what is reachable.
            if (newLeft[(point1 - 1), point2, 1] <= 
                  left[point1, point2, 1]) {
              newLeft[point1, point2, 1] <- left[point1, point2, 1]
              newLeft[point1, point2, 2] <- left[point1, point2, 2]
            } else if (newLeft[(point1 - 1), point2, 1] <=
                       left[point1, point2, 2]) {
              newLeft[point1, point2, 1] <- newLeft[(point1 - 1), point2, 1]
              newLeft[point1, point2, 2] <- left[point1, point2, 2]
            }
          }
        }
      }
    }
  }
  # If the monotone freespace reaches the final point then
  # the leash is successful.
  if (newLeft[length1, (length2 - 1), 2] == 1 |
        newBot[(length1 - 1), length2, 2] == 1) {
    return(TRUE)
  } else {
  # Otherwise the leash is unsuccessful.
    return(FALSE)
  }
}

DTW <- function(traj1, traj2, pointSpacing=-1) {
  # A function for calculating the Dynamic Time Warping value between
  # two trajectories.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   pointSpacing: An integer value of the maximum index difference between 
  #                 trajectory1 and trajectory2 allowed in the calculation. 
  #                 A negative value sets the point spacing to unlimited.
  #
  # Returns:
  #   A floating point value representing the smallest warp path is returned. 
  #   If a problem occurs, then a string containing information about the 
  #   problem is returned.
  
  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # Checking lengths and dimensions.
  if (length1 == 0 | length2 == 0) {
    warning("At least one trajectory contains 0 points.")
    return("At least one trajectory contains 0 points.")
  }
  if (dimensions == 0) {
    warning("The dimension is 0.")
    return("The dimension is 0.")
  }
  # Setting the default point spacing if required.
  if (pointSpacing < 0) {
    pointSpacing <- length1
    if (length1 < length2) {
      pointSpacing <- length2
    }
  }
  # Initialising the matrices required.
  warpPaths <- matrix(-1.0, nrow=length1, ncol=length2)
  dist <- DistanceSq(traj1[1, ], traj2[1, ], dimensions)
  warpPaths[1, 1] <- sqrt(dist)
  # Setting the first row and column of the warp path matrix.
  if (length1 > 1 & pointSpacing > 0) {
    for (point in 2:(min(length1, (pointSpacing + 1)))) {
      dist <- DistanceSq(traj1[point, ], traj2[1, ], dimensions)
      warpPaths[point, 1] <- sqrt(dist) + warpPaths[(point - 1), 1]
    }
  }
  if (length2 > 1 & pointSpacing > 0) {
    for (point in 2:(min(length2, (pointSpacing + 1)))) {
      dist <- DistanceSq(traj1[1, ], traj2[point, ], dimensions)
      warpPaths[1, point] <- sqrt(dist) + warpPaths[1, (point - 1)]
    }
  }
  # Setting the rest of the warp path matrix.
  if (length1 > 1 & length2 > 1 & pointSpacing >= 0) {
    for (point1 in 2:length1) {
      for (point2 in 2: length2) {
        pointDiff <- point1 - point2
        # When the points are within point distance.
        if (abs(pointDiff) <= pointSpacing) {
          # Calculating the square of the distance between the points.
          dist <- DistanceSq(traj1[point1, ], traj2[point2, ], dimensions)
          path <- -1.0
          # When no point spacing is allowed.
          if (pointSpacing == 0) {
            path <- warpPaths[(point1 - 1), (point2 - 1)]
          } else if (pointDiff == pointSpacing) {
            # The furthest distance forward point calculation.
            path <- min(warpPaths[(point1 - 1), (point2 - 1)],
                        warpPaths[(point1 - 1), point2])
          } else if ((-pointDiff) == pointSpacing) {
            # The furthest distance backwards point calculation.
            path <- min(warpPaths[(point1 - 1), (point2 - 1)],
                        warpPaths[point1, (point2 - 1)])
          } else {
            # All of the other points.
            path <- min(warpPaths[(point1 - 1), (point2 - 1)],
                        warpPaths[point1, (point2 - 1)],
                        warpPaths[(point1 - 1), point2])
          }
          # Storing the new best path for these points.
          warpPaths[point1, point2] <- path + sqrt(dist)
        }
      }
    }
  }
  # Returning the best warp path distance.
  return(warpPaths[length1, length2])
}

EditDist <- function(traj1, traj2, pointDistance=20) {
  # A function for calculating the Edit Distance between two trajectories.
  # This takes two trajectories and the maximum distance between points.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #   pointDistance: A floating point number representing the maximum 
  #                  distance in each dimension allowed for points to be 
  #                  considered equivalent.
  #
  # Returns:
  #   An integer representing the minimum number of edits required is 
  #   returned. If a problem occurs, then a string containing information 
  #   about the problem is returned.
  
  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  # Calculating dimension and lengths.
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # When one trajectory has no points.
  if (length1 == 0) {
    warning("At least one trajectory contains 0 points.")
    return(length2)
  }
  if (length2 == 0) {
    warning("At least one trajectory contains 0 points.")
    return(length1)
  }
  # When the dimension is 0.
  if (dimensions == 0) {
    warning("The dimension is 0.")
    return("The dimension is 0.")
  }
  # Initialising the path matrix and setting up the first column and row.
  editPaths <- matrix(-1, nrow=(length1 + 1), ncol=(length2 + 1))
  for (point1 in 1:(length1 + 1)) {
    editPaths[point1, 1] <- point1 - 1
  }
  for (point2 in 2:(length2 + 1)) {
    editPaths[1, point2] <- point2 - 1
  }
  # Setting up the rest of the matrix.
  for (point1 in 2:(length1 + 1)) {
    for (point2 in 2:(length2 + 1)) {
      # Setting the diagonal increment depending on whether the
      # current points are within range.
      diag <- 1
      if (DistanceCheck(traj1[(point1 - 1), ], traj2[(point2 - 1), ],
                        pointDistance, dimensions)) {
        diag <- 0
      }
      # Setting the path to the current two points to the minimum value.
      pathValue <- min((editPaths[(point1 - 1), point2] + 1),
                       (editPaths[point1, (point2 - 1)] + 1),
                       (editPaths[(point1 - 1), (point2 - 1)] + diag))
      editPaths[point1, point2] <- pathValue
    }
  }
  # Returning the final minimum path value as the Edit Distance.
  return(editPaths[(length1 + 1), (length2 + 1)])
}

Dot <- function(vect1, vect2, dimensions=length(vect1)) {
  # A function to calculate dot products between two vectors.
  #
  # Args:
  #   vect1: An n dimensional vector representing the first vector.
  #   vect2: An n dimensional vector representing the second vector.
  #   dimensions: An integer representing the number of dimensions being 
  #               used for the dot product calculation. This defaults to 
  #               the length of the first vector.
  #
  # Returns:
  #   A floating point value is returned, representing the dot product 
  #   between the two vectors.

  total <- 0.0
  if (dimensions == 0) {
    return(total)
  }
  # Calculating the dot product.
  for (d in 1:dimensions) {
    total <- total + vect1[d] * vect2[d]
  }
  return(total)
}

TrajCheck <- function(traj1, traj2) {
  # A function to check that trajectories are in matrix form and have
  # the same dimension.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #
  # Returns:
  #   If there is a problem with one of the checks then a string containing 
  #   information is returned. If all of the checks pass then -1 is returned.
  
  # Checking that the trajectories are matrices.
  if (!(is.matrix(traj1)) | !(is.matrix(traj2))) {
    warning("At least one trajectory is not a matrix.")
    return("At least one trajectory is not a matrix.")
  }
  # Checking trajectories have the same dimension.
  if (dim(traj1)[2] != dim(traj2)[2]) {
    warning("The dimension of the trajectories does not match.")
    return("The dimension of the trajectories does not match.")
  } else {
    return(-1)
  }
}

AveTranslate <- function(traj1, traj2) {
  # A function to calculate a translation vector for trajectory 2 using the
  # average of the two trajectories.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #
  # Returns:
  #   A vector of length n is returned containing the translation in each 
  #   dimension. If there is a problem this returns a string of information 
  #   instead.
  
  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  # Calculating the number of points in each trajectory and the dimensions.
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # If a trajectory has no points or is dimension 0 then there is
  # no translation.
  if (length1 == 0 | length2 == 0 | dimensions == 0) {
    if (length1 == 0 | length2 == 0) {
      warning("At least one trajectory contains 0 points.")
    }
    if (dimensions == 0) {
      warning("The dimension is 0.")
    }
    return(rep(0.0, dimensions))
  }
  # Creating a translation vector to store the translations
  # of each dimension.
  translation <- vector()
  for (d in 1:dimensions) {
    newTrans <- (mean(traj1[, d]) - mean(traj2[, d]))
    translation <- c(translation, newTrans)
  }
  # Returning a vector of the translations.
  return(translation)
}

StartEndTranslate <- function(traj1, traj2) {
  # A function to calculate a variation of trajectory 2 by aligning
  # its start and end points with those of trajectory 1.
  #
  # Args:
  #   traj1: An m x n matrix containing trajectory1. Here m is the number 
  #          of points and n is the dimension of the points.
  #   traj2: A k x n matrix containing trajectory2. Here k is the number 
  #          of points and n is the dimension of the points. The two 
  #          trajectories are not required to have the same number of points.
  #
  # Returns:
  #   An m x n matrix containing the new variation of trajectory2 is 
  #   returned. Here m is the number of points and n is the dimension of 
  #   the points.
  
  # Checking the trajectories.
  trajTest <- TrajCheck(traj1, traj2)
  if (is.character(trajTest)) {
    return(trajTest)
  }
  # Calculating the number of points in each trajectory and the dimensions.
  dimensions <- dim(traj1)[2]
  length1 <- dim(traj1)[1]
  length2 <- dim(traj2)[1]
  # If a trajectory has no points or is dimension 0 then there is
  # no translation.
  if (length1 == 0 | length2 == 0 | dimensions == 0) {
    if (length1 == 0 | length2 == 0) {
      warning("At least one trajectory contains 0 points.")
    }
    if (dimensions == 0) {
      warning("The dimension is 0.")
    }
    return(traj2)
  }
  newTraj <- traj2
  for (d in 1:dimensions) {
    diff1 <- traj1[length1, d] - traj1[1, d]
    diff2 <- traj2[length2, d] - traj2[1, d]
    if (diff2 == 0) {
      warning("Equivalent start and end points in 1 dimension.")
    } else {
      for (point in 1:length2) {
        pointDiff <- traj2[point, d] - traj2[1, d]
        newTraj[point, d] <- (pointDiff / diff2) * diff1 + traj1[1, d]
      }
    }
  }
  return(newTraj)
}
