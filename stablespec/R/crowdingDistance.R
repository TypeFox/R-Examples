crowdingDistance <- function (pop, rnk, rng, varNo) {
  # All credit goes to Ching-Shih Tsou who originally wrote this function.
  # We borrow this function and modified it so as to meet our need.
  #
  # Ching-Shih Tsou (2013). nsga2R: Elitist Non-dominated Sorting
  # Genetic Algorithm based on R. R package version 1.0.
  # https://CRAN.R-project.org/package=nsga2R

  popSize <- nrow(pop)
  objDim <- length(rng)
  cd <- matrix(Inf, nrow = popSize, ncol = objDim)

  for (i in 1:length(rnk)) {

      len <- length(rnk[[i]])
      if (len > 2) {

        #For objective 1
        #to get the index of front i ordered ascendingly; Chi2
        originalIdx <- rnk[[i]][order(pop[rnk[[i]], 1])]

        #compute the crowded distance of objective 1
        cd[originalIdx[2:(len - 1)], 1] <- abs(pop[originalIdx[3:len],
                                              1] - pop[originalIdx[1:(len - 2)],
                                                      1])/rng[1]
        #For objective 2; model complexity (df)
        originalIdx <- rnk[[i]][order(pop[rnk[[i]], 2])]

        #compute the crowded distance of objective 2
        cd[originalIdx[2:(len - 1)], 2] <- abs(pop[originalIdx[3:len],
                                              2] - pop[originalIdx[1:(len - 2)],
                                                      1])/rng[2]
    }
  }
  return(as.matrix(apply(cd, 1, sum)))
}

