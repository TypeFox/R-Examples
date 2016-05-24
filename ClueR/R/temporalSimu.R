#' Temporal data simulation
#' 
#' This function simulates time-series data using 14 pre-defined temporal profile templates. Type 'temporalSimu' to see the details of the templates.
#' 
#' @param seed to seed the simulation. Default is current system time.
#' @param groupSize the number of the temporal profiles to simulate from each template. The total number of profiles will be the number of templates used times the size of each group.
#' @param sdd the standard deviation to be used to generate randomness for each temporal profile.
#' @param numGroups number of templates to be used for generating data.
#' @return a matrix containing simulated time-series dataset.
#' @export
#' @examples
#' # simulate a time-series data with four distinctive profile groups and each group with 
#' # a size of 500 phosphorylation sites
#'
#' simulated.temporal <- temporalSimu(seed=1, groupSize=500, sdd=1, numGroups=4)
#' 
temporalSimu <- function(seed=unclass(Sys.time()), groupSize, sdd, numGroups) {
   set.seed(seed)

   ## create 14 temporal profile templates
   # accute up response
   group1 <- cbind(rnorm(groupSize, 0, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd))
   # delayed up response
   group2 <- cbind(rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 4, sdd))

   # accute down response
   group3 <- cbind(rnorm(groupSize, 4, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd))
   # delayed down response
   group4 <- cbind(rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 0, sdd))
   
   # median up response
   group5 <- cbind(rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd))
   # median down response
   group6 <- cbind(rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd))

   # v shape up then down response
   group7 <- cbind(rnorm(groupSize, 4, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 4, sdd))
   # v shape down then up response
   group8 <- cbind(rnorm(groupSize, 0, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 3, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 0, sdd))
   
   # v shape up then down then flat response
   group9 <- cbind(rnorm(groupSize, 0, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd))
   # v shape down then up then flat response
   group10 <- cbind(rnorm(groupSize, 4, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd), rnorm(groupSize, 4, sdd))
   
   # over shoot response
   group11 <- cbind(rnorm(groupSize, 0, sdd), rnorm(groupSize, 0.2, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 5, sdd), rnorm(groupSize, 5, sdd), rnorm(groupSize, 5, sdd), rnorm(groupSize, 5, sdd))
   # under shoot reponse
   group12 <- cbind(rnorm(groupSize, 5, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 0.2, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 0, sdd))

   # s shape responses
   group13 <- cbind(rnorm(groupSize, 0, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 2, sdd))
   group14 <- cbind(rnorm(groupSize, 2, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 0, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 2, sdd), rnorm(groupSize, 1, sdd), rnorm(groupSize, 0, sdd))

   groups <- list()
   groups[[1]] <- group1
   groups[[2]] <- group2
   groups[[3]] <- group3
   groups[[4]] <- group4
   groups[[5]] <- group5
   groups[[6]] <- group6
   groups[[7]] <- group7
   groups[[8]] <- group8
   groups[[9]] <- group9
   groups[[10]] <- group10
   groups[[11]] <- group11
   groups[[12]] <- group12
   groups[[13]] <- group13
   groups[[14]] <- group14
   
   
   simulated.dat <- c()
   # random sampling from the 14 templates
   ints <- sample(1:14, numGroups, replace=FALSE)

   simulated.dat <- groups[[ints[1]]]
   for (i in 2:numGroups){
		simulated.dat <- rbind(simulated.dat, groups[[ints[i]]])
   }

   # standardize the simulated data
   tmp <- sweep(simulated.dat, 1, apply(simulated.dat, 1, mean), FUN="-")
   simulated.stand <- sweep(tmp, 1, apply(simulated.dat, 1, sd), FUN="/")
   
   rownames(simulated.stand) <- paste("p", 1:nrow(simulated.stand), sep="_")
   colnames(simulated.stand) <- paste("t", 1:ncol(simulated.stand), sep="_")
   return(simulated.stand)
}
