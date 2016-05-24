################################################################################
# 
# HuiFunctions.R
# Version 1.0
# 13/04/2015
#
# Updates:
# 
# Set of functions for Hui model
#
################################################################################

prob_q11 <- function(presences, unit){
  # Function to calulcate q+/+ which is the conditional probability that a 
  # randomly chosen cell adjacent to an occupied cell is also occupied
  #
  # Args:
  #   presences: data frame of presence or presence-absence data along at atlas scale
  #   unit: Resolution of atlas data - distance between cell centres
  coords <- presences[presences$presence == 1, c("lon", "lat")]
  QQ <- matrix(NA, nrow = length(coords))
  for (i in 1:nrow(coords)){ 
    # loop through each of the focal occupied cells calculating the proportion 
    # of neightbouring cells that are occupied for each one
    ## Easting and Northing of occupied focal cell
    E <- coords[i, "lon"] 
    N <- coords[i, "lat"] 
    # four adjacent neighbours
    n1 <- presences[presences$lon == (E + unit) & 
                      presences$lat == N, "presence"]
    n2 <- presences[presences$lon == (E - unit) & 
                      presences$lat == N, "presence"]
    n3 <- presences[presences$lon == E & 
                      presences$lat == (N + unit), "presence"]
    n4 <- presences[presences$lon == E & 
                      presences$lat == (N - unit), "presence"]
    # number of neighbours adjacent to each occupied cell (0-4)
    ncells <- length(n1) + length(n2) + length(n3) + length(n4)
    occ <- sum(n1, n2, n3, n4, na.rm = TRUE)
    QQ[i]<- occ / ncells 
  }
  # q+/+ for a species is mean of all proportions for all occupied focal cells
  return(mean(QQ, na.rm = TRUE))
}

prob_absence <- function(p0_fine, n, q00, p0_coarse) {
  # Function to calculate the probability of absence at the finer scale
  #
  # Args:
  #   p0_fine: probability of absence at fine grain - this is the unknown.
  #   n: ratio between atlas grid size and fine grid size
  #       eg. if atlas scale = 10km2, fine scale = 2km2, n = 5.
  #   q00: the conditional probability that a randomly chosen cell adjacent to
  #        an empty cell is occupied.
  #   p0_coarse: observed probability of absence at coarse grain.
  part1 <- p0_fine * (((p0_fine ^ (-(1 / (-1 + n)))) * 
                    (p0_coarse ^ (1 / (-1 + n))) * 
                    (q00 ^ (-1 / n)))
                 ^ (2 * (-1 + n)))
  part2 <- (p0_fine ^ (1 - (2 / (-1 + n)))) * 
    (p0_coarse ^ (2 / (-1 + n))) * 
    (q00 ^ (-2 / n))
  part3 <- (p0_fine ^ (1 - (2 / (-1 + n)))) * 
    (p0_coarse ^ (2 / (-1 + n))) * 
    (q00 ^ (-2 / n))
  part4 <- (p0_fine ^ 2) * 
    ((1 - 
        (p0_fine ^ (-(1 / (-1 + n)))) *
        (p0_coarse ^ (1 / (-1 + n))) *
        (q00 ^ (-1 / n))) ^ 2) 
  part5 <- (1 - p0_fine)
  return(part1 * ((part2 / (part3 + (part4 / part5))) ^ ((-1 + n) ^ 2)))
}
