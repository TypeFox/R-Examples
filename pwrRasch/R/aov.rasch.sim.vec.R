##########################################################################################################
#
# pwrRasch: Statistical Power Simulation for Testing the Rasch Model
#
# Internal function: Vectorized computation of Three-Way ANOVA for balanced design
#
# Authors: Takuya Yanagida <takuya.yanagida@univie.ac.at>
#  		     Jan Steinfeld <jan.steinfeld@univie.ac.at>
#
##########################################################################################################

aov.rasch.sim.vec <- function(data) {
  
  group <- data$group
  person <- data$person
  item <- data$item
  response <- data[, grep("response", colnames(data)):ncol(data)]
  
  a <- length(unique(group))
  b <- length(unique(person))
  c <- length(unique(item))  
  
  ###
  
  SS.T <- (1/(a*b/a*c))*colSums(response)^2
  
  SS.A <- colSums(rowsum(response, group)^2)*(1/(b/a*c))
  
  SS.B <- colSums(rowsum(response, person)^2)*(1/c)
  
  SS.C <- colSums(rowsum(response, item)^2)*(1/(a*b/a))
  
  ###
  
  SS.interm. <- rbind(rowsum(response[group == 1,], item[group == 1]),
                      rowsum(response[group == 2,], item[group == 2]))
  
  SS.BC <- colSums(SS.interm.)
  
  SS.AC <- colSums(SS.interm.^2)*(1/(b/a))
  
  ###
  
  F.AC <- ((SS.AC - SS.A - SS.C + SS.T)/(c - 1)) / ((SS.BC - SS.B - SS.AC + SS.A)/(a*(b/a - 1)*(c - 1)))
  p.AC <- pf(F.AC, c - 1, a*(b/a - 1)*(c - 1), lower.tail = FALSE)
  
  ###
  
  return(p.AC)
  
}