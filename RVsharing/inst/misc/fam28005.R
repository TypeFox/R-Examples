# Example of Syrian family 28005, shown in Figure 3B of Bureau et al. 2013
library(Bureau)
trio1 <- new("Trio", id = "33", spouse = "46", offspring = list("48"))
trio2 <- new("Trio", id = "42", spouse = "50", offspring = list("66"))
trio3 <- new("Trio", id = "23", spouse = "22", offspring = list(trio1))
trio4 <- new("Trio", id = "18", spouse = "19", offspring = list(trio2))
trio5 <- new("Trio", id = "8", spouse = "9", offspring = list(trio3))
trio6 <- new("Trio", id = "10", spouse = "9", offspring = list(trio4))
trio7 <- new("Trio", id = "1", spouse = "2", offspring = list(trio5,trio6))
geno.vec <- c(0,0,NA,0,NA,NA,0,0,NA,NA,NA,0,NA,0,NA)
names(geno.vec) = c(1,2,8:10,18,19,22,23,33,42,46,48,50,66)

seed = .Random.seed
.Random.seed = seed
nrep=1e4
share.vec <- occur.vec <- logical(nrep)
for(i in 1:nrep ){
       founder <- sample(c(1,2,4,7,8,12,14),1,replace = FALSE)
       geno.vec[founder] <- 1
       geno.vec.sim <- GeneDrop(trio7, geno.vec)
       share.vec[i] <- geno.vec.sim["48"]==1 & geno.vec.sim["66"]==1
       occur.vec[i] <- geno.vec.sim["48"]==1 | geno.vec.sim["66"]==1
       geno.vec[founder] <- 0

 }
sum(share.vec)/sum(occur.vec)

# The exact value is 0.01185771