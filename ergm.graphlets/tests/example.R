#
# Toy Example: For CRAN test
#

# 1 - Load the necessary libraries 
library("ergm.graphlets")

# 2 - Load the EMON dataset from statnet package for the illustration
data("emon3")
data("spi")

# 3- Estimate a simple model with the network
emon.ergm <- ergm(emon.3 ~ graphletCount(0), control = control.ergm(seed=1))

summary(emon.ergm)




# ORIGINAL EXAMPLES IN THE PAPER (commented out due to huge computational cost -- slowing down the CRAN tests)
# NOTE: This procedure also includes the data preprocessing which we simplified with the data("emon3") and data("spi") objects
#
# EXAMPLE 1: Lake Pomona emergent multi-organizational network (EMON)
#
# Note: the execution of this example takes a long time
#

# 5 - Estimate the parameters for the ERGM model
# emon.ergm <- ergm(emon.3 ~ edges + nodefactor("Sponsorship") + nodecov("Command.Rank.Score") + grorbitFactor("Location", c(9:11)), control = control.ergm(seed=1, MCMC.samplesize=50000, MCMC.interval=100000, MCMC.burnin=50000, parallel=60))

# 6 - Get the summary of the estimated model
# summary(emon.ergm)

# 7 - Evalute the goodness-of-fit for the estimated model
# EMONgof <- gof(emon.ergm, GOF = ~ degree + distance + espartners + triadcensus)
## par(mfrow=c(2, 2))
## plot(EMONgof)





#
# Example 2: Protein secondary structure network 
#
# Note: The execution of this example may take a long time depending performance of the tested machine

# 2 - Load the protein structure network of serine protease inhibitor
#spi <- read.table("serineProteaseInhibitor.txt")		# Modify this path to show the correct file location, if it fails to load 
#spi <- as.matrix(spi)
#attr(spi, "n") <- 53
#spi <- as.sociomatrix.sna(spi)
#spi <- network(spi, directed=F)

# 3 - Define the assembly node attribute as a categorical node attribute
#spi.assembly <- 1:53
#spi.assembly[1:25] <- "1"
#spi.assembly[26:53] <- "2"
#spi %v% "Assembly" <- spi.assembly

# 4 - Estimate the parameters of the first ERGM model
#spi.ergm.34678 <- ergm(spi ~ edges + nodematch("Assembly") + triangle("Assembly") + gwdegree(.5, fixed=T) + graphletCount(c(3, 4, 6, 7, 8)),control = control.ergm(seed=1, MCMC.samplesize=500000, MCMC.interval=75000, MCMC.burnin=300000, parallel=60))

# 5 - Get the summary of the first estimated model
#summary(spi.ergm.34678)

# 6 - Estimate the parameters of the second ERGM model
#spi.ergm.all<-ergm(spi ~ edges + nodematch("Assembly") + triangle("Assembly") + gwdegree(.5, fixed=T) + graphletCount(c(6, 7, 8)), control = control.ergm(seed=1, MCMC.samplesize=15000, MCMC.interval=2000, MCMC.burnin=15000))

# 7 - Get the summary of the second estimated model
#summary(spi.ergm.all)

# 8 - Evalute the goodness-of-fit for the second estimated model
#SPIgof <- gof(spi.ergm.all, GOF = ~ degree + distance + espartners + triadcensus)
##par(mfrow=c(2, 2))
##plot(SPIgof)
