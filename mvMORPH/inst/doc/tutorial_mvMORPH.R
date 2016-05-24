## ---- comment=">"--------------------------------------------------------
# Load the package and dependencies (ape, phytools, corpcor, subplex, spam)
library(mvMORPH)
# Use a specified random number seed for reproducibility
set.seed(14)


## ------------------------------------------------------------------------
tree<-pbtree(n=100)

## ---- comment=">"--------------------------------------------------------
# Simulate two selective regimes
state<-as.vector(c(rep("Forest",60),rep("Savannah",40))); names(state)<-tree$tip.label

# Make the tree with mapped states using SIMMAP
tree<-make.simmap(tree, state, model="ER", nsim=1)


## ---- comment=">"--------------------------------------------------------
# Plot the phylogeny with the mapped discrete trait
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")
plotSimmap(tree,col, fsize=0.6, node.numbers=FALSE, lwd=3, pts=FALSE)


## ---- comment=">"--------------------------------------------------------
# 2 Random traits evolving along the phylogeny as a two-optimum OU
set.seed(101)

alpha<-matrix(c(1.1,-0.9,-0.9,1),2)
sigma<-matrix(c(0.35,0.06,0.06,0.35),2)
theta<-c(5.5,5.1,1.2,1.4)
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, mu=theta, 
names_traits=c("limb.length","limb.width")), model="OUM", nsim=1)


## ---- comment=">"--------------------------------------------------------
# Fitting the Ornstein Uhlenbeck on the whole tree 
trait1_OU1<- mvOU(tree, data[,1], model="OU1", diagnostic=FALSE, echo=FALSE)
trait2_OU1<- mvOU(tree, data[,2], model="OU1", diagnostic=FALSE, echo=FALSE)
# Fitting the Ornstein Uhlenbeck with multiple optimums 
trait1_OUM<- mvOU(tree, data[,1], model="OUM", diagnostic=FALSE, echo=FALSE)
trait2_OUM<- mvOU(tree, data[,2], model="OUM", diagnostic=FALSE, echo=FALSE)

# Compare the AIC values between models fit
AIC(trait1_OUM); AIC(trait1_OU1)
AIC(trait2_OUM); AIC(trait2_OU1)

# Now compare with the multivariate fit
OUM<- mvOU(tree, data, model="OUM")
OU1<- mvOU(tree, data, model="OU1")
AIC(OUM); AIC(OU1)


## ---- eval=FALSE---------------------------------------------------------
#  
#  # Simulate 1000 traits under the two optima (OUM) and the unique optimum OU (OU1) process
#  library(parallel)
#  nsim=1000
#  
#  # Dataset simulated with the OUM maximum likelihood estimates
#  data1<-simulate(OUM, nsim=nsim, tree=tree)
#  # Dataset simulated with the MLE of the OU1 model
#  data2<-simulate(OU1, nsim=nsim, tree=tree)
#  
#  # Fit of the models using the parallel package (we will use 2 cores), can take a while...
#  
#  library(parallel)
#  nb_cores=2L
#  oum_data1<- mclapply(1:nsim, function(x){
#      mvOU(tree, data1[[x]], model="OUM", method="sparse", diagnostic=F, echo=F)
#  }, mc.cores = getOption("mc.cores", nb_cores))
#  
#  ou1_data1<- mclapply(1:nsim, function(x){
#      mvOU(tree, data1[[x]], model="OU1", method="sparse", diagnostic=F, echo=F)
#  }, mc.cores = getOption("mc.cores", nb_cores))
#  
#  
#  # Now same simulations on the second dataset
#  oum_data2<- mclapply(1:nsim, function(x){
#      mvOU(tree, data2[[x]], model="OUM", method="sparse", diagnostic=F, echo=F)
#  }, mc.cores = getOption("mc.cores", nb_cores))
#  
#  ou1_data2<- mclapply(1:nsim, function(x){
#      mvOU(tree, data2[[x]], model="OU1", method="sparse", diagnostic=F, echo=F)
#  }, mc.cores = getOption("mc.cores", nb_cores))
#  
#  # Retrieve the results from the simulations
#  OUM_simul<-sapply(1:nsim, function(x){
#      c(oum_data1[[x]]$AICc,ou1_data1[[x]]$AICc)
#  })
#  
#  OU1_simul<-sapply(1:nsim, function(x){
#      c(oum_data2[[x]]$AICc,ou1_data2[[x]]$AICc)
#  })
#  
#  # Now compute the type I error and power (type II)
#  sum(OU1_simul[1,]<OU1_simul[2,])/nsim
#  [1] 0.135
#  sum(OUM_simul[1,]<OUM_simul[2,])/nsim
#  [1] 1
#  
#  

## ---- comment=">"--------------------------------------------------------
# We now try to test for significant "selective" interactions toward the optima
# First: we fit a OUM model without interactions
OUM_1<-mvOU(tree, data, model="OUM", param=list(alpha="diagonal"), 
            echo=FALSE, diagnostic=FALSE)
AIC(OUM_1)

# We then compare the results to the original fit :
AIC(OUM)

# Log-likelihood ratio test
LRT(OUM,OUM_1)

## ---- comment=">"--------------------------------------------------------
stationary(OUM)

# we can use the cov2cor function from the R stats package
# to get the evolutionary correlations between traits
cov2cor(stationary(OUM))

## ----results="hide", message=FALSE, comment=">"--------------------------
# Simulated dataset
set.seed(14)

# Generating a random tree
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states

tree<-make.simmap(tree,sta , model="ER", nsim=1)


###-------------------------Parameters------------------------###
# Define the correlation structure in each group
C1<-matrix(c(1,0.7,0.7,1),2)
C2<-matrix(c(1,0.3,0.3,1),2)

# Define the variance in each group
D1<-diag(c(sqrt(0.2),sqrt(0.1)))
D2<-diag(c(sqrt(0.03),sqrt(0.05)))

# Compute the rate matrices
sigma_1<-D1%*%C1%*%D1
sigma_2<-D2%*%C2%*%D2

# Ancestral states
theta<-c(0,0)

###---------------------Simulate the traits------------------###

data<-mvSIM(tree, param=list(sigma=list(sigma_1,sigma_2), ntraits=2, mu=theta, 
                             names_traits=c("Trait 1","Trait 2")), model="BMM", nsim=1)


## ---- comment=">"--------------------------------------------------------
# Fitting the models

# BM1 - (Equal rate matrix)
model_1<-mvBM(tree, data, model="BM1", diagnostic=FALSE, echo=FALSE)

# BMM - (Proportional rate matrices)
model_2<-mvBM(tree, data, param=list(constraint="proportional")
              , diagnostic=FALSE, echo=FALSE)

# BMM - (Shared eigenvectors between rate matrices)
model_3<-mvBM(tree, data, param=list(constraint="shared")
              , diagnostic=FALSE, echo=FALSE)

# BMM - (Similar correlations between rate matrices)
model_4<-mvBM(tree, data, param=list(constraint="correlation")
              , diagnostic=FALSE, echo=FALSE)

# BMM - (independent rate matrices)
model_5<-mvBM(tree, data, model="BMM", diagnostic=FALSE, echo=FALSE)

# Compare the models with AIC
AIC(model_1)

AIC(model_2)

AIC(model_3)

AIC(model_4)

AIC(model_5)

# Test significance with LRT

LRT(model_5,model_4)

LRT(model_5,model_3)

LRT(model_5,model_2)

LRT(model_5,model_1)


## ---- comment=">"--------------------------------------------------------
# Forest species
cov2cor(model_5$sigma[,,1])

# Savannah species
cov2cor(model_5$sigma[,,2])

