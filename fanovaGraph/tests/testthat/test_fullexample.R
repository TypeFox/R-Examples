#################################################################
# One full application example of the package 'fanovaGraph'
#################################################################

test_that("the full example works as before", {
  
d <- 6     
domain <- c(-1,1)

Bsp <- function(x){
  beta <- c(-0.8, -1.1, 1.1, 1)
  gamma <- c(-0.5, 0.9, 1, -1.1)
  result <- cos(cbind(1, x[,c(1,5,3)]) %*% beta) + 
    sin(cbind(1, x[,c(4,2,6)]) %*% gamma)   # function a
  return(result)
}

### maximin design with package "lhs" 

data(L)

### kriging model with package "DiceKriging"

y <- Bsp(L)
covtype <- "matern5_2"
set.seed(1)
KM <- km( ~ 1, design = data.frame(L), response = y, covtype = covtype,)

### prediction function (the new metamodel!)

krigingMean <- function(Xnew) 
  predict(object = KM, newdata = Xnew, type = "UK", se.compute = FALSE, checkNames=FALSE)$mean

### model validation

plot(KM)
par(mfrow = c(1, 1))
npred <- 1000
xpred <- matrix(runif(d * npred,domain[1],domain[2]), ncol = d) 
yexact <- Bsp(xpred)
yKM <- krigingMean(xpred)
plot(yexact, yKM, asp = 1)
abline(0, 1)

expect_equal(0.13046535, sqrt(mean((yKM- yexact)^2)))

### standard Sobol indices with package "sensitivity"
### and function "krigingMean"

soverall <- fast99(model = krigingMean, factors = d, n = 2000, 
                   q = "qunif", q.arg = list(min = domain[1], max = domain[2]))
expect_equal(0.7772051, mean(soverall$V))

#############################################################
### estimate total interaction indices for the graph edges
#############################################################

### estimation by function "estimateGraph" with fixing method
set.seed(1)
totalInt <- estimateGraph(f.mat=krigingMean, d=d, n.tot=10000,
                          q.arg=list(min=domain[1],max=domain[2]), method="FixFast")
expect_equivalent(c(0.003276626653, 0.009495300733), totalInt$tii[1:2,1])  

set.seed(1)
totalInt <- estimateGraph(f.mat=krigingMean, d=d, n.tot=10000,
                          q.arg=list(min=domain[1],max=domain[2]), method="LiuOwen")
expect_equivalent(c(0.001869581276145, 0.043586803962427), totalInt$tii[1:2,1])

##########################################
### Graph plotting
##########################################

### creating the graph with cut at indices > 0.01 by function "threshold"

graph <- threshold(totalInt, delta = 0.01)

cl <- graph$cliques
expect_equal(list(c(1,3,5),c(2,4,6)),list(as.numeric(cl[[1]]), as.numeric(cl[[2]])))


##############################
### estimating the new model
##############################

## new kernel estimation by function "kmAdditive"
## and prediction by function "predictAdditive"

cl <- list(c(1, 3, 5),c(2, 4, 6))
eps.R <-  1e-04
eps.Var <- 1e-4
nMaxit <- 30
nInitial <- 20

set.seed(1)
parameter <- kmAdditive(x=L, y, n.initial.tries = nInitial, eps.R = eps.R, 
      cl=cl, covtype = covtype, eps.Var = eps.Var, max.it = nMaxit, iso = FALSE)
expect_equal(0.517164345, parameter[[1]]$alpha)
ypred <- predictAdditive(xpred, x=L, y, parameter, covtype = covtype, eps.R = eps.R, cl)
expect_equal(0.05050225006, sqrt(mean((ypred$mean - yexact)^2)))

### isotropic clique
set.seed(1)
parameter <- kmAdditive(x=L, y, n.initial.tries = nInitial, eps.R = eps.R, 
                                cl=cl, covtype = covtype, eps.Var = eps.Var, 
                                max.it = nMaxit, iso = c(FALSE,TRUE))
expect_equal(1.902815638, parameter[[2]]$theta)
ypred <- predictAdditive(xpred, x=L, y, parameter, covtype = covtype, eps.R = eps.R, 
              cl, iso=c(FALSE,TRUE))
expect_equal(0.04277603821, sqrt(mean((ypred$mean- yexact)^2)))

### two isotropic cliques
set.seed(1)
parameter <- kmAdditive(x=L, y, n.initial.tries = nInitial, 
                                eps.R = eps.R, cl=cl, covtype = covtype, eps.Var = eps.Var, 
                                max.it = nMaxit, iso = c(TRUE,TRUE))
expect_equal( 1.904219623, parameter[[2]]$theta)
ypred <- predictAdditive(xpred, x=L, y, parameter, covtype = covtype, eps.R = eps.R, 
              cl, iso=c(TRUE,TRUE))
expect_equal( 0.04246029631, sqrt(mean((ypred$mean- yexact)^2)))

### single clique (km is used for estimation and prediction)
cl <- list(1:6)
set.seed(1)
parameter <- kmAdditive(x=L, y, cl=cl, covtype=covtype)
expect_equal(parameter@covariance@sd2, KM@covariance@sd2)
ypred <- predictAdditive(xpred, x=L, y, parameter, covtype = covtype, cl=cl)
expect_true(all(ypred$mean == yKM))
})

test_that("simulation works", {
cl <- list(c(1, 3, 5),c(2, 4, 6))
d <- 6
parameter <- list(list(alpha=0.53,theta=c(1.7,1.8,1.8)),
                  list(alpha=0.47,theta=c(1.9,1.9,1.7)))
set.seed(1)
xsimu <- matrix(runif(d*3,-1,1),3,d)

expect_equal(c(-0.305388389,  0.792219843,  0.313896841),
       simAdditive(xsimu, mu=0, parameter, "matern5_2", cl, iso=FALSE))    

parameter <- list(list(alpha=0.53,theta=c(1.7,1.8,1.8)),
                  list(alpha=0.47,theta=c(1.1)))
set.seed(1)
xsimu <- matrix(runif(d*3,-1,1),3,d)
expect_equal(c(-0.305388389, 0.812151639, 0.288423824),
     simAdditive(xsimu, mu=0, parameter, "matern5_2", cl, iso=c(FALSE,TRUE)))
})
