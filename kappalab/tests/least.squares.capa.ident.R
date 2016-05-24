library(kappalab)


## the number of alternatives
n.a <- 300

## a randomly generated 5-criteria matrix
C <- matrix(rnorm(5*n.a,10,2),n.a,5)

## the corresponding global scores
g <- numeric(n.a)
mu <- capacity(c(0:29,29,29)/29)
for (i in 1:n.a)
  g[i] <- Choquet.integral(mu,C[i,])

## the full solution 
lsc <- least.squares.capa.ident(5,5,C,g)
a <- lsc$solution
a
mu.sol <- zeta(a)

## the difference between mu and mu.sol
mu@data - mu.sol@data

## the residuals
lsc$residuals

## the mean square error
mean(lsc$residuals^2)

## a 3-additive solution 
lsc <- least.squares.capa.ident(5,3,C,g)
a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals




## a similar example based on the Sipos integral

n.a <- 300
## a randomly generated 5-criteria matrix
C <- matrix(rnorm(5*n.a,0,2),n.a,5)

## the corresponding global scores
g <- numeric(n.a)
mu <- capacity(c(0:29,29,29)/29)
for (i in 1:n.a)
  g[i] <- Sipos.integral(mu,C[i,])

## the full solution 
lsc <- least.squares.capa.ident(5,5,C,g,Integral = "Sipos")
a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals

## a 3-additive solution 
lsc <- least.squares.capa.ident(5,3,C,g,Integral = "Sipos")
a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals




## additional constraints

## a Shapley preorder constraint matrix
## Sh(1) - Sh(2) >= -delta.S
## Sh(2) - Sh(1) >= -delta.S
## Sh(3) - Sh(4) >= -delta.S
## Sh(4) - Sh(3) >= -delta.S
## i.e. criteria 1,2 and criteria 3,4
## should have the same global importances
delta.S <- 0.01    
Asp <- rbind(c(1,2,-delta.S),
             c(2,1,-delta.S),
             c(3,4,-delta.S),
             c(4,3,-delta.S)
            )

## a Shapley interval constraint matrix
## 0.3 <= Sh(1) <= 0.9 
Asi <- rbind(c(1,0.3,0.9))


## an interaction preorder constraint matrix
## such that I(12) = I(45)
delta.I <- 0.01
Aip <- rbind(c(1,2,4,5,-delta.I),
             c(4,5,1,2,-delta.I))

## an interaction interval constraint matrix
## i.e. -0.20 <= I(12) <= -0.15 
delta.I <- 0.01
Aii <- rbind(c(1,2,-0.2,-0.15))

## an inter-additive partition constraint
## criteria 1,2,3 and criteria 4,5 are independent 
Aiap <- c(1,1,1,2,2)





## a more constrained solution

lsc <- least.squares.capa.ident(5,5,C,g,Integral = "Sipos",
                              A.Shapley.preorder = Asp,
                              A.Shapley.interval = Asi,
                              A.interaction.preorder = Aip,
                              A.interaction.interval = Aii,
                              A.inter.additive.partition = Aiap,
                              sigf = 5)

a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals
summary(a)

