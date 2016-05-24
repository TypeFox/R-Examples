### R code from vignette source 'Rnews-paper.Rnw'

###################################################
### code chunk number 1: firstChunk
###################################################
## likelihood pivots
wald.stat <- function(theta, y) { 
    sqrt(2) * (y - theta) }
lik.root <- function(theta, y) { 
    sign(y - theta) * sqrt( 2 * log(1 + (y - theta)^2) ) }
score.stat <- function(theta, y) { 
    ( sqrt(2) * (y - theta) )/( 1 + (y - theta)^2 ) }
rstar <- function(theta, y) { 
    lik.root(theta, y) + 1/lik.root(theta, y) * 
        log( score.stat(theta, y)/lik.root(theta, y) ) }

## significance functions : Figure 1
theta.seq <- seq(-4, 4, length = 100)
par( las = 1, mai = c(0.9, 0.9, 0.2, 0.2) )
plot( theta.seq, pcauchy( q = 1.32 - theta.seq ), type = "l", lwd = 2,
      ylim = c(0,1), xlab = expression(theta), 
      ylab = "significance function", cex.lab = 1.5, cex.axis = 1.5 )
lines( theta.seq, pnorm( wald.stat(theta.seq, 1.32) ), lty = "dotted" )    
lines( theta.seq, pnorm( lik.root(theta.seq, 1.32) ), lty = "dashed" )
lines( theta.seq, pnorm( rstar(theta.seq, 1.32) ), lty = "solid" )
abline( v = 0, lty = "longdash" )


###################################################
### code chunk number 2: secondChunk
###################################################
## likelihood pivots
wald.stat <- function(theta, y) { 
    sqrt(2) * (y - theta) }
lik.root <- function(theta, y) { 
    sign(y - theta) * sqrt( 2 * log(1 + (y - theta)^2) ) }
score.stat <- function(theta, y) { 
    ( sqrt(2) * (y - theta) )/( 1 + (y - theta)^2 ) }
rstar <- function(theta, y) { 
    lik.root(theta, y) + 1/lik.root(theta, y) * 
        log( score.stat(theta, y)/lik.root(theta, y) ) }

## significance functions : Figure 1
theta.seq <- seq(-4, 4, length = 100)
par( las = 1, mai = c(0.9, 0.9, 0.2, 0.2) )
plot( theta.seq, pcauchy( q = 1.32 - theta.seq ), type = "l", lwd = 2,
      ylim = c(0,1), xlab = expression(theta), 
      ylab = "significance function", cex.lab = 1.5, cex.axis = 1.5 )
lines( theta.seq, pnorm( wald.stat(theta.seq, 1.32) ), lty = "dotted" )    
lines( theta.seq, pnorm( lik.root(theta.seq, 1.32) ), lty = "dashed" )
lines( theta.seq, pnorm( rstar(theta.seq, 1.32) ), lty = "solid" )
abline( v = 0, lty = "longdash" )


###################################################
### code chunk number 3: thirdChunk
###################################################
## exact p-value
round( 2 * ( min( tp <- pt(1.32, df = 1), 1 - tp ) ), digits = 3 )


###################################################
### code chunk number 4: fourthChunk
###################################################
## Wald pivot p-value
round( 2 * ( min( tp <- pnorm( wald.stat(0, 1.32) ), 1 - tp ) ), 
       digits = 3 )
## likelihood root p-value
round( 2 * ( min( tp <- pnorm( lik.root(0, 1.32) ), 1 - tp ) ), 
       digits = 3 )
## modified likelihood root p-value
round( 2 * ( min( tp <- pnorm( rstar(0, 1.32) ), 1 - tp ) ), 
       digits = 3 )


###################################################
### code chunk number 5: fifthChunk
###################################################
## simulated data
library(hoa)
set.seed(321)
y <- rt(n = 15, df = 3)
y.rsm <- rsm(y ~ 1, family = student(3))
y.cond <- cond(y.rsm, offset = 1)
summary(y.cond, test = 0)


###################################################
### code chunk number 6: sixthChunk
###################################################
## `airway' data
#library(cond)
head( airway, n = 3 )  


###################################################
### code chunk number 7: seventhChunk
###################################################
## binomial model fit
airway.glm <- glm( formula(airway), family = binomial, data = airway )
summary( airway.glm ) 


###################################################
### code chunk number 8: eighthChunk
###################################################
## higher order inference
airway.cond <- cond( airway.glm, offset = type1 )
summary( airway.cond )    #    produces 95% confidence intervals
plot(airway.cond, which = 4)    # Figure 2


###################################################
### code chunk number 9: ninthChunk
###################################################
plot(airway.cond, which = 4)    # Figure 2


###################################################
### code chunk number 10: tenthChunk
###################################################
## `calcium' data
library(boot) 
head( calcium, n = 3 )


###################################################
### code chunk number 11: eleventhChunk
###################################################
#library(nlreg)

## maximum likelihood fit
calcium.nl <- nlreg( cal ~ b0 * (1 - exp(-b1 * time)), 
                     weights = ~ (1 + time)^g, data = calcium, 
                     start = c(b0 = 4, b1 = 0.1, g = 0) )
summary( calcium.nl )    # yields estimates and standard errors


###################################################
### code chunk number 12: twelfthChunk
###################################################
## pivot profiling for \gamma
calcium.prof <- profile( calcium.nl, offset = g )
summary( calcium.prof )  


###################################################
### code chunk number 13: thirteenthChunk
###################################################
## inference on proportion of maximum
calcium.nl <- nlreg( cal ~ b0 * (1 - exp(- log(1 + exp(psi)) * time / 15)), 
                     data = calcium, start = c(b0 =4.3, psi =2) )
calcium.prof <- profile( calcium.nl, offset = psi )
calcium.sum <- summary( calcium.prof )    
exp(calcium.sum$CI) / (1 + exp(calcium.sum$CI))    # 95% confidence intervals for \pi


###################################################
### code chunk number 14: fourteenthChunk
###################################################
## profile and contour plots : Figure 3
calcium.prof <- profile( calcium.nl )
par( las = 1, mai = c(0.5, 0.5, 0.2, 0.2) )
contour( calcium.prof, alpha = 0.05, cl1 = "black", cl2 = "black", 
         lwd2 = 2 )    


###################################################
### code chunk number 15: fifteenthChunk
###################################################
## profile and contour plots : Figure 3
par( las = 1, mai = c(0.5, 0.5, 0.2, 0.2) )
contour( calcium.prof, alpha = 0.05, cl1 = "black", cl2 = "black",
         lwd2 = 2 )    


