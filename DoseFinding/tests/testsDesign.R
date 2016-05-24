require("DoseFinding")
## Some examples from the JASA paper (for validation)
########################################################################
# Emax model p.1228 l. 5
fMod <- Mods(emax = 25, doses = c(0,150), placEff=0, maxEff=0.4)
fMod$emax[2] <- 0.6666667
doses <- c(0, 18.75, 150)
probs <- 1
deswgts1 <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD",
                      optimizer="Nelder-Mead")
deswgts2 <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD",
                      optimizer="nlminb")
## efficiency compared to standard design (last column)
crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                Delta=0.2, designCrit = "TD")
exp(deswgts1$crit - crt)

# Paper p. 1228 l. 2
fMod <- Mods(emax = 25, doses = c(0,150), placEff=0, maxEff=0.4)
doses <- c(0, 18.75, 150)
probs <- 1
deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
deswgts

########################################################################
#### exponential
# Paper p.1229 2nd line
fMod <- Mods(exponential=85, doses = c(0, 150), placEff=0, maxEff=0.4)
doses <- c(0, 50, 104.52, 150)
probs <- 1
deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD",
                     optimizer="Nelder-Mead")
deswgts
## efficiency compared to standard design (last column)
crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                Delta=0.2, designCrit = "TD")
exp(deswgts$crit - crt)

# Paper p.1229 1st line
fMod <- Mods(exponential=65, doses=c(0, 150), placEff=0, maxEff=0.4)
fMod$exponential[2] <- 0.08264711
doses <- c(0, 101.57, 150)
probs <- 1
deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
deswgts

########################################################################
#### Logistic
#### Paper: p.1230 7th line
fMod <- Mods(logistic=c(50, 10.881), doses = c(0, 150), placEff=0, maxEff=0.4)
doses <- c(0, 37.29, 64.44, 150)
probs <- 1
deswgts <- optDesign(fMod, probs, doses, Delta=0.05, designCrit = "TD")
deswgts
## efficiency compared to standard design (last column)
crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                Delta=0.05, designCrit = "TD")
exp(deswgts$crit - crt)


#### Paper p.1230 line 1
fMod <- Mods(logistic=c(50, 10.881), doses = c(0, 150), placEff=0, maxEff=0.4)
doses <- c(0, 50.22)
probs <- 1
deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
deswgts

########################################################################
#### beta
# Paper p.1230 line 5
fMod <- Mods(betaMod = c(0.33, 2.31), doses = c(0,150), addArgs=list(scal=200),
              placEff=0, maxEff=0.4)
doses <- c(0, 0.49, 25.2, 108.07, 150)
probs <- 1
deswgts <- optDesign(fMod, probs, doses, Delta=0.1,
                     control=list(maxit=1000), designCrit = "TD")
deswgts
## efficiency compared to standard design (last column)
crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                Delta=0.1, designCrit = "TD")
exp(deswgts$crit - crt)

# Paper p. 1230 line 10
fMod <- Mods(betaMod = c(1.39, 1.39), doses=c(0, 150), addArgs=list(scal=200),
              placEff=0, maxEff=0.4)
#doses <- c(0, 10, 25, 50, 100, 150)
doses <- c(0, 27, 94.89, 150)
probs <- 1
deswgts <- optDesign(fMod, probs, doses, Delta=0.1, designCrit = "TD")
deswgts
## efficiency compared to standard design (last column)
crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                Delta=0.1, designCrit = "TD")
exp(deswgts$crit - crt)

# Paper p. 1230 line 1
fMod <- Mods(betaMod = c(0.23, 2.31), doses=c(0,150), addArgs=list(scal=200),
             placEff=0, maxEff=0.4)
doses <- c(0, 0.35, 150)
probs <- 1
deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
deswgts
## efficiency compared to standard design (last column)
crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                Delta=0.2, designCrit = "TD")
exp(deswgts$crit - crt)

########################################################################
#### mixed Paper p. 1233, l. 2 (note the off and probably also the
#### scal parameter were treated as unknown in this example in the paper, 
#### hence the results need not be consistent with paper)
doses <- c(0, 9.9, 49.5, 115.4, 150)
fMod <- Mods(linear = NULL, emax = 25, exponential = 85,
             linlog = NULL, logistic = c(50, 10.8811),
             doses=doses, addArgs=list(off=1),
             placEff=0, maxEff=0.4)
probs <- rep(1/5, 5)
deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
deswgts2 <- optDesign(fMod, probs, doses, Delta=0.2, optimizer = "nlminb",
                      designCrit = "TD")

# Some other examples 
########################################################################
doses <- c(0, 62.5, 125, 250, 500)
fMod <- Mods(emax = c(25, 107.14), linear = NULL,
             logistic = c(150, 45.51), betaMod = c(1,1),
             doses = doses, addArgs=list(scal=1.2*500),
             placEff=60, maxEff=280)
probs <- rep(0.2, length=5)
deswgts <- optDesign(fMod, probs, Delta=200, designCrit = "TD")

########################################################################
#### using already allocated patients
fMod <- Mods(betaMod = c(0.33, 2.31), doses = c(0,150), addArgs=list(scal=200),
             placEff=0, maxEff=0.4)
doses <- c(0, 0.49, 25.2, 108.07, 150)
probs <- 1
# no previously allocated patients
deswgts <- optDesign(fMod, probs, doses=doses, Delta=0.1,
                     control=list(maxit=1000), designCrit = "TD")

# now use previously allocated patients
nold <- c(45, 50, 0, 0, 0)
deswgts2 <- optDesign(fMod, probs, doses=doses, Delta=0.1, n=30,
                      control=list(maxit=1000), nold=nold, designCrit = "TD")
# the overall design
(30*deswgts2$design+nold)/(30+sum(nold))
deswgts$design

########################################################################
#### Dopt Examples
doses <- c(0, 62.5, 125, 250, 500)
fMod <- Mods(emax = c(25, 107.14), logistic = c(150, 45.51),
             linear = NULL, betaMod = c(1,1),
             doses=doses, addArgs=list(scal=500*1.2),
             placEff=60, maxEff=280)
probs <- rep(0.2, 5)
des1 <- optDesign(fMod, probs, doses, Delta = 200, scal = 500*1.2, designCrit = "TD")
des2 <- optDesign(fMod, probs, doses, Delta = 200, scal = 500*1.2, designCrit = "Dopt")
des3 <- optDesign(fMod, probs, doses, Delta = 200, scal = 500*1.2, designCrit = "Dopt&TD")

########################################################################
#### optimizer = "exact" and "solnp"
doses <- c(0, 62.5, 125, 250, 500)
fMod <- Mods(emax = c(25, 107.14), logistic = c(150, 45.51),
             linear = NULL, betaMod = c(1,1),
             doses=doses, addArgs=list(scal=500*1.2),
             placEff=60, maxEff=280)
probs <- rep(0.2, 5)
des41 <- optDesign(fMod, probs, doses=doses, Delta = 200, n = 10,
                   optimizer = "exact", lowbnd = c(0.3,0,0,0,0), designCrit = "TD")
des42 <- optDesign(fMod, probs, doses=doses, Delta = 200, 
                   optimizer = "solnp", designCrit = "TD",
                   lowbnd = c(0.1,0,0,0,0))
des51 <- optDesign(fMod, probs, doses=doses, Delta = 200, n = 10,
                   designCrit = "Dopt", optimizer = "exact",
                   uppbnd = rep(0.5,5))
des52 <- optDesign(fMod, probs, doses=doses, Delta = 200, 
                   designCrit = "Dopt", optimizer = "solnp",
                   uppbnd = rep(0.5,5))
des61 <- optDesign(fMod, probs, doses=doses, Delta = 200, n = 10,
                   optimizer = "exact", designCrit = "Dopt&TD")
des62 <- optDesign(fMod, probs, doses=doses, Delta = 200, 
                   optimizer = "solnp", designCrit = "Dopt&TD")

########################################################################
#### Example from Padmanabhan and Dragalin, Biometrical Journal 52 (2010)
#### p. 836-852
fm <- Mods(sigEmax = c(4, 5), doses = 0:8,
           placEff=0, maxEff=-1.65)
fm$sigEmax <- c(0, -1.70, 4, 5)
## compare to Figure 1, p. 841
desSED <- optDesign(fm, 1, designCrit="Dopt", optimizer = "solnp")
desSEM <- optDesign(fm, 1, Delta = 1.3, designCrit = "TD",
                    optimizer = "solnp")

## designs underlying Table 2, p. 843 (from an e-mail of Vlad)
## I cannot reproduce the displayed efficiencies exactly
## (most probably due to numerical round-off)
##LDoD
## [1,] 0.246 0.141 0.123 0.000 0.000 0.240    0    0 0.250
## [2,] 0.248 0.233 0.061 0.210 0.000 0.000    0    0 0.248
## [3,] 0.246 0.000 0.000 0.223 0.081 0.204    0    0 0.246
## [4,] 0.250 0.247 0.045 0.210 0.000 0.000    0    0 0.248
## [6,] 0.250 0.249 0.192 0.062 0.000 0.000    0    0 0.246
## MEDoD
## [1,] 0.49 0.01 0.00 0.00 0.00 0.00 0.36 0.14    0
## [2,] 0.49 0.02 0.00 0.15 0.35 0.00 0.00 0.00    0
## [3,] 0.23 0.26 0.01 0.00 0.00 0.46 0.04 0.00    0
## [4,] 0.50 0.00 0.49 0.01 0.00 0.00 0.00 0.00    0
## [6,] 0.49 0.01 0.47 0.02 0.00 0.00 0.00 0.00    0
doses <- 0:8
fm <- list()
fm[[1]] <- Mods(sigEmax = c(23.07, 1.18), doses=doses, placEff=0, maxEff=-1.65);fm[[1]]$sigEmax <- c(0, -7.29, 23.07, 1.18)
fm[[2]] <- Mods(sigEmax = c(2, 2.22), doses=doses, placEff=0, maxEff=-1.65);fm[[2]]$sigEmax <- c(-0.08, -1.71, 2, 2.22)
fm[[3]] <- Mods(sigEmax = c(4, 5), doses=doses, placEff=0, maxEff=-1.65);fm[[3]]$sigEmax <- c(0, -1.70, 4, 5)
fm[[4]] <- Mods(sigEmax = c(0.79, 1), doses=doses, placEff=0, maxEff=-1.65);fm[[4]]$sigEmax <- c(0, -1.81, 0.79, 1.00)
fm[[5]] <- Mods(sigEmax = c(0.74, 1.18), doses=doses, placEff=0, maxEff=-1.65);fm[[5]]$sigEmax <- c(-0.03, -1.72, 0.74, 1.18)

desD <- desM <- matrix(ncol = 9, nrow = 5)
for(i in 1:5){
  cc1 <- optDesign(fm[[i]], 1, doses=doses, designCrit = "TD", optimizer = "solnp",
                   Delta = 1.3)
  cc2 <- optDesign(fm[[i]], 1, doses=doses, designCrit="Dopt", optimizer = "solnp")
  desM[i,] <- cc1$design
  desD[i,] <- cc2$design
}
round(desD, 3)
round(desM, 2)

## compare criterion for TD design under model 2
crDrag <- calcCrit(c(0.49,0.02,0,0.15,0.34,0,0,0,0), models=fm[[2]],
                   probs=1, doses=doses, designCrit="TD", Delta=1.3)
crDF <- optDesign(fm[[i]], 1, doses=doses, designCrit = "TD", optimizer = "solnp",
                  Delta = 1.3)$crit
exp(crDF-crDrag) ## design calculated by P and Dragalin only has 88% efficacy?


################################################################################
#### look at standardized Dopt and Dopt&TD criteria
doses <- c(0, 62.5, 125, 250, 500)
fMod1 <- Mods(sigEmax = rbind(c(25, 5), c(107.14, 2)), doses=doses, placEff=60, maxEff=280)
fMod2 <- Mods(sigEmax = rbind(c(25, 5), c(107.14, 2)), linear = NULL,
                  doses=doses, placEff=60, maxEff=280)
w1 <- rep(0.5, 2)
w2 <- rep(1/3, 3)
## des1 and des2 should be exactly the same
des1 <- optDesign(fMod1, w1, doses, designCrit = "Dopt", standDopt = FALSE)
des2 <- optDesign(fMod1, w1, doses, designCrit = "Dopt", standDopt = TRUE)

## des1 and des2 should be different (as linear and emax have
## different number of parameters)
des1 <- optDesign(fMod2, w2, doses, designCrit = "Dopt", standDopt = FALSE,
                  optimizer = "solnp")
des2 <- optDesign(fMod2, w2, doses, designCrit = "Dopt", standDopt = TRUE,
                  optimizer = "solnp")

## same with Dopt&TD criterion
## des1 and des2 will differ (due to different scaling
## of Dopt and TD criteria)
des1 <- optDesign(fMod1, w1, doses, designCrit = "Dopt&TD",
                  Delta = 100, standDopt = FALSE,
                  optimizer = "solnp")
des2 <- optDesign(fMod1, w1, doses, designCrit = "Dopt&TD",
                  Delta = 100, standDopt = TRUE,
                  optimizer = "solnp")


########################################################################
#### optimial design logistic regression
## compare this to Atkinson et al. (2007), p. 400
## theoretically the D-opt design should have weights 0.5,0.5 at points where
## the probability is 0.176 and 1-0.176 (0.3456 and 0.6544 in this case)
doses <- seq(0, 1, length = 21)
fMod <- Mods(linear = NULL, doses=doses, placEff=-5, maxEff = 10)
pp <- 1 # just one model
## by default calculates TD optimal design
mu <- as.numeric(getResp(fMod, doses=doses))
mu <- 1/(1+exp(-mu))
weights <- mu*(1-mu)
des1 <- optDesign(fMod, pp, doses, weights = weights, optimizer = "solnp") 
des2 <- optDesign(fMod, pp, doses, designCrit = "TD", Delta=0.2,
                  optimizer = "solnp", weights = weights)
des3 <- optDesign(fMod, pp, doses, Delta=0.2, designCrit = "Dopt&TD",
                  optimizer = "solnp", weights = weights)

########################################################################
#### code using lower and upper bound (previous to version 0.9-6 this
#### caused problems as the starting value for solnp rep(0.2, 5) was
#### on the boundary, now a feasible starting values is used
doses <- seq(0, 1, length=5)
nold <- rep(0, times=5)
lowbnd <- c(0.2,0.0,0.0,0.0,0.2)
uppbnd <- c(1.0,0.3,1.0,1.0,1.0)
trueModels <- Mods(linear=NULL, doses=doses, placEff = 0, maxEff = 1)
optDesign(models=trueModels, probs=1, doses=doses, designCrit="Dopt",
          lowbnd=lowbnd,uppbnd=uppbnd)

########################################################################
## TD optimal design for beta model (previously instabilities for
## numerical gradients)
mm <- Mods(betaMod=c(1.5,0.8), doses=seq(0,1,by=0.25), placEff=0, maxEff=1)
optDesign(mm, probs=1, designCrit="TD", Delta=0.5)
## Output from GUI
## placEff=0, maxEff=1
## TD-optimalität mit Delta= 0.5
## Model: BetaMod mit delta1=1.5, delta2=0.8
## Dosen      0        0.25     0.5      0.75    1
## Design     0.4895   0.3552   0.1448   0       0.0105
