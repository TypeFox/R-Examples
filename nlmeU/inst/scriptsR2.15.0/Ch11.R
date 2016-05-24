###################################################
### code chunk: Chap11init
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)
date()
packageVersion("nlme")
packageVersion("lattice")
sessionInfo()

library(nlme)
library(lattice)

###################################################
### code chunk number 7: R11.1
###################################################
 tx <- c(0, 10^-2, 0.8)          # Auxilary vector
 cX <-                           # corExp object defined
   corExp(value = c(1, 0.2),     # range rho:(10.16), nugget rho0: (10.18)
          form = ~tx, 
          nugget = TRUE)         # Nugget defined 
 Dtx <- data.frame(tx)         
 (cXi <-                         # corExp object initialized
   Initialize(cX, data = Dtx))
 (getCovariate(cXi))             # tx diffs: 2-1, 3-1, 3-2
 Vrg <- Variogram(cXi)           # Semi-variogram created ...
 plot(Vrg, smooth = FALSE,       # ... and plotted. Fig. 10.1a
      type = "l")
 corFunDt <-                     # Data for correlation function
    data.frame(dist = Vrg$dist,
               corF = 1 - Vrg$variog)
 plot(corFunDt,                  # Corr function plotted with ...
      type = "l", ylim = c(0,1)) # ... traditional graphics ... 
 xyplot(corF ~ dist,             # ... and xyplot().  Fig. 10.1b
    data = corFunDt, type = "l")


###################################################
### code chunk: Syntax for Fig. 10.1a
###################################################

xlab <- "Distance (s)"
myPanel1 <- function(x,y,subscripts,...){
  panel.xyplot(x, y, ...)
  panel.grid(h = -1, v = -1)
  panel.xyplot(x, y, ...) 
  lpoints(0, 0, type = "p", pch = 19, cex = 1.2) 
  lpoints(0, 0.2, type = "p", pch = 1, cex = 1.2) 
}
 
xyplot(variog ~ dist, 
       data = Vrg,
       type = "l",
       ylim = c(-0.1,1.1),
       panel = myPanel1, 
       ylab = "", 
       xlab = xlab)
  
###################################################
### code chunk: Syntax for Fig. 10.1b
###################################################

myPanel2 <- function(x,y,subscripts,...){
  panel.xyplot(x, y, ...)
  panel.grid(h = -1, v = -1)
  panel.xyplot(x, y, ...) # over-write grid
  lpoints(0, 1, type = "p", pch = 19, cex = 1.2) 
  lpoints(0, 0.8, type ="p", pch = 1, cex = 1.2) 
}

xyplot(corF ~ dist,
  data = corFunDt,
  type = "l",
  ylim = c(-0.1, 1.1),
  panel = myPanel2, 
  ylab = "",
  xlab= xlab
)

###################################################
### code chunk number 2: R11.2
###################################################
subj <- rep(1:2, each = 4)               # Two subjects
occ  <- rep(1:4, 2)                      # Four observations each
loc1 <- rep(c(0, 0.2, 0.4, 0.8), 2)      # First coordinate
loc2 <-                                  # Second coordinate 
   c(0, 0.2, 0.4, 0.8, 0, 0.1, 0.2, 0.4) 
df0  <-                                  # Hypothetical data frame
   data.frame(subj, occ, loc1, loc2)
(df  <-                                  # Occ = 3 for subj.2 deleted
   subset(df0, subj != 2 | occ != 3))    


###################################################
### code chunk number 3: R11.3
###################################################
cs <-                                       # Object defined...
   corCompSymm(value = 0.3, form = ~1|subj) 
cs <- Initialize(cs, df)                    # ... initialized
coef(cs, unconstrained = FALSE)             # Constrained coefficient
coef(cs)                       # Unconstrained = log((1/3+.3)/(1-.3))
getCovariate(cs)                            # Positions in series
corMatrix(cs)                               # Corr. matrix displayed

###################################################
### code chunk: R11.4
###################################################
cs1 <- corAR1(0.3, form = ~tx)   # Uninitialized corAR1 struct
coef(cs1, unconstrained = FALSE) # Constrained coefficient 
coef(cs1)                        # Unconstrained = log((1+.3)/(1-.3))
tx  <- 1:4                       # A covariate with values 1, 2, 3, 4
corMatrix(cs1, covariate = tx)   # Corr(Ri) of uninitialized object 
df2 <- data.frame(tx)            # An auxiliary data frame
cs1i <-                          # Initialized corAR1 object
   Initialize(cs1, data = df2) 
corMatrix(cs1i)                  # corAR1 matrix displayed
(chL <-                          # Cholesky factor L= (U') ^(-1)
   corMatrix(cs1i, corr = FALSE))
solve(t(chL) %*% chL)            # Back to Corr(Ri) = U'U =(L'L)^(-1) 

###################################################
### code chunk: R11.5a
###################################################
car <-                            # Not-recommended syntax ...
   corAR1(value = 0.3, form = ~1|subj) 
carI <- Initialize(car, df)       # corAR1 class object initialized
getCovariate(carI)                # Position=order of observations for a subject     
corMatrix(carI)[[1]]              # Correct matrix for the 1st subject
corMatrix(carI)[[2]]              # Incorrect matrix for the 2nd subject

###################################################
### code chunk: R11.5b
###################################################
car1 <- corAR1(value = 0.3, form = ~occ|subj)   # Recommended syntax
car1 <- Initialize(car1, df)      # corAR1 classs object initialized  
getCovariate(car1)                # Correct positions based on the occ variable
corMatrix(car1)[[2]]              # Correct matrix for the 2nd subject


###################################################
### code chunk: R11.6a
###################################################
ceE <- corExp(value = 1, form = ~loc1 + loc2 | subj)# Euclidean metric 
ceE <- Initialize(ceE, df)     
corMatrix(ceE)         # List with corr matrices for both subjects

###################################################
### code chunk: R11.6b
###################################################
ceM <-                                          # Manhattan metric
   corExp(1, ~ loc1 + loc2 | subj, metric = "man") 
ceM <- Initialize(ceM, df)
corMatrix(ceM)[[1]]              # Corr matrix for the 1st subject

###################################################
### code chunk: R11.6c
###################################################
ceEn <-                          # nugget = 0.2
   corExp(c(1, 0.2), ~ loc1 + loc2 | subj, nugget = TRUE) 
ceEn <- Initialize(ceEn, df)
coef(ceEn, unconstrained=FALSE)  # Constrained rho, rho0
corMatrix(ceEn)[[1]]             # Corr matrix for the 1st subject

### SessionInfo 
sessionInfo()     # with packages attached
detach(package:nlme)
detach(package:lattice)
