date()
packageVersion("nlmeU")
message("Code for first part of Chapter 17 (from Panel R17.1 til R17.13) is executed below.") 
message("Code for remaining panels pertaining to pdKronecker class is distributed with nlmeUpdK package (Jun.10, 2013)")

###################################################
### code chunk: Chap17init
###################################################
options(digits = 5, show.signif.stars = FALSE)
sessionInfo()
library(nlme)
library(lattice)

###################################################
### code chunk: R17.1a
###################################################
data(prt, package = "nlmeU")
lme.spec.form1 <- 
   formula(spec.fo ~ (prt.f + occ.f)^2 + sex.f + age.f + 
             sex.f:age.f + bmi) 
prt1 <- subset(prt, fiber.f == "Type 1", select = -fiber.f)
fm17.1 <-                                        # M17.1 (17.1)
   lme(lme.spec.form1,                    
       random = ~occ.f - 1|id,                   # Random effects structure (including D)
       data = prt1) 


###################################################
### code chunk: R17.1b
###################################################
getGroupsFormula(fm17.1)
str(grpF <- getGroups(fm17.1))
nF1 <- xtabs(~grpF)       # Number of type-1 fibers per subject
range(nF1)                # Min, max number of type-1 fibers
nF1[which.min(nF1)]       # Subject with the minimum number of fibers
str(fm17.1$dims)          # Basic dimensions used in the fit


###################################################
### code chunk: R17.2
###################################################
fixed1 <- summary(fm17.1)$tTable           # beta, se(beta), t-test
nms <- rownames(fixed1)                    # beta names
nms[7:8] <- c("fLow:fPos", "fMale:fOld")   # Selected names shortened 
rownames(fixed1) <- nms                    # New names assigned
printCoefmat(fixed1, digits = 3,           # See also Table 17.1 
             has.Pvalue = TRUE, P.values = TRUE)


###################################################
### code chunk: R17.3a
###################################################
getVarCov(fm17.1)             # D: (17.3)
VarCorr(fm17.1)


###################################################
### code chunk: R17.3b
###################################################
Ri <-                         # Ri is a list containing R_i ...
  getVarCov(fm17.1, c("5", "275"),# ... for subjects "5" and "275". 
            type = "conditional")
Ri$"275"                      # R_i for the subject "275" (17.2)
Ri.5   <- Ri$"5"              # R_i for the subject "5" ...
dim(Ri.5)                       # ... with large dimensions ...
(Ri.5d <- diag(Ri.5)[1:6])      # ... its first 6 diagonal elements.
sgma <- summary(fm17.1)$sigma # sigma          
sgma^2                        # sigma^2


###################################################
### code chunk: R17.4a
###################################################
dt5 <-                             # Data with 30 observations 
   subset(prt1,               
          select = c(id, occ.f),     # ... and 2 variables                    
          id == "5")                 # ... for the subject "5".
auxF1 <- function(elv) {
   idx <- 1:min(length(elv), 2)    # Up to two indices per vector 
   elv[idx]                          # ... returned.
}
(i.u5 <-                           # Selected indices printed
   unlist(
      tapply(rownames(dt5),          # ... for the subject "5"
             dt5$occ.f,              # ... by occ.f subgroups 
             FUN = auxF1)))
dt.u5  <- dt5[i.u5, ]              # Raw data for selected indices
(nms.u5 <-                         # Row names constructed
   paste(i.u5, dt.u5$occ.f, sep = ".")) 


###################################################
### code chunk: R17.4b
###################################################
Vi <-                                # Vi is a list containing ...  
  getVarCov(fm17.1, "5",               # ... matrix V_i for subject "5".
            type = "marginal")
Vi.5 <- Vi$"5"                       # Vi.5 is a V_i matrix: (17.4)
Vi.u5 <- Vi.5[i.u5, i.u5]            # A sub-matrix selected, ...
rownames(Vi.u5) <- nms.u5              # ... row/column names changed,
Vi.u5                                  # ... the sub-matrix printed.
cov2cor(Vi.u5)                       # Corr(V_i) 


###################################################
### code chunk: R17.5
###################################################
rnf <- ranef(fm17.1)         # b_i: (13.50)
(vrnf <- var(rnf))           # var(b_i). Compare to D in R17.13a.
plot(rnf)                    # Side-by-side plot (Fig. 17.1a)


library(ellipse)
myPanel <- function(x,y, ...){
  panel.grid(h = -1, v = -1)
  panel.xyplot(x, y)
  ex1 <-                     # Ellipse based on D: (17.3)
     ellipse(getVarCov(fm17.1)) 
  panel.xyplot(ex1[, 1], ex1[, 2], type = "l", lty = 1) 
  ex2 <- ellipse(vrnf)       # Ellipse based on var(b_i). 
  panel.xyplot(ex2[ ,1], ex2[, 2], type = "l", lty = 2)
}


xyplot(rnf[, 2] ~ rnf[, 1],  # Scatterplot b_i1 versus b_i0 (Fig. 17.1b) 
       xlab = "Pre-intervention", 
       ylab = "Post-intervention",
       xlim = c(-40, 40), ylim = c(-40, 40), 
       panel = myPanel)


###################################################
### code chunk: R17.6
###################################################
prt1r <-                                     # Auxiliary data
   within(prt1, 
          {                                  # Pearson residuals                       
            residP1 <- residuals(fm17.1, type = "p")  
            fitted1 <- fitted(fm17.1)
          })
range(prt1r$residP1)                    # Info for y-axis range
xyplot(residP1 ~ fitted1| occ.f,   # Resids vs fitted (Fig. 17.2a)
       data = prt1r, ylim = c(-6, 6), 
       type = c("p", "smooth"),
       grid = TRUE)
qqnorm(prt1r$residP1); qqline(prt1r$residP1) # Q-Q plot (Fig. 17.3a)



###################################################
### code chunk: R17.7
###################################################
fm17.2 <-                             #  M17.2 <- M17.1
   update(fm17.1,
          weights = varPower(form = ~fitted(.)),
          data = prt1) 
intervals(fm17.2)$varStruct           # 95% CI for delta, (17.5)
anova(fm17.1, fm17.2)                 # H0: delta = 0 (M17.1 nested in M17.2)



###################################################
### code chunk: R17.8a
###################################################
lme.spec.form3 <- 
   update(lme.spec.form1,                    # M17.3  <-  M17.1 
          . ~ . + fiber.f + prt.f:fiber.f + occ.f:fiber.f)
fm17.3 <- 
   lme(lme.spec.form3,                       # (17.6)
       random = ~occ.f:fiber.f - 1|id,       # D(17.7)
       data = prt) 



###################################################
### code chunk: R17.8b
###################################################
fixed.D4 <- summary(fm17.3)$tTable         # beta, se(beta), t-test  
rnms <- rownames(fixed.D4)                 # beta names (not shown)
rnms[8:11] <-                              # Selected names shortened
   c("Low:Pos", "Low:Type2", "Pos:Type2", "Male:Old") 
rownames(fixed.D4) <- rnms                 # Short names assigned
printCoefmat(fixed.D4, digits = 3, zap.ind = 5)  


###################################################
### code chunk: R17.9a
###################################################
fm17.3cov <-                                   # D: (17.7) extracted
   getVarCov(fm17.3, type = "random.effect") 
rownames(fm17.3cov)                            # Long names ...
nms. <- c("T1.Pre", "T1.Pos", "T2.Pre", "T2.Pos")# ... abbreviated
dimnames(fm17.3cov) <- list(nms., nms.)          # ... and reassigned.
fm17.3cov                                      # D: (17.7) printed
fm17.3cor <- cov2cor(fm17.3cov)                # Corr(D) ...
print(fm17.3cor, digits = 2,                     # ...printed.
      corr = TRUE, stdevs = FALSE)             


###################################################
### code chunk: R17.9b
###################################################
dim(R.5 <-                                     # Dims of R_i ...
   getVarCov(fm17.3,                    
             type = "conditional")[["5"]])        # ... for subject "5".
diag(R.5)[1:6]                            # First 6 diagonal elements
(sgma <- fm17.3$sigma)                         # sigma
print(sgma^2)                                  # sigma^2


###################################################
### code chunk: R17.10
###################################################
CI <- intervals(fm17.3, which = "var-cov") # 95% CIs for theta_D
interv <- CI$reStruct$id                   
# rownames(interv)                         # Long names (not shown)
thDnms  <- 
   c("sd(T1Pre)", "sd(T1Pos)", "sd(T2Pre)", "sd(T2Pos)",
     "cor(T1Pre,T1Pos)", "cor(T1Pre,T2Pre)", "cor(T1Pre,T2Pos)", 
                         "cor(T1Pos,T2Pre)", "cor(T1Pos,T2Pos)", 
                                             "cor(T2Pre,T2Pos)")   
rownames(interv) <- thDnms                 # Short names assigned
interv                                     # CIs printed

###################################################
### code chunk: R17.11
###################################################
residP3 <- residuals(fm17.3, type =  "p") # Pearson residuals 
xyplot(residP3 ~ fitted(fm17.3)|   # Scatterplots ...
       fiber.f:occ.f,              # ...per type*occasion (Fig. 17.4)
       data = prt,
       type = c("p", "smooth"))
qqnorm(residP3); qqline(residP3)   # Q-Q plot (Fig. 17.5)



###################################################
### code chunk: R17.12a
###################################################
Vx <-                                # Vx is a list ...
   getVarCov(fm17.3, type = "marginal",                    
             individual = "5")       # ... with one component.
Vmtx.5 <- Vx$"5"                     # Vmtx.5 is V_i matrix (17.8)...  
dim(Vmtx.5)                            # ... with large dimensions.
dt5 <-                               # Data with 41 rows ...
   subset(prt,                   
          select = c(id, fiber.f, occ.f), # ... and 3 variables ...
          id == "5")                   # ... for subject "5".


###################################################
### code chunk: R17.12b
###################################################
(i.u5  <- unlist(                    # Selected indices printed.
   tapply(rownames(dt5),             # Indices for subject "5" ...                       
          list(dt5$fiber.f, dt5$occ.f), # ... by fiber.f and occ.f.
          FUN = auxF1))) 
dt.u5  <- dt5[i.u5, ]                # Raw data for selected indices
nms.u5 <- 
   paste(format(i.u5, 2, justify = "right"),   
         abbreviate(dt.u5$fiber.f, 2),     # Row names abbreviated 
         dt.u5$occ.f, sep = ".") 


###################################################
### code chunk: R17.12c
###################################################
Vmtx.u5 <- Vmtx.5[i.u5, i.u5]      # Submatrix of V_i for subject "5"
dimnames(Vmtx.u5) <- list(nms.u5, i.u5) # dimnames assigned
Cmtx.u5 <- cov2cor(Vmtx.u5)        # Submatrix of Corr(V_i)
uptri <- upper.tri(Cmtx.u5)        # Logical matrix
Vmtx.u5[uptri] <- Cmtx.u5[uptri]      
print(Vmtx.u5, digits = 2)         # Submatrix printed


###################################################
### code chunk: R17.13a
###################################################
fm17.3a <- 
   lme(lme.spec.form3,                     # M17.3a
       random = ~1 + fiber.f + occ.f + fiber.f:occ.f|id,
       data = prt) 
print(fm17.3a$sigma, digits = 4)           # sigma
fm17.3acov <-                              # D
   getVarCov(fm17.3a,                
             type = "random.effect", individual = "5")
dimnames(fm17.3acov)[[1]]                  # Row/col D names ...
nms <- c("(Int)", "T2", "Pos", "T2:Pos")   # ... shortened
dimnames(fm17.3acov) <- list(nms,nms)      # ... and assigned.
print(fm17.3acov, digits = 4)              # D printed


###################################################
### code chunk: R17.13b
###################################################
td <-                                      # T_D: (17.12) created...
   matrix(c(1, 0, 0, 0,                   
            1, 0, 1, 0,                     
            1, 1, 0, 0,
            1, 1, 1, 1), 
          nrow = 4, ncol = 4, byrow = TRUE)
mat.D4 <- td %*% fm17.3acov %*% t(td)          # ... and applied.
dimnames(mat.D4) <- list(nms., nms.)       # Row/col names shortened.
print(mat.D4, digits = 5)                  # D:(17.7); see R17.9.



###################################################
### code chunk: fig 17.6 using splom() function
###################################################

D173 <- getVarCov(fm17.3)
D173a <- getVarCov(fm17.3a)
nms    <- c("T1:Pre/(Int)","T2:Pre/T2","T1:Pos/Pos","T2:Pos/T2:Pos","fitName")
dtref1 <- within(ranef(fm17.3),  fitName <- "fm17.3")
names(dtref1)
names(dtref1) <- nms
dtref2 <- within(ranef(fm17.3a), fitName <- "fm17.3a")
names(dtref2)
names(dtref2) <- nms
dtref  <- rbind(dtref1, dtref2)
(lims <- range(dtref[,1:4]))
lims <- c(-40,40)  # user-defined limits for every variable
atx <- -1


myFunL <- function(corx) ltext(-15, 25, substitute(paste(rho, corx), list(corx = corx)), cex = 0.9)

myFunU <- function(corx) ltext(-15, -32, substitute(paste(rho, corx), list(corx = corx)), cex = 0.9)

my.upperPanel <-  ## pairwise.complete.obs 
  function(x, y, subscripts, ...){
  vr <- dtref$fitName == "fm17.3a" 
  subs <- subscripts[vr]         
  x1 <- x[subs]
  y1 <- y[subs]
  panel.grid(h = atx, v = atx, col = "grey", ...)
  panel.xyplot(x1, y1, ...)
  corx <- round(cor(x1, y1, use = "complete.obs"), 2)
  abs.corx <- abs(corx)
  corx <- paste("=", corx, sep = "")
  myFunU(corx)
}

my.lowerPanel <-  ## pairwise.complete.obs 
  function(x, y, subscripts, ...){
  vr <- dtref$fitName == "fm17.3" 
  subs <- subscripts[vr]         
  x1 <- x[subs]
  y1 <- y[subs]
  panel.grid(h=atx,v=atx, col="grey") ##  ...lty="13",...)
  panel.xyplot(x1, y1, ...)
  corx <- round(cor(x1, y1, use = "complete.obs"), 2)
  abs.corx <- abs(corx)
  corx <- paste("=", corx, sep = "")
  print(corx)
  cex.value <- 2
  rx <- expression(paste(rho,corx))
  myFunL(corx)
}

mySuperPanel <- function(z, subscripts, panel.subscripts,...){
  panel.pairs(z, subscripts = subscripts,
              panel.subscripts = panel.subscripts,
              as.matrix = TRUE, 
              upper.panel = "my.upperPanel",
              lower.panel = "my.lowerPanel",
              ## Possible to shorten syntax. See other splom figures
              pscales =list(
                "T1:Pre/(Int)"  = list(limits = lims),
                "T2:Pre/T2"     = list(limits = lims),
                "T1:Pos/Pos"    = list(limits = lims),
                "T2:Pos/T2:Pos" = list(limits = lims)) )
print(names(z))
}

abbrev.names <- c("vis0", "vis4", "vis12", "vis24", "vis52")

splom.form <- formula(~as.matrix(dtref[,1:4])) 
splom(splom.form,
  data = dtref, 
  as.matrix = TRUE,  #### varnames = abbrev.names, 
  xlab = "",
  superpanel = mySuperPanel 
)

sessionInfo()
detach(package:nlme)

