
###################################################
### code chunk: Chap9init
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)
date()
packageVersion("nlmeU")
packageVersion("nlme")
packageVersion("ellipse")
packageVersion("lattice")
sessionInfo()
data(armd, package = "nlmeU")

library(nlme)
lm1.form <-                             # See also R6.1
   formula(visual ~ -1 + visual0 + time.f + treat.f:time.f)
fm6.1 <- gls(lm1.form, data = armd)

###################################################
### code chunk: R9.1a
###################################################
lm1.form <-                             # See also R6.1
   formula(visual ~ -1 + visual0 + time.f + treat.f:time.f)
fm9.1 <-                                # M9.1
   gls(lm1.form,                  
       weights = varIdent(form = ~1|time.f), # Var. function; <delta>-group 
       data = armd)
summary(fm9.1)                  
fm9.1$modelStruct$varStruct             # (9.3): delta1=1, ...  
(intervals(fm9.1, which = "var-cov"))   # 95% CI for delta2,3,4 & sigma

###################################################
### code chunk: R9.1b
###################################################
anova(fm6.1, fm9.1)                     # M6.1 nested in M9.1


###################################################
### code chunk: R9.2a
###################################################
fm9.2 <-                                    # M9.2 <- M9.1
   update(fm9.1,                    
          weights = varPower(form = ~time)) # (9.4); <delta>-group 
fm9.3 <-                                    # M9.3 <- M9.1
   update(fm9.1,                            # (9.5); strata=treat.f 
          weights = varPower(form = ~time|treat.f))
fm9.4 <-                                    # M9.4 <- M9.1
  update(fm9.1, weights = varPower())       # (9.6); <delta,mu>-group
fm9.5 <-                                    # M9.5 <- M9.1
  update(fm9.1,                             
         weights = varPower(fixed = 1))     # (9.7), <mu>-group


###################################################
### code chunk: R9.2b
###################################################
anova(fm9.2, fm9.3)                       # M9.2 nested in M9.3


###################################################
### code chunk: R9.2c
###################################################
anova(fm9.2, fm9.1)                       # M9.2 nested in M9.1


###################################################
### code chunk: R9.2d
###################################################
anova(fm9.5, fm9.4)                       # M9.5 nested in M9.4 

###################################################
### code chunk: R9.2e
###################################################
AIC(fm9.1, fm9.2, fm9.3,                  # Non-nested models
    fm9.4, fm9.5)                         # Smaller AIC corresponds to a better fit


###################################################
### code chunk: R9.3a
###################################################
mSt2 <- fm9.2$modelStruct           # Model structure
vF2 <- mSt2$varStruct               # Variance function:(9.4)
summary(vF2)                        # Summary: delta.
summary(fm9.2)$sigma                # sigma



###################################################
### code chunk: R9.3b
###################################################
mSt3 <- fm9.3$modelStruct           # Model structure
vF3  <- mSt3$varStruct              # Variance function:(9.5)
summary(vF3)                        # Summary: delta1,2
coef(vF3)                           # delta1,2
formula(vF3)                        # Variance function formula
varWeights(vF3)[3:10]               # Weights for two subjects

###################################################
### code chunk: R9.4a
###################################################
library(lattice)
plot(fm9.2,                                   # Fig. 9.1a
     resid(., type = "response") ~ fitted(.)) # Raw vs fitted  
plot(fm9.2,                                   # Raw vs time (not shown)
     resid(., type = "response") ~ time)      # (See Fig. 9.1a)     
bwplot(resid(fm9.2) ~ time.f,                 # Fig. 9.1b         
       pch = "|", data = armd)                # Raw vs time.f.                      

###################################################
### code chunk: R9.4b
###################################################

plot(fm9.2,                                   # Fig. 9.1c
     resid(., type = "pearson" ) ~ fitted(.)) # Pearson vs fitted 
plot(fm9.2,                                   #  vs time (not shown)
     resid(., type = "pearson") ~ time)       # (See Fig. 9.1c)
bwplot(                                       # Fig. Fig. 9.1d
  resid(fm9.2, type = "pearson") ~ time.f,    # Pearson vs time.f                          
  pch = "|", data = armd)                 


###################################################
### code chunk: R9.4c
###################################################
plot(fm9.2,                               # Fig. 9.2a
     sqrt(abs(resid(., type = "response"))) ~ fitted(.),
     type = c("p", "smooth"))
plot(fm9.2,                               # Fig. 9.2b
     sqrt(abs(resid(., type = "pearson"))) ~ fitted(.),
     type = c("p", "smooth"))

###################################################
### code chunk: More elaborate syntax for Figure 9.2
###################################################

ylab <- c("Residuals", "Standardized residuals")
xlab <- c("Fitted values", "Time in weeks")
ylim1 <- c(-50,50)                          # For Figs. 9.1a and 9.1b
ylim2 <- c(-4.2, 4.2)                       # For Figs. 9.1a and 9.1b
xlims <- c("4", "12","24","52wks")          # For Figs. 9.1b and 9.1d

xyPanel <- function(x,y,subscripts,...){   # For Figs. 9.1a and 9.1c
  panel.grid(h = -1, v = -1)  # vertical and horizontal
  panel.xyplot(x, y, ...)     # points over the grid
}

bwPanel <- function(x,y,subscripts,...){   # For Figs. 9.1b and 9.1d
  panel.grid(h = -1, v = 0)
  panel.bwplot(x, y, ...)         # bwplot over the grid
}

### Fig. 9.1a. Raw vs fitted
plot(fm9.2,              
  resid(., type = "response") ~ fitted(.), panel = xyPanel,
  ylim = ylim1, xlab = xlab[1], ylab = ylab[1])  

### Fig. 9.1b Raw vs time 
bwp <- bwplot(                        # Fig. 9.1b         
  resid(fm9.2, type = "response") ~ time.f,              # Raw  vs time.f 
  data = armd, 
  panel= bwPanel,  
  pch = "|"
)                 
update(bwp, xlim = xlims, xlab=xlab[2], ylab=ylab[1], ylim = ylim1)


### Fig. 9.1c. Pearson vs fitted
plot(fm9.2, 
  resid(., type = "pearson") ~ fitted(.), 
  panel = xyPanel,
  ylim = ylim2, xlab = xlab[1], ylab = ylab[2])  

### Fig. 9.1d Pearson vs time 
bwp <- bwplot(                        # Fig. 9.1b         
  resid(fm9.2, type = "pearson") ~ time.f,              # Raw  vs time.f 
  data = armd, 
  panel= bwPanel,  
  pch = "|"
)                 
update(bwp, xlim = xlims, xlab=xlab[2], ylab=ylab[2], ylim = ylim2)

###################################################
### code chunk: Syntax for Fig. 9.2a
###################################################

myPanel <- function(x,y,...){
  panel.grid(h = -1, v = -1)
  panel.xyplot(x, y, ...)
}

plot(fm9.2,
  sqrt(abs(resid(., type = "response"))) ~ fitted(.),
  type = c("p", "smooth"), 
  panel=myPanel,
  xlab = xlab[1],
  ylab = "Raw residuals"
)

###################################################
### code chunk: Syntax for Fig. 9.2b
###################################################

plot(fm9.2,
  sqrt(abs(resid(., type = "pearson"))) ~ fitted(.),
  type = c("p", "smooth"), 
  panel=myPanel,
  xlab = xlab[1],
  ylab = "Pearson residuals"
)


###################################################
### code chunk: Scatterplot matrix: Data preparation
###################################################
residP <- resid(fm9.2, type="p")
dtAux <- subset(armd, select = c(subject, visual, time, time.f, treat.f)) 
require(reshape)
dtP <- data.frame(dtAux, residP)
dtPm <- melt(dtP,
     measure.var=c("residP"),
     id.var = c("subject","time.f"))
dtPc <- cast(subject ~ time.f ~ variable, data = dtPm) # array
dtPc <- data.frame(dtPc) 
names(dtPc) <- c("P4wks","P12wks","P24wks","P52wks")
head(dtPc)
range(dtPc, na.rm=TRUE)

###################################################
### code chunk: Syntax to create Fig. 9.3 using splom() function
###################################################

library(ellipse)
my.upperPanel <-                           ## pairwise.complete.obs 
  function(x, y, subscripts, ...){
  panel.xyplot(x, y, type = "n", ...)      # no plot
  ic <- complete.cases(cbind(x,y))
  mn <- c(mean(x[ic]), mean(y[ic]))
  covx <- var(cbind(x,y), use="complete.obs")
  # print(covx)
  # ex <- ellipse(covx)
  corrx <- cov2cor(covx)
  corx <- round(corrx[1,2], 2)
  abs.corx <- abs(corx)
  # print(corx)
  cex.value <- 3
  ltext(0, 0, corx, cex = abs.corx * cex.value)
}

my.lowerPanel <-                          ## pairwise.complete.obs 
  function(x,y,subscripts,...){
  panel.grid(h = -1, v = -1)
  covx <- var(cbind(x, y), use = "complete.obs")
  # print(covx)
  ex <- ellipse(covx)
  panel.xyplot(ex[ ,1], ex[ ,2], lty = 2, type = "l", ...)
  panel.xyplot(x, y, ...)
}


mySuperPanel <- function(z, subscripts, panel.subscripts,...){
  panel.pairs(z, subscripts = subscripts,
              panel.subscripts = panel.subscripts,
              as.matrix=TRUE, 
              upper.panel = "my.upperPanel",
              lower.panel = "my.lowerPanel",
              prepanel.limits = function(z) return(c(-4,4))
)
}


splom.form <- formula(~cbind(P4wks,P12wks,P24wks,P52wks))
splom.object <- splom(splom.form,
  data=dtPc,             #### subset(armd240,miss.pat =="----"),
  as.matrix=TRUE,  #### varnames = abbrev.names, 
  xlab="",
  superpanel = mySuperPanel
)
### fig. 9.3
print(splom.object)
rm(my.upperPanel,mySuperPanel,splom.object)

#### sessionInfo with packages attached

sessionInfo()

detach(package:nlme)
detach(package:lattice)
detach(package:ellipse)



