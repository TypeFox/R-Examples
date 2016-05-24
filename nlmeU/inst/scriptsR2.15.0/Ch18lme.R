###################################################
### code chunk: Chap18init
###################################################
options(width = 65, digits = 5, show.signif.stars = FALSE)
date()
packageVersion("nlmeU")
packageVersion("nlme")
packageVersion("lattice")
packageVersion("splines")
require(nlme)    
require(lattice) 
sessionInfo()





###################################################
### code chunk: R18.1
###################################################
data(SIIdata, package = "nlmeU")
form1 <- 
   formula(mathgain ~ ses + minority  +    # (18.1)
             mathkind + sex + housepov)       
(fm18.1 <- 
   lme(form1, 
       random = list(schoolid = ~1,        # See Table 14.1, syntax (a)
                     classid = ~1), 
       data = SIIdata,  method = "ML"))
       
update(fm18.1,                             # An alternative syntax
       random = ~1 | schoolid/classid)     # See Table 14.1, syntax (d)



###################################################
### code chunk:  R18.2
###################################################
getGroupsFormula(fm18.1)                 # Grouping formula
str(grpF1 <- getGroups(fm18.1, level=1)) # Grouping factor at level 1 of the data hierarchy
str(grpF2 <- getGroups(fm18.1))          # Grouping factor at level 2 of the data hierarchy
grpF2


###################################################
### code chunk: R18.3a
###################################################
(fxd <- fixef(fm18.1))                   # beta 


###################################################
### code chunk: R18.3b
###################################################
vcov1 <- vcov(fm18.1)                   # vcov (beta)
nms <- abbreviate(names(fxd))           # Abbreviated beta names ...
dimnames(vcov1) <- list(nms, nms)       # ... assigned.
print(vcov1, digits = 2)


###################################################
### code chunk: R18.4
###################################################
## getVarCov(fm18.1)                    # Error: Not implemented
VarCorr(fm18.1)




###################################################
### code chunk: R18.5a
###################################################
anova(fm18.1, type = "marginal")


###################################################
### code chunk: R18.5b
###################################################
anova(fm18.1, Terms = c("housepov"))
anova(fm18.1, Terms = c("sex"))

###################################################
### code chunk: R18.5c
###################################################

## anova(fm18.1, Terms = c("housepov", "sex"))  # Error: Terms must all have the same dnominator DF 


###################################################
### code chunk: R18.6a
###################################################
rsd1 <-                               # Marginal residuals 
    resid(fm18.1, level = 0)     
range(rsd1)                           # Range
outi <- abs(rsd1) > 120               # Selection of outliers
as.numeric(SIIdata$childid[outi])     # Outliers' ids
rsd1[outi]                            # Outliers' values and labels

###################################################
### code chunk: R18.6b
###################################################
myPanel <- function(x,y, subscripts, ...){
  panel.xyplot(x,y,...)
  outi <- abs(y) > 120
  y1   <- y[outi]
  x1   <- x[outi]  
  ltext(x1, y1, names(y1), pos=3)
}

xyplot(rsd1 ~ housepov|sex, SIIdata,  # Fig. 18.1 
       type = c("p","r"),
       panel = myPanel)




###################################################
### code chunk: R18.7a
###################################################
form2 <- update(form1, . ~ . + sex:housepov)   # (18.4)
fm18.2 <- update(fm18.1, form2)                # M18.2 <- M18.1  
summary(fm18.2)                                # Summary




###################################################
### code chunk: R18.7b
###################################################
anova(fm18.2, Terms = "sex:housepov")            # Approximate F-test
anova(fm18.1, fm18.2)                            # M18.1 nested in M18.2 


###################################################
### code chunk: R18.8
###################################################
form3 <- update(form1, . ~ . - sex - housepov) # (18.5)
fm18.3 <- update(fm18.1, form3)                # M18.3 <- M18.1
anova(fm18.1, fm18.3, fm18.2)                  # M18.3 nested in M18.1 nested in M18.2


### Syntax for Fig. 18.2  #################
rsd3 <-                                        # Marginal residuals. Syntax similar to that in Panel R18.6 
    resid(fm18.3, level = 0)      
xyplot(rsd3 ~ mathkind, SIIdata,               # Fig. 18.2a
    type = c("p", "smooth"))
xyplot(rsd3 ~ ses, SIIdata,                    # Fig. 18.2b
    type = c("p", "smooth"))



###################################################
### code chunk: R18.9
###################################################
form4 <-                                         # (18.6)
   formula(mathgain ~ ses + minority + poly(mathkind, 3)) 
fm18.4 <- update(fm18.3, form4)                  # M18.4 <- M18.3
anova(fm18.3, fm18.4)                            # M18.3 nested in M18.4


###################################################
### code chunk: R18.10
###################################################
auxL <-                                       # Auxiliary list
   list(ses = 0,              
        minority = factor(c("Mnrt=No", "Mnrt=Yes")),
        mathkind = seq(290, 625, by = 5))
dim (auxDt <-  expand.grid(auxL))             # Data frame created
names(auxDt)
prd   <- predict(fm18.4, auxDt, level = 0)    # Predicted values 
prd4Dt <- data.frame(auxDt, pred4 = prd)
head(prd4Dt)
xyplot (pred4 ~ mathkind, groups = minority,  # Fig. 18.3a
        data = prd4Dt, type = "l", grid = TRUE)


###################################################
### Code for Fig. 18.3a
###################################################
xyplot (pred4 ~ mathkind, groups = minority, data = prd4Dt,  
        type = "l",  grid = TRUE, 
        key = list(
              lines = list(lty = c(1,2)),
              text = list(c("Mnrt=No", "Mnrt=Yes")),
              columns = 2, cex = 0.9), ylim = c(-60,220)
)

###################################################
### Code for Fig. 18.3b
###################################################

xyplot (pred4 ~ ses, groups = minority, data = prd4Dt,  
        type = "l",  grid = TRUE, 
        key = list(
              lines = list(lty = c(1,2)),
              text = list(c("Mnrt=No", "Mnrt=Yes")),
              columns = 2, cex= 0.9), ylim = c(-60,220)
)

### Code for Fig. 18.4  #################
rsd4 <-                                        # Marginal residuals. Syntax similar to that in Panel R18.6  
    resid(fm18.4, level = 0)      
xyplot(rsd4 ~ mathkind, SIIdata,               # Fig. 18.4a
    type = c("p", "smooth"))

xyplot(rsd4 ~ ses, SIIdata,                    # Fig. 18.4b
    type = c("p", "smooth"))



###################################################
### code chunk: R18.11
###################################################
options(digits = 7) 
require(splines)                 
form5 <-                                        # (18.7)
   formula(mathgain ~ ses + minority + bs(mathkind, df = 4))
fm18.5 <- update(fm18.4, form5)                 # M18.5 <- M18.4
AIC(fm18.3, fm18.4, fm18.5)
detach(package:splines)
options(digits = 5) # Going back to 5


###################################################
### code chunk: R18.12
###################################################
form6 <-                                        # (18.8)            
   formula(mathgain ~ ses + minority + poly(mathkind, 3) +
             ses:minority) 
fm18.6 <- update(fm18.4, form6)                 # M18.6 <- M18.4
anova(fm18.4, fm18.6)                           # M18.4 nested in M18.6


###################################################
### code chunk: R18.13a
###################################################
rsd6 <- resid(fm18.6, level = 0) 
xyplot(rsd6 ~ ses | minority, SIIdata,
       type = c("p", "smooth"))              # Fig. 18.5

###################################################
### code chunk: R18.13b
###################################################
qqnorm(fm18.6)                          # Fig. 18.6a 
qqnorm(fm18.6,                          # Equivalent call
       form = ~resid(., type = "p", level = 2)) 
qqnorm(fm18.6,                          # Fig. 18.6b 
       form = ~resid(., type = "p")     # Residuals... 
                | sex*minority,           # ... by sex and minority.
       id = 0.0005)                     # Outliers identified.

###################################################
### code chunk: R18.13c
###################################################

qqnorm(fm18.6,                          # Plot not shown
       form = ~resid(., type = "p", 
                     level = 1))        # School level


###################################################
### code chunk: R18.14a
###################################################
rsd6p <- resid(fm18.6, type = "p")
keep <- abs(rsd6p) < 3
rsd6x <- rsd6p[keep]
rsd6p[!keep]


###################################################
### code chunk: R18.14b (Fig. 18.7a, outlying residuals omitted)
###################################################
qqDtx <- qqnorm(rsd6x, plot.it = FALSE)

xp1 <-  xyplot(x ~ y, data.frame(qqDtx))         # Draft plot 
update(xp1,                                      # Plot updated 
       ylab = "Quantiles of standard normal", 
       xlab = "Standardized residuals", 
       grid = TRUE)


###################################################
### code chunk: R18.14c (Fig. 18.7b, outlying residuals omitted)
###################################################
qqDtx2 <- cbind(SIIdata[keep, ], qqDtx) 
xp2 <-                              # See R18.14b how to update xp2
   xyplot(x ~ y | sex*minority, data = data.frame(qqDtx2)) 
update(xp2,                                      # Plot updated 
       ylab = NULL, 
       xlab = "Standardized residuals", 
       grid = TRUE)


###################################################
### code chunk: R18.15a  (Fig. 18.8a)
###################################################
ref6 <- ranef(fm18.6)                  # Random effects for classes.
mode(ref6)                             # A list ...
length(ref6)                           # ... with two components.

pref6 <- plot(ref6)          # Default plot for classes; not legible.
pref6lims <- pref6$y.limits            # Y-labels extracted
len  <- length(pref6lims)              # No. of labels
sel  <- seq(1, len, by = 15)           # Select every 15-th label.
pref6lims[-sel] = ""                   # Other labels set to blank.
update(pref6, ylim = pref6lims,        # Assign new Y-labels.
       ylab = "classid %in% schoolid") # Y-axis label


###################################################
### code chunk: R18.15b (Fig. 18.8b)
###################################################
ref61 <- ranef(fm18.6, level = 1)      # Random effects for schools.
plot(ref61)                            # Plot the random effects.

###################################################
### code chunk: R18.16 (Fig. 18.9)
###################################################

qqnorm(fm18.6, ~ranef(., level = 2), # Random effects for classes
       id = 0.2,                     # Fig. 18.9a
       xlab = "Random effects for classes")
qqnorm(fm18.6, ~ranef(., level=1),   # Random effects for schools
       id = 0.2,                     # Fig. 18.9b
       xlab = "Random effects for schools")


#### sessionInfo  ####

sessionInfo()                         # with packages attached    
detach(package:nlme)   
detach(package:lattice)


