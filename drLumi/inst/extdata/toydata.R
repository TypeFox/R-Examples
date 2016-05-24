
####################################
### Data, functions and plots
### for the drLumi paper
### date: 2014 - 01 - 07
####################################

### Load packages and ToyData
rm(list=ls())
library(gridExtra)
library(ggplot2)
library(drLumi)



######################
#### Generate data ###
######################

# 4 parameter logistic regression 
# based on VEGF analyte from plate_01. 
# l5p[[1]]$VEGF$model

### Test 1
############
par.ll4 <- c(-0.82485, 1.57269, 3.95003, 1.68015, 1)
sde <- 0.01538

dil <- c(7600,7600,3800,1900,950,475,237.5,118.75,59.38,29.69,14.84,7.42,3.71,1.86,0.93,0.46,0.23)

log10dil <- log10(dil)
mm.ll4 <- par.ll4[2] + ((par.ll4[3]-par.ll4[2])/ (1 + 10^(par.ll4[1]*(log10dil-par.ll4[4])))^par.ll4[5])

set.seed(1234)
std.out <- rnorm(length(log10dil), mm.ll4 , sd = sde) 
test1 <- as.data.frame(cbind(10^std.out, 10^log10dil))
names(test1) <- c("median","ec")
test1$plate <- "plate_ll4"
test1$well <- c(paste0("A",1:12),paste0("B",1:5))
test1$analyte <- "test1"
test1$sample <- paste0("standard",1:17)

# Test2: missing
##################
test2 <-   test1
test2[c(1,2,4,6,9,16,17),1] <- NA
test2$plate <- "plate_ll4"
test2$well <- c(paste0("A",1:12),paste0("B",1:5))
test2$analyte <- "test2"
test2$sample <- paste0("standard",1:17)

# Test4: outliers
##################
set.seed(1234)
std.out <- rnorm(length(log10dil), mm.ll4 , sd = c(sde,rep(sde,3), sde*10^2, rep(sde,3), sde*10^2, rep(sde,6)))
test4 <- as.data.frame(cbind(10^std.out, 10^log10dil))
names(test4) <- c("median","ec")
test4$plate <- "plate_ll4"
test4$well <- c(paste0("A",1:12),paste0("B",1:5))
test4$analyte <- "test4"
test4$sample <- paste0("standard",1:17)


# Test3: outliers and missing
###############################
test3 <- test4
test3[c(1,2,4,6,9,16,17),1] <- NA
test3$plate <- "plate_ll4"
test3$well <- c(paste0("A",1:12),paste0("B",1:5))
test3$analyte <- "test3"
test3$sample <- paste0("standard",1:17)

# All data
dat <- rbind(test1, test2, test3, test4)


# Background data 
######################
analytes <- c("test1","test2","test3","test4") 
values <- c(40,45)
all.background <- expand.grid(analytes, values)
all.background$well <- c(rep("C6",4), rep("C7",4))
names(all.background) <- c("analyte", "median", "well")
all.background$ec <- 0




#####################
###### Functions ####
#####################

### No flag data
## Fit functions
mod.ig <- scluminex("toydata", dat, all.background, lfct="SSl4", bkg="ignore")
mod.in <- scluminex("toydata", dat, all.background, lfct="SSl4", bkg="include")
mod.sub <- scluminex("toydata", dat, all.background, lfct="SSl4", bkg="subtract")
mod.cons <- scluminex("toydata", dat, all.background, lfct="SSl4", bkg="constraint")


## Plot: standard curve
q1 <- plot(mod.ig, size.legend=4) + ggtitle("Ignore")
q2 <- plot(mod.in, size.legend=4) + ggtitle("Include")
q3 <- plot(mod.sub, size.legend=4) + ggtitle("Subtract")
q4 <- plot(mod.cons, size.legend=4) + ggtitle("Constraint")
grid.arrange(q1,q2,q3,q4,ncol=2)

## Plot: residuals
q1 <- plot(mod.ig, "res", size.text=5, out.limit=2) + ggtitle("Ignore")
q2 <- plot(mod.in, "res", size.text=5, out.limit=2) + ggtitle("Include")
q3 <- plot(mod.sub, "res", size.text=5, out.limit=2) + ggtitle("Subtract")
q4 <- plot(mod.cons, "res", size.text=5, out.limit=2) + ggtitle("Constraint")
grid.arrange(q1,q2,q3,q4,ncol=2)

## Identify outliers
out.ig <- get_outliers(mod.ig, out.limit=2)
out.in <- get_outliers(mod.in, out.limit=2)
out.sub <- get_outliers(mod.sub, out.limit=2)
out.cons <- get_outliers(mod.cons, out.limit=2)


### Flag data
###############
## Flag outliers
flag.ig <- merge(dat, out.ig, by=c("well","analyte"), all.x=TRUE)
flag.in <- merge(dat, out.in, by=c("well","analyte"), all.x=TRUE)
flag.sub <- merge(dat, out.sub, by=c("well","analyte"), all.x=TRUE)
flag.cons <- merge(dat, out.cons, by=c("well","analyte"), all.x=TRUE)

# Fit models with flag data
mod.ig.flag <- scluminex("toydata_flag", flag.ig, all.background, lfct="SSl4",bkg="ignore")
mod.in.flag <- scluminex("toydata_flag", flag.in, all.background, lfct="SSl4",bkg="include")
mod.sub.flag <- scluminex("toydata_flag", flag.sub, all.background, lfct="SSl4", bkg="subtract")
mod.cons.flag <- scluminex("toydata_flag", flag.cons, all.background, lfct="SSl4", bkg="constraint")


## Plot: standard curve
q1 <- plot(mod.ig.flag, size.legend=4) + ggtitle("Ignore flagged")
q2 <- plot(mod.in.flag, size.legend=4) + ggtitle("Include flagged")
q3 <- plot(mod.sub.flag, size.legend=4) + ggtitle("Subtract flagged")
q4 <- plot(mod.cons.flag, size.legend=4) + ggtitle("Constraint flagged")
grid.arrange(q1,q2,q3,q4,ncol=2)

## Plot: residuals
q1 <- plot(mod.ig.flag, "res", size.text=5, out.limit=2) + ggtitle("Ignore flagged")
q2 <- plot(mod.in.flag, "res", size.text=5, out.limit=2) + ggtitle("Include flagged")
q3 <- plot(mod.sub.flag, "res", size.text=5, out.limit=2) + ggtitle("Subtract flagged")
q4 <- plot(mod.cons.flag, "res", size.text=5, out.limit=2) + ggtitle("Constraint flagged")
grid.arrange(q1,q2,q3,q4,ncol=2)

#################################
## Comparison: No flag vs. Flag
#################################
## Summary of Data (I)
summary(mod.ig)
summary(mod.ig.flag)

summary(mod.in)
summary(mod.in.flag)

summary(mod.sub)
summary(mod.sub.flag)

summary(mod.cons)
summary(mod.cons.flag)

## Summary of Data (II)
as.data.frame(summary(mod.ig))
as.data.frame(summary(mod.ig.flag))

as.data.frame(summary(mod.in))
as.data.frame(summary(mod.in.flag))

as.data.frame(summary(mod.sub))
as.data.frame(summary(mod.sub.flag))

as.data.frame(summary(mod.cons))
as.data.frame(summary(mod.cons.flag))

#########################################################
#### Limits of Quantification for flag and no flag data
#########################################################

# Ignore
der <- summary(loq_derivatives(mod.ig))
der.flag <- summary(loq_derivatives(mod.ig.flag))
int <- summary(loq_interval(mod.ig, high.asymp=3, low.asymp=2))
int.flag <- summary(loq_interval(mod.ig.flag, high.asymp=3 , low.asymp=2))
cv <- summary(loq_cv(mod.ig, max.cv=0.3))
cv.flag <- summary(loq_cv(mod.ig.flag, max.cv=0.3))

# Subtract
der <- summary(loq_derivatives(mod.sub))
der.flag <- summary(loq_derivatives(mod.sub.flag))
int <- summary(loq_interval(mod.sub, high.asymp=3, low.asymp=2))
int.flag <- summary(loq_interval(mod.sub.flag, high.asymp=3 , low.asymp=2))
cv <- summary(loq_cv(mod.sub, max.cv=0.3))
cv.flag <- summary(loq_cv(mod.sub.flag, max.cv=0.3))

# Constraint
der <- summary(loq_derivatives(mod.cons))
der.flag <- summary(loq_derivatives(mod.cons.flag))
int <- summary(loq_interval(mod.cons, high.asymp=3, low.asymp=2))
int.flag <- summary(loq_interval(mod.cons.flag, high.asymp=3 , low.asymp=2))
cv <- summary(loq_cv(mod.cons, max.cv=0.3))
cv.flag <- summary(loq_cv(mod.cons.flag, max.cv=0.3))

# Include
der <- summary(loq_derivatives(mod.in))
der.flag <- summary(loq_derivatives(mod.in.flag))
int <- summary(loq_interval(mod.in, high.asymp=3, low.asymp=2))
int.flag <- summary(loq_interval(mod.in.flag, high.asymp=3 , low.asymp=2))
cv <- summary(loq_cv(mod.in, max.cv=0.3))
cv.flag <- summary(loq_cv(mod.in.flag, max.cv=0.3))

#################################
### Outliers: no flag vs. flag
#################################

# No flag
ignore <- get_outliers(mod.ig, out.limit=2)
ignore$bkg_method <- "ignore"
subtract <- get_outliers(mod.sub, out.limit=2)
subtract$bkg_method <- "subtract"
include <- get_outliers(mod.in, out.limit=2)
include$bkg_method <- "include"
constraint <- get_outliers(mod.cons, out.limit=2)
constraint$bkg_method <- "constraint"

noflag <- rbind(ignore, subtract, include,constraint)

# Flag
ignore <- get_outliers(mod.ig.flag, out.limit=2)
ignore$bkg_method <- "ignore"
subtract <- get_outliers(mod.sub.flag, out.limit=2)
subtract$bkg_method <- "subtract"
include <- get_outliers(mod.in.flag, out.limit=2)
include$bkg_method <- "include"
constraint <- get_outliers(mod.cons.flag, out.limit=2)
constraint$bkg_method <- "constraint"

flag <- rbind(ignore, subtract, include,constraint)

noflag
flag




