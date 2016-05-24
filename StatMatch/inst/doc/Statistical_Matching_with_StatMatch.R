### R code from vignette source 'Statistical_Matching_with_StatMatch.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Statistical_Matching_with_StatMatch.Rnw:104-108
###################################################
options(useFancyQuotes="UTF-8")
#options(useFancyQuotes=FALSE)
options(width=66)
options(warn=-1)


###################################################
### code chunk number 2: Statistical_Matching_with_StatMatch.Rnw:112-118
###################################################
library(StatMatch) #loads pkg StatMatch
data(samp.A) # sample A in SM examples
str(samp.A)

data(samp.B) # sample B in the SM examples
str(samp.B)


###################################################
### code chunk number 3: Statistical_Matching_with_StatMatch.Rnw:123-128
###################################################
X.vars <- intersect(names(samp.A), names(samp.B))
X.vars

setdiff(names(samp.A), names(samp.B)) # available just in A
setdiff(names(samp.B), names(samp.A)) # available just in B


###################################################
### code chunk number 4: Statistical_Matching_with_StatMatch.Rnw:147-150
###################################################
require(Hmisc)
spearman2(n.income~area5+urb+hsize+age+sex+marital+edu7, 
          p=2, data=samp.A)


###################################################
### code chunk number 5: Statistical_Matching_with_StatMatch.Rnw:157-158
###################################################
pw.assoc(labour5~area5+urb+hsize5+c.age+sex+marital+edu7, data=samp.B)


###################################################
### code chunk number 6: Statistical_Matching_with_StatMatch.Rnw:198-208
###################################################
# choiche of the matching variables based on uncertainty
xx <- xtabs(~c.age+sex+marital+edu7, data=samp.A)
xy <- xtabs(~c.age+sex+marital+edu7+c.neti, data=samp.A)
xz <- xtabs(~c.age+sex+marital+edu7+labour5, data=samp.B)

out.fbw <-  Fbwidths.by.x(tab.x=xx, tab.xy=xy, tab.xz=xz)

# sort output according to average width
sort.av <- out.fbw$sum.unc[order(out.fbw$sum.unc$av.width),]
head(sort.av) # best 6 combinations of the Xs


###################################################
### code chunk number 7: Statistical_Matching_with_StatMatch.Rnw:228-232
###################################################
group.v <- c("area5","sex")
X.mtc <- "age" 
out.nnd <- NND.hotdeck(data.rec=samp.A, data.don=samp.B,
                       match.vars=X.mtc, don.class=group.v)


###################################################
### code chunk number 8: Statistical_Matching_with_StatMatch.Rnw:237-239
###################################################
summary(out.nnd$dist.rd) # summary distances rec-don
summary(out.nnd$noad) # summary available donors at min. dist.


###################################################
### code chunk number 9: Statistical_Matching_with_StatMatch.Rnw:244-250
###################################################
head(out.nnd$mtc.ids)
fA.nnd <- create.fused(data.rec=samp.A, data.don=samp.B,
                       mtc.ids=out.nnd$mtc.ids,
                       z.vars="labour5")

head(fA.nnd) #first 6 obs.


###################################################
### code chunk number 10: Statistical_Matching_with_StatMatch.Rnw:257-266
###################################################
group.v <- c("sex","area5")
X.mtc <- "age"
out.nnd.c <- NND.hotdeck(data.rec=samp.A, data.don=samp.B, 
                         match.vars=X.mtc, don.class=group.v, 
                         dist.fun="Manhattan", constrained=TRUE, 
                         constr.alg="Hungarian")
fA.nnd.c <- create.fused(data.rec=samp.A, data.don=samp.B,
                         mtc.ids=out.nnd.c$mtc.ids,
                         z.vars="labour5")


###################################################
### code chunk number 11: Statistical_Matching_with_StatMatch.Rnw:271-274
###################################################
#comparing distances
sum(out.nnd$dist.rd) # unconstrained
sum(out.nnd.c$dist.rd) # constrained


###################################################
### code chunk number 12: Statistical_Matching_with_StatMatch.Rnw:279-289
###################################################
# estimating marginal distribution of labour5
tt0 <- xtabs(~labour5, data=samp.B) # reference distr.
tt <- xtabs(~labour5, data=fA.nnd)  # synt unconstr.
ttc <- xtabs(~labour5, data=fA.nnd.c) #synt. constr.
#
# comparing marginal distributions
cp1 <- comp.prop(p1=tt, p2=tt0, n1=nrow(fA.nnd), n2=NULL, ref=TRUE)
cp2 <- comp.prop(p1=ttc, p2=tt0, n1=nrow(fA.nnd), n2=NULL, ref=TRUE)
cp1$meas
cp2$meas


###################################################
### code chunk number 13: Statistical_Matching_with_StatMatch.Rnw:299-306
###################################################
# random hot deck in classes formed crossing "area5" and "sex"
group.v <- c("area5","sex")
rnd.1 <- RANDwNND.hotdeck(data.rec=samp.A, data.don=samp.B, 
                          match.vars=NULL, don.class=group.v)
fA.rnd <- create.fused(data.rec=samp.A, data.don=samp.B,
                       mtc.ids=rnd.1$mtc.ids, 
                       z.vars="labour5")


###################################################
### code chunk number 14: Statistical_Matching_with_StatMatch.Rnw:313-324
###################################################
# random choice of a donor among the closest k=20 wrt age
# sharing the same values of "area5" and "sex"
group.v <- c("area5","sex")
X.mtc <- "age"
rnd.2 <- RANDwNND.hotdeck(data.rec=samp.A, data.don=samp.B, 
                          match.vars=X.mtc, don.class=group.v, 
                          dist.fun="Manhattan", 
                          cut.don="exact", k=20)
fA.knnd <- create.fused(data.rec=samp.A, data.don=samp.B,
                        mtc.ids=rnd.2$mtc.ids, 
                        z.vars="labour5")


###################################################
### code chunk number 15: Statistical_Matching_with_StatMatch.Rnw:329-330
###################################################
head(rnd.2$sum.dist)


###################################################
### code chunk number 16: Statistical_Matching_with_StatMatch.Rnw:348-357
###################################################
# distance computed on the percentage points of ecdf of "age"
rnk.1 <- rankNND.hotdeck(data.rec=samp.A, data.don=samp.B, 
                         var.rec="age", var.don="age")
#create the synthetic data set
fA.rnk <- create.fused(data.rec=samp.A, data.don=samp.B,
                       mtc.ids=rnk.1$mtc.ids, 
                       z.vars="labour5", 
                       dup.x=TRUE, match.vars="age")
head(fA.rnk)


###################################################
### code chunk number 17: Statistical_Matching_with_StatMatch.Rnw:362-372
###################################################
# distance computed on the percentage points of ecdf of "age"
# computed separately by "sex"
rnk.2 <- rankNND.hotdeck(data.rec=samp.A, data.don=samp.B, var.rec="age",
                         var.don="age", don.class="sex",
                         constrained=TRUE, constr.alg="Hungarian")
fA.grnk <- create.fused(data.rec=samp.A, data.don=samp.B,
                        mtc.ids=rnk.2$mtc.ids, 
                        z.vars="labour5",
                        dup.x=TRUE, match.vars="age")
head(fA.grnk)


###################################################
### code chunk number 18: Statistical_Matching_with_StatMatch.Rnw:397-424
###################################################
# step 0) introduce missing values in iris
data(iris, package="datasets")
set.seed(1324)
miss <- rbinom(150, 1, 0.30) #generates randomly missing
iris.miss <- iris
iris.miss$Petal.Length[miss==1] <- NA
summary(iris.miss$Petal.L)
#
# step 1) separate units in two data sets
rec <- subset(iris.miss, is.na(Petal.Length), select=-Petal.Length)
don <- subset(iris.miss, !is.na(Petal.Length))
#
# step 2) search for closest donors
X.mtc <- c("Sepal.Length", "Sepal.Width", "Petal.Width")
nnd <- NND.hotdeck(data.rec=rec, data.don=don,
                         match.vars=X.mtc, don.class="Species",
                         dist.fun="Manhattan")
# fills rec
imp.rec <- create.fused(data.rec=rec, data.don=don,
                        mtc.ids=nnd$mtc.ids, z.vars="Petal.Length")
imp.rec$imp.PL <- 1 # flag for imputed
#
# step 3) re-aggregate data sets
don$imp.PL <- 0
imp.iris <- rbind(imp.rec, don)
#summary stat of imputed and non imputed Petal.Length
tapply(imp.iris$Petal.Length, imp.iris$imp.PL, summary)


###################################################
### code chunk number 19: Statistical_Matching_with_StatMatch.Rnw:440-457
###################################################
# uses iris data set
iris.A <- iris[101:150, 1:3]
iris.B <- iris[1:100, c(1:2,4)]

X.mtc <- c("Sepal.Length","Sepal.Width") # matching variables

# parameters estimated using ML
mix.1 <- mixed.mtc(data.rec=iris.A, data.don=iris.B, match.vars=X.mtc,
                    y.rec="Petal.Length", z.don="Petal.Width", 
                    method="ML", rho.yz=0, 
                    micro=TRUE, constr.alg="Hungarian")

mix.1$mu #estimated means
mix.1$cor #estimated cor. matrix

head(mix.1$filled.rec) # A filled in with Z
cor(mix.1$filled.rec)


###################################################
### code chunk number 20: Statistical_Matching_with_StatMatch.Rnw:464-471
###################################################
# parameters estimated using ML and rho_YZ|X=0.85
mix.2 <- mixed.mtc(data.rec=iris.A, data.don=iris.B, match.vars=X.mtc,
                    y.rec="Petal.Length", z.don="Petal.Width", 
                    method="ML", rho.yz=0.85, 
                    micro=TRUE, constr.alg="Hungarian")
mix.2$cor
head(mix.2$filled.rec)


###################################################
### code chunk number 21: Statistical_Matching_with_StatMatch.Rnw:476-482
###################################################
mix.3 <- mixed.mtc(data.rec=iris.A, data.don=iris.B, match.vars=X.mtc,
                    y.rec="Petal.Length", z.don="Petal.Width", 
                    method="MS", rho.yz=0.75, 
                    micro=TRUE, constr.alg="Hungarian")

mix.3$rho.yz


###################################################
### code chunk number 22: Statistical_Matching_with_StatMatch.Rnw:497-519
###################################################
# summary info on the weights
sum(samp.A$ww) # estimated pop size from A
sum(samp.B$ww) # estimated pop size from B
summary(samp.A$ww)
summary(samp.B$ww)

# NND constrained hot deck
group.v <- c("sex","area5")
out.nnd <- NND.hotdeck(data.rec=samp.A, data.don=samp.B,
                       match.vars="age", don.class=group.v,
                       dist.fun="Manhattan",
                       constrained=TRUE, constr.alg="Hungarian")

fA.nnd.m <- create.fused(data.rec=samp.A, data.don=samp.B,
                         mtc.ids=out.nnd$mtc.ids,
                         z.vars="labour5")

# estimating distribution of labour5 using weights
t1 <- xtabs(ww~labour5, data=fA.nnd.m) # imputed in A
t2 <- xtabs(ww~labour5, data=samp.B) # ref. estimate in B
c1 <- comp.prop(p1=t1, p2=t2, n1=nrow(fA.nnd.m), ref=TRUE)
c1$meas


###################################################
### code chunk number 23: Statistical_Matching_with_StatMatch.Rnw:526-539
###################################################
group.v <- c("sex","area5")
rnd.2 <- RANDwNND.hotdeck(data.rec=samp.A, data.don=samp.B, 
                          match.vars=NULL, don.class=group.v, 
                          weight.don="ww")
fA.wrnd <- create.fused(data.rec=samp.A, data.don=samp.B, 
                        mtc.ids=rnd.2$mtc.ids,
                        z.vars="labour5")

# comparing marginal distribution of labour5 using weights
tt.0w <- xtabs(ww~labour5, data=samp.B)
tt.fw <- xtabs(ww~labour5, data=fA.wrnd)
c1 <- comp.prop(p1=tt.fw, p2=tt.0w, n1=nrow(fA.wrnd), ref=TRUE)
c1$meas


###################################################
### code chunk number 24: Statistical_Matching_with_StatMatch.Rnw:550-567
###################################################
rnk.w <- rankNND.hotdeck(data.rec=samp.A, data.don=samp.B, 
                         don.class="area5", var.rec="age", 
                         var.don="age", weight.rec="ww",
                         weight.don="ww", constrained=TRUE,
                         constr.alg="Hungarian")
#
#create the synthetic data set
fA.wrnk <- create.fused(data.rec=samp.A, data.don=samp.B,
                        mtc.ids=rnk.w$mtc.ids, 
                        z.vars="labour5", 
                        dup.x=TRUE, match.vars="age")

# comparing marginal distribution of labour5 using weights
tt.0w <- xtabs(ww~labour5, data=samp.B)
tt.fw <- xtabs(ww~labour5, data=fA.wrnk)
c1 <- comp.prop(p1=tt.fw, p2=tt.0w, n1=nrow(fA.wrnk), ref=TRUE)
c1$meas


###################################################
### code chunk number 25: Statistical_Matching_with_StatMatch.Rnw:590-613
###################################################
tt.A <- xtabs(ww~sex+c.age, data=samp.A)
tt.B <- xtabs(ww~sex+c.age, data=samp.B)
(prop.table(tt.A)-prop.table(tt.B))*100
comp.prop(p1=tt.A, p2=tt.B, n1=nrow(samp.A),
          n2=nrow(samp.B), ref=FALSE)

library(survey, warn.conflicts=FALSE) # loads survey
# creates svydesign objects
svy.samp.A <- svydesign(~1, weights=~ww, data=samp.A)
svy.samp.B <- svydesign(~1, weights=~ww, data=samp.B)
#
# harmonizes wrt to joint distr. of gender vs. c.age
out.hz <- harmonize.x(svy.A=svy.samp.A, svy.B=svy.samp.B,
                      form.x=~c.age:sex-1)
#
summary(out.hz$weights.A) # new calibrated weights for A
summary(out.hz$weights.B) # new calibrated weights for B

tt.A <- xtabs(out.hz$weights.A~sex+c.age, data=samp.A)
tt.B <- xtabs(out.hz$weights.B~sex+c.age, data=samp.B)
c1 <- comp.prop(p1=tt.A, p2=tt.B, n1=nrow(samp.A),
                n2=nrow(samp.B), ref=FALSE)
c1$meas


###################################################
### code chunk number 26: Statistical_Matching_with_StatMatch.Rnw:629-635
###################################################
# estimating c.netI vs. labour5 under the CI assumption
out <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                    svy.C=NULL, y.lab="c.neti", z.lab="labour5",
                    form.x=~c.age:sex-1)
#
addmargins(t(out$yz.CIA))  # table estimated under the CIA


###################################################
### code chunk number 27: Statistical_Matching_with_StatMatch.Rnw:640-653
###################################################
data(samp.C, package="StatMatch")
str(samp.C)

#
svy.samp.C <- svydesign(~1, weights=~ww, data=samp.C) 

#
# incomplete two-way estimation
out.inc <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                        svy.C=svy.samp.C, y.lab="c.neti", 
												z.lab="labour5", form.x=~c.age:sex-1, 
												estimation="incomplete")
addmargins(t(out.inc$yz.est))           


###################################################
### code chunk number 28: Statistical_Matching_with_StatMatch.Rnw:658-680
###################################################
new.ww <- weights(out.inc$cal.C) #new cal. weights for C 
#
# marginal distributions of c.neti
m.work.cA <- xtabs(out.hz$weights.A~c.neti, data=samp.A)
m.work.cC <- xtabs(new.ww~c.neti, data=samp.C)
m.work.cA-m.work.cC
#
# marginal distributions of labour5
m.cnetI.cB <- xtabs(out.hz$weights.B~labour5, data=samp.B)
m.cnetI.cC <- xtabs(new.ww~labour5, data=samp.C)
m.cnetI.cB-m.cnetI.cC 

# joint distribution of the matching variables
tt.A <- xtabs(out.hz$weights.A~sex+c.age, data=samp.A)
tt.B <- xtabs(out.hz$weights.B~sex+c.age, data=samp.B)
tt.C <- xtabs(new.ww~sex+c.age, data=samp.C)
c1 <- comp.prop(p1=tt.A, p2=tt.B, n1=nrow(samp.A),
                n2=nrow(samp.B), ref=FALSE)
c2 <- comp.prop(p1=tt.C, p2=tt.A, n1=nrow(samp.C),
                n2=nrow(samp.A), ref=FALSE)
c1$meas
c2$meas


###################################################
### code chunk number 29: Statistical_Matching_with_StatMatch.Rnw:685-692
###################################################
# synthetic two-way estimation
out.synt <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                         svy.C=svy.samp.C, y.lab="c.neti", 
												 z.lab="labour5", form.x=~c.age:sex-1, 
												 estimation="synthetic")
#
addmargins(t(out.synt$yz.est))           


###################################################
### code chunk number 30: Statistical_Matching_with_StatMatch.Rnw:701-714
###################################################
# predicting prob of labour5 in A under the CI assumption
out <- comb.samples(svy.A=out.hz$cal.A, svy.B=out.hz$cal.B,
                    svy.C=NULL, y.lab="c.neti", z.lab="labour5",
                    form.x=~c.age:sex-1, micro=TRUE)
head(out$Z.A)
sum(out$Z.A<0) # negative est. prob.
sum(out$Z.A>1) # est. prob. >1

# compare marginal distributions of Z
t.zA <- colSums(out$Z.A * out.hz$weights.A)
t.zB <- xtabs(out.hz$weights.B ~ samp.B$labour5)
c1 <- comp.prop(p1=t.zA, p2=t.zB, n1=nrow(samp.A), ref=TRUE)  
c1$meas


###################################################
### code chunk number 31: Statistical_Matching_with_StatMatch.Rnw:719-738
###################################################
# predicting categories of labour5 in A
# randomized prediction with prob proportional to estimated prob.
pps1 <- function(x) sample(x=1:length(x), size=1, prob=x)
pred.zA <- apply(out$Z.A, 1, pps1)
samp.A$labour5 <- factor(pred.zA, levels=1:nlevels(samp.B$labour5), 
                       labels=as.character(levels(samp.B$labour5)), 
                       ordered=T)

# comparing marginal distributions of Z
t.zA <- xtabs(out.hz$weights.A ~ samp.A$labour5)
c1 <- comp.prop(p1=t.zA, p2=t.zB, n1=nrow(samp.A), ref=TRUE)  
c1$meas

# comparing joint distributions of X vs. Z
t.xzA <- xtabs(out.hz$weights.A~c.age+sex+labour5, data=samp.A)
t.xzB <- xtabs(out.hz$weights.B~c.age+sex+labour5, data=samp.B)
out.comp <- comp.prop(p1=t.xzA, p2=t.xzB, n1=nrow(samp.A), ref=TRUE)  
out.comp$meas
out.comp$chi.sq


###################################################
### code chunk number 32: Statistical_Matching_with_StatMatch.Rnw:779-790
###################################################
#comparing joint distribution of the X_M variables in A and in B
t.xA <- xtabs(ww~c.age+sex, data=samp.A)
t.xB <- xtabs(ww~c.age+sex, data=samp.B)
comp.prop(p1=t.xA, p2=t.xB, n1=nrow(samp.A), n2=nrow(samp.B), ref=FALSE)
#
#computing tables needed by Frechet.bounds.cat
t.xy <- xtabs(ww~c.age+sex+c.neti, data=samp.A)
t.xz <- xtabs(ww~c.age+sex+labour5, data=samp.B)
out.fb <- Frechet.bounds.cat(tab.x=t.xA, tab.xy=t.xy, tab.xz=t.xz, 
                             print.f="data.frame")
out.fb


###################################################
### code chunk number 33: Statistical_Matching_with_StatMatch.Rnw:798-807
###################################################
# continuous variables
samp.A$log.netI <- log(ifelse(samp.A$n.income>0, samp.A$n.income, 0) + 1)
lab <- as.integer(samp.B$labour5)
samp.B$work <- factor(ifelse(lab<3, 1, 2)) # binary variable working status

X.mtc <- c("age", "sex")
mix.3 <- mixed.mtc(data.rec=samp.A, data.don=samp.B, match.vars=X.mtc,
                   y.rec="log.netI", z.don="work", 
                   method="MS")


