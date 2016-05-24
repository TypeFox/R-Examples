##############################################
#                                            #
#                 Chapter 4                  #
#                                            #
##############################################

 
library(WWGbook)

head(classroom)

attach(classroom)


###############################
#      Descriptive Analyses   #
###############################


 
#### Descriptive Summary #1

# Level 1 Descriptive Statistics

 
level1 <- data.frame(sex,minority,mathkind,mathgain,ses)
summary(level1)
dim(level1)

# Level 2 Descriptive Statistics

level2 <- aggregate(classroom,list(classid = classroom$classid),mean)
summary(level2$yearstea)
dim(level2)

# Level 3 Descriptive Statistics

level3 <- aggregate(classroom,list(schoolid = classroom$schoolid),mean)
summary(level3$housepov)
dim(level3)


###############################
#      Descriptive Analyses   #
###############################


 
#### Descriptive Summary #1

# Level 1 Descriptive Statistics

 
level1 <- data.frame(sex,minority,mathkind,mathgain,ses)
summary(level1)
dim(level1)

# Level 2 Descriptive Statistics

level2 <- aggregate(classroom,list(classid = classroom$classid),mean)
summary(level2$yearstea)
dim(level2)

# Level 3 Descriptive Statistics

level3 <- aggregate(classroom,list(schoolid = classroom$schoolid),mean)
summary(level3$housepov)
dim(level3)

 
#### Descriptive Summary #2
 
classroom.nomiss <- subset(classroom, !is.na(mathknow))

# Level 1 Descriptive Statistics

level1 <- data.frame(classroom.nomiss$sex,classroom.nomiss$minority,classroom.nomiss$mathkind,classroom.nomiss$mathgain,classroom.nomiss$ses)
summary(level1)
dim(level1)

# Level 2 Descriptive Statistics

level2.agg <- aggregate(classroom.nomiss,list(classid = classroom.nomiss$classid),mean)
level2 <- data.frame(level2.agg$yearstea,level2.agg$mathknow,level2.agg$mathprep)
summary(level2)
dim(level2)

# Level 3 Descriptive Statistics

level3.agg <- aggregate(classroom.nomiss,list(schoolid = classroom.nomiss$schoolid),mean)
summary(level3.agg$housepov)
dim(level3.agg)

 
#### Boxplots for Figure 4.2. 
 
classroom.first8 <- classroom[classroom$schoolid <= 8,]
par(mfrow=c(4,2))
for (i in 1:8) 
{boxplot(classroom.first8$mathgain[classroom.first8$schoolid==i] ~ classroom.first8$classid[classroom.first8$schoolid==i],ylab="Mathgain", xlab="Classid" )}

###############################
#      Models Fitted          #
###############################

library(nlme)

# Model 4.1.
model4.1.fit <- lme(mathgain ~ 1, random = ~1 | schoolid/classid, classroom, method = "REML")
summary(model4.1.fit)
random.effects(model4.1.fit)

# Model 4.1A.
model4.1A.fit <- lme(mathgain ~ 1, random = ~1 | schoolid, classroom, method = "REML")
anova(model4.1.fit, model4.1A.fit)

# Model 4.2.
model4.2.fit <- lme(mathgain ~ mathkind + sex + minority + ses, random = ~1 | schoolid/classid, classroom, 
   na.action = "na.omit", method = "REML")
summary(model4.2.fit)
ranef(model4.2.fit)

# Model 4.1: ML estimation.
model4.1.ml.fit <- lme(mathgain ~ 1, 
random = ~1 | schoolid/classid, classroom, method = "ML")

# Model 4.2: ML estimation.
model4.2.ml.fit <- lme(mathgain ~ mathkind + sex + minority + ses, random = ~1 | schoolid/classid, classroom, 
na.action = "na.omit", method = "ML")

anova(model4.1.ml.fit, model4.2.ml.fit)

# Model 4.3.
model4.3.fit <- update(model4.2.fit, 
fixed = ~ mathkind + sex + minority + ses + yearstea + mathprep + mathknow)

summary(model4.3.fit)

# Model 4.4.
model4.4.fit <- update(model4.2.fit, 
fixed = ~ mathkind + sex + minority + ses + housepov)

summary(model4.4.fit)


intervals(model4.4.fit)

###getVarCov(model4.4.fit, individual="1", type="marginal")

############################
#  Model Diagnostics       #
############################

library(lattice)
trellis.device(color=F)

#### Figure 4.4
#pred <- predict(model4.2.fit)
#classroom.1 <- data.frame(classroom,pred)
#postscript(name="threshold_all.jpf",onfile=F, print.it=F)
x <- mathkind
y <- mathgain
classroom.1<-data.frame(classroom,x,y)
fig.4panel <- function(x,y,subscripts,groups) {
                panel.superpose(x,y,subscripts,groups)
                for(g in as.character(unique(groups[subscripts]))) {
                        sub <- groups[subscripts] == g
                        fit <- lm(y~x,subset=sub)
                        ord <- order(x[sub])
                        lines(x[sub][ord],fitted(fit)[ord]) 
                }
        }

xyplot(y ~ x, groups=minority,data=classroom.1) #panel=fig.4panel)


#### Figure 4.5
qqnorm(model4.2.fit, ~ranef(.,level=2))

###Figure 4.6
qqnorm(model4.2.fit, ~ranef(.,level=1))

#### Figure 4.7
qqnorm(model4.2.fit, ~resid(.),plot.it = TRUE ) 
 

#### Figure 4.8 
plot(model4.2.fit, resid(., type="p") ~ fitted(.), abline=0)


detach(classroom) 

 

