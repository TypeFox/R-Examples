#####################################
#### Data Analysis for Coral Data ###
#####################################

data(corals)

#################
## ABUNDANCE  ###
#################

fit <- glm(abundance~tr+yr,data=corals,family="Gamma")
anova(fit,test="F")
drop1(fit,test="F")
par(mfrow=c(2,2),mar=c(3,4,1.5,1.5),mgp=c(2,0.7,0))
plot(fit)
par(mfrow=c(1,1))

#################
## DIVERSITY  ###
#################

y <- y2p(corals[,-c(1:6)])
fit0 <- mdm(y~1,data=corals,maxit=mx)
fit1 <- mdm(y~tr,data=corals,maxit=mx)
fit2 <- mdm(y~tr+yr,data=corals,maxit=mx)
fit3 <- mdm(y~factor(1:nrow(corals)),data=corals,maxit=mx)

sapply(list(fit0,fit1,fit2,fit3),dev2div)
anova(fit0,fit1,fit2,fit3)

###################
## Entropy plots ##
###################

e0 <- fitted(fit0)
e1 <- fitted(fit1)
e2 <- fitted(fit2)
e3 <- fitted(fit3)

x11(hei=12,wid=7)
par(mar=c(4,11,1,1),las=1,mgp=c(2.5,0.75,0))
labs <- c(rep("A",25),rep("P",9),rep("O",8),rep("P",7),rep("O",10),rep("A",3),"O")
entropy.plot(list(e0,e3,e1,e2),labs=paste(dimnames(e0)[[2]],labs,sep=" "),ord=T)
