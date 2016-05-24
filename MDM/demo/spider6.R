####################################
### Spider Example -- 6 Species ####
####################################

## Read Spider data ##
data(spider6)

########################
## ABUNDANCE ANALYSIS ##
########################

spider6$Abundance <- apply(spider6[,1:6],1,sum)
fit0 <- glm(Abundance~1,data=spider6,family=quasipoisson)
fit1 <- glm(Abundance~Herbs,data=spider6,family=quasipoisson)
fit2 <- glm(Abundance~Water,data=spider6,family=quasipoisson)
fit3 <- glm(Abundance~Herbs+Water,data=spider6,family=quasipoisson)
anova(fit0,fit2,fit3,test="F")

########################
## DIVERSITY ANALYSIS ##
########################

## scale rows to sum p's = 1 ##
y <- as.matrix(spider6[,1:6])
y <- y/apply(y,1,sum)
## fit models ##
mx <- 1000
fit0 <- mdm(y~1,data=spider6,maxit=mx)
fit1 <- mdm(y~Water,data=spider6,maxit=mx)
fit2 <- mdm(y~Water+Herbs,data=spider6,maxit=mx)
fit3 <- mdm(y~Water+Herbs+I(Water^2),data=spider6,maxit=mx)
fit4 <- mdm(y~Water+Herbs+I(Water^2)+I(Herbs^2),data=spider6,maxit=mx)
fit00 <- mdm(y~Site,data=spider6,maxit=mx)
anova(fit0,fit1,fit2,fit3,fit4,fit00)
