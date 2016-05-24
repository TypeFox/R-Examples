cat("\ntest sp+nsp ranefs:")
data(blackcap)
set.seed(123)
somegrp <- cbind(blackcap,grp=sample(2,14,replace=TRUE))
#somegrp <- cbind(blackcap,grp=c(rep(1,7),rep(2,7))) ## to test cov mat with nonzero var of grp effect
fitobject <- corrHLfit(migStatus ~ 1 +  (1|grp) +Matern(1|latitude+longitude),data=somegrp,
                       ranFix=list(nu=4,rho=0.4,phi=0.05))
res <- get_predVar(fitobject,newdata=somegrp[1:5,])
expect_equal(res,c("Gibraltar"=0.04458288, "CapeVerde"=0.04779076, "SouthernFrance"=0.03575583, "LaPalma"=0.02247332, "Madeira"=0.04352935),tolerance=1e-5)

grouped <- cbind(blackcap,grp=c(rep(1,7),rep(2,7))) 
fitobject <- corrHLfit(migStatus ~ 1 +  (1|grp) +Matern(1|latitude+longitude),
                       data=grouped,  ranFix=list(nu=4,rho=0.4,phi=0.05))
p1 <- predict(fitobject)
p2 <- predict(fitobject,newdata=grouped)
expect_equal(max(abs(p1-p2)),0,tolerance=1e-8)

