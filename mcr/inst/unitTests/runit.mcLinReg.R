cat("\n\nmcLinReg.R method comparison test cases\n\n")

## check computation of analytical CIs

test.mc.linreg.call <- function() 
{
    data(creatinine)
    set.seed(19061978)
    smpl <- sample(110, 30)
    
    checkException(mcr:::mc.linreg(as.numeric(NA), creatinine[smpl,1]))
    checkException(mcr:::mc.linreg(numeric(0), numeric(0)))
    checkException(mcr:::mc.linreg(creatinine[smple[-1], 1], creatinine[smpl,2]))
    
    res <- mcr:::mc.linreg(creatinine[smpl,1], creatinine[smpl,2])
    
    lm.fit <- lm(plasma.crea~serum.crea, creatinine[smpl,])
    Coef <- coef(lm.fit)
    Smry <- summary(lm.fit)$coefficients
    checkEquals( res, list(b0=as.numeric(Coef[1]), b1=as.numeric(Coef[2]),  
                           se.b0=Smry[1,2], se.b1=Smry[2,2], 
                           xw=mean(creatinine[smpl,1]), 
                           weight=rep(1,dim(creatinine[smpl,])[1])) )
}

test.mc.wlinreg.call <- function()
{
    data(creatinine)
    set.seed(19061978)
    smpl <- sample(110, 30)
    
    checkException(mcr:::mc.wlinreg(as.numeric(NA), creatinine[smpl,1]))
    checkException(mcr:::mc.wlinreg(numeric(0), numeric(0)))
    checkException(mcr:::mc.wlinreg(creatinine[smple[-1], 1], creatinine[smpl,2]))
    
    res <- mcr:::mc.wlinreg(creatinine[smpl,1], creatinine[smpl,2])
    
    lm.fit <- lm(plasma.crea~serum.crea, creatinine[smpl,], weights=1/creatinine[smpl,1]^2)

    w <- creatinine[smpl,1]^-2

    Coef <- coef(lm.fit)
    Smry <- summary(lm.fit)$coefficients
    checkEquals( res, list(b0=as.numeric(Coef[1]), b1=as.numeric(Coef[2]), 
                           se.b0=Smry[1,2],se.b1=Smry[2,2], 
                           xw=sum(creatinine[smpl,1]*1/creatinine[smpl,1]^2)/sum(1/creatinine[smpl,1]^2), 
                           weight=w) )
}