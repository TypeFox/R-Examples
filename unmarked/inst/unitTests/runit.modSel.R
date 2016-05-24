test.fitList <- function() {  
    y <- matrix(rep(1, 10), 5, 2)
    umf <- unmarkedFrameOccu(y = y, siteCovs=data.frame(x=-2:2), 
        obsCovs= data.frame(z=-5:4))    
    obsCovs(umf)[3, 1] <- NA
    fm1 <- occu(~ 1 ~ 1, data = umf)
    fm2 <- occu(~ 1 ~ x, data = umf)
    
    fits1.1 <- fitList(m1=fm1, m2=fm2)
    fits1.2 <- fitList(fm1, fm2)
    fits2.1 <- fitList(fits = list(m1=fm1, m2=fm2))
    fits2.2 <- fitList(fits = list(fm1, fm2))
        
    checkIdentical(fits1.1, fits2.1)
    
    checkException(fitList(fm1, fm2, fits=list(fm1, fm2)))  # ... and fits    

    siteCovs(umf) <- data.frame(x=-3:1)
    fm2 <- occu(~ 1 ~ x, data = umf)    
    checkException(fitList(fm1, fm2))   # Different umf used

    fm3 <- occu(~ z ~ 1, data = umf)
    checkException(fitList(fm1, fm3))   # Missing value problem
    }





test.modSel <- function() {  
    y <- matrix(rep(1, 10), 5, 2)
    umf <- unmarkedFrameOccu(y = y, siteCovs=data.frame(x=-2:2), 
        obsCovs= data.frame(z=-5:4))    
    fm1 <- occu(~ 1 ~ 1, data = umf)
    fm2 <- occu(~ 1 ~ x, data = umf)
    
    fits <- fitList(m1=fm1, m2=fm2)
    ms1 <- modSel(fits)
    
    checkTrue(all(is.na(ms1@Full$Rsq)))
    checkEqualsNumeric(sum(ms1@Full$AICwt), 1)
    checkEqualsNumeric(ms1@Full$delta[1L], 0)
    
    checkException(modSel(fits, nullmod=fm2))

    ms2 <- modSel(fits, nullmod='m1')    
    
    checkIdentical(
        ms1@Full[,-which(colnames(ms1@Full)=="Rsq")],
        ms1@Full[,-which(colnames(ms2@Full)=="Rsq")]        
        )
    
    # Fake hessian problem    
    fm1@opt$hessian[] <- NA
    fm1@estimates@estimates$state@covMat[] <- NA
    fits2 <- fitList(m1=fm1, m2=fm2)
    ms3 <- modSel(fits2)
    checkEquals(coef(ms1), coef(ms3))
    
    
    }

