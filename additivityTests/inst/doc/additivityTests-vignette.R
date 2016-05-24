
## ----'generating data with no interaction'-------------------------------
    require(additivityTests)
    set.seed(123)
    subjects = rnorm(10)
    treatments = rnorm(10)
    noise = rnorm(100)/100
    Y = matrix(rep(subjects,10), 10, 10) + matrix(rep(treatments, each=10), 10, 10) + noise


## ----'testing for interaction'-------------------------------------------
    tukey.test(Y)
    mandel.test(Y)
    lbi.test(Y)
    tusell.test(Y)
    johnson.graybill.test(Y)
    mandel.test(Y)
    mtukey.test(Y, correction=2, Nboot=1000)


## ----'adding interaction and testing again'------------------------------
    Y[1:5,] = Y[1:5,] + 10*rep(treatments, each=5)
    
    tukey.test(Y)
    mandel.test(Y)
    lbi.test(Y)
    tusell.test(Y)
    johnson.graybill.test(Y)
    mandel.test(Y)
    mtukey.test(Y, correction=2, Nboot=1000)


