context('package')

test_that('all', {
    library(testthat)
    ##linear
    data <- t(t(matrix(rnorm(1000),200)) + 1:5)  
    mod1 <- OLScurve(~ time, data = data)        
    expect_is(mod1, 'OLScurve')          
    plt <- plot(mod1)
    expect_is(plt, 'trellis')              
    
    ##quadratic
    data <- t(t(matrix(rnorm(1000),200)) + (0:4)^2)  
    mod2 <- OLScurve(~ time + I(time^2), data = data)        
    expect_is(mod2, 'OLScurve')          
    plt <- plot(mod2)
    expect_is(plt, 'trellis')          
    
    ##sqrt
    data <- t(t(matrix(rnorm(1000),200)) + 20*sqrt(5:1))   
    mod3 <- OLScurve(~ sqrt(time), data = data)    
    expect_is(mod3, 'OLScurve')          
    plt <- plot(mod3)
    expect_is(plt, 'trellis')          
    
    ##exponential
    data <- t(t(matrix(rnorm(1000,0,5),200)) + exp(0:4))  
    mod4 <- OLScurve(~ exp(time), data = data)    
    expect_is(mod4, 'OLScurve')          
    plt <- plot(mod4)
    expect_is(plt, 'trellis')          
    
    ##combination
    data <- t(t(matrix(rnorm(1000),200)) + 20*sqrt(1:5))  
    mod5 <- OLScurve(~ time + sqrt(time), data = data)        
    expect_is(mod5, 'OLScurve')          
    plt <- plot(mod5)
    expect_is(plt, 'trellis')          
    
    ##piecewise (global linear trend with linear shift at time point 4)
    data <- t(t(matrix(rnorm(1000),200)) + (0:4)^2) 
    time <- data.frame(time1 = c(0,1,2,3,4), time2 = c(0,0,0,1,2))
    mod6 <- OLScurve(~ time1 + time2, data, time=time)    
    expect_is(mod6, 'OLScurve')          
    plt <- plot(mod6)
    expect_is(plt, 'trellis')          
    
    ##two group analysis with linear trajectories
    data1 <- t(t(matrix(rnorm(500),100)) + 1:5) 
    data2 <- t(t(matrix(rnorm(500),100)) + 9:5)
    data <- rbind(data1,data2) 
    group <- c(rep('male',100),rep('female',100)) 
    
    mod <- OLScurve(~ time, data)    
    expect_is(mod, 'OLScurve')          
    plt <- plot(mod,group)
    expect_is(plt, 'trellis') 
    
    plt <- parplot(mod)
    expect_is(plt, 'trellis') 
    plt <- parplot(mod, type = 'boxplot')
    expect_is(plt, 'trellis') 
    plt <- parplot(mod, type = 'splom')
    expect_is(plt, 'trellis') 
     
    plt <- parplot(mod, group=group)
    expect_is(plt, 'trellis') 
    plt <- parplot(mod, type='boxplot', group=group)
    expect_is(plt, 'trellis') 
    plt <- parplot(mod, type='splom', group=group)
    expect_is(plt, 'trellis') 
})

