###hop:9-15
summary(glht(chol.model, mcp(trt = 
    rbind(
        "1time - 2times" = c(1,-1,0,0,0),
        "(1 or 2 times) - 4times" = c(0.5,0.5,-1,0,0),
        "new - old" = c(2,2,2,-3,-3)/6,
        "drugD - drugE" = c(0,0,0,1,-1))
    )))
