M <- cbind(                                    # model matrix
        "C1" = rep(c(-1,-1,1,1),each=4)/2,     # C1
        "C2" = rep(c(-1,1,-1,1),each=4)/2,     # C2
        "C3" = rep(c(1,-1,-1,1),each=4)/4      # C3
        )
taste.lm2 <- lm(score~M,data=tastetest)
###hop:3-9
summary(taste.lm2)
