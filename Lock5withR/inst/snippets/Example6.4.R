p.hat <- 0.28; p.hat            
SE <- sqrt(p.hat * (1 - p.hat) / 800); SE    # est. SE
p.hat - 1.96 * SE                            # lower end of CI 
p.hat + 1.96 * SE                            # upper end of CI
confint(prop.test(224, 800))                 # 224 = 0.28 * 800

