p.hat <- 52/100; p.hat            
SE <- sqrt(p.hat * (1 - p.hat) / 100); SE    # est. SE
p.hat - 1.96 * SE                            # lower end of CI 
p.hat + 1.96 * SE                            # upper end of CI

