# P(55 <= X <= 65)
pbinom(65,100,0.6) - pbinom(54,100,0.6)        
# without continuity correction:
pnorm(65,60,sqrt(100*0.6*0.4)) -  pnorm(55,60,sqrt(100*0.6*0.4))
# with continuity correction:
pnorm(65.5,60,sqrt(100*0.6*0.4)) -  pnorm(54.5,60,sqrt(100*0.6*0.4))
