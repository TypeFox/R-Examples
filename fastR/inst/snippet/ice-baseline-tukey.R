require(multcomp) 
confint(glht(base.lmadd, mcp(Treatment="Tukey")),level=0.9)
