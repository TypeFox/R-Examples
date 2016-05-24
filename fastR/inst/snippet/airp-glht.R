airp.cint <- confint(glht(airp.model,mcp(location="Tukey")))
airp.cint  
plot(TukeyHSD(airp.aov));  plot(airp.cint)   # plots
