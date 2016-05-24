ice.trt2 <- lm(T1930 - B1930 ~ Treatment, ice,
                 subset=Location=='intramuscular')
###hop:3-9
summary(ice.trt2)

###hop:3-13
confint(glht(ice.trt2, mcp(Treatment='Tukey')),level=0.90)
