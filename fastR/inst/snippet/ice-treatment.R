ice.trt <- lm(T1930 - B1930 ~ Treatment * Location, ice)
anova(ice.trt)
