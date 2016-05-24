coag.aov <- aov(coag~diet,coagulation); coag.aov
TukeyHSD(coag.aov)
