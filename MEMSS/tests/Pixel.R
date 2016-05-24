library(MEMSS)
options(show.signif.stars = FALSE, digits = 3)
m1 <- lmer(pixel ~ day + I(day^2) + (1|Dog:Side) + (day|Dog),
           Pixel)
print(m1, corr = FALSE)
unclass(ranef(m1))
