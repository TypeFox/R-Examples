# The data set in Devore6 is slightly different, but the Devore6 package
# has been archived, so we'll substitute the similar one from Devore 7 here.
require(Devore7); data(xmp13.13, package='Devore7')
concrete <- data.frame(
    limestone=xmp13.13$x1,
    water=xmp13.13$x2,
    strength=xmp13.13$strength)
concrete.lm1 <- lm(strength ~ limestone + water, concrete)
###hop:3-6
concrete.lm1
