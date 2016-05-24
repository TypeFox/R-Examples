## HE plots and candisc HE plots for MASS::painters data

library(MASS)
library(heplots)

data(painters)

# longer labels to identify the schools
school <- c("Renaissance", "Mannerist", "Sciento", "Venetian",
		"Lombard", "16th C", "17th C", "French")
levels(painters$School) <- school

# how do the schools differ according to the aesthetic qualities?
painters.mod <- lm(cbind(Composition,Drawing,Colour,Expression)~School, data=painters)
heplot(painters.mod)
pairs(painters.mod)


library(candisc)
# How many dimensions of differences?
painters.can <- candisc(painters.mod)
painters.can
summary(painters.can)
heplot(painters.can)

# There seem to be 3 significant dimensions. View it in 3D
heplot3d(painters.can, col=c("pink", "brown"))
						
						