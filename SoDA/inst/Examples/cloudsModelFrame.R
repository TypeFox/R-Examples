(if(file.access("Examples/clouds.rda") == 0) load("Examples/clouds.rda")
else require(HSAUR))
formula <- rainfall ~ seeding *   (sne + cloudcover + prewetness + echomotion) + time
mf <- model.frame(formula, clouds)
class(mf)
names(attributes(mf))
