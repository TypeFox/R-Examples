library(gnm)
set.seed(1)

count <- with(voting, percentage/100 * total)
yvar <- cbind(count, voting$total - count)

classMobility <- gnm(yvar ~ Dref(origin, destination),
                     family = binomial, data = voting)

print(classMobility$deviance, digits = 10)
print(classMobility$df)

upward <- with(voting, origin != 1 & destination == 1)
downward <- with(voting, origin == 1 & destination != 1)

socialMobility <- gnm(yvar ~ Dref(origin, destination,
                                  delta = ~ 1 + downward + upward),
                      family = binomial, data = voting)

print(socialMobility$deviance, digits = 10)
print(socialMobility$df)
