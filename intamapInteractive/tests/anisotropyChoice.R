library(intamapInteractive)
library(gstat)
data(walker)
object=createIntamapObject(observations=walker)
object=anisotropyChoice(object)

print(summary(object$clusters$index))
print(object$anisPar)



object=doSegmentation(object)

print(summary(object$clusters$index))
