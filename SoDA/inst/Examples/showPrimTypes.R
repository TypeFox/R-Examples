allObjects = objects("package:base", all=TRUE)
primitives = sapply(allObjects, function(x)is.primitive(get(x)))
primFuns  = allObjects[primitives]
functionTypes <- split(primFuns,
               sapply(primFuns, function(x)typeof(get(x))))
names(functionTypes)
sapply(functionTypes, length)
functionTypes$special
