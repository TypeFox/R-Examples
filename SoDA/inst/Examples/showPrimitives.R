allObjects <- objects("package:base", all=TRUE)
primitives <- sapply(allObjects, 
   function(x)is.primitive(get(x)))
primFuns  <- allObjects[primitives]; primFuns
