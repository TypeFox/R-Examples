filt.gen <- setRefClass("Filter", properties(fields = list(cutoff = "numeric",
                                                           weight = "numeric"),
                                  prototype = list(cutoff = 0, weight = 1)),
                                  contains = "PropertySet")
obj <- filt.gen$new()
obj
obj$properties()
as.list(obj)
obj$changed$connect(function(name) print(name))
obj$cutoffChanged$connect(function() print(paste("change to", obj$cutoff)))
obj$cutoff <- 0
obj$cutoff <- 2
obj$weight <- 3


## use setPropertySet, the same thing as above
filt.gen <- setPropertySet("Filter", fields = list(cutoff = "numeric",
                                         weight = "numeric"),
                           prototype = list(cutoff = 0, weight = 1))

obj <- filt.gen$new()
obj
obj$properties()
as.list(obj)
obj$changed$connect(function(name) print(name))
obj$cutoffChanged$connect(function() print(paste("change to", obj$cutoff)))
obj$cutoff <- 0
obj$cutoff <- 2
obj$weight <- 3

