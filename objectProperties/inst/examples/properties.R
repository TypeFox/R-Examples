## we could pass prototype as in S4
GPars.gen <- setRefClass("GraphicProperties",
                fields = properties(fields = list(size = "numeric",
                                                  color = "character"),
                                    prototype = list(size =1,
                                                     color = "red")))

obj <- GPars.gen$new()
## since it's not PropertySet, no global signal
## let's register individual signal
obj$sizeChanged$connect(function(){
  print("size changed")
})
## emit signal
obj$size <- 3
## no signal
obj$color <- "black"

