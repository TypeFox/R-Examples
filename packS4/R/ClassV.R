##################################################
###                   ClassV.R
setClass("ClassV",representation(v1="numeric",v2="numeric"))
setMethod("privateA","ClassV",function(object){object@v1^3})
setMethod("publicA","ClassV",function(object){sqrt(object@v2^3)})
setMethod("plot","ClassV",function(x,y){barplot(c(x@v1,x@v2^2))})
