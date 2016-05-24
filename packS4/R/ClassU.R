##################################################
###                    ClassU.R
setClass("ClassU",representation(u1="numeric",u2="numeric"))
setMethod("privateA","ClassU",function(object){object@u1^2})
setMethod("plot","ClassU",function(x,y){barplot(c(x@u1,x@u2))})
