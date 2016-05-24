setClass("BAYESRETURN", representation(

alpha="numeric",
beta="numeric",
var_alpha="numeric",
a_inc  = "numeric",
fst = "numeric",
mfst="numeric",
P="numeric"
))

setMethod("show", "BAYESRETURN",
 function(object){

out <- data.frame(Slots=c("alpha","beta","a_inc","fst"),Description=c("locus effect","population effect","alpha included","FST values"))
print(out)

})

setGeneric("getalpha", function(object) standardGeneric("getalpha"))
 setMethod("getalpha", "BAYESRETURN",
 function(object){
 return(object@alpha)
 })

