if(!isGeneric("Map")){
    setGeneric("Map", function(f, ...) standardGeneric("Map"))
}
if(!isGeneric("Domain")){
    setGeneric("Domain", function(object) standardGeneric("Domain"))
}
if(!isGeneric("Range")){ 
    setGeneric("Range", function(object) standardGeneric("Range"))
}
if(!isGeneric("Dim")){
    setGeneric("Dim", function(object) standardGeneric("Dim"))
}
if(!isGeneric("nrow")){
    setGeneric("nrow", function(x) standardGeneric("nrow"))
}
if(!isGeneric("ncol")){
    setGeneric("ncol", function(x) standardGeneric("ncol"))
}
if(!isGeneric("Map<-")){
    setGeneric("Map<-", function(object, value) standardGeneric("Map<-"))
}
if(!isGeneric("Domain<-")){
    setGeneric("Domain<-", function(object, value) standardGeneric("Domain<-"))
}
if(!isGeneric("Range<-")){
    setGeneric("Range<-", function(object, value) standardGeneric("Range<-"))
}
if(!isGeneric("Dim<-")){
    setGeneric("Dim<-", function(object, value) standardGeneric("Dim<-"))
}
if(!isGeneric("compatibleDomains")){ 
    setGeneric("compatibleDomains", function(e1, e2) standardGeneric("compatibleDomains"))
}
if(!isGeneric("evalRandVar")){
    setGeneric("evalRandVar", function(RandVar, x, distr) standardGeneric("evalRandVar"))
}
if(!isGeneric("imageDistr")){ 
    setGeneric("imageDistr", function(RandVar, distr) standardGeneric("imageDistr"))
}
if(!isGeneric("t")){ 
    setGeneric("t", function(x) standardGeneric("t"))
}
if(!isGeneric("%m%")){ 
    setGeneric("%m%", function(x, y) standardGeneric("%m%"))
}
if(!isGeneric("numberOfMaps")){
    setGeneric("numberOfMaps", function(object) standardGeneric("numberOfMaps"))
}
