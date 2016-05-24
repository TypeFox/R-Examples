
## Access methods
if(!isGeneric("radDistr")) 
setGeneric("radDistr", function(object) standardGeneric("radDistr"))

if(!isGeneric("sigma")) 
   setGeneric("sigma", function(object) standardGeneric("sigma"))


# Replacement methods
if(!isGeneric("radDistr<-")) 
    setGeneric("radDistr<-", function(object,value) standardGeneric("radDistr<-"))

## wrappers
if(!isGeneric("r.rd")) 
   setGeneric("r.rd", function(object) standardGeneric("r.rd"))
if(!isGeneric("d.rd")) 
   setGeneric("d.rd", function(object) standardGeneric("d.rd"))
if(!isGeneric("p.rd")) 
   setGeneric("p.rd", function(object) standardGeneric("p.rd"))
if(!isGeneric("q.rd")) 
   setGeneric("q.rd", function(object) standardGeneric("q.rd"))
if(!isGeneric("plot.rd")) 
   setGeneric("plot.rd", function(x, ...) standardGeneric("plot.rd"))


