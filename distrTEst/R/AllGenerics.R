## Access methods
if(!isGeneric("call.ev")) 
   setGeneric("call.ev", function(object) standardGeneric("call.ev"))
if(!isGeneric("result")) 
   setGeneric("result", function(object) standardGeneric("result"))
if(!isGeneric("estimator")) 
   setGeneric("estimator", function(object) standardGeneric("estimator"))
if(!isGeneric("Elist")) 
    setGeneric("Elist", function(object) standardGeneric("Elist"))

if(!isGeneric("name")) 
    setGeneric("name", function(object) standardGeneric("name"))
if(!isGeneric("name<-")) 
    setGeneric("name<-", function(object, value) standardGeneric("name<-"))
     
## general methods
if(!isGeneric("evaluate")) 
    setGeneric("evaluate", 
       function(object, estimator, ...) standardGeneric("evaluate"))
    
