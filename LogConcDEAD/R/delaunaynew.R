##Note this has identical functionality to the function "delaunayn" in the package "geometry"
"delaunaynew"<-
function (p, options = "QJ") 
.Call("delaunayn", as.matrix(p), as.character(options),PACKAGE="LogConcDEAD")
