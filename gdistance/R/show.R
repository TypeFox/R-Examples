setMethod ("show" , "TransitionLayer",
           function(object) {
             callNextMethod(object)
             cat("values      :", matrixValues(object), "\n")
             cat("matrix class:", class(transitionMatrix(object)), "\n")
           }
           )

setMethod ("show" , "TransitionStack",
           function(object) {
             callNextMethod(object)
             cat("nlayers      :", nlayers(object), "\n")              
             #show something about the layers if layers <=5
             #cat("values      :", matrixValues(x), "\n")
             #cat("matrix class:", class(transitionMatrix(x)), "\n")
           }
           )

