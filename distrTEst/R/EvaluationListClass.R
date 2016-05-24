
EvaluationList <- function(..., name0 = "a list of \"Evaluation\" objects")
        return(new("EvaluationList", Elist=list(...), name=name0))


setMethod("Elist", "EvaluationList", function(object) object@Elist)
setMethod("name", "EvaluationList", function(object) object@name)
setMethod("Data", "EvaluationList", function(object) (object@Elist[[1]])@Data)

setReplaceMethod("name", "EvaluationList", 
   function(object, value){ object@name <- value; object}) ### new 1.8

