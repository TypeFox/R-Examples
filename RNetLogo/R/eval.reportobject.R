eval.reportobject <- function(resobj) {
  type_ <- .jcall(resobj, "S", "getType")
  
  eval.nestedlist <- function(obj) {
    retlist <- lapply(obj, function(x) {eval.reportobject(x)})
    return (retlist)
  }
  
  value <- switch(type_, 
                  "Boolean" =     .jcall(resobj, "Z", "getResultAsBoolean"), 
                  "String" =      .jcall(resobj, "S", "getResultAsString"), 
                  "Double" =      .jcall(resobj, "D", "getResultAsDouble"),
                  "DoubleList" =  .jcall(resobj, "[D", "getResultAsDoubleArray"),
                  "BoolList" =    .jcall(resobj, "[Z", "getResultAsBooleanArray"),
                  "StringList" =  .jcall(resobj, "[S", "getResultAsStringArray"),
                  "NestedList" =  eval.nestedlist(.jcall(resobj, "[Ljava/lang/Object;", "getResultAsObjectArray"))
                  )
  
  if (is.null(value)) {
    print("unknown datatype")    
  }
  
  return (value)
}