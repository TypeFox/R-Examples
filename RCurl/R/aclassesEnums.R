# EnumValue is typically for individual values.
# Bitwise can support a vector.
# These are for values. The definition is separate.
setClass("SymbolicConstant", representation = c(names = "character"), contains = "integer")
setClass("EnumValue", contains = "SymbolicConstant", prototype = as.integer(NA))
setClass("BitwiseValue", contains = "SymbolicConstant", prototype = as.integer(NA))



##################################
# IGNORE
# These 3 and the next classes look like versions of the same idea.
if(FALSE) {
      setClass("SymbolicConstantsDefinition", representation(name = "character"))
      setClass("EnumValueDefinition", contains = c("EnumValue", "SymbolicConstantsDefinition"))
      setClass("BitwiseValueDefinition", contains = c("BitwiseValue", "SymbolicConstantsDefinition"))


      BitwiseValueDefinition =
        function(values, name = names(values), class = "BitwiseValueDefinition")
          {
            ans = new("BitwiseValueDefinition", wxStretchValues, name = "wxStretch")
            names(ans) = name
            ans
          }
      setMethod("[", "BitwiseValueDefinition",
                function(x, i, j, ..., drop = TRUE) {
                  k = class(x[[1]])
                  ans = callNextMethod()
                  BitwiseValue(ans, class = k)
                })
}
# IGNORE
##################################


setClass("EnumDef",
          representation(EnumName = "character"),
          contains = "integer")

EnumDef =
function(name, values, symbolicNames = names(values))
{
  values = as(values, "integer")
  ans = new("EnumDef", structure(values, names = symbolicNames), EnumName = name)
  names(ans) = symbolicNames
  ans
}


 # General types that can be used in an S4 method signature  in setMethod
 # to cover an enumeration or bitwise in its various forms of specification,
 #  i.e. as a number, a string or as an actual EnumValue or BitwiseValue
setClassUnion("EnumerationValue", c("numeric", "integer", "character", "EnumValue"))



# Display with the name.
setMethod("show", "SymbolicConstant",
           function(object)
                   # avoid the bitlist coercion
              show(structure(as(unclass(object), "numeric"), names = names(object)))
          )

# display as a 1 row matrix with the name of the enum type as the row name.
# Do we need this or can we have the generic SymbolicConstant.
tmp = function(object)
         show(matrix(as(object, "integer"), 1, , dimnames = list(paste(object@EnumName, ":", sep = ""), names(object))))

#setMethod("show", "EnumValue", tmp) # do we want this.
setMethod("show", "EnumDef", tmp)


setMethod("show", "SymbolicConstant", function(object)
             show(matrix(as(object, "integer"), 1, , dimnames = list(class(object), names(object)))))



setMethod("[", "EnumDef",
            function(x, i, j, ..., drop = TRUE) {
              vals = get(paste(x@EnumName, "Values", sep = ""))
              asEnumValue( unclass(x)[i], vals)
            })


makeSymbolicVariables =
  #
  # For element in the def, create a corresponding R variable with that name
  #  which contains that value.
  # e.g.
  #    c(a = 1, b = 2)
  # for a class MyEnum, we would end up with
  #   variables named a and b with values 1 and 2 respectively
  # and each would be of class MyEnum.
  #
  # The target class must have been defined before this.
  #
function(def, className = class(def), where = globalenv())
{
  invisible(sapply(names(def),
          function(i) {
            if(is(def, "BitwiseValue"))
              el = BitwiseValue(def[i], i, class =  className)
            else
              el = def[i]
            assign(i, el, where)
          }))
}  

cumBitOr = bitlist =
function(...)
{
  x = unlist(list(...))
  if(length(x) == 1)
    return(x)

  ans = x[1]
  for(i in 2:length(x)) {
    ans = bitOr(ans, x[i])
  }
  ans
}




asEnumValue =
  #
  # if fixCloseMatches is TRUE, we continue on if we can find a match
  # for all of the possible values that were specified slightly incorrectly.
  #
function(val, values, class = values@EnumName, fromString = NA,
         fixCloseMatches = TRUE, prefix = character(), S3 = is.null(getClassDef(class)))
{
   # handle multiple entries.
  if(length(val) > 1) {
    tmp = sapply(val, asEnumValue, values, class, fromString, fixCloseMatches, prefix, USE.NAMES = FALSE)
      # if we have multiple values and they relate to a BitwiseValue enumeration,
      # collapse them into a single value.
          # class was getClass(class)
    if(extends(class, "BitwiseValue")) #XXX augment for S3
       return(bitlist(tmp))
    if(S3)
      class(tmp) = c(class, "EnumValue")
    else
      tmp = as(tmp, class)
    return(tmp)
  }
  
  if(is.na(fromString))
     fromString = is(val, "character")

  if(fromString) {
     i = pmatch(val, names(values))  # allowing pmatch, but should type it explicitly in code.
       # deal with lowercase matches for covenience
     if(is.na(i))
       i = pmatch(val, tolower(names(values)))

       # and if still not there, remove the prefix.
     if(is.na(i) && length(prefix))
       i = pmatch(val, tolower(gsub(paste("^", prefix, sep = ""), "",  names(values))))
  } else
     i = match(val, values)

  if(any(is.na(i))) {
      i = raiseEnumError(val, values, fromString, fixCloseMatches, index = i)
  }

  if(S3) {
    ans = structure(unclass(values)[i], names = names(values)[i])
    class(ans) =  c(class, "EnumValue")
    return(ans)
  } else {
    ans = new(class, unclass(values)[i])
    names(ans) = names(values)[i]
  }
  ans
}  


GenericEnumValue =
function(name, val, class = "EnumValue", S3 = FALSE)
{
  if(S3) {
    ans = structure(val, names = name, class = unique(c(class, "EnumValue")))
  } else {
    ans = new(class, val)
    names(ans) = name
  }
  
  ans
}  


raiseEnumError =
function(val, values, fromString = is(fromString, "character"), fixCloseMatches = TRUE,
          index = match(val, if(fromString) names(values) else values))
{
        # see if we can find values that were close to the ones the user gave us incorrectly.
    if(fromString) {
      possibles = names(values)[m <- agrep(val[is.na(index)], names(values))]
    } else {
      possibles = values[ m <- agrep(as.character(val), as.character(values)) ] 
    }
    if(length(possibles)) {
      txt = paste("\n\tPerhaps you meant",  if(length(possibles) > 1) "one of", paste(possibles, collapse = ", "))
      txt = paste("No such value(s) ", val[is.na(index)], " in ", paste(names(values), collapse = ", "), txt, sep = "")

      msg = list(message = txt, call = NULL,
                 possibleValues = possibles,
                 class = class)

      if(fixCloseMatches && all(!is.na(m))) {
        class(msg) =  c("EnumCoercionWarning", "warning", "condition")                
        warning(msg)
        index[is.na(index)] = m
        return(index)
      } else {
        class(msg) =  c("EnumCoercionError", "error", "condition")
        stop(msg)
      }
    }
    else 
      stop("No such value(s) ", val[is.na(index)], " in ", paste(names(values), collapse = ", "))
}
