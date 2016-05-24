asBitwiseValue =
function(val, defValues, className, prefix = NA)
{
  if(is.character(val)) {
    i = pmatch(val, names(defValues))
    if(!is.na(prefix) && any(is.na(i))) {
      i[is.na(i)] = pmatch(paste(prefix, val[i], sep = ""), names(defValues))
    }
  } else {
        # check for a single scalar value
    i = match(val, defValues)
        # we didn't get a match, check for a combination
    if(any(is.na(i))) {
       tmp = sapply(defValues, function(x) bitAnd(val, x))
       i = which(tmp != 0)
       if(length(i) == 0)
          stop("can't match ", val, " to elements of bitwise enumeration: ", paste(names(defValues), collapse = ", "))
    }
  }

  if(any(is.na(i)))
    stop("unmatched name(s) ", val[is.na(i)])

  value = bitlist(defValues[i])

    # if className = NA, the we return the number. This is useful
    # to avoid the expense of creating a new instance only to
    # take the value to pass to C/C++ code.
  if(is.na(className))
    return(value)

    # Create a new instance of the given class.
  id = paste(names(defValues)[i], collapse = " | ")
  ans = new(className, value)
  names(ans) = id
  ans
}



BitwiseValue =
  #
  # create one or more bitwise values.
  #
function(val, name = names(val),
         class = if(is(val, "BitwiseValue")) class(val) else "BitwiseValue",
         asVector = TRUE, S3 = FALSE,
         defValues = getEnumValues(class, name))
{
  if(length(val) > 1 &&  !asVector) {
      lapply(seq(along = val), function(i) BitwiseValue(val[i], name[i], class))
  } else {
    if(FALSE && S3) {
      ans = structure(val, names = name, class = unique(c(class, "BitwiseValue")))
    } else {
      if(is(val, "character")) {
         if(missing(name))
           name = val

         i = match(val, names(defValues))
         if(any(is.na(i)))
           raiseEnumError(val, defValues, TRUE, FALSE, index = i)
         
         val = defValues[i]
      }

      if(!is.integer(val) && is.na(as.integer(val))) {
         # the value is too big for an integer in R.
        warning("BitwiseValue ", val, " doesn't fit into integer value")
      }
          
      ans =  new(class, val)
      
      names(ans) = name
    }
    ans
  }
}

getEnumValues =
function(class = NA, ids = character())
{
  if(!is.na(class)) {
     var = sprintf("%sValues", class)
     if(exists(var))
       return(get(var))
  }

  structure(sapply(ids, get), names = ids)
}


setAs("BitwiseValue", "numeric",
       function(from) {
          ans = bitlist(from)
          names(ans) = paste(names(from), collapse = " | ")
          ans
       })


setMethod("|", c("BitwiseValue", "BitwiseValue"),
          function(e1, e2) {
            if(class(e1) != class(e2) &&
                 !( is(e1, class(e2)) || is(e2, class(e1))))
              warning("OR'ing two BitwiseValue objects of different class")

            v = bitlist(e1, e2)
                 # if these are of different classes, we should decide what class the result
                 # should be.
            BitwiseValue(v, paste(names(e1), names(e2), sep = " | "), class = class(e1))
          })


setMethod("&", c("BitwiseValue", "BitwiseValue"),
          function(e1, e2) {
            if(class(e1) != class(e2) &&
                 !( is(e1, class(e2)) || is(e2, class(e1))))
              warning("OR'ing two BitwiseValue objects of different class")

            v = bitAnd(e1, e2)
                 # if these are of different classes, we should decide what class the result
                 # should be.
            BitwiseValue(v, paste(names(e1), names(e2), sep = " & "), class = class(e1))
          })




 # Should this return a vector with the same number of elements as x
 # a single value?
 # Version of 
 # Perhaps introduce a BitwiseValueVector to represent more than one element.
if(FALSE) {
setMethod("c", c("BitwiseValue"),
           function(x, ..., recursive = FALSE) {
             val = callNextMethod()
             val = val[ - grep("\\.\\.\\.[0-9]+", names(val)) ]
             names(val)[1] = gsub("^x\\.", "", names(val)[1])
             val = structure(unlist(val), names = names(val))
             BitwiseValue(val, class = class(x))
           })
}

setMethod("c", c("BitwiseValue"),
           function(x, ..., recursive = FALSE) {
             els = c(as(x, "numeric"), sapply(list(...), as, "numeric"))
             defValues = get(paste(class(x), "Values", sep = ""))
             return(asBitwiseValue(unlist(els), defValues, class(x)))
            })

