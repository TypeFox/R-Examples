emptyNamedList = structure(list(), names = character())

trim =
function (x) 
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)

dQuote =
function(x)
{
  if(length(x) == 0)
     character(0)
  else
      paste('"', x, '"', sep = "")
}

isContainer =
function(x, asIs, .level)
  (is.na(asIs) && .level == 1) || (!is.na(asIs) && asIs) ||
       .level == 1L || length(x) > 1  || length(names(x)) > 0 || is.list(x)

setGeneric("toJSON",
  function(x, container = isContainer(x, asIs, .level),
           collapse = "\n", ...,
           .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0,
            .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity")   {

  container; .withNames  # force these values.
  
  ans <- standardGeneric("toJSON")

  if(pretty)
     jsonPretty(ans)
  else
     ans
  })



setMethod("toJSON", "NULL",
           function(x, container = isContainer(x, asIs, .level), collapse = "\n", ...,
                     .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null",
                      .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {
             if(container) "[ null ] " else "null"
           })

setMethod("toJSON", "array",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", ...,
                       .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0,
                      .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {

             d = dim(x)
             txt = apply(x, length(d), toJSON, collapse = collapse, ..., .level = .level + 1L, .withNames = .withNames, .na = .na, .escapeEscapes = .escapeEscapes, pretty = pretty, asIs = asIs, .inf = .inf)
              paste(c("[", paste(txt, collapse = ", "), "]"), collapse = collapse)
           })


setMethod("toJSON", "function",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", ...,
                       .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0,
                      .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {
             toJSON(paste(deparse(x), collapse = collapse), container, collapse, ..., .level = .level, .withNames = .withNames, .na = .na, .escapeEscapes = .escapeEscapes, pretty = pretty, asIs = asIs, .inf = .inf)
           })



setMethod("toJSON", "ANY",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", ...,
                       .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0,
                      .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {

             if(isS4(x)) {
               paste("{", paste(dQuote(slotNames(x)), sapply(slotNames(x),
                                                              function(id)
                                                                 toJSON(slot(x, id), ..., .level = .level + 1L,
                                                                        .na = .na, .escapeEscapes = .escapeEscapes, asIs = asIs,
                                                                         .inf = .inf)),
                                 sep = ": ", collapse = ","),
                     "}", collapse = collapse)
             } else {
#cat(class(x), "\n")
               if(is.language(x)) {
                  return(toJSON(as.list(x), asIs = asIs, .inf = .inf, .na = .na))
                  stop("No method for converting ", class(x), " to JSON")
               }

               toJSON(unclass(x), container, collapse, ..., .level = .level + 1L,
                         .withNames = .withNames, .na = .na, .escapeEscapes = .escapeEscapes, asIs = asIs)
               # stop("No method for converting ", class(x), " to JSON")
             }
             
           })



setMethod("toJSON", "integer",
           function(x, container =  isContainer(x, asIs, .level),
                      collapse = "\n  ", ..., .level = 1L,
                      .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null",
                       .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity")
          {
            if(any(is.infinite(x)))
              warning("non-fininte values in integer vector may not be approriately represented in JSON")
            
             if(any(nas <- is.na(x)))
                 x[nas] = .na

             if(container) {
                if(.withNames)
                   paste(sprintf("{%s", collapse), paste(dQuote(names(x)), x, sep = ": ", collapse = sprintf(",%s", collapse)), sprintf("%s}", collapse))
                else
                   paste("[", paste(x, collapse = ", "), "]")
              } else
                 as.character(x)               
           })

setOldClass("hexmode")

setMethod("toJSON", "hexmode",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n   ", ...,
                     .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0,
                       .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {

             
             tmp = paste("0x", format(x), sep = "")
             if(any(nas <- is.na(x)))
                 tmp[nas] = .na             
             
             if(container) {
                if(.withNames)
                   paste(sprintf("{%s", collapse), paste(dQuote(names(x)), tmp, sep = ": ", collapse = sprintf(",%s", collapse)), sprintf("%s}", collapse))
                else               
                paste("[", paste(tmp, collapse = ", "), "]")
             } else
                tmp
           })


setMethod("toJSON", "factor",
           function(x, container =  isContainer(x, asIs, .level),
                     collapse = "\n", ..., .level = 1L,
                     .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null", pretty = FALSE, asIs = NA, .inf = " Infinity") {

             toJSON(as.character(x), container, collapse, ..., .level = .level, .na = .na, .escapeEscapes = .escapeEscapes, asIs = asIs)
           })

setMethod("toJSON", "logical",
           function(x, container =  isContainer(x, asIs, .level),
                     collapse = "\n", ..., .level = 1L,
                     .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null", pretty = FALSE, asIs = NA, .inf = " Infinity") {
             tmp = ifelse(x, "true", "false")
             if(any(nas <- is.na(tmp)))
                 tmp[nas] = .na             

             if(container) {
                if(.withNames)
                   paste(sprintf("{%s", collapse), paste(dQuote(names(x)), tmp, sep = ": ", collapse = sprintf(",%s", collapse)), sprintf("%s}", collapse))
                else               
                   paste("[", paste(tmp, collapse = ", "), "]")
             } else
                tmp
           })

setMethod("toJSON", "numeric",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", digits = 5, ...,
                      .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0,
                        .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA,
                         .inf = " Infinity") {

            if(any(is.infinite(x)))
              warning("non-fininte values in numeric vector may not be approriately represented in JSON")             

             tmp = formatC(x, digits = digits)
              
             if(any(nas <- is.na(x)))
                 tmp[nas] = .na
             if(any(inf <- is.infinite(x))) 
                 tmp[inf] = sprintf(" %s%s", ifelse(x[inf] < 0, "-", ""), .inf)
             
             if(container) {
                if(.withNames)
                   paste(sprintf("{%s", collapse),
                         paste(dQuote(names(x)), tmp, sep = ": ", collapse = sprintf(",%s", collapse)),
                         sprintf("%s}", collapse))
                else
                   paste("[", paste(tmp, collapse = ", "), "]")
             } else
               tmp
           })


setMethod("toJSON", "character",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", digits = 5, ..., .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {
# Don't do this: !             tmp = gsub("\\\n", "\\\\n", x)

#             if(length(x) == 0)    return("[ ]")

             tmp = x

             tmp = gsub('(\\\\)', '\\1\\1', tmp)
             if(.escapeEscapes) {
               tmp = gsub("\\t", "\\\\t", tmp)
               tmp = gsub("\\n", "\\\\n", tmp)
               tmp = gsub("\b", "\\\\b", tmp)               
               tmp = gsub("\\r", "\\\\r", tmp)
               tmp = gsub("\\f", "\\\\f", tmp)                              
             }
             tmp = gsub('"', '\\\\"', tmp)
             tmp = dQuote(tmp)

             if(any(nas <- is.na(x)))
                 tmp[nas] = .na                          
             
             if(container) {
                if(.withNames)
                   paste(sprintf("{%s", collapse),
                          paste(dQuote(names(x)), tmp, sep = ": ", collapse = sprintf(",%s", collapse)),
                         sprintf("%s}", collapse))
                else               
                   paste("[", paste(tmp, collapse = ", "), "]")
             } else if(length(x) == 0)
                      "[ ]"
               else
                      tmp
           })



# Symbols.
# names can't be NA
setMethod("toJSON", "name",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", ..., .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {
             as.character(x)
           })

setMethod("toJSON", "name",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", ..., .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA) {
             
               sprintf('"%s"', as.character(x))
           })


setOldClass("AsIs")
setMethod("toJSON", "AsIs",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", ..., .level=1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null", .escapeEscapes = TRUE, asIs = NA, .inf = " Infinity") {
              toJSON(structure(x, class = class(x)[-1]), container = TRUE, collapse = collapse, ..., .level = .level + 1L, .withNames = .withNames, .na = .na, .escapeEscapes = .escapeEscapes, asIs = asIs, .inf = .inf)
           })



setMethod("toJSON", "matrix",
           function(x, container =  isContainer(x, asIs, .level), collapse = "\n", ...,
                    .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity")
           {
             tmp = paste(apply(x, 1, toJSON, .na = .na, ..., .escapeEscapes = .escapeEscapes, .inf = .inf), collapse = sprintf(",%s", collapse))
             if(!container)
               return(tmp)

              if(.withNames)
                paste("{", paste(dQuote(names(x)), tmp, sep = ": "), "}")                
              else
                paste("[", tmp, "]")
           })

setMethod("toJSON", "list",
           function(x, container = isContainer(x, asIs, .level), collapse = "\n", ..., .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null", .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {
                # Degenerate case.

             if(length(x) == 0) {
                          # x = structure(list(), names = character()) gives {}
                return(if(is.null(names(x))) "[]" else "{}")
             }

             els = lapply(x, toJSON, ..., .level = .level + 1L, .na = .na, .escapeEscapes = .escapeEscapes, asIs = asIs, .inf = .inf)

             if(all(sapply(els, is.name)))
               names(els) = NULL
             if(missing(container) && is.na(asIs))
                container = TRUE
               
             if(!container)
               return(els)

#             els = unlist(els)
   w = sapply(els, length) == 0
   els[w] = "[]" # or "" or "null"
             
             if(.withNames)
                paste(sprintf("{%s", collapse),
                      paste(dQuote(names(x)), els, sep = ": ", collapse = sprintf(",%s", collapse)),
                      sprintf("%s}", collapse))
             else
                 paste(sprintf("[%s", collapse), paste(els, collapse = sprintf(",%s", collapse)), sprintf("%s]", collapse))
           })


if(TRUE)
setMethod("toJSON", "data.frame",
           function(x, container = isContainer(x, asIs, .level), collapse = "\n", ..., byrow = FALSE, colNames = FALSE, .level = 1L,
                    .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null",
                    .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity")
 {
  if(byrow) {
     tmp = lapply(1:nrow(x), function(i) { row = as.list(x[i, ])
                                           if(colNames)
                                             row                                             
                                           else
                                             unname(row)
                                         })
     toJSON(tmp, container, collapse, ..., .level = .level, .withNames = FALSE,
            .na = .na, .escapeEscapes = .escapeEscapes, pretty = pretty, asIs = asIs, .inf = .inf)
  } else
              # Data frame columns should always be treated as containers. From Joe Cheng
      toJSON(lapply(as.list(x),
                    function(col) 
                       I(col)
                   ),
              container = container, collapse = collapse, ...,
              .level = .level, .withNames = .withNames, .na = .na,
              .escapeEscapes = .escapeEscapes, pretty = pretty, asIs = asIs, .inf = .inf)
           })


jsonPretty =
function(txt)
{
   txt = paste(as.character(txt), collapse = "\n")
   enc = mapEncoding(Encoding(txt))
   .Call("R_jsonPrettyPrint", txt, enc, package = "RJSONIO")
}




setMethod("toJSON", "environment",
           function(x, container = isContainer(x, asIs, .level), collapse = "\n  ", ...,
                     .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null",
                      .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {
      toJSON(as.list(x), container, collapse, .level = .level, .withNames = .withNames, .escapeEscapes = .escapeEscapes, asIs = asIs, .inf = .inf)
           })


setMethod("toJSON", "function",
           function(x, container = isContainer(x, asIs, .level), collapse = "\n  ", ...,
                     .level = 1L, .withNames = length(x) > 0 && length(names(x)) > 0, .na = "null",
                      .escapeEscapes = TRUE, pretty = FALSE, asIs = NA, .inf = " Infinity") {
             warning("converting an R function to JSON as null. To change this, define a method for toJSON() for a 'function' object.")
             "null"
           })
