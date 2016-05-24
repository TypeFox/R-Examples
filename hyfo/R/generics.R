##### Prepared for future use, when hyfo becomes a class.


# # change to hyfo
# setGeneric('as.hyfo', function(x) {
#   standardGeneric('as.hyfo')
# })
# 
# setMethod('as.hyfo', signature('list'),
#           function(x) {
#             
#             if (!is.null(x$Members)) {
#               hyfo <- new("hyfo.multiMember", varName = x$Variable$varName, xyCoords = x$xyCoords, Dates = x$Dates, Data = x$Data,
#                           Members = x$Members)
#             } else {
#               hyfo <- new("hyfo", varName = x$Variable$varName, xyCoords = x$xyCoords, Dates = x$Dates, Data = x$Data)
# 
#             }
#             return(hyfo)            
#             
#           })
# 
