# # redefition of print() for JunctionTree objects
# #' @rdname print.JunctionTree-methods
# #' @aliases print.Inference
# setMethod("print.JunctionTree",
#           "InferenceEngine",
#           function(x, ...)
#           {
#             str <- "\nJunction Tree "
#             str <- paste(str, "with ", sep = '')
#             str <- paste(str, x@num.nodes, sep = '')
#             str <- paste(str, " cliques", sep = '')
#             message(str)
#             
#             if (x@num.nodes > 0)
#             {          
#               str <- ""
#               clnames <- NULL
#               for (i in 1:x@num.nodes)
#               {
#                 clname <- ''
#                 clname <- paste(clname, '(', sep = '')
#                 clname <- paste(clname, paste(sort(x@cliques[[i]]), sep="", collapse=','))
#                 clname <- paste(clname, ')', sep = '')
#                 clnames[[i]] <- clname
#                 str <- paste(str, clname, sep = ' ')
#               }
#               colnames(x@junction.tree) <- clnames
#               rownames(x@junction.tree) <- clnames
#               
#               message(str, '\n\nAdjacency matrix:')
#               print(x@junction.tree)
#               message("")
#             }
#             
#           })
