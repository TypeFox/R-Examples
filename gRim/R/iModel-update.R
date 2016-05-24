
update.iModel <- function(object, items, fit=TRUE, details=0, ...){

##   cat("update.iModel\n")
##   print(class(object))
  
  glist           <- modify_glist(object$glist, items)
  object$glist    <- glist

  switch(class(object)[1],
         "dModel"={
           upd <- .dModel_finalize(glist, object$varNames)    
           object[names(upd)] <- upd
         },
         "cModel"={
           upd <- .cModel_finalize(glist, object$varNames)    
           object[names(upd)] <- upd
         },
         "mModel"={
           upd <- .mModel_finalize(glist, object$varNames, object$datainfo)    
           object[names(upd)] <- upd
         }
         )
  
  if (fit){    
    object <- fit(object)
  }  
  object
}


































































