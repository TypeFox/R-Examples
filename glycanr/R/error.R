# derived traits error handling
# Useful if some glycans are missing
dteh <- function(expr, mess=NULL){
    tryCatch(
             expr
             , error = function(e) {
                 if(grepl("object '[^ ]*' not found", e$message)){
                     if(is.null(mess)){
                         warning("not all derived traits could be calculated!")
                     } else {
                         warning(mess)
                     }
                 } else {
                     stop(e)
                 }
             }
                 
    )
}
