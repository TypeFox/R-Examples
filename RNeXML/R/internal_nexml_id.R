
nexml_env = new.env(hash=TRUE)

# If no prefix is given, will use a UUID
# Generates an id number by appending a counter to the prefix
# Will keep track of the counter for each prefix for that session.  
#' @import uuid
nexml_id <- function(prefix = "",
                     use_uuid = getOption("uuid", FALSE)){
  if(use_uuid){
      uid <- paste0("uuid-", UUIDgenerate())
  } else {
    if((prefix %in% ls(envir=nexml_env)))
      id_counter <- get(prefix, envir=nexml_env)
    else {
      assign(prefix, 1, envir=nexml_env)
      id_counter <- 1
    } 

    uid <- paste0(prefix, id_counter)
    id_counter <- id_counter + 1
    assign(prefix, id_counter, envir=nexml_env)
  }
  uid
}

#' reset id counter
#' 
#' reset the id counter
#' @export
reset_id_counter <- function(){
  rm(list=ls(envir=nexml_env), envir=nexml_env)
}

# use an environment to store counter 
