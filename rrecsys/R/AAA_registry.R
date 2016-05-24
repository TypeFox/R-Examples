# setup registry for recommender methods
rrecsysRegistry <- registry(registry_class = "rrecsys_registry", entry_class = "rrecsys_entries")

rrecsysRegistry$set_field("alg", type = "character", is_key = TRUE, index_FUN = match_partial_ignorecase)

rrecsysRegistry$set_field("fun", type = "function")

rrecsysRegistry$set_field("description", type = "character")

rrecsysRegistry$set_field("reference", type = "character")

rrecsysRegistry$set_field("parameters", type = "list")



print.rrecsys_registry <- function(x, ...) with(x, {
    writeLines(sprintf("Regirsty defined for rrecsys with %i entries as follows:\n", length(x)))
    a <- x$get_field_entries("alg")
    for (i in 1:length(a)) {
        writeLines(sprintf("Algorithm: %s", x$get_entry(alg = a[i])$alg))
        writeLines(sprintf("Reference: %s", x$get_entry(alg = a[i])$reference))
        
        if (is.list(x$get_entry(alg = a[i])$parameters)) {
            writeLines("Parameter and default values:")
            print(as.data.frame(x$get_entry(alg = a[i])$parameters))
        } else {
            writeLines("No parameter.")
        }
        writeLines("")
    }
    
})



print.rrecsys_entries <- function(x, ...) with(x, {
    writeLines(sprintf("Algorithm: %s", description))
    writeLines(sprintf("Reference: %s", reference))
  
  
    if (is.list(parameters)) {
        writeLines("Parameter and default values:")
        print(as.data.frame(parameters))
    } else writeLines("No parameter.")
}) 
