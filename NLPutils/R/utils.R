.document_id <-
function(a)
{
    ## Determine the first non-missing id of a document annotation in a,
    ## or the next available id if there if there is no such id.
    id <- a$id[a$type == "document"]
    if(length(id))
        id <- id[!is.na(id)]
    if(length(id))
        id[1L]
    else
        next_id(a$id)
}
