variable.names.hdlm <-
function(object, full = FALSE, ...)
{
    if(object$rank[[1]]) colnames(model.matrix(object))
    else character()
}

