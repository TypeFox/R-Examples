filter.words <-
function (documents, to.remove) 
{
    lapply(documents, function(x) x[, !(x[1, ] %in% to.remove), 
        drop = FALSE])
}
