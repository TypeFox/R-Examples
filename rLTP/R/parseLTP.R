parseLTP = function(result)
{
    result = strsplit(result,'\n\n')[[1]]
    result = sapply(result,strsplit,'\n')
    names(result) = NULL
    result = unlist(result)
    result = strsplit(result,' ')
    result = unlist(result)
    result
}
