"blockTreeToList" <-
function (tree) 
{
    result <- list(tree$block)
    if (!is.null((tree$sub.blocks))) 
        for (i in 1:length(tree$sub.blocks)) result <- append(result, 
            blockTreeToList(tree$sub.blocks[[i]]))
    class(result) <- "dg.BlockList"
    return(result)
}
