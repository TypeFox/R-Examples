links.as.edgelist <-
function (links) 
{
    cbind(rep(0:(length(links) - 1), sapply(links, length)), 
        unlist(links))
}
