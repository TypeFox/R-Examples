arrayInd <-
function (ind, .dim, .dimnames = NULL, useNames = FALSE) 
{
    rank = length(.dim)
    call = if (length(ind)) {
        .C("arrayInd", as.integer(ind - 1), as.integer(.dim), 
            length(ind), rank, integer(rank * length(ind)), package = "rje")[[5]]
    }
    else integer(0)
    out = matrix(call, nrow = length(ind), ncol = rank)
    if (useNames) 
        dimnames(out) = list(.dimnames[[1L]][ind], if (rank == 
            2L) c("row", "col") else paste0("dim", seq_len(rank)))
    out
}
