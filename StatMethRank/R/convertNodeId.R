convertNodeId <- function(id)
{
    id = as.numeric(id)
    tmp = id
    ans = 1
    while (tmp > 0.001)
    {
        tmp = tmp * 10
        LR = round(tmp)
        ans = ans * 2 + LR - 1
        tmp = tmp - LR
    }
    return(ans)
}