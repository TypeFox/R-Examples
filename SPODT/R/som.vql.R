som.vql <-
function(vql, data.som)
{
    return(by(data.som, vql, sum))
}
