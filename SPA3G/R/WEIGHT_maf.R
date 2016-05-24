WEIGHT_maf <-
function(G)
{
qs <- apply(G, 1, sum)/nrow(G)
return(1/sqrt(qs))
}
