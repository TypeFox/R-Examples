random.cluster.C <-
function(info.to.C)
{
  .Call("random_cluster_C",
        as.numeric(as.double(info.to.C))
        )
}
