simuler.spodt <-
function(d.sim, data, vqt, vql, weight, graft,
           level.max, min.parent, min.child, rtwo.min)
{
    data$z <- d.sim
    
    spodt.sim <- sp1(data, vqt, vql,
                       weight, graft,
                       level.max, min.parent, min.child, rtwo.min)

    return(spodt.sim@R2)
}
