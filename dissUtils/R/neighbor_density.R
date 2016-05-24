neighbor.density <- function(neigh.dists, D, k, N){

    return(ifelse(neigh.dists > 0,
                  exp(log(k)
                      - ( log(N)
                         + log(unit.hypersphere.volume(D))
                         + D * log(neigh.dists)
                         )
                      ),
                  0)
           );
}
