`MOD2VEC` <-
function(MOD)
  {
    VEC = NULL
    for(i in 1:length(MOD$MOD))
      {
        VEC = c(VEC, as.vector(MOD$MOD[[i]]))
      }

    attr(VEC, 'x') = MOD$x
    attr(VEC, 'y') = MOD$y
    attr(VEC, 'D') = MOD$D

    invisible(VEC)
  }

