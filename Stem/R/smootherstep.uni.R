`smootherstep.uni` <-
function(m,C,Gmatx,Wtx,mx,Cx) {
    Rx <- Gmatx * C * Gmatx + Wtx# Gmatx = G_{t+1}
    B  <- C * Gmatx / Rx
    ms <- m + B*(mx - Gmatx *m) # mx=ms_{t+1}
    Cs <- C + B*(Cx-Rx)*B
    list(ms=ms,Cs=Cs)
  }

