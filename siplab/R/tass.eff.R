tass.eff <-
function(dx, dy, marks, par=list(b=3.52*0.975, c=6.1, smark=1)) {
# TASS-like efficiency function  1 - c[exp(R/b) - 1] / s
    with(as.list(par),
        pmax(0, 1 - c * (exp(sqrt(dx^2 + dy^2) / b) - 1) / marks[[smark]])
    )
}
