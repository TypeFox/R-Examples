tass.inf <-
function(dx, dy, marks, par=list(b=3.52*0.975, c=6.1, smark=1)) {
# TASS influence function  s - c[exp(R/b) - 1]
    with(as.list(par),
        pmax(0, marks[[smark]] - c * (exp(sqrt(dx^2 + dy^2) / b) - 1))
    )
}
