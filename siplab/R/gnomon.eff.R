gnomon.eff <-
function(dx, dy, marks, par=list(a=1, b=4, smark=1)) {
# Gnomonic efficiency functions  1 - b R^a /s
    with(as.list(par),
        pmax(0, 1 - b * (dx^2 + dy^2)^(a/2) / marks[[smark]])
    )
}
