gnomon.inf <-
function(dx, dy, marks, par=list(a=1, b=4, smark=1)) {
# Gnomonic influence functions  s - b R^a
    with(as.list(par),
        pmax(0, marks[[smark]] - b * (dx^2 + dy^2)^(a/2))
    )
}
