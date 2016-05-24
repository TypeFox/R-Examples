gates.eff <-
function(dx, dy, marks, par=list(a=1, b=4, smark=1)) {
# Gates et al efficiency functions  [1 - (b R / s)^a]^(1/a) 
    with(as.list(par),
         (pmax(0, 1 - (b * sqrt(dx^2 + dy^2) / marks[[smark]])^a))^(1/a)
    )
}
