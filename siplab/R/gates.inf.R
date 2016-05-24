gates.inf <-
function(dx, dy, marks, par=list(a=1, b=4, smark=1)) {
# Gates et al influence functions  [(s/b)^a - R^a]^(1/a) 
    with(as.list(par),
         (pmax(0, (marks[[smark]] / b)^a - (dx^2 + dy^2)^(a/2)))^(1/a)
    )
}
