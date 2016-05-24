zoi.inf <-
function(dx, dy, marks, par=list(k=0.2, smark=1)) {
# ZOI influence function  1 if R < k s, 0 otherwise
    with(as.list(par),
         ifelse(dx^2 + dy^2 < (k * marks[[smark]])^2,
                rep(1, length(dx)), rep(0, length(dx)))
    )
}
