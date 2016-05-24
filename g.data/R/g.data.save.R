## Save objects in position "pos" to a delayed-data package:
g.data.save <- function(dir=attr(env, "path"), obj=ls(env, all.names=TRUE), pos=2, rm.obj=NULL) {
    if (is.character(pos)) pos <- match(pos, search())
    if (is.na(pos)) stop("pos not found")
    env <- pos.to.env(pos)
    if (isTRUE(attr(env, "readonly"))) stop("Read-Only!")
    if (!file.exists(dir)) dir.create(dir)
    if (length(rm.obj)) {rm(list=rm.obj, pos=pos); file.remove(g.data.mash(dir, rm.obj))}
    is.promise <- function(i) is.call(eval(parse(text=paste("substitute(", i, ", env)"))))
    for (i in obj) if (!is.promise(i)) save(list=i, file=g.data.mash(dir, i), envir=env)
}
