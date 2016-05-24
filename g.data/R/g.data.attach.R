## Attach (or virtually create) a delayed-data package ("DDP"):
g.data.attach <- function(dir, pos=2, warn=TRUE, readonly=FALSE) {
    env <- attach(NULL, pos, basename(dir))
    attr(env, "path")     <- dir
    attr(env, "readonly") <- readonly
    if (!file.exists(dir)) {if (warn) warning("New DDP: ", dir); return(invisible())}
    for (fn in dir(dir, pattern="\\.RData$", all.files=TRUE, full.names=TRUE))
      eval(substitute(delayedAssign(OB, get(load(FN))), list(OB=g.data.unmash(fn), FN=fn)), env)
}
