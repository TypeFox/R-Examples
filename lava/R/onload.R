'.onLoad' <- function(libname, pkgname="lava") {
    addhook("heavytail.init.hook","init.hooks")
    addhook("glm.estimate.hook","estimate.hooks")
    addhook("cluster.post.hook","post.hooks")
    addhook("ordinal.sim.hook","sim.hooks")
    addhook("color.ordinal","color.hooks")

}

'.onAttach' <- function(libname, pkgname="lava") {
    desc <- utils::packageDescription(pkgname)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}
