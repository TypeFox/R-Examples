'.onAttach' <- function(lib, pkg="lava.tobit") {
  addhook("lava.tobit.estimate.hook","estimate.hooks")
  addhook("lava.tobit.color.hook","color.hooks")
  addhook("lava.tobit.sim.hook","sim.hooks")
  addhook("lava.tobit.init.hook","init.hooks")

  lava.options(tobitAlgorithm=mvtnorm::GenzBretz(abseps=1e-5),
               tobitseed=1, threshold=1)
  desc <- utils::packageDescription(pkg)
  packageStartupMessage("\nLoading '", desc$Package, "' package...\n",
                        "Version    : ", desc$Version, "\n")
}
