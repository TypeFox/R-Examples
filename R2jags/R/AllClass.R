setOldClass("jags")
setOldClass("bugs")
setOldClass("mcmc.list")

setClass("rjags",
     representation(
            model = "jags",
            BUGSoutput = "bugs")
)

setClass("rjags.parallel",
     representation(
            BUGSoutput = "bugs"),
     contains = "rjags"
)
