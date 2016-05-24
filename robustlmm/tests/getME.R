## testing getME function
require(robustlmm)
formatNum <- function(x, ...)
     format(round(x, 8), trim = TRUE, drop0trailing = TRUE, ...)
options(str = strOptions(formatNum = formatNum))

sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])
rfm <- rlmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
             rho.e = smoothPsi, rho.b=smoothPsi,
             rho.sigma.e = psi2propII(smoothPsi, k=2.28),
             rho.sigma.b = psi2propII(smoothPsi, k=2.28),
             doFit=FALSE)

(nmME <- eval(formals(robustlmm:::getME.rlmerMod)$name))
for (nm in nmME) {
    cat("\nName:", nm, "\n")
    str(getME(rfm, name=nm))
}
str(g.all <- getME(rfm, "ALL"),
    max.level = 2)
