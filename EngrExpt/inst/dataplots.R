library(EngrExpt)
#lattice.options(default.theme = standard.theme())
trellis.device(pdf, file = "dataplots.pdf", height = 6, width = 10)

ls.str("package:EngrExpt")

densityplot(~ absorb, absorb, xlab = "Average oil absorption of silica samples")
qqmath(~ absorb, absorb, ylab = "Average oil absorption of silica samples",
       xlab = "Standard normal quantiles", aspect = 1)
xyplot(adhesion ~ ph, adhesion, ylab = "Adhesion of a lens coating",
       xlab = "pH", type = c("g", "p", "a"))
dotplot(as.factor(ph) ~ adhesion, adhesion, ylab = "pH", type =
       c("p","a"), xlab = "Adhesion of a lens coating")
dotplot(as.factor(ph) ~ adhesion, adhesion2, groups = cat,
        type = c("p","a"), ylab = "pH", xlab = "Adhesion of a lens coating",
        auto.key = list(space = "right", lines = TRUE, title = "Catalyst"))
dotplot(as.factor(charge) ~ gloss|as.factor(flowrate), applicat, groups = distance,
        strip = FALSE, type = c("p", "a"), layout = c(1,2), ylab = "Charge",
        strip.left = TRUE, auto.key = list(space = "right", lines = TRUE, title = "Distance"))
qqmath(~ yield, assay, groups = process, auto.key = list(space = "right", title = "Process"),
       type = c("g", "p"), aspect = 1, xlab = "Standard normal quantiles")
dotplot(reorder(as.factor(type), lifetime) ~ lifetime, battery, xlab = "Lifetime", ylab = "Type", type = c("p", "a"))
densityplot(~ thickness, ccthickn)
qqmath(~ thickness, ccthickn, aspect = 1, xlab = "Standard normal quantiles", type = c("g","p"),
       panel=function(...){panel.qqmathline(..., alpha = 0.5, lty = 2); panel.qqmath(...)})
qqmath(~ time, cement, groups = type, aspect = 1, xlab = "Standard normal quantiles",
#       type = c("g","p"),
#       panel=function(...){panel.qqmathline(..., alpha = 0.5, lty = 2); panel.qqmath(...)},
       auto.key = list(space = "right", title = "Type"))
dotplot(as.factor(temp) ~ yield | as.factor(time), chemreac, groups = cat,
        strip = FALSE, strip.left = TRUE, type = c("p", "a"), layout = c(1,2),
        auto.key = list(space = "right", title = "Catalyst", lines = TRUE))

