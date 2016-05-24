## ---- echo = FALSE-------------------------------------------------------
library(titrationCurves)

## ------------------------------------------------------------------------
wb_sa(eqpt = TRUE, main = "Titration of WB w/ SA")
wb_sa(pka = 7, col = "blue", overlay = TRUE)

## ---- fig.keep='none'----------------------------------------------------
wb1 = wb_sa(pka = 8)
wb2 = wb_sa(pka = 6)
head(wb1)

## ------------------------------------------------------------------------
plot(wb1, ylim = c(0, 12), xlim = c(0, 80), type = "l", col = "blue",
     lwd = 2, xlab = "volume of titrant in mL")
lines(wb2, col = "green", lwd = 2)
abline(v = 50, col = "red", lty = 2)
legend(x = "topright", legend = c("pKa = 8", "pKa = 6"), 
       col = c("blue", "green"), lty = 1, lwd = 2)

## ------------------------------------------------------------------------
metal_edta(eqpt = TRUE)
metal_edta(logkf = 6, col = "blue", overlay = TRUE)

## ------------------------------------------------------------------------
redox_titration(eqpt = TRUE)
redox_titration(pot.analyte = 0.5, pot.titrant = 1.5, col = "blue", overlay = TRUE)

## ------------------------------------------------------------------------
ppt_mixture(eqpt = TRUE)
ppt_mixture(pksp1 = 12, pksp2 = 8, col = "blue", overlay = TRUE)

## ---- fig.width=6--------------------------------------------------------
wbd = derivative(wb1)

## ------------------------------------------------------------------------
str(wbd)

## ------------------------------------------------------------------------
plot(wbd$first_deriv, xlim = c(48, 52), col = "blue", type = "l", lwd = 2,
     xlab = "volume of titrant in mL", ylab = "first derivative")
abline(v = 50, col = "red", lty = 2)

## ------------------------------------------------------------------------
triwa_sb(conc.acid = 0.0400, conc.base = 0.120, pka1 = 3.128, pka2 = 4.761, 
         pka3 = 6.396, col = "blue", eqpt = TRUE)

## ------------------------------------------------------------------------
wa_sb(pka = 8, pkw = 20, col = "blue", eqpt = TRUE)
wa_sb(pka = 8, col = "green", overlay = TRUE)
legend(x = "topleft", legend = c("non-aqueous", "aqueous"), 
       col = c("blue", "green"), lty = 1, lwd = 2)

## ------------------------------------------------------------------------
metal_edta(col = "blue", eqpt = TRUE)
metal_edta(ph =  7, col = "green", overlay = TRUE)
legend(x = "topleft", legend = c("pH = 10", "pH = 7"), 
       col = c("blue", "green"), lty = 1, lwd =2)

## ------------------------------------------------------------------------
metal_edta(conc.metal = 0.0500, conc.edta = 0.025, vol.metal = 25.0, 
           alpha.metal = 0.00415, logkf = 18.80, col = "blue", eqpt = TRUE)
metal_edta(conc.metal = 0.0500, conc.edta = 0.0250, vol.metal = 25.0, 
           alpha.metal = 4.63e-10, logkf = 18.80, col = "green", overlay = TRUE)
legend(x = "topleft", 
       legend = c(expression(paste("0.0010 M N", H[3])), 
                  expression(paste("0.10 M N", H[3]))), 
       col = c("blue", "green"), lty = 1, lwd = 2)

## ------------------------------------------------------------------------
redox_titration(pot.analyte = 0.154, elec.analyte = 2, pot.titrant = 1.72, 
                col = "blue", eqpt = TRUE)

## ------------------------------------------------------------------------
redox_titration(pot.analyte = 0.771, pot.titrant = 1.51, elec.titrant = 5, 
                col = "black", eqpt = TRUE)
redox_titration(pot.analyte = 0.771, pot.titrant = 1.415, elec.titrant = 5, 
                col = "blue", overlay =TRUE)
redox_titration(pot.analyte = 0.771, pot.titrant = 1.321, elec.titrant = 5, 
                col = "green", overlay = TRUE)
legend(x = "topleft", legend = c("pH = 0", "pH = 1", "ph = 2"),
       col = c("black", "blue", "green"), lty = 1, lwd =  2)

## ---- fig.keep= 'none'---------------------------------------------------
p.a = ppt_analyte(eqpt = TRUE)
p.t = ppt_titrant(overlay = TRUE)

## ------------------------------------------------------------------------
plot(p.a, col = "blue", type = "l", lwd = 2, xlim = c(0,50), ylim = c(0,15), 
     xlab = "volume of titrant (mL)", ylab = "pAg or pI")
lines(p.t, col = "green", lwd = 2)
legend(x = "left", legend = c("pAg", "pI"), col = c("blue", "green"),
       lty = 1, lwd = 2)

## ------------------------------------------------------------------------
ppt_mixture(col = "blue", eqpt = TRUE)

