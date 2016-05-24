### R code from vignette source 'marelac.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("marelac")
options(prompt = "> ")
options(width=75)


###################################################
### code chunk number 2: marelac.Rnw:133-137
###################################################
data.frame(cbind(acronym = names(Constants),
             matrix(ncol = 3, byrow = TRUE, data = unlist(Constants),
             dimnames=list(NULL, c("value", "units", "description")))))



###################################################
### code chunk number 3: marelac.Rnw:145-148
###################################################
data.frame(cbind(acronym = names(Oceans),
             matrix(ncol = 3, byrow = TRUE, data = unlist(Oceans),
             dimnames = list(NULL, c("value", "units", "description")))))


###################################################
### code chunk number 4: marelac.Rnw:183-187
###################################################
SURF <- outer(X = Bathymetry$x,
              Y = Bathymetry$y,
              FUN <- function(X, Y) earth_surf(Y, X))
sum(SURF)


###################################################
### code chunk number 5: marelac.Rnw:192-195
###################################################
 sum(SURF*(Bathymetry$z < 0))

- sum(SURF*Bathymetry$z*(Bathymetry$z < 0))


###################################################
### code chunk number 6: marelac.Rnw:202-210
###################################################
SurfDepth <- vector()

dseq <- seq(-7500, -250, by = 250)

for (i in 2:length(dseq)) {
  ii <- which (Bathymetry$z > dseq[i-1] & Bathymetry$z <= dseq[i])
  SurfDepth[i-1]<-sum(SURF[ii])
}


###################################################
### code chunk number 7: ocean
###################################################
plot(dseq[-1], SurfDepth, xlab="depth, m", log = "y",
     ylab = "m2", main = "Surface at ocean depths")


###################################################
### code chunk number 8: ocean
###################################################
plot(dseq[-1], SurfDepth, xlab="depth, m", log = "y",
     ylab = "m2", main = "Surface at ocean depths")


###################################################
### code chunk number 9: marelac.Rnw:239-242
###################################################
AtomicWeight
AtomicWeight[8, ]
(W_H2O<- with (atomicweight, 2 * H + O))


###################################################
### code chunk number 10: marelac.Rnw:250-253
###################################################
atmComp("O2")
atmComp()
sum(atmComp())    #!


###################################################
### code chunk number 11: cor
###################################################
plot(-90:90, coriolis(-90:90), xlab = "latitude, dg North",
     ylab = "/s" , main = "Coriolis factor", type = "l", lwd = 2)


###################################################
### code chunk number 12: figcor
###################################################
plot(-90:90, coriolis(-90:90), xlab = "latitude, dg North",
     ylab = "/s" , main = "Coriolis factor", type = "l", lwd = 2)


###################################################
### code chunk number 13: marelac.Rnw:288-290
###################################################
diffcoeff(S = 15, t = 15)*24*3600*1e4  # cm2/day
diffcoeff(t = 10)$O2


###################################################
### code chunk number 14: marelac.Rnw:294-295
###################################################
difftemp <- diffcoeff(t = 0:30)[ ,1:13]


###################################################
### code chunk number 15: diff
###################################################
matplot(0:30, difftemp, xlab = "temperature", ylab = " m2/sec",
        main = "Molecular/ionic diffusion", type = "l")
legend("topleft", ncol = 2, cex = 0.8, title = "mean", col = 1:13, lty = 1:13, 
   legend = cbind(names(difftemp), format(colMeans(difftemp),digits = 4)))


###################################################
### code chunk number 16: figdif
###################################################
matplot(0:30, difftemp, xlab = "temperature", ylab = " m2/sec",
        main = "Molecular/ionic diffusion", type = "l")
legend("topleft", ncol = 2, cex = 0.8, title = "mean", col = 1:13, lty = 1:13, 
   legend = cbind(names(difftemp), format(colMeans(difftemp),digits = 4)))


###################################################
### code chunk number 17: visc
###################################################
plot(0:30, viscosity(S = 35, t = 0:30, P = 1), xlab = "temperature",
      ylab = "g/m/s", main = "shear viscosity of water",type = "l")
lines(0:30, viscosity(S = 0, t = 0:30, P = 1), col = "red")
lines(0:30, viscosity(S = 35, t = 0:30, P = 100), col = "blue")
legend("topright", col = c("black", "red", "blue"), lty = 1,
        legend = c("S=35, P=1", "S=0, P=1", "S=35, P=100"))


###################################################
### code chunk number 18: figshear
###################################################
plot(0:30, viscosity(S = 35, t = 0:30, P = 1), xlab = "temperature",
      ylab = "g/m/s", main = "shear viscosity of water",type = "l")
lines(0:30, viscosity(S = 0, t = 0:30, P = 1), col = "red")
lines(0:30, viscosity(S = 35, t = 0:30, P = 100), col = "blue")
legend("topright", col = c("black", "red", "blue"), lty = 1,
        legend = c("S=35, P=1", "S=0, P=1", "S=35, P=100"))


###################################################
### code chunk number 19: marelac.Rnw:348-350
###################################################
gas_O2sat(t = 20)
t <- seq(0, 30, 0.1)


###################################################
### code chunk number 20: marelac.Rnw:355-356
###################################################
gas_O2sat(S=35, t=20)*1000/molweight("O2")


###################################################
### code chunk number 21: o2sat
###################################################
plot(t, gas_O2sat(t = t), type = "l", ylim = c(0, 15), lwd = 2,
     main = "Oxygen saturation", ylab = "mg/l", xlab = "temperature")
lines(t, gas_O2sat(S = 0, t = t, method = "Weiss"), col = "green",
      lwd = 2, lty = "dashed")
lines(t, gas_O2sat(S = 35, t = t, method = "Weiss"), col = "red", lwd = 2)
legend("topright", c("S=35", "S=0"), col = c("red","green"), 
       lty = c(1, 2), lwd = 2)


###################################################
### code chunk number 22: figo2sat
###################################################
plot(t, gas_O2sat(t = t), type = "l", ylim = c(0, 15), lwd = 2,
     main = "Oxygen saturation", ylab = "mg/l", xlab = "temperature")
lines(t, gas_O2sat(S = 0, t = t, method = "Weiss"), col = "green",
      lwd = 2, lty = "dashed")
lines(t, gas_O2sat(S = 35, t = t, method = "Weiss"), col = "red", lwd = 2)
legend("topright", c("S=35", "S=0"), col = c("red","green"), 
       lty = c(1, 2), lwd = 2)


###################################################
### code chunk number 23: marelac.Rnw:386-389
###################################################
gas_satconc(species = "O2")
Temp <- seq(from = 0, to = 30, by = 0.1)
Sal  <- seq(from = 0, to = 35, by = 0.1)


###################################################
### code chunk number 24: solub
###################################################
#
mf  <-par(mfrow = c(2, 2))
#
gs  <-gas_solubility(t = Temp)
species   <- c("CCl4", "CO2", "N2O", "Rn", "CCl2F2")
matplot(Temp, gs[, species], type = "l", lty = 1, lwd = 2, xlab = "temperature",
     ylab = "mmol/m3", main = "solubility (S=35)")
legend("topright", col = 1:5, lwd = 2, legend = species)
#
species2 <- c("Kr", "CH4", "Ar", "O2", "N2", "Ne")
matplot(Temp, gs[, species2], type = "l", lty = 1, lwd = 2, xlab = "temperature",
     ylab = "mmol/m3", main = "solubility (S=35)")
legend("topright", col = 1:6, lwd = 2, legend = species2)
#

species <- c("N2", "CO2", "O2", "CH4", "N2O")
gsat  <-gas_satconc(t = Temp, species = species)
matplot(Temp, gsat, type = "l", xlab = "temperature", log = "y", lty = 1,
     ylab = "mmol/m3", main = "Saturated conc (S=35)", lwd = 2)
legend("right", col = 1:5, lwd = 2, legend = species)
#
gsat  <-gas_satconc(S = Sal, species = species)
matplot(Sal, gsat, type = "l", xlab = "salinity", log = "y", lty = 1,
     ylab = "mmol/m3", main = "Saturated conc (T=20)", lwd = 2)
legend("right", col = 1:5, lwd = 2, legend = species)
#
par("mfrow" = mf)


###################################################
### code chunk number 25: figsolub
###################################################
#
mf  <-par(mfrow = c(2, 2))
#
gs  <-gas_solubility(t = Temp)
species   <- c("CCl4", "CO2", "N2O", "Rn", "CCl2F2")
matplot(Temp, gs[, species], type = "l", lty = 1, lwd = 2, xlab = "temperature",
     ylab = "mmol/m3", main = "solubility (S=35)")
legend("topright", col = 1:5, lwd = 2, legend = species)
#
species2 <- c("Kr", "CH4", "Ar", "O2", "N2", "Ne")
matplot(Temp, gs[, species2], type = "l", lty = 1, lwd = 2, xlab = "temperature",
     ylab = "mmol/m3", main = "solubility (S=35)")
legend("topright", col = 1:6, lwd = 2, legend = species2)
#

species <- c("N2", "CO2", "O2", "CH4", "N2O")
gsat  <-gas_satconc(t = Temp, species = species)
matplot(Temp, gsat, type = "l", xlab = "temperature", log = "y", lty = 1,
     ylab = "mmol/m3", main = "Saturated conc (S=35)", lwd = 2)
legend("right", col = 1:5, lwd = 2, legend = species)
#
gsat  <-gas_satconc(S = Sal, species = species)
matplot(Sal, gsat, type = "l", xlab = "salinity", log = "y", lty = 1,
     ylab = "mmol/m3", main = "Saturated conc (T=20)", lwd = 2)
legend("right", col = 1:5, lwd = 2, legend = species)
#
par("mfrow" = mf)


###################################################
### code chunk number 26: vapor
###################################################
plot(0:30, vapor(t = 0:30), xlab = "Temperature, dgC", ylab = "pH2O/P",
     type = "l")


###################################################
### code chunk number 27: figvapor
###################################################
plot(0:30, vapor(t = 0:30), xlab = "Temperature, dgC", ylab = "pH2O/P",
     type = "l")


###################################################
### code chunk number 28: marelac.Rnw:462-465
###################################################
gas_schmidt(species = "CO2", t = 20)

useq <- 0:15


###################################################
### code chunk number 29: transfer
###################################################
plot(useq, gas_transfer(u10 = useq, species = "O2"), type = "l", 
     lwd = 2, xlab = "u10,m/s", ylab = "m/s", 
     main = "O2 gas transfer velocity", ylim = c(0, 3e-4))
lines(useq, gas_transfer(u10 = useq, species = "O2", method = "Nightingale"),
      lwd = 2, lty = 2)
lines(useq, gas_transfer(u10 = useq, species = "O2", method = "Wanninkhof1"), 
      lwd = 2, lty = 3)
lines(useq, gas_transfer(u10 = useq, species = "O2", method = "Wanninkhof2"),
      lwd = 2, lty = 4)

legend("topleft", lty = 1:4, lwd = 2, legend = c("Liss and Merlivat 1986",
   "Nightingale et al. 2000", "Wanninkhof 1992", "Wanninkhof and McGills 1999"))


###################################################
### code chunk number 30: figtransfer
###################################################
plot(useq, gas_transfer(u10 = useq, species = "O2"), type = "l", 
     lwd = 2, xlab = "u10,m/s", ylab = "m/s", 
     main = "O2 gas transfer velocity", ylim = c(0, 3e-4))
lines(useq, gas_transfer(u10 = useq, species = "O2", method = "Nightingale"),
      lwd = 2, lty = 2)
lines(useq, gas_transfer(u10 = useq, species = "O2", method = "Wanninkhof1"), 
      lwd = 2, lty = 3)
lines(useq, gas_transfer(u10 = useq, species = "O2", method = "Wanninkhof2"),
      lwd = 2, lty = 4)

legend("topleft", lty = 1:4, lwd = 2, legend = c("Liss and Merlivat 1986",
   "Nightingale et al. 2000", "Wanninkhof 1992", "Wanninkhof and McGills 1999"))


###################################################
### code chunk number 31: marelac.Rnw:500-501
###################################################
sw_conserv(S = seq(0, 35, by = 5))


###################################################
### code chunk number 32: marelac.Rnw:521-523
###################################################
convert_AStoPS(S = 35)
convert_PStoAS(S = 35)


###################################################
### code chunk number 33: marelac.Rnw:550-559
###################################################
convert_PStoAS(S = 35, lat = -10, lon = 0)
convert_PStoAS(S = 35, lat = 0, lon = 0)
convert_PStoAS(S = 35, lat = 10, lon = 0)
convert_PStoAS(S = 35, lat = -10, lon = 0)

convert_PStoAS(S = 35, lat = -10, DSi = 1:10, Ocean = "Pacific")

dsal <- t(sw_sfac$del_sa[1, , ])
dsal [dsal < -90] <- NA


###################################################
### code chunk number 34: sfac
###################################################
image(sw_sfac$longs, sw_sfac$lats, dsal, col = femmecol(100),
      asp = TRUE, xlab = "dg", ylab = "dg",
      main = "salinity conversion - p = 0 bar")
contour(sw_sfac$longs, sw_sfac$lats, dsal, asp = TRUE, add = TRUE)


###################################################
### code chunk number 35: sfac
###################################################
image(sw_sfac$longs, sw_sfac$lats, dsal, col = femmecol(100),
      asp = TRUE, xlab = "dg", ylab = "dg",
      main = "salinity conversion - p = 0 bar")
contour(sw_sfac$longs, sw_sfac$lats, dsal, asp = TRUE, add = TRUE)


###################################################
### code chunk number 36: sfac2
###################################################
ii <- c(6, 21, 24, 43)
par(mfrow = c(2, 2))
for ( i in ii)
{
dsal <- t(sw_sfac$del_sa[ ,i, ])
dsal [dsal < -90] <- 0
image(sw_sfac$longs, sw_sfac$p, dsal, col = c("black", femmecol(100)),
      xlab = "longitude, dg", ylab = "depth, m", zlim = c(0, 0.018),
      main = sw_sfac$lat[i], ylim = c(6000, 0))
contour(sw_sfac$longs, sw_sfac$p, dsal, asp = TRUE, add = TRUE)
}


###################################################
### code chunk number 37: sfac2
###################################################
ii <- c(6, 21, 24, 43)
par(mfrow = c(2, 2))
for ( i in ii)
{
dsal <- t(sw_sfac$del_sa[ ,i, ])
dsal [dsal < -90] <- 0
image(sw_sfac$longs, sw_sfac$p, dsal, col = c("black", femmecol(100)),
      xlab = "longitude, dg", ylab = "depth, m", zlim = c(0, 0.018),
      main = sw_sfac$lat[i], ylim = c(6000, 0))
contour(sw_sfac$longs, sw_sfac$p, dsal, asp = TRUE, add = TRUE)
}


###################################################
### code chunk number 38: marelac.Rnw:620-622
###################################################
sw_cp(S = 40,t = 1:20)
sw_cp(S = 40,t = 1:20, method = "UNESCO")


###################################################
### code chunk number 39: marelac.Rnw:628-641
###################################################
t  <-  25.5
p  <- 1023/10  # pressure in bar
S  <- 35.7
sw_alpha(S, t, p)               -0.0003098378393192645
sw_beta(S, t, p)                -0.0007257297978386655
sw_cp(S,t, p)                   -3974.42541259729
sw_tpot(S, t, p)                -25.2720983155409
sw_dens(S, t, p)                -1027.95249315662
sw_enthalpy(S, t, p)            -110776.712408975
sw_entropy(S, t, p)             -352.81879771528
sw_kappa(S, t, p)               -4.033862685464779e-6
sw_kappa_t(S, t, p)             -4.104037946151349e-6
sw_svel(S, t, p)                -1552.93372863425


###################################################
### code chunk number 40: sw
###################################################
plotST <- function(fun, title)
{
  Sal  <-  seq(0, 40, by = 0.5)
  Temp <- seq(-5, 40, by = 0.5)

  Val <- outer(X = Sal, Y = Temp, FUN = function(X, Y) fun(S = X, t = Y))
  contour(Sal, Temp, Val, xlab = "Salinity", ylab = "temperature",
          main = title, nlevel = 20)
}

par (mfrow = c(3, 2))
par(mar = c(4, 4, 3, 2))
plotST(sw_gibbs, "Gibbs function")
plotST(sw_cp, "Heat capacity")
plotST(sw_entropy, "Entropy")
plotST(sw_enthalpy, "Enthalpy")
plotST(sw_dens, "Density")
plotST(sw_svel, "Sound velocity")


###################################################
### code chunk number 41: sw
###################################################
plotST <- function(fun, title)
{
  Sal  <-  seq(0, 40, by = 0.5)
  Temp <- seq(-5, 40, by = 0.5)

  Val <- outer(X = Sal, Y = Temp, FUN = function(X, Y) fun(S = X, t = Y))
  contour(Sal, Temp, Val, xlab = "Salinity", ylab = "temperature",
          main = title, nlevel = 20)
}

par (mfrow = c(3, 2))
par(mar = c(4, 4, 3, 2))
plotST(sw_gibbs, "Gibbs function")
plotST(sw_cp, "Heat capacity")
plotST(sw_entropy, "Entropy")
plotST(sw_enthalpy, "Enthalpy")
plotST(sw_dens, "Density")
plotST(sw_svel, "Sound velocity")


###################################################
### code chunk number 42: sw2
###################################################
par (mfrow = c(3, 2))
par(mar = c(4, 4, 3, 2))
plotST(sw_kappa, "Isentropic compressibility")
plotST(sw_kappa_t, "Isothermal compressibility")
plotST(sw_alpha, "Thermal expansion coefficient")
plotST(sw_beta, "Haline contraction coefficient")
plotST(sw_adtgrad, "Adiabatic temperature gradient")
par (mfrow = c(1, 1))


###################################################
### code chunk number 43: sw2
###################################################
par (mfrow = c(3, 2))
par(mar = c(4, 4, 3, 2))
plotST(sw_kappa, "Isentropic compressibility")
plotST(sw_kappa_t, "Isothermal compressibility")
plotST(sw_alpha, "Thermal expansion coefficient")
plotST(sw_beta, "Haline contraction coefficient")
plotST(sw_adtgrad, "Adiabatic temperature gradient")
par (mfrow = c(1, 1))


###################################################
### code chunk number 44: sw3
###################################################
par (mfrow = c(2, 2))
par(mar = c(4, 4, 3, 2))
plotST(function(S, t) sw_dens(S, t, method = "UNESCO") - sw_dens(S, t),
  "Density UNESCO - Gibbs")
plotST(function(S, t) sw_cp(S, t, method = "UNESCO") - sw_cp(S, t),
  "Heat capacity UNESCO - Gibbs")
plotST(function(S, t) sw_svel(S, t, method = "UNESCO") - sw_svel(S, t),
  "Sound velocity UNESCO - Gibbs")
par (mfrow = c(1, 1))


###################################################
### code chunk number 45: sw3
###################################################
par (mfrow = c(2, 2))
par(mar = c(4, 4, 3, 2))
plotST(function(S, t) sw_dens(S, t, method = "UNESCO") - sw_dens(S, t),
  "Density UNESCO - Gibbs")
plotST(function(S, t) sw_cp(S, t, method = "UNESCO") - sw_cp(S, t),
  "Heat capacity UNESCO - Gibbs")
plotST(function(S, t) sw_svel(S, t, method = "UNESCO") - sw_svel(S, t),
  "Sound velocity UNESCO - Gibbs")
par (mfrow = c(1, 1))


###################################################
### code chunk number 46: marelac.Rnw:741-746
###################################################
1/molweight("CO3")
1/molweight("HCO3")
1/molweight(c("C2H5OH", "CO2", "H2O"))

molweight(c("SiOH4", "NaHCO3", "C6H12O6", "Ca(HCO3)2", "Pb(NO3)2", "(NH4)2SO4"))


###################################################
### code chunk number 47: marelac.Rnw:752-757
###################################################
#species <- colnames(gs)  ## thpe: does not work any more, because 1D return value is vector
species = c("He", "Ne", "N2", "O2",
  "Ar", "Kr", "Rn", "CH4", "CO2", "N2O", "CCl2F2", "CCl3F", "SF6", "CCl4")
gs <- gas_solubility(species = species)
mw <- molweight(species)


###################################################
### code chunk number 48: mwgs
###################################################
plot(mw, gs, type = "n", xlab = "molecular weight", 
     ylab = "solubility", log = "y")
text(mw, gs, species)


###################################################
### code chunk number 49: mwgs
###################################################
plot(mw, gs, type = "n", xlab = "molecular weight", 
     ylab = "solubility", log = "y")
text(mw, gs, species)


###################################################
### code chunk number 50: marelac.Rnw:780-782
###################################################
molvol(species = "ideal")
molvol(species = "ideal", t = 1:10)


###################################################
### code chunk number 51: marelac.Rnw:784-787
###################################################
1/molvol(species = "O2", t = 0)*1000
1/molvol(species = "O2", q = 1:6, t = 0)
1/molvol(t = 1:5, species = c("CO2", "O2", "N2O"))


###################################################
### code chunk number 52: marelac.Rnw:810-811
###################################################
redfield(1, "P")


###################################################
### code chunk number 53: marelac.Rnw:817-818
###################################################
redfield(1, "N")


###################################################
### code chunk number 54: marelac.Rnw:825-826
###################################################
redfield(2, "P", "mass")


###################################################
### code chunk number 55: marelac.Rnw:831-833
###################################################
x <- redfield(1, "P", "mass")
x / sum(x)


###################################################
### code chunk number 56: marelac.Rnw:838-841
###################################################
stumm <- c(C = 106, H = 180, O = 45, N = 16, P = 1)
x <- redfield(1, "P", "mass", ratio = stumm)
x / sum(x)


###################################################
### code chunk number 57: marelac.Rnw:854-855
###################################################
convert_p(1, "atm")


###################################################
### code chunk number 58: marelac.Rnw:863-864
###################################################
convert_T(1, "C")


###################################################
### code chunk number 59: marelac.Rnw:872-875
###################################################
convert_StoCl(S = 35)
convert_RtoS(R = 1)
convert_StoR(S = 35)


