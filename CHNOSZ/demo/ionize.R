## ionize.aa(): Contour plots of net charge and ionization properties of LYSC_CHICK
aa <- ip2aa("LYSC_CHICK")
pH <- seq(0, 14, 0.2)
T <- seq(0, 200, 2)
val <- expand.grid(pH=pH, T=T)
par(mfrow=c(2, 2))
for(X in c("Z", "A", "Cp", "V")) {
  Y <- ionize.aa(aa, property=X,  pH=val$pH, T=val$T)
  contour(pH, T, matrix(Y[, 1], ncol=length(T)),
    xlab="pH", ylab=axis.label("T"))
  title(main=axis.label(X))
}
par(mfrow=c(1, 1))
pu <- par("usr")
text(mean(pu[1:2]), sum(pu[3:4])*0.45, 
  "additive properties of ionization of LYSC_CHICK")
