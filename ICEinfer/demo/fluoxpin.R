require(ICEinfer)
# input the fluoxpin data of Sacristan et al. (2000).
data(fluoxpin)
# Effectiveness = respond, Cost = cost, trtm = flxpin where
# flxpin = 1 ==> fluoxetine plus pindolol and flxpin = 0 ==> fluoxetine alone

cat("\n Display of Lambda => Shadow Price Summary Statistics...\n")
ICEscale(fluoxpin, flxpin, respond, cost)
ICEscale(fluoxpin, flxpin, respond, cost, lambda=100000)

cat("\nBootstrap ICE Uncertainty calculations can be lengthy...\n")
fpunc <<- ICEuncrt(fluoxpin, flxpin, respond, cost, R = 10000, lambda=100000)
fpunc

cat("\nDisplay the Bootstrap ICE Uncertainty Distribution...\n")
plot(fpunc)

fpwdg <- ICEwedge(fpunc)
fpwdg
opar <- par(ask = dev.interactive(orNone = TRUE))
cat("\nClick within graphics window to display the Bootstrap 95% Confidence Wedge...\n")
plot(fpwdg)

cat("\nComputing VAGR Acceptability and ALICE Curves...\n")
fpacc <- ICEalice(fpwdg)
plot(fpacc)

cat("\nColor Interior of Confidence Wedge with LINEAR Economic Preferences...\n")
fpcol <- ICEcolor(fpwdg, gamma=1)
plot(fpcol)

cat("\nIncrease Lambda and Recolor Confidence Wedge with NON-Linear Preferences...\n")
fpcol <- ICEcolor(fpwdg, lfact=10)
plot(fpcol)

cat("\nDecrease Lambda and Recolor Confidence Wedge with LINEAR Preferences...\n")
fpcol <- ICEcolor(fpwdg, lfact=10, gamma=1)
plot(fpcol)
par(opar)