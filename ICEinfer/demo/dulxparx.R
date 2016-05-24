require(ICEinfer)
# input the dulxparx data of Obenchain et al. (2000).
data(dulxparx)
# Effectiveness = idb, Cost = ru, trtm = dulx where
# dulx = 1 ==> Duloxetine treatment and dulx = 0 ==> Paroxetine treatment
#
# Display of Lambda => Shadow Price Summary Statistics...
ICEscale(dulxparx, dulx, idb, ru)
ICEscale(dulxparx, dulx, idb, ru, lambda=0.26)

# Bootstrap ICE Uncertainty calculations can be time consuming...
dpunc <- ICEuncrt(dulxparx, dulx, idb, ru, R = 10000, lambda=0.26)
dpunc

# Display the Bootstrap ICE Uncertainty Distribution...
plot(dpunc)

dpwdg <- ICEwedge(dpunc)
dpwdg
opar <- par(ask = dev.interactive(orNone = TRUE))
# Click within graphics window to display the Bootstrap 95% Confidence Wedge...
plot(dpwdg)

# Computing VAGR Acceptability and ALICE Curves...
dpacc <- ICEalice(dpwdg)
plot(dpacc)

# Color Interior of Confidence Wedge with LINEAR Economic Preferences...
dpcol <- ICEcolor(dpwdg, gamma=1)
plot(dpcol)

# Increase Lambda and Recolor Confidence Wedge with NON-Linear Preferences...
dpcol <- ICEcolor(dpwdg, lfact=10)
plot(dpcol)

# Decrease Lambda and Recolor Confidence Wedge with LINEAR Preferences...
dpcol <- ICEcolor(dpwdg, lfact=10, gamma=1)
plot(dpcol)
par(opar)