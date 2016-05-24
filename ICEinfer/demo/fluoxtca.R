require(ICEinfer)
# input the fluoxtca data of Obenchain et al. (1997).
data(fluoxtca)
# Effectiveness = stable, Cost = cost, Trtm = fluox where
# fluox = 1 ==> fluoxetine treatment and fluox = 0 ==> TCA or HCA
#
# Display of Lambda => Shadow Price Summary Statistics...

ICEscale(fluoxtca, fluox, stable, cost)
ICEscale(fluoxtca, fluox, stable, cost, lambda=10000)

# Bootstrap ICE Uncertainty calculations can be lengthy...
ftunc <- ICEuncrt(fluoxtca, fluox, stable, cost, R = 10000, lambda=10000)
ftunc

# Display the Bootstrap ICE Uncertainty Distribution...\n")
plot(ftunc)

ftwdg <- ICEwedge(ftunc)
ftwdg
opar <- par(ask = dev.interactive(orNone = TRUE))
# Click within graphics window to display the Bootstrap 95% Confidence Wedge...
plot(ftwdg)

# Computing VAGR Acceptability and ALICE Curves...\n")
ftacc <- ICEalice(ftwdg)
plot(ftacc)

# Color Interior of Confidence Wedge with LINEAR Economic Preferences...
ftcol <- ICEcolor(ftwdg, gamma=1)
plot(ftcol)

# Increase Lambda and Recolor Confidence Wedge with NON-Linear Preferences...
ftcol2 <- ICEcolor(ftwdg, lfact=10)
plot(ftcol2)

# Decrease Lambda and Recolor Confidence Wedge with LINEAR Preferences...
ftcol3 <- ICEcolor(ftwdg, lfact=10, gamma=1)
plot(ftcol3)
par(opar)