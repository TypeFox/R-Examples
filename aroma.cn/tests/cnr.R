library("aroma.cn")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build up a tumor CN profile
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pT <- cnr(1,1000, 2) +
     cnr(400,500) +
     cnr(600,800) +
     cnr(600,700) +
     cnr(100,200) - cnr(850,900)
print(pT)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate CN signals from this profile
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cn <- simulateRawCopyNumbers(pT, n=2000, sd=1/2)
print(cn)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, col="#aaaaaa", ylim=c(0,5))
drawLevels(pT, col="#ff0000", lwd=3)
