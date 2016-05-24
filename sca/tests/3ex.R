library(sca)
data(hearlossC)
sc.hear <- sca(hearlossC)
sc.hear
sca(hearlossC, cluster = "single")
sca(hearlossC, qmin= 6, corblocks = 0.4)

data(reflexesC)
sc.refl <- sca(reflexesC)
sc.refl
sca(reflexesC, cluster = "complete")
sca(reflexesC, qmin = 4, b = 3, corblocks = 0, cluster = "single")
sca(reflexesC, qmin = 4)

data(pitpropC)
sc.pitp <- sca(pitpropC)
sc.pitp

sca(pitpropC, cluster = "single")
sca(pitpropC,        d = 3, qmin= 0)
sca(pitpropC, b = 1, d = 3, qmin= 0, corblocks=0)


