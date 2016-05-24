if (interactive())
   old.prompt <- devAskNewPage(TRUE)

##################################################
#   Some examples from the OpenRepGrid package   #
##################################################

#### Biplots ####
plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="")
text(1, 0, "Biplots", cex=5)

#### Standard biplot ####
biplot2d(feixas2004)

#### Slater's INGRID biplot ####
biplotSlater2d(feixas2004)

#### Creating an ESA biplot ####
biplotEsa2d(feixas2004)

#### Pseudo 3D biplot #####
biplotPseudo3d(feixas2004)

#### 3D biplot #####
x <- scan(n=1)
biplot3d(boeker)

#### Bertin displays ####
plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="")
text(1, 0, "Bertin \ndisplay", cex=5)

#### Bertin display ####
bertin(feixas2004)

#### Colored display ####
bertin(feixas2004, colors=c("white", "darkred"))

#### Clustered Bertin display ####
x <- scan(n=1)
dev.off()
dev.new()

plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="")
text(1, 0, "Clustered \nBertin display", cex=5)

#### Clustered Bertin display ####
bertinCluster(feixas2004)

###############################################
if (interactive())
  devAskNewPage(old.prompt)
  