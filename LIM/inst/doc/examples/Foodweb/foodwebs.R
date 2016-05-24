###############################################
# To read the files in this directory:
# First set the working directory
# setwd("    LIM/examples/Foodweb")
###############################################
require(LIM)
Inverse <- Read("rigaautumn.input")
donali <- Setup(Inverse)
rigaweb <- Flowmatrix(donali)
donweb <- Lsei(donali)
Xranges (donali)
Varranges(donali)
Variables(donali, donweb$X)
plotweb(rigaweb)

LIMTakapoto   <- Setup("takapoto.input", verbose = FALSE)
LIMRigaAutumn <- Setup("rigaautumn.input")
LIMRigaSpring <- Setup("rigaspring.input", verbose = FALSE)
LIMRigaSummer <- Setup("rigasummer.input", verbose = FALSE)
LIMBrouageMudflat <-Setup("BrouageMudflat.input", verbose = FALSE)
LIMEverglades <- Setup("Everglades.input", verbose = FALSE)
LIMScheldtIntertidal <- Setup("ScheldtIntertidal.input", verbose = FALSE)
LIMCaliforniaSediment <- Setup("CaliforniaSediment.input", verbose = FALSE)

ll <- Read("LIM_example.input")
lim <- Setup(ll)
Ldei(lim)$X
# all in one
Lsei("LIM_example.input", verbose = FALSE)$X
Ldei("LIM_example.input", verbose = FALSE)$X

