
#########  load up a data set (volcano image)
data(volcano)
#########  assign information to the data set
attr(volcano, 'dx') =10
attr(volcano, 'dy') =10

###########   set the buttons

mybutts = c("DONE", "REFRESH", "rainbow", "topo", "terrain", "CONT",
"XSEC","PS" )

###  in the following change demo=FALSE to get interactive behavior
cat("Click on buttons for interaction", sep="\n")
cat("Click twice on image and ", sep="\n")

XSECDEM(volcano, mybutts, demo=FALSE)

