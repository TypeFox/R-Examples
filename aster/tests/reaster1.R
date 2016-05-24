
 library(aster)

 data(radish)

 options(digits=4) # avoid rounding differences

 pred <- c(0,1,2)
 fam <- c(1,3,2)

 rout <- reaster(resp ~ varb + fit : (Site * Region), ~ 0 + fit : Block,
     pred, fam, varb, id, root, data = radish)
 summary(rout)
 summary(rout, stand = FALSE)

