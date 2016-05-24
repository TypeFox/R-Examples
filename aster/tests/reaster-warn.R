
 library(aster)

 data(radish)

 options(digits=4) # avoid rounding differences

 pred <- c(0,1,2)
 fam <- c(1,3,2)
 famlist <- fam.default()
 famlist[[2]] <- fam.negative.binomial(3.1415926526)

 rout <- try(reaster(resp ~ varb + fit : (Site * Region),
     list(block = ~ 0 + fit : Block, pop = ~ 0 + fit : Pop),
     pred, fam, varb, id, root, data = radish, famlist = famlist))

