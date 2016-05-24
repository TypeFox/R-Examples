"CIrb" <-
function (x, LEVEL=.95, homogenous=TRUE)
 {
rb <- rbar(x)
 noma <- 1 - LEVEL
if (!homogenous)
{serb <- SERHOM(x)}
else {serb <- SERHET(x)}
 zs <- -qnorm(noma/2)
 merb <- zs*serb
 lcl <- rb - merb
 ucl <- rb + merb
 mat <- list(lcl, ucl)
 return(as.numeric(mat))
}

