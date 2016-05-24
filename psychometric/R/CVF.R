"CVF" <-
function(x)
 {
 ma <- aRxx(x)[[1]]
 va <- aRxx(x)[[2]]
 mb <- bRyy(x)[[1]]
 vb <- bRyy(x)[[2]]
 mc <- cRR(x)[[1]]
 vc <- cRR(x)[[2]]
 cv <- va/ma^2 + vb/mb^2 + vc/mc^2
 return(cv)
 }

