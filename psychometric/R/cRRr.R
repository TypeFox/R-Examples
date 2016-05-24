"cRRr" <-
function (rr, sdy, sdyu)
 {
 rxy <- (rr*(sdyu/sdy))/sqrt(1+rr^2*((sdyu^2/sdy^2)-1))
 return(data.frame(unrestricted = rxy))
 }

