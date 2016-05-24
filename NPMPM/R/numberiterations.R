numberiterations <-
function(iterations,errorb,cfu_average){
if (iterations<100){
return(TRUE)} else {
errorbound <- errorb*cfu_average[iterations]
calcerror <- abs(cfu_average[iterations] - cfu_average[iterations-1])
return(calcerror > errorbound)
} # END ELSE
} # END FUNCTION

