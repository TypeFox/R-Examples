stretchbetaval <- function(mn,std,minb,maxb,fx)
{
   if (std == 0) {bb <- mn} # with no variation, then the value = mean
   else
   {
      	# convert the stretched beta parameters to corresponding
       	# ones for a {0,1} beta
       	mnbeta <- (mn - minb) /(maxb - minb)
       	sdbeta <- std /(maxb-minb)
        # next, check for un-doable parameter combos
       	if (sdbeta < (mnbeta * (1 - mnbeta))^0.5)
        {
            bvalue <- betaval(mnbeta,sdbeta,fx) # find beta value
            bb <- bvalue * (maxb - minb) + minb # convert to stretched value
         
        }
       	else
        {
            maxsd <- ((mnbeta *( 1 - mnbeta))^0.5) * (maxb-minb)
            bb<-paste("WARNING: The std is too high.  The maximum std possible is:", round(maxsd,3)) 
        }
    }
    bb
} 

