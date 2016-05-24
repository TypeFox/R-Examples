saturation <-
function(sigma_inf=NULL,coeff=NULL,pas=NULL)
{
if (is.null(sigma_inf)) sigma_inf <- 0.4                 ## default values
if (is.null(coeff)) coeff <- 1
if (is.null(pas)) pas <- 0.1

TL<- seq(from=2, to=7, by=pas)
sigma<-sigma_inf*(1-exp(-coeff*(TL-2)))
return(sigma)
}

