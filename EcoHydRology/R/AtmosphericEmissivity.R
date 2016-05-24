AtmosphericEmissivity <- function (airtemp, cloudiness, vp=NULL, opt="linear") 
{
    if (opt == "Brutsaert") {
		if(is.null(vp)){ print("To use Brutsaert's 1975 Clear-Sky Emissivity eqn, enter vapor pressure in kPa")
			} else return((1.24*((vp*10)/(T+273.2))^(1/7)) * (1 - 0.84 * cloudiness) + 
				0.84 * cloudiness)
	} else {
	return((0.72 + 0.005 * airtemp) * (1 - 0.84 * cloudiness) + 
        0.84 * cloudiness)
	}
}

