# This function is adapted from gstat:::vgm.panel.xyplot
autokrige.vgm.panel = function(x, y, model, subscripts,...)
{
    # Call original gstat panel
    vgm.panel.xyplot(x, y, model = model, subscripts = TRUE,...)

    # Add the text.
    no_digits = function(a)
    # Determine the number of digits that is nice for the size
    # of a var_model parameter
    {
        if(a > 10)
            return(0)
        else
        {
            if(a < 1)
                return(2)
            else return(1)
        }
    }
    
    nugget = sum(model[1,"psill"])
    sill = sum(model[,"psill"])
    range = sum(model[,"range"])
    txt = paste("Model: ",
              as.character(model[2,"model"]),
              "\nNugget: ",
              round(nugget, digits = no_digits(nugget)),
              "\nSill: ",
              round(sill, digits = no_digits(sill)),
              "\nRange: ",
              round(range, digits = no_digits(range)),
			  sep = "")
    if (model[2,"model"] %in% c("Mat", "Ste")) {
		kappa = model[2,"kappa"]
		txt = paste(txt,
              "\nKappa: ",
              round(kappa, digits = no_digits(kappa)), 
			  sep = "")
	}


    #ltext(max(x) - (0.06 * max(x)) , max(y) - 0.83 *max(y) ,txt,font = 2, cex = 0.7, adj = c(0,.5))
    ltext(max(x), 0.02 * max(y), txt, font = 2, cex = 0.7, adj = c(1,0), col = grey(.3))
}

