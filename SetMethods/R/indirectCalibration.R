indirectCalibration <-
function(x, x_cal, binom = TRUE){
	
    if(binom == TRUE){
		options(warn = -1)
    	Z_bin <- glm(x_cal ~ x, family = binomial)
   		b_bin <- coef(Z_bin)
    	y_hat_bin <- 1/(1 + (exp(-(b_bin[1] + b_bin[2] * x))))
    	return(structure(y_hat_bin, class = "SetMethod"))
		options(warn = 1)
    }
    else{
        N <- length(x_cal)
        x_cal_t <- ((x_cal*(N-1) + 1/2)/N)
        Z_bet <- betareg(x_cal_t ~ x, link = "logit", phi = TRUE)
        y_hat_bet <- b_bet <- coef(Z_bet)
        y_hat_bet <- 1/(1 + (exp(-(b_bet[1] + b_bet[2] * x))))    
        return(structure(y_hat_bet, class = "SetMethod"))
    }
}
