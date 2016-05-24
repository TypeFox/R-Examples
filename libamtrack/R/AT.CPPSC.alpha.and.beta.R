AT.CPPSC.alpha.and.beta <- function( E.MeV.u,
                                     particle.no,
                                     fluence.cm2.or.dose.Gy,
                                     material.no,
                                     RDD.model,
                                     RDD.parameters,
                                     ER.model,
                                     gamma.model,
                                     gamma.parameters,
                                     N2                       = 20,
					                 fluence.factor           = 1.0,
					                 write.output             = FALSE,
                                     shrink.tails             = TRUE,
                                     shrink.tails.under       = 1e-30,
                                     adjust.N2                = TRUE)
{
    df                     <- expand.grid( fluence.factor  = c(seq(0.1, 2, by = .25), seq(2.5, 10, by = 1)),
                                           D.check.Gy      = 0,
                                           response.HCP    = 0,
                                           response.gamma  = 0)

    for (i in 1:nrow(df)){
        CPPSC.res              <- AT.run.CPPSC.method(    E.MeV.u                     = E.MeV.u,
	    								    particle.no                 = particle.no,
		    							    fluence.cm2.or.dose.Gy      = df$fluence.factor[i] * fluence.cm2.or.dose.Gy,
			    						    material.no                 = material.no,
									    rdd.model                   = RDD.model,
									    rdd.parameters              = RDD.parameters,
									    er.model                    = ER.model,
									    gamma.model                 = gamma.model,
									    gamma.parameters            = gamma.parameters,
									    N2                          = N2,
									    fluence.factor              = 1.0,
									    write.output                = write.output,
									    shrink.tails                = shrink.tails,
									    shrink.tails.under          = shrink.tails.under,
									    adjust.N2                   = adjust.N2,
									    lethal.events.mode          = TRUE)
        df$D.check.Gy[i]       <- CPPSC.res[which(names(CPPSC.res) == "d.check")]
        df$response.HCP[i]     <- exp(-CPPSC.res[which(names(CPPSC.res) == "S.HCP")])
        df$response.gamma[i]   <- exp(-CPPSC.res[which(names(CPPSC.res) == "S.gamma")])
    }

    HCP.surv <- nls(-log(response.HCP) ~ alpha * D.check.Gy + beta * D.check.Gy^2,
                    data = df,
                    start = list(alpha = 10, beta = 0.0001))
    
    alpha    <- coef(HCP.surv)[1]
    beta     <- coef(HCP.surv)[2]
    df$response.HCP.LQ         <- exp(-alpha * df$D.check.Gy - beta * df$D.check.Gy^2)
    df$residuals               <- (df$response.HCP.LQ - df$response.HCP)/df$response.HCP 
    return(list( alpha = alpha, beta = beta, data = df))
}
