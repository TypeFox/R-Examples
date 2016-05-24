plot_acet <- function(acet, boot = FALSE, heri = FALSE)
{
	if(!(class(acet) %in% c('AtCtEt_model', 'AtCtEp_mc_model','AtEtp_mc_model','AtCtEtp_mc_model')))
	{
		stop('The first parameter must be an acet object.')
	}

	if(class(acet)=='AtCtEt_model')
	{
		if(heri == FALSE)
		{
			plot_AtCtEt(acet, boot)
		}else{
			plot_AtCtEt_h(acet, boot)
		}
	}
	
	
	#if(class(acet)=='AtCtEp_mc_model')
	#{
	#	plot_AtCtEp(acet)
	#}

	if(class(acet)=='AtCtEtp_mc_model')
	{
		if(heri==FALSE)
		{
	    plot_AtCtEtp(acet)
		}else{
		  plot_AtCtEt_h(acet, boot)
		}
	}

	#if(class(acet)=='AtEtp_mc_model')
	#{
	#	plot_AtEtp(acet)
	#}

}