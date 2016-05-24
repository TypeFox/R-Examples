# Package Code 'fishmove'
# 
# Author: Johannes Radinger
###############################################################################


# Plotting dispersal kernel 
pdk <-function(fishmove,p=0.67,...){
	
	#validation of arguments
	if(!is.numeric(p)) stop("p is not numeric")
	if(0 > p | p > 1) stop("p must be a value between 0 and 1")
	#if(length(fishmove$pred.fishmove["fit","sigma_stat",,,,])>1) stop("for plotting: please provide only single variables in fishmove")
	
	
	#dimesion to plot
	if(length(fishmove$pred.fishmove["fit","sigma_stat",,,,])>1 && !hasArg(dim)){
		warning("Multiple values supplied in fishmove. Only plotted for first values.",call. = FALSE)
	} 
	
	if(hasArg(dim)){
		print("Specifying dimensions not implemented yet")
		#fishmove.array <- fishmove$pred.fishmove[dim]
	} else{
		fishmove.array <- fishmove$pred.fishmove[,,1,1,1,1]
	}
	
	
	# General equation of heterogeneous movement kernel (two superimposed normal distributions)
	eq <- function(x,sigma_stat,sigma_mob,p) {
		(dnorm(x, sd = sigma_stat)*p)+(dnorm(x, sd = sigma_mob)*(1-p))
	}
	
	# Equations for mean, upper and lower movement kernel
	eq_fit <- function(x) {
		eq(x,fishmove.array["fit","sigma_stat"],fishmove.array["fit","sigma_mob"],p)
	}
	eq_lwr <- function(x) {
		eq(x,fishmove.array["lwr","sigma_stat"],fishmove.array["lwr","sigma_mob"],p)
	}
	eq_upr <- function(x) {
		eq(x,fishmove.array["upr","sigma_stat"],fishmove.array["upr","sigma_mob"],p)
	}
	
	# Creating dataframe for plotting, indicating min and max of x axis
	x <- seq(-1.5*fishmove.array["fit","sigma_mob"],1.5*fishmove.array["fit","sigma_mob"],1)
	y <- eq(x,fishmove.array["fit","sigma_stat"],fishmove.array["fit","sigma_mob"],p)
	tmp <- data.frame(x=x, y=y)
	
	# Plot commands to ggplot2
	pdk <- ggplot(tmp, aes(x=x, y=y))+
			stat_function(fun=eq_lwr,aes(colour="lower/upper bound"),n=500)+
			stat_function(fun=eq_upr,aes(colour="lower/upper bound"),n=300)+
			stat_function(fun=eq_fit,aes(colour="fitted mean"),size=0.8,n=300)+
			scale_y_continuous("Probability",limits=c(0,eq(0,fishmove.array["fit","sigma_stat"],fishmove.array["fit","sigma_mob"],p)))+
			scale_x_continuous("Movement Distance (m)")+
			scale_colour_manual(name="",values=c("fitted mean"="green4","lower/upper bound"="grey65"))+
			theme_bw()+
			theme(	legend.position = c(0.95,0.85),
					legend.justification ="right",
					legend.title = element_blank(),
					legend.background = element_rect(colour = 'grey', fill = 'white'),
					panel.grid.minor = element_blank())
	
	pdk #plot final density plot
}
