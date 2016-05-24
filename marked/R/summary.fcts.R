#' Various utility parameter summary functions
#' 
#' Several functions have been added to compute mean values and boxplots of values. Currently
#' they have only been defined for Phi and p for cjs and they are not generic.
#' 
#' @aliases Phi.mean p.mean Phi.boxplot p.boxplot 
#' @usage Phi.mean(x,age=0,time=NULL,age.bins=NULL,age.levels=NULL)
#' 
#'        p.mean(x,age=0,time=NULL,age.bins=NULL,age.levels=NULL)
#' 
#'        Phi.boxplot(x,age=0,time=NULL,sex=NULL)
#' 
#'        p.boxplot(x,age=0,time=NULL,sex=NULL)
#' 
#' @param x dataframe of reals contained in model
#' @param age at which Phi or p should be shown across time
#' @param time at which Phi or p should be shown across ages
#' @param sex for which Phi or p should be shown across ages
#' @param age.bins bins for age in which values are summarized 
#' @param age.levels labels for age.bins
#' @export Phi.mean
#' @export p.mean
#' @export Phi.boxplot
#' @export p.boxplot
#' @return matrix of labelled values for Phi.mean and p.mean or boxplot object 
#' @author Jeff Laake
#' @keywords utility
Phi.mean=function(x,age=0,time=NULL,age.bins=NULL,age.levels=NULL)
{	
	if(is.null(x$sex))x$sex="All"
	if(is.null(time))
	{	
		with(x[x$Time>=x$Cohort&x$Age%in%age,],tapply(Phi,list(sex,Phi.time),mean))
	} else
	{
		x$age=cut(as.numeric(x$Phi.age),age.bins,right=FALSE)
		levels(x$age)=age.levels
		with(x[x$Phi.time%in%time,],tapply(Phi,list(sex,age),mean))
	}
}
p.mean=function(x,age=0,time=NULL,age.bins=NULL,age.levels=NULL)
{	
	if(is.null(x$sex))x$sex="All"
	if(is.null(time))
	{
		with(x[x$Time>=x$Cohort&x$Age%in%age,],tapply(p,list(sex,p.time),mean))
	} else
	{	
		if(is.null(age.bins))
			with(x[x$p.time%in%time,],tapply(p,list(sex,p.age),mean))
		else
		{
			x$age=cut(as.numeric(x$p.age),age.bins,right=FALSE)
			levels(x$age)=age.levels
			with(x[x$p.time%in%time,],tapply(p,list(sex,age),mean))
		}
	}
}


Phi.boxplot=function(x,age=0,time=NULL,sex=NULL){
	if(!is.null(sex))
	{
		if(is.null(time))
		{
			boxplot(Phi~Phi.time,data=x[x$Age==age&x$sex==sex,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age,"and sex ",sex))
		} else
		{
			boxplot(Phi~Phi.age,data=x[x$Phi.time%in%time&x$sex==sex,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time,"and sex ",sex))
			
		}
	}else
	{
		if(is.null(time))
		{
			boxplot(Phi~Phi.time,data=x[x$Age==age,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age))
		} else
		{
			boxplot(Phi~Phi.age,data=x[x$Phi.time==time,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time))
		}
	}
}

p.boxplot=function(x,age=0,time=NULL,sex=NULL){
	if(!is.null(sex))
	{
		if(is.null(time))
		{
			boxplot(p~p.time,data=x[x$Age==age&x$sex==sex,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age,"and sex ",sex))
		} else
		{
			boxplot(p~p.age,data=x[x$p.time==time&x$sex==sex,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time,"and sex ",sex))
			
		}
	}else
	{
		if(is.null(time))
		{
			boxplot(p~p.time,data=x[x$Age==age,],ylim=c(0,1),xlab="Cohort",ylab=paste("Survival for age",age))
		} else
		{
			boxplot(p~p.age,data=x[x$p.time==time,],ylim=c(0,1),xlab="Age",ylab=paste("Survival for time",time))
		}
	}
}





