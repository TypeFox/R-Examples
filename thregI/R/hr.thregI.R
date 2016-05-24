### The following is the definition of the hr.thregI function
hr <- function(object, var, timevalue, scenario)
    UseMethod("hr")
"hr.thregI" <-
function (object,var,timevalue,scenario)
{
	para <- match.call(expand.dots = FALSE)
  indx <- match(c("object", "var", "timevalue", "scenario"), names(para), nomatch=0)
	if (indx[1] ==0) stop("An object argument is required")
	if (indx[2] ==0) stop("A var argument is required")
  if (!inherits(object, 'thregI')) stop("Primary argument must be a thregI object")

	m_timevalue<-match(c("timevalue"), names(para), 0)
 #if "timevalue" option is not specified, use the study time
	if (!m_timevalue) timevalue<- model.extract(object$mf, "response")[,1]
	if(!is.numeric(timevalue)) {
	                            stop(paste("'timevalue' option must specify a numerical value!"))
	                           }
	m_scenario<-match(c("scenario"), names(para), 0)
	m_hr<-match(c("var"), names(para), 0)
  scenario_value<-NULL
  scenario_covariate<-NULL
  if(m_scenario!=0) {
	 cur_scenario_string<-para[[m_scenario]]
	 while(length(cur_scenario_string)>=2) {
	                                        if(length(cur_scenario_string)==3)
	                                        {
			                                     current_scenario_covariate_string<-cur_scenario_string[[3]]
			                                     if(length(current_scenario_covariate_string)==2)
			                                     {
			                                      current_scenario_covariate<-current_scenario_covariate_string[[1]]
				                                    current_covariate_value<-current_scenario_covariate_string[[2]]
				                                    scenario_covariate<-c(current_scenario_covariate,scenario_covariate)
				                                    scenario_value<-c(current_covariate_value,scenario_value)
				                                   }
			                                     else {
			                                           stop("wrong scenario specification")
				                                        }
				                                    cur_scenario_string<-cur_scenario_string[[2]]
			                                    }
		                                    	else if(length(cur_scenario_string)==2)
		                                    	{
				                                   current_scenario_covariate_string<-cur_scenario_string
				                                   current_scenario_covariate<-current_scenario_covariate_string[[1]]
				                                   current_covariate_value<-current_scenario_covariate_string[[2]]

				                                   scenario_covariate<-c(current_scenario_covariate,scenario_covariate)
				                                   scenario_value<-c(current_covariate_value,scenario_value)

				                                   current_scenario_covariate_string<-NULL
				                                   cur_scenario_string<-NULL
				                                   current_scenario_covariate<-NULL
				                                   current_covariate_value<-NULL
			                                    }
		                                     }

  }
  #add intercept term
	scenario_covariate<-c("(Intercept)",scenario_covariate)
	scenario_value<-c(1,scenario_value)
	names(scenario_value)<-scenario_covariate
  m_hr<-match(c("var"), names(para), 0)
	hr_var<-as.character(para[[m_hr]])
  if (!is.factor(object$mf[[hr_var]])) stop("The variable for the hazard ratio calculation should be a factor variable")
  if (length(levels(object$mf[[hr_var]]))<2) stop("The variable for the hazard ratio calculation should have more than one level")
  if (!any(names(object$mf)==hr_var)) stop("The variable ",as.character(para[[m_hr]]) , " for the hazard ratio calculation should be included in the model")
  scenario_value_hr<-rbind(rep(0,length(levels(object$mf[[hr_var]]))-1),diag(rep(1,length(levels(object$mf[[hr_var]]))-1)))
	#names(scenario_value_hr)<-paste(hr_var,levels(object$mf[[hr_var]])[-1],sep="")
  colnames(scenario_value_hr) <- paste(hr_var,levels(object$mf[[hr_var]])[-1],sep="")
  if (any(matrix(rep(names(scenario_value),times=length(dimnames(scenario_value_hr)[[2]])),nrow=length(dimnames(scenario_value_hr)[[2]]),byrow=TRUE)==dimnames(scenario_value_hr)[[2]])) stop("Don't specify the scenario value of the dummy variable for the hazard ratio calculation!")
  scenario_value<-matrix(rep(scenario_value,times=length(levels(object$mf[[hr_var]]))),nrow=length(levels(object$mf[[hr_var]])),byrow=TRUE)
  dimnames(scenario_value)[[2]] <- as.list(scenario_covariate)
  scenario_value<-cbind(scenario_value,scenario_value_hr)
  #judge if any covariate for lny0 lacks of scenario value
	judgematrix_lny0<-matrix(rep(object$lny0,times=dim(scenario_value)[2]),ncol=dim(scenario_value)[2])==matrix(rep(dimnames(scenario_value)[[2]],times=length(object$lny0)),nrow=length(object$lny0),byrow=TRUE)
  if (any(!apply(judgematrix_lny0,1,any)))
               stop(paste("scenario value for",  object$lny0[!apply(judgematrix_lny0,1,any)][1], "is required!"))
	#choose the right position of scenario values to combine
	positionpick_lny0=judgematrix_lny0%*%c(1:dim(scenario_value)[2])
  lny0<-scenario_value[,positionpick_lny0]%*%object$coef[1:length(object$lny0)]
	y0<-exp(lny0)
  #judge if any covariate for mu lacks of scenario value
	judgematrix_mu<-matrix(rep(object$mu,times=dim(scenario_value)[2]),ncol=dim(scenario_value)[2])==matrix(rep(dimnames(scenario_value)[[2]],times=length(object$mu)),nrow=length(object$mu),byrow=TRUE)
	if (any(!apply(judgematrix_mu,1,any))) stop(paste("scenario value for",  object$mu[!apply(judgematrix_mu,1,any)][1], "is required!"))
	#choose the right position of scenario values to combine
	positionpick_mu=judgematrix_mu%*%c(1:dim(scenario_value)[2])
  mu<-scenario_value[,positionpick_mu]%*%object$coef[(length(object$lny0)+1):(length(object$lny0)+length(object$mu))]
  dim(y0)<-NULL
	dim(lny0)<-NULL
	dim(mu)<-NULL
  y0<-c(rep(y0,length(timevalue)))
	lny0<-c(rep(lny0,length(timevalue)))
	mu<-c(rep(mu,length(timevalue)))
	lenth_timevalue<-length(timevalue)
	#replicate values for matrix calculation
	timevalue_rep<-rep(timevalue,each=length(dimnames(scenario_value_hr)[[2]])+1)
  f<-exp((lny0-.5*(log(2*pi*(timevalue_rep^3))+(y0+mu*timevalue_rep)^2/timevalue_rep)))
  S<-exp(log(pnorm((mu*timevalue_rep+y0)/sqrt(timevalue_rep))-exp(-2*y0*mu)*pnorm((mu*timevalue_rep-y0)/sqrt(timevalue_rep))))
	h<-f/S
	dim(h)<-c(length(dimnames(scenario_value_hr)[[2]])+1,lenth_timevalue)
	#list the hazard ratios
  if (dim(h)[[1]]>2) {
	                    hr<-t(h[-1,]/h[1,])
                     }
  else if (dim(h)[[1]]==2) {
                            hr<-cbind(h[-1,]/h[1,])
                           }
	colnames(hr)<-dimnames(scenario_value_hr)[[2]]
	hr<-cbind(timevalue,hr)
 	#rownames(hr)<-c(timevalue)
	#table_hr<-cbind(hr)
	#dimnames(table_hr)<-list(names(hr),"Haz. Ratio")
	#prmatrix(table_hr)
	hr
}
