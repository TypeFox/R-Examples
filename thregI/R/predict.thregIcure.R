### The following is the definition of the predict.thregIcure function
"predict.thregIcure" <-
function (object,timevalue,scenario,...)
{
	para <- match.call(expand.dots = FALSE)

	m_timevalue<-match(c("timevalue"), names(para), 0)



	m_scenario<-match(c("scenario"), names(para), 0)

        if (!inherits(object, 'thregIcure'))
           stop("Primary argument must be a thregIcure object")


        #if "timevalue" option is not specified, use the study time
	if (!m_timevalue) timevalue<- model.extract(object$mf, "response")[,1]

	if(!is.numeric(timevalue)) {
		stop(paste("'timevalue' option must specify a numerical value!"))
	}

	if (m_scenario)
	{
	        #if scenario is specified

		cur_scenario_string<-para[[m_scenario]]
		#extract covariate names and values from scenario
		scenario_value<-NULL
		scenario_covariate<-NULL
		while(length(cur_scenario_string)>=2) {

			if(length(cur_scenario_string)==3) {
				current_scenario_covariate_string<-cur_scenario_string[[3]]
				if(length(current_scenario_covariate_string)==2) {
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
			else if(length(cur_scenario_string)==2) {
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

		#add intercept term for calculation
		scenario_covariate<-c("(Intercept)",scenario_covariate)
		scenario_value<-c(1,scenario_value)

		names(scenario_value)<-scenario_covariate




		#judge if any covariate for lny0 lacks of scenario value
		judgematrix_lny0<-matrix(rep(object$lny0,times=length(scenario_value)),ncol=length(scenario_value))==matrix(rep(names(scenario_value),times=length(object$lny0)),nrow=length(object$lny0),byrow=TRUE)


		if (any(!apply(judgematrix_lny0,1,any)))
		       stop(paste("scenario value for",  object$lny0[!apply(judgematrix_lny0,1,any)][1], "is required!"))

		#choose the right position of scenario values to combine
		positionpick_lny0=judgematrix_lny0%*%c(1:length(scenario_value))

		lny0<-scenario_value[positionpick_lny0]%*%object$coef[1:length(object$lny0)]
		y0<-exp(lny0)


		#judge if any covariate for mu lacks of scenario value
		judgematrix_mu<-matrix(rep(object$mu,times=length(scenario_value)),ncol=length(scenario_value))==matrix(rep(names(scenario_value),times=length(object$mu)),nrow=length(object$mu),byrow=TRUE)
		if (any(!apply(judgematrix_mu,1,any)))
		       stop(paste("scenario value for",  object$mu[!apply(judgematrix_mu,1,any)][1], "is required!"))
 	#choose the right position of scenario values to combine
		positionpick_mu=judgematrix_mu%*%c(1:length(scenario_value))
		mu<-scenario_value[positionpick_mu]%*%object$coef[(length(object$lny0)+1):(length(object$lny0)+length(object$mu))]

		#.............................................
		#judge if any covariate for lamda lacks of scenario value
		judgematrix_lamda<-matrix(rep(object$lamda,times=length(scenario_value)),ncol=length(scenario_value))==matrix(rep(names(scenario_value),times=length(object$lamda)),nrow=length(object$lamda),byrow=TRUE)
		if (any(!apply(judgematrix_lamda,1,any)))
		      stop(paste("scenario value for",  object$lamda[!apply(judgematrix_lamda,1,any)][1], "is required!"))
		#choose the right position of scenario values to combine
		positionpick_lamda=judgematrix_lamda%*%c(1:length(scenario_value))
		lamda<-scenario_value[positionpick_lamda]%*%object$coef[(length(object$lny0)+length(object$mu)+1):(length(object$lny0)+length(object$mu)+length(object$lamda))]
		p0<-exp(lamda)/(1+exp(lamda))
		#........................


		dim(y0)<-NULL
		dim(lny0)<-NULL
		dim(mu)<-NULL
		#...............
		dim(lamda)<-NULL
		dim(p0)<-NULL
		#..................

		y0<-c(rep(y0,length(timevalue)))
		lny0<-c(rep(lny0,length(timevalue)))
		mu<-c(rep(mu,length(timevalue)))
		#...........................
		lamda<-c(rep(lamda,length(timevalue)))
		p0<-c(rep(p0,length(timevalue)))
		#................................

	}
	else
	{
	        #if scenario is not specified, use the covariate value of each subject as scenario
                if (m_timevalue & length(timevalue)>1) stop("Please specify only one time value when no scenario is specified!")


        	f1<-formula(object$call[2])
	        f1<-Formula(f1)
		f_lny0 <-formula(f1,lhs=0,rhs=1)
		f_mu <-formula(f1,lhs=0,rhs=2)
                x_lny0<-model.matrix(f_lny0,object$mf)
                x_mu<-model.matrix(f_mu,object$mf)

                #.................................................
                f_lamda<-formula(f1,lhs=0,rhs=3)
                x_lamda<-model.matrix(f_lamda,object$mf)
                #.................................................

		lny0<-x_lny0%*%object$coef[1:length(object$lny0)]
		y0<-exp(lny0)
		mu<-x_mu%*%object$coef[(length(object$lny0)+1):(length(object$lny0)+length(object$mu))]

		lny0<-function(para_lny0){x_lny0%*%para_lny0}
		mu<-function(para_mu){x_mu%*%para_mu}
		##################################################
		lamda<-x_lamda%*%object$coef[(length(object$lny0)+length(object$mu)+1):(length(object$lny0)+length(object$mu)+length(object$lamda))]
		#................................................


		dim(y0)<-NULL
		dim(lny0)<-NULL
		dim(mu)<-NULL
		#...............
		dim(lamda)<-NULL
		dim(p0)<-NULL
		#..................

	}
	lenth_timevalue<-length(timevalue)

	f<-p0*exp((lny0-.5*(log(2*pi*(timevalue^3))+(y0+mu*timevalue)^2/timevalue)))
	S<-p0*exp(log(pnorm((mu*timevalue+y0)/sqrt(timevalue))-exp(-2*y0*mu)*pnorm((mu*timevalue-y0)/sqrt(timevalue))))+(1-p0)
	h<-f/S

        table<-cbind(timevalue,y0,mu,lamda,p0,f,S,h)

	#rownames(table)<-c(timevalue)
        table

}
