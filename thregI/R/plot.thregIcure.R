### The following is the definition of the plot.thregIcure function
"plot.thregIcure" <-
function (x,var,scenario,graph,nolegend=0,nocolor=0,...)
{
	para <- match.call(expand.dots = FALSE)

        indx <- match(c("object", "var", "timevalue", "scenario"),
		  names(para), nomatch=0)

	if (missing(x)) stop("An x argument is required")
	if (missing(var)) stop("A var argument is required")
	if (missing(graph)) stop("A graph argument is required")

        timevalue<- sort(model.extract(x$mf, "response")[,1])
	names(timevalue)<-NULL
	m_scenario<-match(c("scenario"), names(para), 0)
	m_hr<-match(c("var"), names(para), 0)
	m_graph<-match(c("graph"), names(para), 0)
        scenario_value<-NULL
        scenario_covariate<-NULL
        if(m_scenario!=0) {
		cur_scenario_string<-para[[m_scenario]]
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

        }
        #add intercept term
	scenario_covariate<-c("(Intercept)",scenario_covariate)
	scenario_value<-c(1,scenario_value)

	names(scenario_value)<-scenario_covariate




        m_hr<-match(c("var"), names(para), 0)
	hr_var<-as.character(para[[m_hr]])
        if (!is.factor(x$mf[[hr_var]])) stop("The variable for the hazard ratio calculation should be a factor variable")

        if (length(levels(x$mf[[hr_var]]))<2) stop("The variable for the hazard ratio calculation should have more than one level")


        if (!any(names(x$mf)==hr_var)) stop("The variable ",as.character(para[[m_hr]]) , " for the hazard ratio calculation should be included in the model")




        scenario_value_hr<-rbind(rep(0,length(levels(x$mf[[hr_var]]))-1),diag(rep(1,length(levels(x$mf[[hr_var]]))-1)))
	#names(scenario_value_hr)<-paste(hr_var,levels(x$mf[[hr_var]])[-1],sep="")
        colnames(scenario_value_hr) <- paste(hr_var,levels(x$mf[[hr_var]])[-1],sep="")
        if (any(matrix(rep(names(scenario_value),times=length(dimnames(scenario_value_hr)[[2]])),nrow=length(dimnames(scenario_value_hr)[[2]]),byrow=TRUE)==dimnames(scenario_value_hr)[[2]])) stop("Don't specify the scenario value of the dummy variable for the hazard ratio calculation!")


        scenario_value<-matrix(rep(scenario_value,times=length(levels(x$mf[[hr_var]]))),nrow=length(levels(x$mf[[hr_var]])),byrow=TRUE)

        dimnames(scenario_value)[[2]] <- as.list(scenario_covariate)


	scenario_value<-cbind(scenario_value,scenario_value_hr)



        #judge if any covariate for lny0 lacks of scenario value
	judgematrix_lny0<-matrix(rep(x$lny0,times=dim(scenario_value)[2]),ncol=dim(scenario_value)[2])==matrix(rep(dimnames(scenario_value)[[2]],times=length(x$lny0)),nrow=length(x$lny0),byrow=TRUE)


	if (any(!apply(judgematrix_lny0,1,any)))
               stop(paste("scenario value for",  x$lny0[!apply(judgematrix_lny0,1,any)][1], "is required!"))

	#choose the right position of scenario values to combine
	positionpick_lny0=judgematrix_lny0%*%c(1:dim(scenario_value)[2])

        lny0<-scenario_value[,positionpick_lny0]%*%x$coef[1:length(x$lny0)]
	y0<-exp(lny0)


        #judge if any covariate for mu lacks of scenario value
	judgematrix_mu<-matrix(rep(x$mu,times=dim(scenario_value)[2]),ncol=dim(scenario_value)[2])==matrix(rep(dimnames(scenario_value)[[2]],times=length(x$mu)),nrow=length(x$mu),byrow=TRUE)


	if (any(!apply(judgematrix_mu,1,any)))
               stop(paste("scenario value for",  x$mu[!apply(judgematrix_mu,1,any)][1], "is required!"))

	#choose the right position of scenario values to combine
	positionpick_mu=judgematrix_mu%*%c(1:dim(scenario_value)[2])

  mu<-scenario_value[,positionpick_mu]%*%x$coef[(length(x$lny0)+1):(length(x$lny0)+length(x$mu))]

  #.............................................
  #judge if any covariate for lamda lacks of scenario value
  judgematrix_lamda<-matrix(rep(x$lamda,times=dim(scenario_value)[2]),ncol=dim(scenario_value)[2])==matrix(rep(dimnames(scenario_value)[[2]],times=length(x$lamda)),nrow=length(x$lamda),byrow=TRUE)
  if (any(!apply(judgematrix_lamda,1,any))) stop(paste("scenario value for",  x$lamda[!apply(judgematrix_lamda,1,any)][1], "is required!"))
  #choose the right position of scenario values to combine
  positionpick_lamda=judgematrix_lamda%*%c(1:dim(scenario_value)[2])
  lamda<-scenario_value[,positionpick_lamda]%*%x$coef[(length(x$lny0)+length(x$mu)+1):(length(x$lny0)+length(x$mu)+length(x$lamda))]
  p<-exp(lamda)/(1+exp(lamda))
  #........................


	dim(y0)<-NULL
	dim(lny0)<-NULL
	dim(mu)<-NULL
	#...............
	dim(lamda)<-NULL
	dim(p)<-NULL
	#..................

	y0<-c(rep(y0,length(timevalue)))
	lny0<-c(rep(lny0,length(timevalue)))
	mu<-c(rep(mu,length(timevalue)))
	#...........................
	lamda<-c(rep(lamda,length(timevalue)))
	p<-c(rep(p,length(timevalue)))
	#................................
	lenth_timevalue<-length(timevalue)

	#replicate values for matrix calculation
	timevalue_rep<-rep(timevalue,each=length(dimnames(scenario_value_hr)[[2]])+1)

	f<-p*exp((lny0-.5*(log(2*pi*(timevalue_rep^3))+(y0+mu*timevalue_rep)^2/timevalue_rep)))
	S<-p*exp(log(pnorm((mu*timevalue_rep+y0)/sqrt(timevalue_rep))-exp(-2*y0*mu)*pnorm((mu*timevalue_rep-y0)/sqrt(timevalue_rep))))+(1-p)
	h<-f/S

  if(m_graph!=0) 	graph_type<-para[[m_graph]]
	if(graph_type=="ds"){
	        dim(f)<-c(length(dimnames(scenario_value_hr)[[2]])+1,lenth_timevalue)
        	dt_graph<-cbind(timevalue,t(f))
		y_label="Estimated f(t)"
		legend_pos="topright"
	}
	else if(graph_type=="sv"){
	        dim(S)<-c(length(dimnames(scenario_value_hr)[[2]])+1,lenth_timevalue)
        	dt_graph<-cbind(timevalue,t(S))
		y_label="Estimated S(t)"
		legend_pos="topright"
	}
	else if(graph_type=="hz"){
	        dim(h)<-c(length(dimnames(scenario_value_hr)[[2]])+1,lenth_timevalue)
        	dt_graph<-cbind(timevalue,t(h))
		y_label="Estimated h(t)"
		legend_pos="topright"
		#legend_pos="topleft"
	}
	if(nocolor==0) {

		matplot(dt_graph[order(timevalue),1],dt_graph[order(timevalue),2:(length(dimnames(scenario_value_hr)[[2]])+2)],type = "l",lty=1:(length(dimnames(scenario_value_hr)[[2]])+1),xlab="time",ylab=y_label)
	}
	else {
		matplot(dt_graph[order(timevalue),1],dt_graph[order(timevalue),2:(length(dimnames(scenario_value_hr)[[2]])+2)],type = "l",col=1,lty=1:(length(dimnames(scenario_value_hr)[[2]])+1),xlab="time",ylab=y_label)
	}


	if(nolegend==0) {
		if(nocolor==0) {
		legend(legend_pos,paste(hr_var,"=",as.character(levels(x$mf[[hr_var]]))),col=1:(length(levels(x$mf[[hr_var]]))),lty=1:(length(levels(x$mf[[hr_var]]))),bty="n")
		}
		else {
		legend(legend_pos,paste(hr_var,"=",as.character(levels(x$mf[[hr_var]]))),col=1,lty=1:(length(levels(x$mf[[hr_var]]))),bty="n")
		}
	}
}
