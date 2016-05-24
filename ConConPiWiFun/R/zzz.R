# TODO: Add comment
# 
# Author: robin
###############################################################################

loadModule("mod_cplfunction", TRUE)
loadModule("mod_cpqfunction", TRUE)
evalqOnLoad({
			setMethod( "plot", signature(x="Rcpp_cplfunction",y="ANY") , function(x,y,...) {      
						firstBreakVal=x$FirstBreakVal_;
						tmp=x$get_BreakPoints_();
						if (tmp$Breakpoints[1]==-Inf){
							tmp$Breakpoints[1]=tmp$Breakpoints[2]-1;
							firstBreakVal=   firstBreakVal -tmp$Slopes[1];           
						}
						if (tmp$Breakpoints[length(tmp$Breakpoints)]==Inf){
							tmp$Breakpoints[length(tmp$Breakpoints)]=tmp$Breakpoints[length(tmp$Breakpoints)-1]+1
						}
						xx=tmp$Breakpoints; 
						yy=array(NA,length(xx));
						yy[1]=firstBreakVal;
						for (i in 2:length(yy)){
							yy[i]=yy[i-1]+tmp$Slopes[i-1]*(tmp$Breakpoints[i]-tmp$Breakpoints[i-1]);
						}
						plot(xx,yy,type='l',...);
					} )

			setMethod( "plot", signature(x="Rcpp_cpqfunction",y="ANY") , function(x,y,...) {
						
						firstBreakVal=x$FirstBreakVal_;
						tmp=x$get_BreakPoints_();
						if (tmp$Breakpoints[1]==-Inf){
							tmp$Breakpoints[1]=tmp$Breakpoints[2]-1;
							firstBreakVal= firstBreakVal - (tmp$Slopes1[1]/2-tmp$Slopes0[1])*(2*tmp$Breakpoints[2]-1)-tmp$Slopes0[1]
						}
						if (tmp$Breakpoints[length(tmp$Breakpoints)]==Inf){
							tmp$Breakpoints[length(tmp$Breakpoints)]=tmp$Breakpoints[length(tmp$Breakpoints)-1]+1
						}
						xx=tmp$Breakpoints;
						last_break_val=firstBreakVal
						
						x_list<-NULL
						y_list<-NULL
						for (i in 1:(length(xx)-1)){
							x_coor<-seq(xx[i],xx[i+1],length.out=50)
							a<-(tmp$Slopes1[i]-tmp$Slopes0[i])/2
							b<-tmp$Slopes0[i]
							c<-last_break_val-a*xx[i]**2-b*xx[i]
							last_break_val=a*xx[i+1]**2+b*xx[i+1]+c
							y_coor<-a*x_coor**2+b*x_coor+c
							x_list<-c(x_list,x_coor)
							y_list<-c(y_list,y_coor)
						}
						plot(x_list,y_list,type='l',...)
					} )
			
			setMethod( "show", "Rcpp_cplfunction" , function(object) {      
						cat('\n')
						cat('Value of f at first non infinite break: ',object$FirstBreakVal_,'\n');
						print(object$get_BreakPoints_());
					} )
			
			setMethod( "show", "Rcpp_cplfunction2" , function(object) {      
						cat('\n')
						cat('Value of f at first non infinite break: ',object$FirstBreakVal_,'\n');
						print(object$get_BreakPoints_());
					} )
			
			setMethod( "show", "Rcpp_cpqfunction" , function(object) {      
						cat('\n')
						cat('Value of f at first non infinite break: ',object$FirstBreakVal_,'\n');
						print(object$get_BreakPoints_());
					} )
			
		})