StopConditions<-function(t,y,pars,...){
	yroot<-rep(1,6)
	yroot[1]<- y['M']
	yroot[2]<-y['X']-pars['Xlim']
	yroot[3]<- y['I']
	if(pars['pH0']!=0){
		yroot[4]<-y['HA']-1e-8
#		yroot[5]<-y['M']-y['HA']
	}
	return(yroot)
}
