"fitCBB" <-
function(x,y,theta0,delta0,copulamodel=c("pCBB1", "pCBB2","pCBB3","pCBB4","pCBB5","pCBB6","pCBB7","pCMax","pCMin"),m,step,deltamin,thetamin,test=c("wilcox.test","t.test"),empcumulative=TRUE,cumulative1,cumulative2,parameters1,parameters2)
{if(empcumulative==TRUE){cumulative1<-pempirical;cumulative2<-pempirical};
if(empcumulative==FALSE){cumulative1<-cumulative1;cumulative2<-cumulative2};
if(missing(m)){m<-15};
if(missing(step)){step<-0.01};
if(missing(deltamin) && copulamodel=="pCBB1"){deltamin<-1};
if(missing(deltamin) && copulamodel=="pCBB2"){deltamin<-0.05};
if(missing(deltamin) && copulamodel=="pCBB3"){deltamin<-1};
if(missing(deltamin) && copulamodel=="pCBB4"){deltamin<-0.05};
if(missing(deltamin) && copulamodel=="pCBB5"){deltamin<-0.05};
if(missing(deltamin) && copulamodel=="pCBB6"){deltamin<-1};
if(missing(deltamin) && copulamodel=="pCBB7"){deltamin<-0.05};
if(missing(deltamin) && copulamodel=="pCMax"){deltamin<-1};
if(missing(deltamin) && copulamodel=="pCMin"){deltamin<-1};
#
if(missing(thetamin) && copulamodel=="pCBB1"){thetamin<-0.05};
if(missing(thetamin) && copulamodel=="pCBB2"){thetamin<-0.05};
if(missing(thetamin) && copulamodel=="pCBB3"){thetamin<-0.05};
if(missing(thetamin) && copulamodel=="pCBB4"){thetamin<-0.05};
if(missing(thetamin) && copulamodel=="pCBB5"){thetamin<-1};
if(missing(thetamin) && copulamodel=="pCBB6"){thetamin<-1};
if(missing(thetamin) && copulamodel=="pCBB7"){thetamin<-1};
if(missing(thetamin) && copulamodel=="pCMax"){thetamin<-1};
if(missing(thetamin) && copulamodel=="pCMin"){thetamin<-1};
#
n<-min(length(x),length(y));
Uni<-matrix(nrow=n,ncol=4);
Uni[1:n,1]<-x;
Uni[1:n,2]<-y;
Uni[1:n,3]<-cumulativemarg(cumulative1,x,parameters1);
Uni[1:n,4]<-cumulativemarg(cumulative2,y,parameters2);
lambdaL1<-FE2(Uni[,3],Uni[,4],0.02,0.02)/0.02;
lambdaU1<-SOB2(Uni[,3],Uni[,4],0.98,0.98)/(1-0.98);
	ini<-c(lambdaL1,lambdaU1);
	inipar<-fitlambdas(ini[1],ini[2])
if(copulamodel=="pCBB1" && inipar[1]$BB1.model=="TRUE BB1" && missing(delta0)){delta0<-inipar[2]$BB1.delta};
if(copulamodel=="pCBB1" && inipar[1]$BB1.model=="FALSE BB1" && missing(delta0)){delta0<-1};
if(copulamodel=="pCBB2" && inipar[4]$BB2.model=="TRUE BB2" && missing(delta0)){delta0<-inipar[5]$BB2.delta};
if(copulamodel=="pCBB2" && inipar[4]$BB2.model=="FALSE BB2" && missing(delta0)){delta0<-0.05};
if(copulamodel=="pCBB3" && inipar[7]$BB3.model=="TRUE BB3" && missing(delta0)){delta0<-inipar[8]$BB3.delta};
if(copulamodel=="pCBB3" && inipar[7]$BB3.model=="FALSE BB3" && missing(delta0)){delta0<-1};
if(copulamodel=="pCBB4" && inipar[10]$BB4.model=="TRUE BB4" && missing(delta0)){delta0<-inipar[11]$BB4.delta};
if(copulamodel=="pCBB4" && inipar[10]$BB4.model=="FALSE BB4" && missing(delta0)){delta0<-0.05};
if(copulamodel=="pCBB5" && inipar[13]$BB5.model=="TRUE BB5" && missing(delta0)){delta0<-inipar[14]$BB5.delta};
if(copulamodel=="pCBB5" && inipar[13]$BB5.model=="FALSE BB5" && missing(delta0)){delta0<-0.05};
if(copulamodel=="pCBB6" && inipar[16]$BB6.model=="TRUE BB6" && missing(delta0)){delta0<-sqrt(inipar[17]$BB6.deltaxtheta)};
if(copulamodel=="pCBB6" && inipar[16]$BB6.model=="FALSE BB6" && missing(delta0)){delta0<-1};
if(copulamodel=="pCBB7" && inipar[18]$BB7.model=="TRUE BB7" && missing(delta0)){delta0<-inipar[19]$BB7.delta};
if(copulamodel=="pCBB7" && inipar[18]$BB7.model=="FALSE BB7" && missing(delta0)){delta0<-0.05};
if(copulamodel=="pCMax" && missing(delta0)){delta0<-1};
if(copulamodel=="pCMin" && missing(delta0)){delta0<-1};
#
if(copulamodel=="pCBB1" && inipar[1]$BB1.model=="TRUE BB1" && missing(theta0)){theta0<-inipar[3]$BB1.theta};
if(copulamodel=="pCBB1" && inipar[1]$BB1.model=="FALSE BB1" && missing(theta0)){theta0<-0.05};
if(copulamodel=="pCBB2" && inipar[4]$BB2.model=="TRUE BB2" && missing(theta0)){theta0<-inipar[6]$BB2.theta};
if(copulamodel=="pCBB2" && inipar[4]$BB2.model=="FALSE BB2" && missing(theta0)){theta0<-0.05};
if(copulamodel=="pCBB3" && inipar[7]$BB3.model=="TRUE BB3" && missing(theta0)){theta0<-inipar[9]$BB3.theta};
if(copulamodel=="pCBB3" && inipar[7]$BB3.model=="FALSE BB3" && missing(theta0)){theta0<-0.05};
if(copulamodel=="pCBB4" && inipar[10]$BB4.model=="TRUE BB4" && missing(theta0)){theta0<-inipar[12]$BB4.theta};
if(copulamodel=="pCBB4" && inipar[10]$BB4.model=="FALSE BB4" && missing(theta0)){theta0<-0.05};
if(copulamodel=="pCBB5" && inipar[13]$BB5.model=="TRUE BB5" && missing(theta0)){theta0<-inipar[15]$BB5.theta};
if(copulamodel=="pCBB5" && inipar[13]$BB5.model=="FALSE BB5" && missing(theta0)){theta0<-1};
if(copulamodel=="pCBB6" && inipar[16]$BB6.model=="TRUE BB6" && missing(theta0)){theta0<-sqrt(inipar[17]$BB6.deltaxtheta)};
if(copulamodel=="pCBB6" && inipar[16]$BB6.model=="FALSE BB6" && missing(theta0)){theta0<-1};
if(copulamodel=="pCBB7" && inipar[18]$BB7.model=="TRUE BB7" && missing(theta0)){theta0<-inipar[20]$BB7.theta};
if(copulamodel=="pCBB7" && inipar[18]$BB7.model=="FALSE BB7" && missing(theta0)){theta0<-1};
if(copulamodel=="pCMax" && missing(theta0)){theta0<-1};
if(copulamodel=="pCMin" && missing(theta0)){theta0<-1};
#
#
#
	deltavec<-matrix(nrow=(2*m),ncol=1);
	deltavec[m]<-delta0;
	for(i in 1:(m))
		deltavec[m+i]<-deltavec[m+i-1]+step;
		for(i in 1:(m-1))
			deltavec[m-i]<-max(deltavec[m-i+1]-step,deltamin);
			diracdeuno<-function(w,deltamin)
			{
				if(w==deltamin){diracdeuno<-1};
					if(w!=deltamin){diracdeuno<-0};
						resu<-diracdeuno; 			}
			SecdiracUno<-matrix(nrow=m,ncol=1);
			for(i in 1:m)
			{
			SecdiracUno[i]<-diracdeuno(deltavec[i],deltamin)
			}
	mun<-sum(SecdiracUno[1:m]);
	if(mun!=0){deltaneto<-matrix(nrow=(2*m-mun+1),ncol=1);
	deltaneto[1:(2*m-mun+1)]<-deltavec[mun:(2*m)];
	dimdelta<-2*m-mun+1}
if(mun==0){deltaneto<-deltavec;dimdelta<-2*m}
#
thetavec<-matrix(nrow=(2*m),ncol=1);
	thetavec[m]<-theta0;
	for(i in 1:(m))
	{
		thetavec[(m+i)]<-thetavec[(m+i-1)]+step;
	}
		for(i in 1:(m-1))
	{
		thetavec[(m-i)]<-max(thetavec[(m-i+1)]-step,thetamin);
	}
		SecdiracEpsilon<-matrix(nrow=m,ncol=1);
			for(i in 1:m)
			{
			SecdiracEpsilon[i]<-diracdeuno(thetavec[i]-thetamin+deltamin,deltamin)
			}
	mepsilon<-sum(SecdiracEpsilon[1:m]);
	if(mepsilon!=0){thetaneto<-matrix(nrow=(2*m-mepsilon+1),ncol=1);
	thetaneto[1:(2*m-mepsilon+1)]<-thetavec[mepsilon:(2*m)];
	dimtheta<-2*m-mepsilon+1}
if(mepsilon==0){thetaneto<-thetavec;dimtheta<-2*m}
	pvec<-matrix(nrow=dimtheta,ncol=dimdelta);
#
	Cumulative<-matrix(nrow=n,ncol=1);
	for(i in 1:n)
	{Cumulative[i]<-FE2(Uni[1:n,3],Uni[1:n,4],Uni[i,3],Uni[i,4])}
		for(j in 1:dimdelta)
	for(i in 1:dimtheta)
	{
		pvec[i,j]<-unlist(ftest(Cumulative,fcopulamodel(thetaneto[i],deltaneto[j],Uni[1:n,3],Uni[1:n,4],model=copulamodel),test)$p.value);
	}
	posicionpmax<-order(pvec)[dimdelta*dimtheta];
	x<-as.integer(posicionpmax/dimtheta);
	if(x!=(posicionpmax/dimtheta)){pvalor<-pvec[(posicionpmax-dimtheta*x),(x+1)];delta<-deltaneto[x+1]; theta<-thetaneto[(posicionpmax-dimtheta*x)]};
	if(x==(posicionpmax/dimtheta)){pvalor<-pvec[(posicionpmax-dimtheta*(x-1)),(x)];delta<-deltaneto[x]; theta<-thetaneto[(posicionpmax-dimtheta*(x-1))]}
	mejorcopula<-fcopulamodel(theta,delta,Uni[1:n,3],Uni[1:n,4],model=copulamodel);
		resu<-list(Empirical=matrix(c(Cumulative),nrow=n),Copula=matrix(c(mejorcopula),nrow=n),fit=c(p.value=pvalor,delta=delta,theta=theta),thetai=matrix(c(thetaneto),nrow=dimtheta,ncol=1),deltaj=matrix(c(deltaneto),nrow=dimdelta,ncol=1),pthetaideltaj=matrix(c(pvec),nrow=dimtheta,ncol=dimdelta))
}

