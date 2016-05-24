"OptimCBB" <-
function(x,y,m,step,test=c("wilcox.test","t.test"),empcumulative=TRUE,cumulative1,cumulative2,parameters1,parameters2)
{if(empcumulative==TRUE){cumulative1<-pempirical;cumulative2<-pempirical};
if(empcumulative==FALSE){cumulative1<-cumulative1;cumulative2<-cumulative2};
if(missing(m)){m<-15};
if(missing(step)){step<-0.01};
n<-min(length(x),length(y));
Uni<-matrix(nrow=n,ncol=4);
Uni[1:n,1]<-x;
Uni[1:n,2]<-y;
Uni[1:n,3]<-cumulativemarg(cumulative1,x,parameters1);
Uni[1:n,4]<-cumulativemarg(cumulative2,y,parameters2);
lambdaLE0<-FE2(Uni[,3],Uni[,4],0.02,0.02)/0.02;
lambdaUE0<-SOB2(Uni[,3],Uni[,4],0.98,0.98)/(1-0.98);
	inicio<-c(lambdaLE0,lambdaUE0);
	M<-fitlambdas(inicio[1],inicio[2])
	#
   if(M[1]$BB1.model=="TRUE BB1"){
   optimBB1<-unlist(fitCBB(x,y,M[3]$BB1.theta,M[2]$BB1.delta,"pCBB1",m,step,1,0.05,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
   if(M[1]$BB1.model=="FALSE BB1"){optimBB1<-"FALSE BB1"}
	#
	if(M[4]$BB2.model=="TRUE BB2"){
	optimBB2<-unlist(fitCBB(x,y,0.05,0.05,"pCBB2",m,step,0.05,0.05,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
	if(M[4]$BB2.model=="FALSE BB2"){optimBB2<-"FALSE BB2"}
	#
	if(M[7]$BB3.model=="TRUE BB3"){
		if(M[9]$BB3.theta=="No Calculado"){M[9]$BB3.theta<-0.05;
	optimBB3<-unlist(fitCBB(x,y,M[9]$BB3.theta,M[8]$BB3.delta,"pCBB3",m,step,1,0.05,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])};
	   if(M[9]$BB3.theta!="No Calculado")	{optimBB3<-unlist(fitCBB(x,y,M[9]$BB3.theta,M[8]$BB3.delta,"pCBB3",m,step,1,0.05,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
    }
	if(M[7]$BB3.model=="FALSE BB3"){optimBB3<-"FALSE BB3"}
	#
	if(M[10]$BB4.model=="TRUE BB4"){
	optimBB4<-unlist(fitCBB(x,y,M[12]$BB4.theta,M[11]$BB4.delta,"pCBB4",m,step,0.05,0.05,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
	if(M[10]$BB4.model=="FALSE BB4"){optimBB4<-"FALSE BB4"}
	#
	if(M[13]$BB5.model=="TRUE BB5"){
	optimBB5<-unlist(fitCBB(x,y,M[15]$BB5.theta,M[14]$BB5.delta,"pCBB5",m,step,0.05,1,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
	if(M[13]$BB5.model=="FALSE BB5"){optimBB5<-"FALSE BB5"}
	#
	if(M[16]$BB6.model=="TRUE BB6"){
	optimBB6<-unlist(fitCBB(x,y,sqrt(M[17]$BB6.deltaxtheta),sqrt(M[17]$BB6.deltaxtheta),"pCBB6",m,step,1,1,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
	if(M[16]$BB6.model=="FALSE BB6"){optimBB6<-"FALSE BB6"}
	#
	if(M[18]$BB7.model=="TRUE BB7"){
	optimBB7<-unlist(fitCBB(x,y,M[20]$BB7.theta,M[19]$BB7.delta,"pCBB7",m,step,0.05,1,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
	if(M[18]$BB7.model=="FALSE BB7"){optimBB7<-"FALSE BB7"}
	#
	if(M[21]$CMin.model=="TRUE CMin"){optimCMin<-unlist(fitCBB(x,y,1,1,"pCMin",1,0,0,0,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
	if(M[21]$CMin.model=="FALSE CMin"){optimCMin<-"FALSE CMin"}
	#
	if(M[22]$CMax.model=="TRUE CMax"){optimCMax<-unlist(fitCBB(x,y,1,1,"pCMax",1,0,0,0,test,empcumulative,cumulative1,cumulative2,parameters1,parameters2)[3])}
	if(M[22]$CMax.model=="FALSE CMax"){optimCMax<-"FALSE CMax"}
	Cumulative<-matrix(nrow=n,ncol=1);
	for(i in 1:n)
	{Cumulative[i]<-FE2(Uni[1:n,3],Uni[1:n,4],Uni[i,3],Uni[i,4])};
	#
	nocharacter<-function(character){testcharacter<-is.numeric(character); if(testcharacter==FALSE){character<-"NA"} else {character<-character}};
	optimpvalor<-matrix(as.numeric(c(nocharacter(optimBB1[[1]]),nocharacter(optimBB2[[1]]),nocharacter(optimBB3[[1]]),nocharacter(optimBB4[[1]]),nocharacter(optimBB5[[1]]),nocharacter(optimBB6[[1]]),nocharacter(optimBB7[[1]]),nocharacter(optimCMin[[1]]),nocharacter(optimCMax[[1]]))))
	optimdelta<-matrix(as.numeric(c(optimBB1[2],optimBB2[2],optimBB3[2],optimBB4[2],optimBB5[2],optimBB6[2],optimBB7[2],optimCMin[2],optimCMax[2])))
	optimtheta<-matrix(as.numeric(c(optimBB1[3],optimBB2[3],optimBB3[3],optimBB4[3],optimBB5[3],optimBB6[3],optimBB7[3],optimCMin[3],optimCMax[3])))
	#
	#
	pvalororden<-order(optimpvalor,na.last=F);
	pvaloroptimo<-optimpvalor[pvalororden][9];
	deltaoptimo<-optimdelta[pvalororden][9];
	thetaoptimo<-optimtheta[pvalororden][9];
	Familia<-matrix(c("BB1","BB2","BB3","BB4","BB5","BB6","BB7","CMin","CMax"),nrow=9);
	Familia<-Familia[pvalororden][9]
	if(Familia=="BB1"){mejorcopula<-pCBB1(thetaoptimo,deltaoptimo,Uni[1:n,3],Uni[1:n,4])};
 	if(Familia=="BB2"){mejorcopula<-pCBB2(thetaoptimo,deltaoptimo,Uni[1:n,3],Uni[1:n,4])};
 	if(Familia=="BB3"){mejorcopula<-pCBB3(thetaoptimo,deltaoptimo,Uni[1:n,3],Uni[1:n,4])};
 	if(Familia=="BB4"){mejorcopula<-pCBB4(thetaoptimo,deltaoptimo,Uni[1:n,3],Uni[1:n,4])};
 	if(Familia=="BB5"){mejorcopula<-pCBB5(thetaoptimo,deltaoptimo,Uni[1:n,3],Uni[1:n,4])};
 	if(Familia=="BB6"){mejorcopula<-pCBB6(thetaoptimo,deltaoptimo,Uni[1:n,3],Uni[1:n,4])};
 	if(Familia=="BB7"){mejorcopula<-pCBB7(thetaoptimo,deltaoptimo,Uni[1:n,3],Uni[1:n,4])};
 	if(Familia=="CMin"){mejorcopula<-pCMin(1,1,Uni[1:n,3],Uni[1:n,4])};
	if(Familia=="CMax"){mejorcopula<-pCMax(1,1,Uni[1:n,3],Uni[1:n,4])};
		MSE<-mean((Cumulative-as.matrix(mejorcopula))^2)
		#
	resu<-list(Empirical=Cumulative,Copula=matrix(c(mejorcopula),nrow=n),OptimumFit=c(Family=Familia,p.value=pvaloroptimo,delta=deltaoptimo,theta=thetaoptimo,MSE=MSE),
	Initial.BB1=c(BB1=M[1]$BB1.model,BB1.delta=M[2]$BB1.delta,BB1.theta=M[3]$BB1.theta),Final.BB1=c(optimBB1=optimBB1),
	Initial.BB2=c(BB2=M[4]$BB2.model,BB2.delta=M[5]$BB2.delta,BB2.theta=M[6]$BB2.theta),Final.BB2=c(optimBB2=optimBB2),
	Initial.BB3=c(BB3=M[7]$BB3.model,BB3.delta=M[8]$BB3.delta,BB3.theta=M[9]$BB3.theta),Final.BB3=c(optimBB3=optimBB3),
	Initial.BB4=c(BB4=M[10]$BB4.model,BB4.delta=M[11]$BB4.delta,BB4.theta=M[12]$BB4.theta),Final.BB4=c(optimBB4=optimBB4),
	Initial.BB5=c(BB5=M[13]$BB5.model,BB5.delta=M[14]$BB5.delta,BB5.theta=M[15]$BB5.theta),Final.BB5=c(optimBB5=optimBB5),
	Initial.BB6=c(BB6=M[16]$BB6.model,BB6.deltaxtheta=M[17]$BB6.deltaxtheta),Final.BB6=c(optimBB6=optimBB6),
	Initial.BB7=c(BB7=M[18]$BB7.model,BB7.delta=M[19]$BB7.delta,BB7.theta=M[20]$BB7.theta),Final.BB7=c(optimBB7=optimBB7),
	Initial.CMin=c(CMin=M[21]$CMin.model),Final.CMin=c(optimCMin=optimCMin),
	Initial.CMax=c(CMax=M[22]$CMax.model),Final.CMax=c( optimCMax=optimCMax))
}

