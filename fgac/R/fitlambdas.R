"fitlambdas" <-
function(lambdaLE,lambdaUE)
{
	if(lambdaUE!="NaN" && lambdaLE!="NaN"){if(lambdaUE>1 | lambdaUE<0){lambdaUE<-NaN}
	if(lambdaLE>1 | lambdaLE<0){lambdaLE<-NaN}}
	#para BB1 son las dos cond. siguientes
	#lambdaUE en [0,1)
	#lambdaLE en (0,1)
	BB1lambdas<-function(lambdaLE,lambdaUE)
	{
		if(lambdaUE=="NaN" | lambdaLE=="NaN"){lambdaUE<-1}
		if(lambdaUE==1 | lambdaLE==0 |lambdaLE==1){lambdaUE<-NaN}
		delta<-log(2)/log(2-lambdaUE);
		theta<--log(2-lambdaUE)/log(lambdaLE);
		modelo<-"TRUE BB1";
		if(delta=="NaN" | theta=="NaN") modelo<-"FALSE BB1";
	resulBB1<-list(model=modelo,delta=delta,theta=theta)
	}
#
#la siguiente familia considera posible el valor NaN para lambdaUE
#mientras que el unico valor de lambdaLE es 1  
#
	BB2lambdas<-function(lambdaLE,lambdaUE)
		{
			if(lambdaLE=="NaN")lambdaLE<-5;
			if(lambdaLE==1 && lambdaUE=="NaN"){modelo<-" TRUE BB2";delta<-"No calculado";theta<-"No calculado"};
			if(lambdaLE!=1 | lambdaUE!="NaN"){modelo<-"FALSE BB2";delta<-NaN;theta<-NaN};
			resulBB2<-list(model=modelo,delta=delta,theta=theta)
		}
		#modelo con delta y theta invertidos(papeles)respecto a texto Joe(1997)
		#
#para BB3 son las dos cond. siguientes
	#lambdaUE en [0,1)
	#lambdaLE en (0,1]
#
#
  BB3lambdas<-function(lambdaLE,lambdaUE)
	{
		if(lambdaLE=="NaN"){lambdaUE<-NaN}
		if(lambdaUE=="NaN"){theta<-NaN; delta<-NaN; modelo<-"FALSE BB3"};
			if(lambdaUE != "NaN" && lambdaUE<1){delta<-log(2)/log(2-lambdaUE);
					if(lambdaLE>0 && lambdaLE <1 && delta==1){theta<--log(2)/log(lambdaLE); modelo<-"TRUE BB3"};
					if(lambdaLE==0 && delta==1){theta<-NaN; modelo<-"FALSE BB3"};
					if(lambdaLE ==1 && delta==1){theta<-NaN; modelo<-"FALSE BB3"};
						if(lambdaLE==1 && delta>1){theta<-"No Calculado"; modelo<-"TRUE BB3"}
						if(lambdaLE!=1 && delta>1){theta<-NaN; modelo<-"FALSE BB3"}}
							if(lambdaUE != "NaN" && lambdaUE>=1){delta<-NaN;theta<-NaN;modelo<-"FALSE BB3"}
										resulBB3<-list(model=modelo,delta=delta,theta=theta)
	}
#
#para BB4 son las dos cond. siguientes
	#lambdaUE en (0,1)
	#lambdaLE en (0,1)
#
	BB4lambdas<-function(lambdaLE,lambdaUE)
	{
		if(lambdaLE=="NaN"){lambdaUE<-NaN}
		if(lambdaUE=="NaN"){theta<-NaN;delta<-NaN;modelo<-"FALSE BB4"};
			if(lambdaUE !="NaN"){
				if(lambdaUE==0 | lambdaUE==1 | lambdaLE==0 | lambdaLE==1){theta<-NaN;delta<-NaN; modelo<-"FALSE BB4"};
				if(lambdaUE>0 && lambdaUE<1 && lambdaLE>0 && lambdaLE <1)
				{delta<--log(2)/log(lambdaUE);
					theta<--log(2-2^(-1/delta))/log(lambdaLE);
					modelo<-"TRUE BB4"}};
		resulBB4<-list(model=modelo,delta=delta,theta=theta)
	}
#
#el suguiente ajuste sera realizado fijando theta=1, esto corresponde a la familiaB6
#en H Joe(1997)pp 152
# condiciones: lambdaUE en (0,1) y lambdaLE=0
#
	BB5lambdas<-function(lambdaLE,lambdaUE)
	{
		if(lambdaLE=="NaN"){lambdaUE<-NaN}
		if(lambdaUE=="NaN"){delta<-NaN; theta<-NaN; modelo<-"FALSE BB5"}
			if(lambdaUE!="NaN" && lambdaUE>=1){lambdaLE<-2;delta<-NaN; theta<-NaN; modelo<-"FALSE BB5"}
			if(lambdaUE!="NaN" && lambdaUE<=0){lambdaLE<-2;delta<-NaN; theta<-NaN; modelo<-"FALSE BB5"}
			if(lambdaUE!="NaN" && lambdaUE<1 && lambdaLE==0){delta<--log(2)/log(lambdaUE); theta<-1; modelo<-"TRUE BB5"}
			if(lambdaUE!="NaN" && lambdaUE<1 && lambdaLE!=0){delta<-NaN; theta<-NaN; modelo<-"FALSE BB5"}
	resulBB5<-list(model=modelo,delta=delta,theta=theta)
	}
#
#para BB6 son las dos cond. siguientes
	#lambdaUE en [0,1)
	#lambdaLE =0
#
	BB6lambdas<-function(lambdaLE,lambdaUE)
	{
		if(lambdaLE=="NaN"){lambdaUE<-NaN}
		if(lambdaUE=="NaN")lambdaLE<-1;
	 	if(lambdaLE==0 && lambdaUE!=1) {deltaxtheta<-log(2)/log(2-lambdaUE);modelo<-"TRUE BB6"} else {deltaxtheta<-NaN; modelo<-"FALSE BB6"};
		if(lambdaLE!=0 | lambdaUE==1){deltaxtheta<-NaN;modelo<-"FALSE BB6"};
				resulBB6<-list(model=modelo,deltaxtheta=deltaxtheta)
	}
#
#para BB7 son las dos cond. siguientes
	#lambdaUE en [0,1)
	#lambdaLE en (0,1)
#
	BB7lambdas<-function(lambdaLE,lambdaUE)
	{
		if(lambdaLE=="NaN"){lambdaUE<-NaN}
		if(lambdaUE=="NaN"){theta<-NaN;delta<-NaN; modelo<-"FALSE BB7"};
			if(lambdaUE!="NaN"){if(lambdaLE>0 && lambdaLE<1 && lambdaUE!=1){delta<--log(2)/log(lambdaLE);
		theta<-log(2)/log(2-lambdaUE);modelo<- "TRUE BB7"};
		if(lambdaLE==0 | lambdaLE==1 | lambdaUE==1){delta<-NaN;theta<-NaN;modelo<-"FALSE BB7"}}
						resulBB7<-list(model=modelo,delta=delta,theta=theta)
	}
#  
 CMinlambdas<-function(lambdaLE,lambdaUE)
{ 
	if(lambdaLE=="NaN"){modelo<-"FALSE CMin"}
	if(lambdaUE=="NaN"){modelo<-"FALSE CMin"}
		if(lambdaUE!="NaN" && lambdaLE!="NaN"){
	if(lambdaLE==0 && lambdaUE==0){modelo<-"TRUE CMin"} else {modelo<-"FALSE CMin"}}
		resu<-list(model=modelo)
}
CMaxlambdas<-function(lambdaLE,lambdaUE)
{ 
	if(lambdaLE=="NaN"){modelo<-"FALSE CMax"}
	if(lambdaUE=="NaN"){modelo<-"FALSE CMax"}
		if(lambdaUE!="NaN" && lambdaLE!="NaN"){
	if(lambdaLE==1 && lambdaUE==1){modelo<-"TRUE CMax"} else {modelo<-"FALSE CMax"}}
		resu<-list(model=modelo)
}
resulfitlambdas<-c(BB1=BB1lambdas(lambdaLE,lambdaUE),BB2=BB2lambdas(lambdaLE,lambdaUE),BB3=BB3lambdas(lambdaLE,lambdaUE),
BB4=BB4lambdas(lambdaLE,lambdaUE),BB5=BB5lambdas(lambdaLE,lambdaUE),BB6=BB6lambdas(lambdaLE,lambdaUE),
BB7=BB7lambdas(lambdaLE,lambdaUE),CMin=CMinlambdas(lambdaLE,lambdaUE),CMax=CMaxlambdas(lambdaLE,lambdaUE))	
}

