ranktes<-function(r, n, index = "spearman", approx = "exact", CC = FALSE, type = "two-sided", print = TRUE){
# Computes the p-value of an observed rank correlation
	 if (!is.numeric(r)) {cat("correlation: ",r,"\n")
	 		stop("Non-numeric argument to mathematical function")}
	 if (!is.numeric(n)) {cat("number of ranks: ",n,"\n")
	 		stop("Non-numeric argument to mathematical function")}
	 if (!is.logical(CC)){type<-CC;CC<-FALSE}
	 if (!(approx=="exact") | index=="r4")	{CC<-FALSE}
	 if (n<5){stop("For n<5 rank correlation tests of independence are not always reliable.")}
	 if (abs(r)>1) {cat("Observed value: ",r,"\n")
	 					 stop("The value of the rank correlation should range between 1 to -1")}
	 cifer<-9;options(digits=9)
	 ain<-c("spearman","kendall","gini","r4","fy","filliben")
	 apx<-c("gaussian","student","vggfr","exact")
	 alter<-c("two-sided", "greater","less")
	 index<-tolower(index);approx<-tolower(approx);type=tolower(type)
	 index<-match.arg(index, ain, several.ok = TRUE) 
	 approx<-match.arg(approx, apx, several.ok = TRUE)
	 type<-match.arg(type, alter, several.ok = TRUE)
	 (opr <- mpfr_default_prec())
	 stopifnot(opr == (oprec <- mpfr_default_prec(300)), 300  == mpfr_default_prec())
	 mpfr_default_prec(opr)
	 ksw<-FALSE
	 if ((index=="spearman" & n>26) | (index=="kendall" & n>60) | (index=="gini" & n>24)) {ksw<-TRUE}
	 if ((index=="r4" & n>15) | (index=="fy" & n>15) | (index=="filliben" & n>15)) {ksw<-TRUE}
	 if (approx=="exact" & ksw){
		    cat("The exact p-value is not available yet. An approximate p-value is computed \n")
			if (index=="r4") {approx<-"vggfr"}
			if (index=="fy" | index=="filliben") {approx<-"gaussian"}
			if (index=="spearman" |index=="kendall") {approx<-"gaussian"}
			if (index=="gini") {approx<-"student"}
			}
#
	C1<-mpfr(25,320);C2<-mpfr(38,320);C3<-mpfr(35,320);C4<-mpfr(72,320)
	C5<-mpfr(5,320);C6<-mpfr(9,320);C7<-mpfr(100,320);C8<-mpfr(328,320);C9<-mpfr(127,320)
	C0<-mpfr(997,320);C11<-mpfr(372,320);C12<-mpfr(1350,320);C13<-mpfr(111,320)
	C14<-mpfr(153,320);C15<-mpfr(366,320);C16<-mpfr(301,320);C17<-mpfr(456,320)
	C18<-mpfr(912,320);C19<-mpfr(1248,320);C20<-mpfr(107,320);C21<-mpfr(4,320)
	C22<-mpfr(76,320);C23<-mpfr(182,320);C24<-mpfr(307,320);C25<-mpfr(315,320)
	C26<-mpfr(342,320);C27<-mpfr(420,320);C28<-mpfr(315,320);C29<-mpfr(105,320)
	C32<-mpfr(1.00762,320);C33<-mpfr(2.01524,320);C34<-mpfr(1,320);C35<-mpfr(10,320)
	C37<-mpfr(0.5,320);C38<-mpfr(3,320);C39<-mpfr(2,320);C40<-mpfr(6,320)
	C41<-mpfr(1.5,320);C42<-mpfr(8,320);C43<-mpfr(0.6,320)
	nu<-mpfr(n,320);ccf<-mpfr(0,320)
	kn<-n%%2;kkn<-mpfr(kn,320)
	nm1<-nu-1;nm2<-nu-2;nm3<-nu-3;np1<-nu+1;np2<-nu+2
	ccf<-0
	if (CC & index=="spearman") {
			S1<-trunc((nu*nu*nu-nu)/C40)-(1-r)*(nu*nu*nu-nu)/C40
			if (S1<0){sig<-1} else {sig<- -1}
			ccf<- sig*C40/(nu*nu*nu-nu)} 
	if (CC & index=="gini") {
			S1<-trunc(0.25*r*(nu*nu-kkn))
			if (S1<0){sig<- 1} else {sig<- -1}
			ccf<- sig*C39/(nu^2-kkn)
			} 
	 if (CC & index=="kendall") {
			S1<-trunc(C37*r*nu*nm1)
			if (S1<0){sig<- 1} else {sig<- -1}
			ccf<- sig*C39/(nu*nm1)}
	 if (CC & index=="r4"| index=="fy" | index=="filliben") {ccf<-0} 
	 rc<-mpfr(r,320)-ccf;r<-asNumeric(rc) 
	 Medun<-rep(0,n)
	 if (index=="filliben"){k<-0
	  		for (i in 1:n){k<-k+1;Medun[k]<-qbeta(0.5,i,n+1-i)}
			Medun<-as.double(Medun)}
	if (approx=="exact") {rc<- -abs(rc)}
#
    Vian<-function(x,L1,L2){
		La1<-lgamma(1/L1);Lb1<-lgamma(C34+L2)
		Lc1<-lgamma(1/L1+C34+L2);d1<-exp(La1+Lb1-Lc1)
		L1*(C34-abs(x)^L1)^L2/(C39*d1)}
#
		if (approx=="gaussian"){
				eappr<-"Gaussian"
				if ("spearman" == index) {zx<-rc*sqrt(nm1);estat<-"Spearman's rho"}
		 		if ("gini" == index) {zx<-rc*sqrt(C41*nu);estat<-"Gini's cograduation coefficient"}
			 	if ("kendall" == index) {zx<-rc*sqrt((C6*nu*nm1)/(C21*nu+C35));estat<-"Kendall's tau"}	
		 		if ("r4" == index) {zx<-rc*sqrt(C32)*sqrt(nm1);estat<-"r4"}
		 		if ("fy"==index){estat<-"Fisher-Yates coefficient";zx<-rc*sqrt(nm1)}
		 		if ("filliben"==index){estat<-"Filliben's rank correlation";zx<-rc*sqrt(nm1)}					
				if (type=="two-sided") {Pv<-2*pnorm(-as.numeric(abs(zx)), mean = 0, sd=1)}
		 		if (!(type=="two-sided")) {Pv<-pnorm(-as.numeric(abs(zx)), mean = 0, sd=1)}
		 		if((r<=0  & type=="greater") | (r>0 & type=="less")) {Pv<-1-Pv}
				r<-trunc(r*100000)/100000
				c1<-trunc(as.numeric(Pv)*1000000000)/1000000000
				if (print){
							cat("Statistic:",estat,"-- Approximation: ",eappr," Alternative :",type," Contin. Adjust. :",CC,"\n")
				 		    cat("n=",n,"--- Observed value:",noquote(sprintf("%.5f",r)),"Est. p-value=",
				 		    noquote(sprintf("%.6f",c1)),"\n")}
				outrank<-list(n=n,Statistic=estat,Value=r,Approx=eappr,Tails=type, Cpv=c1,Lpv=c1)
				return(outrank)}
#	   
		if (approx =="student"){eappr<-"t-Student"
				if ("spearman"== index){estat<-"Spearman's rho"
							zx<-nm2/(C34-rc^2);rcp<-abs(rc)*sqrt(zx);Lamx<-nm2}
				if ("gini" == index) {estat<-"Gini's cograduation coefficient"					
							if (kn==0) {mu2n<-C39*(nu^2+C39)/(C38*nm1*nu*nu)}
							else {mu2n<-C39*(nu^2+C38)/(C38*np1*nm1*nm1)}
							mstar<-(C34/mu2n-C34)*C37
							zx<-C39*mstar/(C34-rc^2);rcp<-abs(rc)*sqrt(zx)
				    		Lamx<-C39*trunc(mstar+C37)}				    						    					    		
				if ("kendall" == index){estat<-"Kendall's tau"
								   z1<-(C21*nu+C35)/(C6*nu*nm1);Lm<-C37*(C34/z1-C34)
			 					   rcp<- abs(rc)*sqrt(C39*Lm/(C34-rc^2));Lamx<-C39*trunc(Lm)}
			 	if ("r4" == index){estat<-"r4"
						    ws1<-(nu-C32)/C33;zx<-C39*ws1/(C34-rc^2);rcp<-abs(rc)*sqrt(zx)
							Lamx<-trunc(C39*ws1)}									   				   										   
		 		if ("fy"== index){estat<-"Fisher_Yates coefficient"
							zx<-nm2/(C34-rc^2);rcp<- abs(rc)*sqrt(zx);Lamx<-nm2}
			    if ("filliben"== index){estat<-"Filliben's rank correlation"
							zx<-nm2/(C34-rc^2);rcp<- abs(rc)*sqrt(zx);Lamx<-nm2}
				Lamx<-asNumeric(Lamx);rcp<-asNumeric(rcp)
				if (type=="two-sided") {Pv<-2*pt(-rcp,Lamx)}
			 	if (!(type=="two-sided")) {Pv<-pt(-abs(rcp),Lamx)}
			 	if((r<=0  & type=="greater") | (r>0 & type=="less")) {Pv<-1-Pv}
				r<-trunc(r*100000)/100000
				c1<-trunc(Pv*1000000000)/1000000000
				if (print){
							cat("Statistic:",estat,"-- Approximation: ",eappr," Alternative :",type," Contin. Adjust. :",CC,"\n")
				 		    cat("n=",n,"--- Observed value:",noquote(sprintf("%.5f",r)),"Est. p-value=",
				 		    noquote(sprintf("%.6f",c1)),"\n")}
				outrank<-list(n=n,Statistic=estat,Value=r,Approx=eappr,Tails=type, Cpv=c1,Lpv=c1)
				return(outrank)}
# 
			if (approx=="vggfr") {				
				eappr<-"Vianelli-GGFR"
				if ("spearman"== index){estat<-"Spearman's rho";icoef<-1
													  mu2n<-C34/nm1
													  mu4n=C38*(C1*nu^3-C2*nu^2-C3*nu+C4)
			   	 									  mu4n<-mu4n/(C1*nu*np1*nm1^3)}	  			   	 									  
				if ("gini"== index){estat<-"Gini's cograduation coefficient";icoef<-2
							if (kn==0) {mu2n<-C39*(nu^2+C39)/(C38*nm1*nu*nu)}
									 else {mu2n<-C39*(nu^2+C38)/(C38*np1*nm1*nm1)}
							if (kn==0) {ws3<-C29*(nu^7)*nm1*nm3
											 ws1<-C3*nu^7-C13*nu^6+C14*nu^5-C15*nu^4+C16*nu^3-C17*nu^2-C18*nu+C19}
								else  {ws3<-C29*(np1^3)*nu*(nm1^4)*nm2
									 	  ws1<-C3*nu^7-C22*nu^6+C23*nu^5-C24*nu^4+C28*nu^3-C26*nu^2-C27*nu-C28}
							mu4n<-C21*ws1/ws3}						
				if  ("kendall"== index) {estat<-"Kendall's tau";icoef<-3
 	 											mu2n<-C39*(C39*nu+C5)/(C6*nu*nm1)
						 						mu4n<-C7*nu^4+C8*nu^3-C9*nu^2-C0*nu-C11
												mu4n<-mu4n/(C12*(nu*nm1/2)^3)}
				if ("r4"== index){estat<-"r4";icoef<-4
												mu2n<-C32/nm1
											    Dc1<-mpfr(1.0949159471,320);Dc2<-mpfr(38.7820781157,320)
											    Dc3<-mpfr(208.8267798530,320);Dc4<-mpfr(396.3338168921,320)
												mu4n<- Dc1/sqrt(nm1)+Dc2/nm1/nm1-Dc3/nm1/nm1/nm1+Dc4/nm1/nm1/nm1/nm1}	
			   if ("fy"== index){estat<-"Fisher-Yates coefficient";icoef<-5
											    mu2n<-C34/nm1;EvosG<-rep(0,n)
												for (i in 1:n){EvosG[i]<-evNormOrdStatsScalar(r=i, n=n,approximate=FALSE)}
												k2a<-sum(EvosG^2);k2b<-sum(EvosG^4);k2<-k2a/(n-1)
												k4<-n*((n+1)*k2b-(3*(n-1)/n)*k2a^2)/(n-1)/(n-2)/(n-3)
												k4a<-3*(n-1)/(n+1)+((n-2)*(n-3)/(n*(n^2-1)))*(k4/k2^2)^2
												k4a<-k4a/(n-1)/(n-1);mu4n<-mpfr(k4a,320)}
				if ("filliben"== index){estat<-"Filliben's rank correlation";icoef<-6
											    mu2n<-C34/nm1;EvosG<-rep(0,n);ix<-n%%2
												if(ix==0){n2<-n/2} else {n2<-(n+1)/2}
  												for (i in 1:n2){
  													EvosG[i]<-qnorm(Medun[i],mean=0,sd=1,lower.tail=TRUE,log.p=FALSE)
  												    EvosG[n-i+1]=-EvosG[i]}
												k2a<-sum(EvosG^2);k2b<-sum(EvosG^4);k2<-k2a/(n-1)
												k4<-n*( (n+1)*k2b - (3*(n-1)/n)*k2a^2)/(n-1)/(n-2)/(n-3)
												k4a<-3*(n-1)/(n+1)+((n-2)*(n-3)/(n*(n^2-1)))*(k4/k2^2)^2
												k4a<-k4a/(n-1)/(n-1);mu4n<-mpfr(k4a,320)}											
			 		a<-Timc(n,asNumeric(mu2n),asNumeric(mu4n),icoef)
			 		Lam<-mpfr(a[[1]],320)
		 			Pv<-integrateR(Vian,L1=Lam[1],L2=Lam[2],-1,-abs(rc),rel.tol=1e-09)$value
		 			if (type=="two-sided") {Pv<-C39*Pv}
		 			if((r<=0  & type=="greater") | (r>0 & type=="less")) {Pv<-1-Pv}
		 			Pv<-formatMpfr(Pv,digits=cifer);r<-trunc(r*100000)/100000
		 			c1<-trunc(asNumeric(Pv)*1000000000)/1000000000
		 			La1<-formatMpfr(Lam[1],digits=cifer);La2<-formatMpfr(Lam[2],digits=cifer)
		 			Lz<-as.numeric(c(noquote(La1),noquote(La2)))
		 			Lz<-trunc(as.numeric(Lz)*1000000000)/1000000000}
#
				if (approx=="vggfr"){
						if(print){
							cat("Statistic:",estat,"-- Approximation: ",eappr," Alternative :",type," Contin. Adjust. :",CC,"\n")
							cat("n=",n,"Observed value:",noquote(sprintf("%.5f",r))," Est. p-value=",
							noquote(sprintf("%.6f",c1)),"\n")}
					   outrank<-list(n=n,Statistic=estat,Value=r,Approx=eappr,Tails=type, Cpv=c1,Lpv=c1, Lambda=Lz)
					   return(outrank)}	
#
	if (approx=="exact") {
		eappr<-"Exact_level"
		if  ("kendall" == index) {estat<-"Kendall's tau"
		 				Cs1<-2/(n^2-n)
		 				fname<-paste("Kend",n,".txt.zip",sep="")
						gname<-paste("Kend",n,".txt",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))			 				
						B[,1]<-mpfr(noquote(XX[,1]),320);B[,2]<-mpfr(noquote(XX[,2]),320)
						B[,1]<-Cs1*B[,1];B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}
		 if  ("gini" == index) {estat<-"Gini's cograduation coefficient"
		 				kn<-n%%2;Cs1<-2/(n^2-kn)
		 				fname<-paste("Gini",n,".txt.zip",sep="")
						gname<-paste("Gini",n,".txt",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)	 				
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))	
		        		B[,1]<-Cs1*B[,1];B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}
	  	if ("spearman" == index & n<=26) {estat<-"Spearman's rho"
		 				Cs1<-(C40/(nu^3-nu));Cs2<-nu*(nu+C34)*(C39*nu+C34)/C38
		 				fname<-paste("Spear",n,".txt.zip",sep="")
						gname<-paste("Spear",n,".txt",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))	
						B[,1]<-C34-Cs1*(Cs2-C39*B[,1])
						B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}	       			
 		 if ("r4" == index) {estat<-"r4"
		 				fname<-paste("R4",n,".csv.zip",sep="")
						gname<-paste("R4",n,".csv",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
						a<-noquote(XX[,1]);b<-noquote(XX[,2])
						B<-cbind(mpfr(a,320),mpfr(b,320))
						B[,1]<-B[,1]/10000;ws1<-sum(as.numeric(B[,2]))
		    			B[,2]<-B[,2]/ws1;B[,2]<-cumsum(B[,2])}
		 if ("fy" == index) {estat<-"Fisher-Yates"
		 				fname<-paste("Rg",n,".txt.zip",sep="")
		 				gname<-paste("Rg",n,".txt",sep="")
		 				XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
					    XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
						a<-noquote(XX[,1]);b<-noquote(XX[,2])
						B<-cbind(mpfr(a,320),mpfr(b,320))
					    B[,1]<-B[,1]/10000;ws1<-sum(as.numeric(B[,2]))
					    B[,2]<-B[,2]/ws1;B[,2]<-cumsum(B[,2])}						
 		 if ("filliben" == index) {estat<-"Filliben's rank correlation"
		 				fname<-paste("Fil",n,".txt.zip",sep="")
						gname<-paste("Fil",n,".txt",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))
						B[,1]<-B[,1]/10000;ws1<-sum(as.numeric(B[,2]))
						B[,2]<-B[,2]/ws1;B[,2]<-cumsum(B[,2])}							
# 			
		 if (type=="two-sided") {
		 		j1<-which(B[,1]<=rc)[length(which(B[,1]<=rc))];j2<-which(B[,1]>=rc)[1]
		    	ws1<-C39*B[j1,2];ws2<-C39*B[j2,2]
		    	if (ws2>1) {ws2<-1}
		    	Cpv<-ws1;Lpv<-ws2
		    	Cpv<-formatMpfr(Cpv,digits=cifer);Lpv<-formatMpfr(Lpv,digits=cifer)
		    	}
		if (!(type=="two-sided")){
				j1<-which(B[,1]<=rc)[length(which(B[,1]<=rc))]
				j2<-which(B[,1]>rc)[1]
		 		if (r>=0) {Cpv<-formatMpfr(1-B[j2,2],digits=cifer);Lpv<-formatMpfr(1-B[j1,2],digits=cifer)
		 					   ws1<-formatMpfr(B[j1,2],digits=cifer);ws2<-formatMpfr(B[j2,2],digits=cifer)} 
		 		if (r<0)   {Cpv<-formatMpfr(B[j1,2],digits=cifer);Lpv<-formatMpfr(B[j2,2],digits=cifer)
		 					   ws1<-formatMpfr(1-B[j2,2],digits=cifer);ws2<-formatMpfr(1-B[j1,2],digits=cifer)}
		 		if ((r>0 & type=="greater") | (r<0 & type=="greater") ) {Cpv<-ws1;Lpv<-ws2}
		 		}
		 }
		 r<-trunc(r*100000)/100000
 	 	 c1<-trunc(asNumeric(Cpv)*1000000000)/1000000000
 	 	 c2<-trunc(asNumeric(Lpv)*1000000000)/1000000000
 	 	 if (print){
 	 				cat("Statistic:",estat,"value: ",noquote(sprintf("%.5f",r)),"\n","Approximation: ",eappr,"  Alternative :",
 	 				type, "Corr. Cont.",CC,"\n")
			        cat(" Conservative p-value:",noquote(sprintf("%.6f",c1))," Liberal p-value:",
			        noquote(sprintf("%.6f",c2)),"\n")}
  		outrank<-list(n=n,Statistic=estat,Value=r,approx=eappr,tails=type,Cpv=c1,Lpv=c2)
return(outrank)}