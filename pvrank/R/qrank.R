qrank<-function(prob, n, index="spearman", approx="vggfr", print=FALSE, lower.tail=TRUE){
	if (!is.numeric(prob)| !is.numeric(n)) {stop("Non-numeric argument to mathematical function")}
	if (!(prob>=0 & prob<=1)) {stop("The nominal significance level must satisfy 0<=prob<=1)");print(prob)}
	nn<-floor(n);cifer<-9
	if(n<0) {stop("A negative number of ranks is indicated")}
	if(!(nn==n)) {warning("A non integer number of ranks is indicated. A truncation is necessary.");n<-nn}
	if (n<5) {stop("The number of ranks must be at least 5");print(n)}
	Lowt<-"(lower tail)."
	if (!lower.tail) {Lowt<-"(upper tail)."}
    ain<-c("spearman","kendall","gini","r4","fy","filliben")
	index<-tolower(index)
	index<-match.arg(index, ain, several.ok = TRUE)
	apx<-c("gaussian","student","vggfr","exact");approx<-tolower(approx)
	approx<-match.arg(approx, apx, several.ok = TRUE)	
 	ksw<-FALSE
	if ((index=="spearman" & n>26) | (index=="kendall" & n>60) | (index=="gini" & n>24)) {ksw<-TRUE}
	if ((index=="r4" & n>15) | (index=="fy" & n>15) | (index=="filliben" & n>15)) {ksw<-TRUE}
	if (approx=="exact" & ksw){
		    cat("The exact p-value is not available. An approximate p-value is computed \n")
			if (index=="r4" | index=="spearman") {approx<-"vggfr"}
			if (index=="kendall" | index=="fy" | index=="filliben") {approx<-"gaussian"}
			if (index=="gini") {approx<-"student"}
			}
(opr <- mpfr_default_prec())
	stopifnot(opr == (oprec <- mpfr_default_prec(300)), 300  == mpfr_default_prec())
	mpfr_default_prec(opr)
	C1<-mpfr(25,320);C2<-mpfr(38,320);C3<-mpfr(35,320);C4<-mpfr(72,320)
	C5<-mpfr(5,320);C6<-mpfr(9,320);C7<-mpfr(100,320);C8<-mpfr(328,320);C9<-mpfr(127,320)
	C0<-mpfr(997,320);C11<-mpfr(372,320);C12<-mpfr(1350,320);C13<-mpfr(111,320)
	C14<-mpfr(153,320);C15<-mpfr(366,320);C16<-mpfr(301,320);C17<-mpfr(456,320)
	C18<-mpfr(912,320);C19<-mpfr(1248,320);C20<-mpfr(107,320);C21<-mpfr(4,320)
	C22<-mpfr(76,320);C23<-mpfr(182,320);C24<-mpfr(307,320);C25<-mpfr(315,320)
	C26<-mpfr(342,320);C27<-mpfr(420,320);C28<-mpfr(315,320);C29<-mpfr(105,320)
	C30<-mpfr(6,320);C32<-mpfr(1.00762,320);C33<-mpfr(2.01524,320);C34<-mpfr(1,320);C35<-mpfr(10,320)
	C36<-mpfr(3,320);C37<-mpfr(2,320);C38<-mpfr(8,320);C39<-mpfr(0.5,320);C40<-mpfr(1.5,320)	
#
	nu<-mpfr(n,320);nm1<-nu-C34;nm2<-nu-C37;nm3<-nu-C36
	kn<-n%%2;kkn<-mpfr(kn,320)
	if (approx=="vggfr"){eappr<-"VGGFR"
		    a<-ranktes(0.0,n,index,"vggfr",F,"less",F)
			L1<-a$Lambda[1];L2<-a$Lambda[2]
			dVian<-function(x,L1,L2){L1*(1-abs(x)^L1)^L2/(2*beta(1/L1,1+L2))}
			fVian<-function(x,L1,L2){integrate(dVian,L1=L1,L2=L2,-1,x)$value}
			gVian<-function(x,prob1,L1,L2){abs(prob1-fVian(x,L1,L2))}			
			if(prob>0.5){prob1<-1-prob} else {prob1<-prob}
			rc<-optimize(gVian,c(-1,0),tol=0.000000001,prob1,L1,L2)$minimum	
			if (prob>0.5){rc<- -rc}
			}
	 if (approx=="student"){eappr<-"t-Student"
	 		 if (prob<1e-323){prob<-1e-323}
			 if (index=="spearman"){Lamx<-nm2
			 		Tx<-qt(prob, asNumeric(Lamx));a<-Tx^2/Lamx;rc<-sqrt(a/(C34+a))}			 	 		  	 		  			 	 		  	
 			 if (index=="gini"){z1<-C36*nm1*(nu^2-kkn);z2<-C37*(nu^2+C37+kkn)
 			 		Lm<-C39*(z1/z2 -C34);Lamx<-C37*ceiling(Lm+C39)
		  	    	Tx<-qt(prob,asNumeric(Lamx))
		  	    	a<-Tx^2/Lamx;rc<-sqrt(a/(C34+a))}			  	    							   		  	    							   		  	    				
			 if (index=="kendall"){z1<-(C21*nu+C35)/(C6*nu*nm1);Lm<-C39*(C34/z1-C34)
			 	    Lamx<-C37*ceiling(Lm+C39);Tx<-qt(prob, asNumeric(Lamx))
			 		a<-Tx^2/Lamx;rc<-sqrt(a/(C34+a))}			 								 			 				
		     if (index=="r4") {Lamx<-ceiling((nu-C33)/C32)
		     		Tx<-qt(prob, asNumeric(Lamx))
		     		a<-Tx^2/Lamx;rc<-sqrt(a/(C34+a))}
		     if (index=="fy"){Lamx<-nm2
			 		Tx<-qt(prob, asNumeric(Lamx)); a<-Tx^2/Lamx;rc<-sqrt(a/(C34+a))}
			 if (index=="filliben"){Lamx<-nm2
			 		Tx<-qt(prob, asNumeric(Lamx)); a<-Tx^2/Lamx;rc<-sqrt(a/(C34+a))}
			if (prob>0.5){rc<- -rc}	
			 }	
	if (approx=="gaussian"){eappr<-"Gaussian"
		     if (prob<1e-323){prob<-1e-323}
			 zx<-qnorm(prob, mean = 0, sd = 1)
			 if (index=="spearman"){rc<-zx/sqrt(nm1)}
			 if (index=="gini"){rc<-zx/sqrt(C40*nu)}
			 if (index=="kendall"){rc<-zx/sqrt((C6*nu*nm1)/(C21*nu+C35))}
			 if (index=="r4"){rc<-sqrt(C32)*zx/sqrt(nm1)}
			 if (index=="fy"|index=="filliben"){rc<-zx/sqrt(nm1)}
			 if (prob>0.5){rc<- -rc}
			 }
	if (!(approx=="exact")){
		prob<-floor(prob*100000)/100000
		c1<-ceiling(asNumeric(rc)*10000000000)/10000000000
		if(rc<0) {c1<-max(c1,-1)}
		if(rc>0) {c1<-min(c1,1)}
		if (print){
					cat("Statistic: ",index," n:",n, " Nominal significance level:",noquote(sprintf("%.6f",asNumeric(prob))),Lowt,"\n")
					cat("Rank correlation approximate quantile:",noquote(sprintf("%.5f",asNumeric(c1))),"\n")
					if (!(approx=="vggfr")){cat("Approx.: ",eappr,"\n")}
					else {cat("Approx. : VGGFR with L=(",noquote(sprintf("%.6f",asNumeric(L1))),",",
						     noquote(sprintf("%.6f",asNumeric(L2))),")","\n")}
					}
		outrank<-list(n=n, Statistic=index, Level=prob, Cq=c1)
		return(outrank)}
##
	eappr<-"Exact"
	nu<-mpfr(n,320);nm1<-nu-1;nm2<-nu-2;nm3<-nu-3;np1=nu+1
	if  ("kendall" == index) {estat<-"Kendall's tau"
		 				Cs1<-C37/(nu^2-nu)
		 				fname<-paste("Kend",n,".txt.zip",sep="")
						gname<-paste("Kend",n,".txt",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol=2, byrow=TRUE)
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))			 				
						B[,1]<-mpfr(noquote(XX[,1]),320);B[,2]<-mpfr(noquote(XX[,2]),320)
						B[,1]<-Cs1*B[,1];B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}
	 if  ("gini" == index) {estat<-"Gini's cograduation coefficient"
		 				Cs1<-C37/(n^2-kkn)
		 				fname<-paste("Gini",n,".txt.zip",sep="")
						gname<-paste("Gini",n,".txt",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)	 				
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))	
		        		B[,1]<-Cs1*B[,1];B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}
	  if ("spearman" == index & n<=26) {estat<-"Spearman's rank correlation"
		 				Cs1<-(C30/(nu^3-nu));Cs2<-nu*(nu+1)*(C37*nu+C34)/C36
		 				fname<-paste("Spear",n,".txt.zip",sep="")
						gname<-paste("Spear",n,".txt",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))	
						B[,1]<-1-Cs1*(Cs2-2*B[,1])
						B[,2]<-cumsum(B[,2]);B[,2]<-B[,2]/factorialMpfr(n)}	       			
 		 if ("r4" == index) {estat<-"r4"
		 				fname<-paste("R4",n,".csv.zip",sep="")
						gname<-paste("R4",n,".csv",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))
						B[,1]<-B[,1]/100000;ws1<-sum(asNumeric(B[,2]))
		    			B[,2]<-B[,2]/ws1;B[,2]<-cumsum(B[,2])}
		 if ("fy" == index) {estat<-"Fisher-Yates coefficient"
		 				fname<-paste("Rg",n,".txt.zip",sep="")
		 				gname<-paste("Rg",n,".txt",sep="")
		 				XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
					    XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
					    B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))
					    B[,1]<-B[,1]/10000;ws1<-sum(asNumeric(B[,2]))
					    B[,2]<-B[,2]/ws1;B[,2]<-cumsum(B[,2])}						
 		 if ("filliben" == index) {estat<-"Filliben's rank correlation"
		 				fname<-paste("Fil",n,".txt.zip",sep="")
						gname<-paste("Fil",n,".txt",sep="")
						XX <- read.table(unz(system.file(package="pvrank","extdata", fname), gname), header=F)
						XX<-as.matrix(XX,ncol = 2, byrow = TRUE)
						B<-cbind(mpfr(XX[,1],320),mpfr(XX[,2],320))
						B[,1]<-B[,1]/10000;ws1<-sum(asNumeric(B[,2]))
						B[,2]<-B[,2]/ws1;B[,2]<-cumsum(B[,2])}					
#  	
		j1<-which(B[,2]<=prob)[length(which(B[,2]<=prob))]
		j2<-which(B[,2]>=prob)[1]	       					
		Cpv<-formatMpfr(B[j1,1],digits=cifer);Lpv<-formatMpfr(B[j2,1],digits=cifer)
		Cqv<-formatMpfr(B[j1,2],digits=cifer);Lqv<-formatMpfr(B[j2,2],digits=cifer)
#
	if(!lower.tail){Cqv<-1-asNumeric(Cqv);Lqv<-1-asNumeric(Lqv)}
	prob<-floor(prob*1000000)/1000000
	c1<-ceiling(asNumeric(Cpv)*1000000000)/1000000000
	c2<-ceiling(asNumeric(Cqv)*1000000000)/1000000000
	c3<-ceiling(asNumeric(Lpv)*1000000000)/1000000000
	c4<-ceiling(asNumeric(Lqv)*1000000000)/1000000000
	if (print){
		    cat("Statistic:",index,"n:",n, "Nominal significance level:",
		    noquote(sprintf("%.6f",asNumeric(prob))),Lowt,"Approx.:",eappr,"\n")
		   	cat("Rank correlation conservative quantile:",noquote(sprintf("%.4f",c1))," Real level:",noquote(sprintf("%.6f",c2)),"\n")
		    cat("Rank correlation liberal quantile          :",noquote(sprintf("%.4f",c3)), " Real level:",noquote(sprintf("%.6f",c4)),"\n")
		    }		    
	outrank<-list(n=n,Statistic=index,Level=prob,Cq=c1,Cv=c2,Lq=c3,Lv=c4)
	return(outrank)}