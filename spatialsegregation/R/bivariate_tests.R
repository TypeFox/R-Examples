# Tests for segregation

##### For two species patterns:

dixon<-function(X, prepR=0) {
	dbg<-FALSE
	if(!is.factor(X$marks))warning("Marks of X are not in factor form. Transforming.")
	X$marks<-as.factor(X$marks)
	splist<-levels(X$marks)
	if(length(splist)!=2) stop("Test defined only for 2-type pattern.")
	A<-splist[1]
	B<-splist[2]
	# start computing:
	Na<-sum(X$marks==A)
	Nb<-sum(X$marks==B)
	N<-Na+Nb
	#
	Paa<-Na*(Na-1)/(N*(N-1))
	Paaa<-Paa*(Na-2)/(N-2)
	Paaaa<-Paaa*(Na-3)/(N-3)
	Pbb<-Nb*(Nb-1)/(N*(N-1))
	Pbbb<-Pbb*(Nb-2)/(N-2)
	Pbbbb<-Pbbb*(Nb-3)/(N-3)
	Paabb<-Na*(Na-1)*Nb*(Nb-1)/(N*(N-1)*(N-2)*(N-3))
	#
	# that's the simple stuff, now the neighbours:
	pp<-sg.modify.pp(X)
	pp$area<-as.numeric(area.owin(X$window))
	
	# compute neighbours
	res<-.External("graph_c", 
			pp, #pp 
			as.integer(1), # knn 
			as.numeric(1), # k=1
			as.numeric(prepR), #
			as.integer(0),     # no toroidal
			as.integer(rep(1, N)), # inclusion vector
			as.integer(dbg),		
			PACKAGE="spatialsegregation")
	E<-unlist(res)
	
	# compute frequencies Naa, Nab, Nba, Nbb i.e. how many A->A, A->B, B->A, B->B
	a<-which(X$marks==A)
	b<-which(X$marks==B)
	Naa<-sum(E[a]%in%a)
	Nbb<-sum(E[b]%in%b)
	Nab<-sum(E[a]%in%b)
	Nba<-sum(E[b]%in%a)
	na<-Naa+Nba
	nb<-Nab+Nbb
	# expected frequencies
	ENaa<-Na*(Na-1)/(N-1)
	ENab<-Na*Nb/(N-1)
	ENba<-Nb*Na/(N-1)
	ENbb<-Nb*(Nb-1)/(N-1)
	# compute R = 2*(n. of pairs of reflexive nearest neighbours)
	#           = (n. of points i whose n-neighbour is i)
	R<-0
	R<-sum((E[E]==(1:N)))	
	# compute Q = 2*(N2+3N3+6N4+10N3+15N6), where
	# Ni = number of points of which are the nearest neighbour of i other points
    Ns<-tabulate(E)
	Ni<-NULL;for(i in 2:6) Ni[i-1]<-sum(Ns==i)
	Q<-2*sum(Ni*c(1,3,6,10,15))
	# Compute Var Naa
	VNaa<-(N+R)*Paa + (2*N-2*R+Q)*Paaa + (N^2 -3*N-Q+R)*Paaa-N^2*Paa^2
	VNbb<-(N+R)*Pbb + (2*N-2*R+Q)*Pbbb + (N^2 -3*N-Q+R)*Pbbb-N^2*Pbb^2
	CNaaNbb<-(N^2-3*N-Q+R)*Paabb-N^2*Paa*Pbb
	# That all the pieces we need. 
    
    #Then we compute the test statistics:
	# z scores
	zaa <- (Naa - ENaa)/sqrt(VNaa)
	zbb <- (Nbb - ENbb)/sqrt(VNbb)
	
	# combined C~Chi^2
	r<-CNaaNbb/sqrt(VNaa*VNbb)
	C<-(zaa^2 + zbb^2 -2*r*zaa*zbb)/(1-r^2)
	
	# measure of segregation
	Saa<-log((Naa/Nab)/((Na-1)/Nb) )
	Sbb<-log((Nbb/Nba)/((Nb-1)/Na) )
	
	# symmetry
	zs<-(Nba-Nab)/sqrt(VNaa+VNbb-2*CNaaNbb)
	
	
	# summary prints
	# table:
	ntable<-cbind(A=c(A=Naa,B=Nab),B=c(A=Nba,B=Nbb))
	etable<-cbind(A=c(A=ENaa,B=ENab),B=c(A=ENba,B=ENbb))
	ptable<-cbind(A=c(A=Na*na/N, B=Na*nb/N), B=c(A=Nb*na/N,B=Nb*nb/N))
	
	# return results
	fx<-list(observed=ntable,
		 dixontable=etable, 
	     pieloutable=ptable,
		 zscores=data.frame(zaa=zaa, zbb=zbb, C=C), 
	     seg.measures=data.frame(Saa,Sbb),
		 symmetry=data.frame(zs),
		 A=A, B=B
 )
 	class(fx)<-"segtest"
	fx
}

# print method
print.segtest<-function(x, ...)
{
	
	cat("Dixon's test of spatial segregation using nearest neighbour contingency table\n******\n")
	cat("Table of counts, using labels A=",x$A,", B=",x$B,":\n",sep="")
	cat("Pairs\tObs.\t Pielou expect.\t Random labeling expect.\n")
	cat("A->A\t",x$observed[1,1],"\t",x$pieloutable[1,1],"\t",x$dixontable[1,1],"\n")
	cat("A->B\t",x$observed[1,2],"\t",x$pieloutable[1,2],"\t",x$dixontable[1,2],"\n")
	cat("B->A\t",x$observed[2,1],"\t",x$pieloutable[2,1],"\t",x$dixontable[2,1],"\n")
	cat("B->B\t",x$observed[2,2],"\t",x$pieloutable[2,2],"\t",x$dixontable[2,2],"\n")
	
	cat("\nSegregation measures:\n S(A->A)=",x$seg.measures$Saa, "\n S(B->B)=",x$seg.measures$Sbb,"\n",sep="")
	
	cat("\nSegregation test:\n")
	cat("Combined\n  C=",x$zscores$C,"\t Chi2(2) p=",pchisq(x$zscores$C,df=2),"\n")
	cat("Individual\n  z(A->A)=",x$zscores$zaa,"\t N(0,1) p~=",pnorm(abs(x$zscores$zaa)),"\n")
	cat("  z(B->B)=",x$zscores$zbb,"\t N(0,1) p~=",pnorm(abs(x$zscores$zbb)),"\n------\n")
}


