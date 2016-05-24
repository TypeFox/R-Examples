make.isos <-
function(
	isotopes,
	use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K"),
	use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2)
	
){

	############################################################################
	if(length(use_isotopes)!=length(use_charges)){ stop("length(use_isotopes) not equal to length(use_charges) !") }
	if(any(is.na(match(use_isotopes,isotopes[,2])))){ stop("unexpected entry in use_isotopes = some isotope not found among isotopes argument")}
	############################################################################
	# on (1) isos & (5) elements ###############################################	
	elements<-c()
	for( i in 1:length(use_isotopes) ){
		elements<-c( elements,
			unique(as.character(isotopes[as.character(isotopes[,2])==use_isotopes[i],1][1]))
		)
	}
	#elements<-unique(elements)
    get_iso<-rep(FALSE,length(isotopes[,1]))
    for(i in 1:length(elements)){get_iso[isotopes[,1]==elements[i]]<-TRUE}
    isos<-isotopes[get_iso,]
    isos<-isos[isos[,4]!=0,]
	############################################################################
	# on (2) list of isotope masses ############################################
	use_charges2<-abs(use_charges)
    isomat<-data.frame(matrix(ncol=7,nrow=length(use_isotopes),0));
    isomat[,5]<-1;
    colnames(isomat)<-c("name","dmass","dabund","how_often","#atoms","/C","z");
	for(i in 1:length(use_isotopes)){
		if(length(isotopes[isotopes[,1]==elements[i],1])==1){
			stop(paste(elements[i]," has only one isotope ... abort",sep=""))
		}
		isosub<-isotopes[isotopes[,1]==elements[i],];
        isomat[i,1]<-use_isotopes[i];
        isomat[i,2]<-c(abs(isosub[isosub[,2]==use_isotopes[i],3]-min(isosub[,3])) / use_charges2[i]);
		if(isomat[i,2]==0){
			stop(paste(use_isotopes[i]," is the lightest (monoisotopic) isotope of ",elements[i]," ... abort",sep=""))
		}
        isomat[i,3]<-c((1/isosub[isosub[,3]==min(isosub[,3]),4])*isosub[isosub[,2]==use_isotopes[i],4]);
        isomat[i,6]<-c(isosub[isosub[,2]==use_isotopes[i],5]);	
		isomat[i,7]<-use_charges2[i];
	}
    isomat<-isomat[order(isomat[,2]),];
	############################################################################
    iso<-list(isos,isomat,unique(use_charges),length(use_isotopes),unique(elements));
    names(iso)<-c("list of isotopes","list of isotope masses","charges","number of isotope m/z","elements");
    return(iso)
	############################################################################
	
}
