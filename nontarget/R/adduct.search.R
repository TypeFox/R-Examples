adduct.search <-
function(
    peaklist,
	adducts,
    rttol=0,
	mztol=2,
	ppm=TRUE,
    use_adducts=c("M+H","M+K","M+Na"),
	ion_mode="positive"
){

    ############################################################################
    # (0) check inputs #########################################################
    if(ion_mode!="positive"&ion_mode!="negative"){stop("ion mode: positive or negative?")}
    use_adducts<-unique(use_adducts)
	if(!is.character(use_adducts)){stop("invalid use_adducts - not a vector of character strings")}
	if(length(use_adducts)<2){stop("at least two entries in use_adducts required")}
	for(i in 1:length(use_adducts)){if(any(adducts[,1]==use_adducts[i])!=TRUE&any(adducts[adducts[,1]==use_adducts[i],6]==ion_mode)){stop(paste("Adduct ",use_adducts[i]," not in adducts!",sep=""))}};
    for(i in 1:length(use_adducts)){if((adducts[adducts[,1]==use_adducts[i],6]!=ion_mode)){stop(paste(use_adducts[i]," not in ion mode ",ion_mode,sep=""))}};
	if(!is.data.frame(peaklist)){stop("peaklist must be a data.frame")}
	if(length(peaklist[1,])>3){stop("peaklist with > 3 columns not allowed")}
	if(!length(peaklist[,1])>1){stop("peaklist with one entry - doesn`t make sense ...")}
	if(!is.numeric(peaklist[,1])||!is.numeric(peaklist[,2])||!is.numeric(peaklist[,3])){stop("peaklist columns not numeric")}
    ############################################################################
    cat("\n(1) Combine adducts...");
    these<-match(use_adducts,adducts[,1]);
    these<-these[is.na(these)==FALSE];
    add<-adducts[these,];
    add<-add[add[,6]==ion_mode,];
    if(length(add[,1])<1){stop("No selected use_adducts among adducts ... abort!")};
    add2<-data.frame(0,0,0,0,0,0,0);
    add3_a<-c();
	add3_b<-c();
    names(add2)<-c("charge1","mult1","mass1","charge2","mult2","mass2","count")
    that<-c(2);
    for(i in 1:(length(add[,1])-1)){
		for(j in (i+1):length(add[,1])){
			add2<-rbind(add2,rep(0,7));
			add2[that,1]<-abs(add[i,3]);
			add2[that,2]<-add[i,4];
			add2[that,3]<-add[i,5];
			add2[that,4]<-abs(add[j,3]);
			add2[that,5]<-add[j,4];
			add2[that,6]<-add[j,5];
			that<-c(that+1);
			add3_a<-c(add3_a,paste(add[i,1],add[j,1],sep="<->"));
			add3_b<-c(add3_b,paste(add[j,1],add[i,1],sep="<->"));				
		}; #j
    }; #i
    add2<-add2[-1,];
    cat("done.");
    ############################################################################
	cat("\n(2) Build peaklist kd-tree, screen, ... \n");
	inter<-as.numeric(interactive())
	pBar <- txtProgressBar( min = 0, max = length(peaklist[,1]), style = 3 )
	peakTree<-.Call("kdtree4", 
		as.matrix(peaklist[,c(1,3)]),
		as.integer(inter),
		pBar,
		PACKAGE="nontarget"
	);
	close(pBar)
	peakTree<-peakTree[,1:4,drop=FALSE];
	if(ppm=="TRUE"){ppm2<-1}else{ppm2<-0}
	pBar <- txtProgressBar(min = 0, max = length(peaklist[,1]), style = 3 )
 	relat<-.Call("adduct_search",
		as.matrix(peaklist[,c(1,3)]),  	# peaklist
		as.matrix(peakTree),			# peaks - search tree
		as.matrix(add2),				# adduct table
		as.numeric(mztol), 				# precision measurement mass
		as.integer(ppm2),				# precision measurement - mass in ppm?
		as.numeric(rttol),				# precision measurement RT
		as.integer(inter),
		pBar,	
		PACKAGE="nontarget"
	);
	close(pBar)
	# check for self-referencing ###############################################	
	if(any(relat[,1]==relat[,2])){
		cat("remove self-references ...");
		relat<-relat[relat[,1]!=relat[,2],];
	};
	if(length(relat)==0){stop("\n No matches found \n ")}
	# form groups ##############################################################
	relat<-rbind(
		cbind(relat,rep(1,length(relat[,1]))),
		cbind(relat[,c(2,1,3),drop = FALSE],rep(2,length(relat[,1])))
	);
	relat<-relat[order(relat[,1],decreasing=FALSE),];	
	groups<-.Call("metagroup",
		as.integer(relat[,1]),
		as.integer(relat[,2]),
		PACKAGE="nontarget" 
	);
    cat("done.");
    ############################################################################
	cat("\n(3) Create output ... ");
    # (3.1) Peakwise ###########################################################
	getit1<-rep("none",length(peaklist[,1]));     # (1) which adduct removed <-> added?
    getit2<-rep("0",length(peaklist[,1]));        # (2) to which peak?
    getit3<-rep("0",length(peaklist[,1]));        # (3) within [1] large or [2] small mass tolerance?
	getit4<-rep("0",length(peaklist[,1]));        # (4) which group?
	for(i in 1:length(relat[,1])){
		if(relat[i,4]==1){
			getit1[relat[i,1]]<-paste(getit1[relat[i,1]],add3_a[relat[i,3]],sep="//")
		}else{
			getit1[relat[i,1]]<-paste(getit1[relat[i,1]],add3_b[relat[i,3]],sep="//")		
		}
		getit2[relat[i,1]]<-paste(getit2[relat[i,1]],as.character(relat[i,2]),sep="/")
		getit3[relat[i,1]]<-paste(getit3[relat[i,1]],"large",sep="/")
		getit4[relat[i,1]]<-paste(getit4[relat[i,1]],as.character(groups[i]),sep="/")
	}
	for(i in 1:length(getit1)){
		if(getit1[i]!="none"){
			getit1[i]<-substr(getit1[i],7,nchar(getit1[i]))
			getit2[i]<-substr(getit2[i],3,nchar(getit2[i]))
			getit3[i]<-substr(getit3[i],3,nchar(getit3[i]))
			getit4[i]<-strsplit(getit4[i],"/")[[1]][2]
		}
	}	
    ID<-seq(1:length(peaklist[,1]));	
    list_adducts<-data.frame(peaklist[,1:3],ID,getit4,getit2,getit1,getit3,stringsAsFactors=FALSE);
    names(list_adducts)<-c("m/z","int","ret","peak ID","group ID","to ID","adduct(s)","mass tolerance");
	# (3.2) Groupwise ##########################################################
    group1<-paste("/",as.character(1:max(groups)),"/",sep="") 	# groupnumber?
    group2<-rep("",max(groups)); 								# which peaks?
    group3<-rep("",max(groups)); 								# which use_adducts?
	group4<-rep("FALSE",max(groups)); 							# any(ambiguous adduct)?
	for(i in 1:length(getit1)){
		if(getit1[i]!="none"){
			this<-strsplit(getit1[i],"//")[[1]]
			this<-strsplit(this,"<->")
			that<-c() # adduct uniquely assigned?
			for(j in 1:length(this)){
				that<-c(that,this[[j]][1])
			}
			that<-unique(that)
			if(length(that)>1){
				group4[as.numeric(getit4[i])]<-"TRUE";
			}
			for(j in 1:length(that)){
				group2[as.numeric(getit4[i])]<-paste(group2[as.numeric(getit4[i])],as.character(i),sep=",")
				group3[as.numeric(getit4[i])]<-paste(group3[as.numeric(getit4[i])],that[j],sep="/")
			}
		}
	}
	group2<-substr(group2,2,nchar(group2))
	group3<-substr(group3,2,nchar(group3))
    grouping<-data.frame(group1,group2,group3,group4,stringsAsFactors=FALSE);
    names(grouping)<-c("group ID","peak IDs","adducts","ambig.?");
	overlaps<-sum(as.numeric(as.logical(group4)));
    # (3.3) counts of adduct matches ###########################################
    hits<-data.frame(use_adducts,rep(0,length(use_adducts)),stringsAsFactors=FALSE);
    names(hits)<-c("names","counts");
    for(i in 1:length(group1)){
		this1<-strsplit(as.character(group3[i]),"/")[[1]];
		for(j in 1:length(this1)){
			hits[hits[,1]==this1[j],2]<-hits[hits[,1]==this1[j],2]+1;
		};
    };
	# list #####################################################################
    parameters<-data.frame(rttol,mztol,ppm,ion_mode,stringsAsFactors=FALSE);	
    adduct<-list(list_adducts,parameters,grouping,hits,overlaps);
    names(adduct)<-c("adducts","Parameters","Peaks in adduct groups","Adduct counts","Overlaps");
    cat("done.\n\n");
    ############################################################################
    return(adduct);
    ############################################################################

}








