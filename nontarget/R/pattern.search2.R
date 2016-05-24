pattern.search2<-function(
	peaklist,
	quantiz,
	mztol=2, 
	ppm=TRUE,
	inttol=0.5,
	rttol=0.3,
	use_isotopes=c("13C","37Cl","15N","81Br","34S","18O"),
	use_charges=c(1,2),
	use_marker=TRUE,
	quick=FALSE,
	isotopes
){
	
	########################################################################################
	# read quantiz & checks ################################################################
	size_deltamass<-quantiz[[5]][[1]]
	size_mass<-quantiz[[5]][[2]]
	size_intens<-quantiz[[5]][[3]]
	adductmass_LB<-quantiz[[5]][[7]]
	adductmass_UB<-quantiz[[5]][[8]]	
	max_d_delmz<-quantiz[[5]][[9]]
	max_d_ratio<-quantiz[[5]][[10]]
	max_d_mass<-quantiz[[5]][[11]]
	isotope_key<-quantiz[[3]]
	charge_key<-quantiz[[4]]
	size_deltamass<-(size_deltamass+max_d_delmz)
	size_mass<-(size_mass+max_d_mass)
	size_intens<-(size_intens+max_d_ratio)
	mass_slots<-quantiz[[7]]		
	cat("\n(1) Check inputs ...");
	if(mztol<0){warning("mztol should be >=0!")};
    if(inttol>1 || inttol<0){ stop("inttol must be >0 and <=1") };
	if(!is.data.frame(peaklist)){stop("peaklist must be a data.frame")}
	if(length(peaklist[1,])>3){stop("peaklist with > 3 columns not allowed")}
	if(!length(peaklist[,1])>1){stop("peaklist with one entry - doesn`t make sense ...")}
	if(!is.numeric(peaklist[,1]) || !is.numeric(peaklist[,2]) || !is.numeric(peaklist[,3]) ){stop("peaklist columns not numeric")}
	if(ppm=="TRUE"){ppm2<-1}else{ppm2<-0}
	if(use_isotopes[1]!="FALSE"){
		if(any(is.na(match(use_isotopes,isotope_key)))){
			paste("invalid use_isotopes, available only: ",sep="");print(isotope_key);stop();
		}
	}
	if(use_charges[1]!="FALSE"){
		use_charges2<-abs(use_charges)
		if(any(is.na(match(use_charges2,charge_key)))){
			paste("invalid use_charges, available only: ",sep="");print(charge_key);stop();
		}
	}
	cat(" done.");
	########################################################################################
	# prescreen for relevant peak-to-peak mass differences: mass + int slots ############### 
	cat("\n(2) Build peaklist kd-tree, screen slots, query quantized data: \n");
	pBar <- txtProgressBar( min = 0, max = length(peaklist[,1]), style = 3 )
	inter<-as.numeric(interactive())
	peakTree<-.Call("kdtree4", 
			as.matrix(peaklist[,1:3]),
			as.integer(inter),
			pBar,
			PACKAGE="nontarget"
	);
	close(pBar);
	peakTree<-peakTree[,1:4,drop=FALSE];
	cat("\n screen ... ");
	mass_slots<-quantiz[[7]]
	int_slots<-(10^quantiz[[8]])
	pBar <- txtProgressBar( min = 0, max = length(peaklist[,1]), style = 3 )
	relat<-.Call("peak_search", 
			as.matrix(peaklist[,1:3]),
			as.matrix(peakTree),		# peaks - search tree
			as.matrix(mass_slots),		# prefilter on mass
			as.matrix(int_slots),		# prefilter on intensity ratios
			as.numeric(mztol), 			# precision measurement mass
			as.numeric(ppm2),			# precision measurement - mass in ppm?
			as.numeric(inttol),			# precision measurement, %percent, NOT fraction
			as.numeric(rttol),			# precision measurement RT
			as.integer(inter),
			pBar,
			PACKAGE="nontarget"
	)	
	close(pBar)
	if(length(relat)<1){stop("\n No matches found \n ")}
	########################################################################################
	# find matches in quantized data #######################################################
	done<-matrix(ncol=length(charge_key),nrow=length(isotope_key),FALSE) # w/ marker necesssary?
	colnames(done)<-charge_key
	rownames(done)<-isotope_key
	search_bounds<-rep(0,6)
	marker_bounds<-matrix(ncol=2,nrow=3,0)
	marker_bounds[2,1]<-min(peaklist[,2])
	marker_bounds[2,2]<-max(peaklist[,2])
	bound_int<-log10((1+inttol)/(1-inttol)) # epsilon_4
	from_peak<-c()
	to_peak<-c()
	isotope<-c()
	charge<-c()
	retr_1<-0; # number of queries 	
	pBar <- txtProgressBar( min = 0, max = length(relat[,1]), style = 3 )
	for(j in 1:length(relat[,1])){
		done[,]<-FALSE;
		got<-FALSE;
		# del_mz bounds: extend by mass precision
		del_mass<-(peaklist[relat[j,2],1]-peaklist[relat[j,1],1])
		if(ppm==TRUE){
			search_bounds[1]<-(del_mass-(2*mztol*peaklist[relat[j,1],1]/1E6))
			search_bounds[2]<-(del_mass+(2*mztol*peaklist[relat[j,1],1]/1E6))
		}else{
			search_bounds[1]<-(del_mass-(2*mztol))
			search_bounds[2]<-(del_mass+(2*mztol))
		}
		# mass bounds: extend by neglected adduct masses
		search_bounds[3]<-(peaklist[relat[j,1],1]+adductmass_LB)
		search_bounds[4]<-(peaklist[relat[j,1],1]+adductmass_UB)
		# intensity bounds: extend by intensity tolerance
		log_int<-log10(peaklist[relat[j,1],2]/peaklist[relat[j,2],2])
		search_bounds[5]<-(log_int-bound_int)	# Lower bound
		search_bounds[6]<-(log_int+bound_int) 	# Upper bound
		# iterate over all quantizations	
		for(i in 1:length(quantiz[[6]])){
			# use this isotope & charge ? #########################################
			do<-TRUE;
			if(use_isotopes[1]!=FALSE){
				if(!any(use_isotopes==isotope_key[as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][1])])){
					do<-FALSE;
				}
			}
			if(use_charges[1]!=FALSE){
				if(!any(use_charges2==as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][2]))){
					do<-FALSE;
				}
			}
			if(do==FALSE){next}
			# without marker peak #################################################
			if(strsplit(names(quantiz[[6]])[i],"_")[[1]][3]=="wo"){
				found <- .Call("search_boxtree", 
						quantiz[[6]][[i]][,1:6],
						quantiz[[6]][[i]][,16:20],
						as.numeric(search_bounds),
						as.integer(0),
						PACKAGE="nontarget"
				)
				retr_1<-c(retr_1+1)					
				if(found==-2){
					done[as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][1]),as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][2])]<-TRUE;
					from_peak<-c(from_peak,relat[j,1])
					to_peak<-c(to_peak,relat[j,2])
					isotope<-c(isotope,as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][1]))
					charge<-c(charge,as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][2]))
					got<-TRUE;
				}				
			}
			if(got & quick) break;
			if(got) next;
			# with marker peak ####################################################
			if(strsplit(names(quantiz[[6]])[i],"_")[[1]][3]=="w"){
				if(done[as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][1]),as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][2])]==FALSE){  # if not hit w/o marker - just check
					if(use_marker!="TRUE"){ # just check for intersection
						found <- .Call("search_boxtree", 
								quantiz[[6]][[i]][,1:6],
								quantiz[[6]][[i]][,16:20],
								as.numeric(search_bounds),
								as.integer(0),
								PACKAGE="nontarget"
						)	
						if(found==-2){
							done[as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][1]),as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][2])]<-TRUE;
							from_peak<-c(from_peak,relat[j,1])
							to_peak<-c(to_peak,relat[j,2])
							isotope<-c(isotope,as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][1]))
							charge<-c(charge,as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][2]))
							got<-TRUE;
						}
						retr_1<-c(retr_1+1)							
					}else{ # check explicitly for marker peak 
						found <- .Call("search_boxtree", 
								quantiz[[6]][[i]][,1:6],
								quantiz[[6]][[i]][,16:20],
								as.numeric(search_bounds),
								as.integer(1), # return full findings
								PACKAGE="nontarget"
						)	
						retr_1<-c(retr_1+1)	
						if(length(found)>0){
							for( k in 1:length(found) ){
								marker_delmass<-c((peaklist[relat[j,2],1]-quantiz[[6]][[i]][found[k],8]),(peaklist[relat[j,2],1]-quantiz[[6]][[i]][found[k],7]))
								# marker m/z bounds
								if(ppm==TRUE){
									marker_bounds[1,1]<-(min(marker_delmass)-(2*mztol*peaklist[relat[j,1],1]/1E6))
									marker_bounds[1,2]<-(max(marker_delmass)+(2*mztol*peaklist[relat[j,1],1]/1E6))
								}else{
									marker_bounds[1,1]<-(min(marker_delmass)-(2*mztol))
									marker_bounds[1,2]<-(max(marker_delmass)+(2*mztol))
								}
								# marker intensity bounds. set above: marker_bounds[4]; reset lower bound
								max_int<-max(peaklist[relat[j,],2])
								marker_bounds[2,1]<-(max_int*(1-inttol))
								# marker RT bounds
								marker_bounds[3,1]<-(peaklist[relat[j,1],3]-rttol)
								marker_bounds[3,2]<-(peaklist[relat[j,1],3]+rttol)
								found_m<-.Call("search_kdtree", 
										as.matrix(peaklist[,1:3]),
										as.matrix(peakTree),		# peaks - search tree
										as.matrix(marker_bounds),
										PACKAGE="nontarget"
								)	
								if(length(found_m)>1){
									done[as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][1]),as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][2])]<-TRUE;
									from_peak<-c(from_peak,relat[j,1])
									to_peak<-c(to_peak,relat[j,2])
									isotope<-c(isotope,as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][1]))
									charge<-c(charge,as.numeric(strsplit(names(quantiz[[6]])[i],"_")[[1]][2]))
									got<-TRUE
									break; # on this isotope & charge
								}
							}
						}
					}
				}
			}
			#######################################################################		
			if(inter==1){
				setTxtProgressBar(pBar,j,title = NULL, label = NULL)			
			}
			if(got & quick) break;
		}				
	}
	close(pBar)
	if(length(from_peak)==0){stop("\n No matches found \n")}
	########################################################################################
	# for grouping, include reverse relation - omit the reverse thereafter... ##############
	# ... as opposed to adduct.search(), where relations go for- AND backward ##############
	# sort relations by increasing ID of first entry #######################################
	use<-c(rep(1,length(to_peak)),rep(2,length(to_peak)))
	isotope<-c(isotope,isotope)
	charge<-c(charge,charge)
	to_peak2<-c(to_peak,from_peak)
	from_peak2<-c(from_peak,to_peak)
	use<-use[order(from_peak2,decreasing=FALSE)]
	isotope<-isotope[order(from_peak2,decreasing=FALSE)]
	charge<-charge[order(from_peak2,decreasing=FALSE)]
	to_peak2<-to_peak2[order(from_peak2,decreasing=FALSE)]
	from_peak2<-from_peak2[order(from_peak2,decreasing=FALSE)]
	groups<-.Call("metagroup",
		as.integer(from_peak2),
		as.integer(to_peak2),
		PACKAGE="nontarget" 
	)
	to_peak<-to_peak2[use==1]
	from_peak<-from_peak2[use==1]
	charge<-charge[use==1]
	isotope<-isotope[use==1]
	groups<-groups[use==1]
	cat(paste("\n  ",length(to_peak)," of ",length(relat[,1])," candidate linkages accepted.","\n",sep=""));
	########################################################################################
	# generate output ######################################################################
	cat("(3) Create output ...");
	pattern<-list(0)
	# (1) Dataframe of peaks and their isotope relations within isotope pattern groups #####
    alls<-length(peaklist[,1]);
    ID<-seq(1:alls);			# (4) peak ID
    getit1<-rep("0",alls);   	# (5) to which peak?
    getit2<-rep("none",alls);   # (6) which isotope?
    getit3<-rep("0",alls);      # (7) mass tolerance
    getit4<-rep("0",alls);      # (8) charge level
	getit5<-rep("0",alls);		# (9) group ID
	getit6<-rep("0",alls);		# (10) interaction level
	for(i in 1:length(from_peak)){
		getit1[from_peak[i]]<-paste(getit1[from_peak[i]],as.character(to_peak[i]),sep="/")
		getit2[from_peak[i]]<-paste(getit2[from_peak[i]],quantiz[[3]][isotope[i]],sep="/")
		getit3[from_peak[i]]<-paste(getit3[from_peak[i]],"large",sep="/")
		getit4[from_peak[i]]<-paste(getit4[from_peak[i]],as.character(charge[i]),sep="/")
		getit5[from_peak[i]]<-paste(getit5[from_peak[i]],as.character(groups[i]),sep="/")
		getit5[to_peak[i]]<-paste(getit5[to_peak[i]],as.character(groups[i]),sep="/")		
		#getit6[relat[i,1]]<-paste(getit6[relat[i,1]],"-",sep="/")
	}
	for(i in 1:alls){	  
		if(getit1[i]!="0"){getit1[i]<-substr(getit1[i],3,nchar(getit1[i]))};
		if(getit2[i]!="none"){getit2[i]<-substr(getit2[i],6,nchar(getit2[i]))};
		if(getit3[i]!="0"){getit3[i]<-substr(getit3[i],3,nchar(getit3[i]))};
		if(getit4[i]!="0"){getit4[i]<-substr(getit4[i],3,nchar(getit4[i]))};	  
		if(getit5[i]!="0"){
			getit5[i]<-substr(getit5[i],3,nchar(getit5[i]))
			those<-strsplit(getit5[i],"/",fixed = TRUE)[[1]]
			those<-unique(those);
			getit5[i]<-those[1];
			if(length(those)>1){
				for(j in 2:length(those)){
					getit5[i]<-paste(getit5[i],those[j],sep="/");
				}
			}
		};
		if(getit6[i]!="0"){getit6[i]<-substr(getit6[i],3,nchar(getit6[i]))};
	}
	grouped_peaks<-data.frame(peaklist,ID,getit5,getit6,getit1,getit2,getit3,getit4,stringsAsFactors=FALSE)
    names(grouped_peaks)<-c(names(peaklist),"peak ID","group ID","interaction level","to ID",
		"isotope(s)","mass tolerance","charge level")
	pattern[[1]]<-grouped_peaks
	# (2) Parameters #######################################################################
	parameters<-data.frame(-rttol,rttol,mztol,0,ppm,inttol,0,0,adductmass_LB,adductmass_UB,
		size_deltamass,size_mass,size_intens,stringsAsFactors=FALSE)
	names(parameters)<-c("rttol","rttol","mztol","mzfrac","ppm","inttol","cutint",
		"deter","adductmass_LB","adductmass_UB","size_deltamass","size_mass","size_intens")
	pattern[[2]]<-parameters
	# (3) Peaks in pattern groups ##########################################################
	groupID<-unique(groups);
	groupID<-groupID[order(groupID)];
	peakIDs<-c();
	charge_group<-c();
	charge_count<-unique(charge)
	counted<-rep(0,length(charge_count))
	for(i in 1:length(groupID)){
		those<-c(from_peak[groups==groupID[i]],to_peak[groups==groupID[i]])
		those<-unique(those)
		get1<-as.character(those[1])
		for(j in 2:length(those)){
			get1<-paste(get1,",",those[j],sep="")
		}
		peakIDs<-c(peakIDs,get1);
		charges<-charge[groups==groupID[i]]
		charges<-unique(charges)
		counted[charge_count==charges[1]]<-(counted[charge_count==charges[1]]+1)
		charge_level<-as.character(charges[1])
		if(length(charges)>1){
			for(j in 2:length(charges)){	
				counted[charge_count==charges[1]]<-(counted[charge_count==charges[j]]+1)
				charge_level<-paste(charges[j],"/",charge_level,sep="")
			}
		}
		charge_group<-c(charge_group,charge_level);
	}
	for(i in 1:length(groupID)){
		groupID[i]<-paste("/",as.character(groupID[i]),"/",sep="")
	}
	grouping<-data.frame(groupID,peakIDs,charge_group,stringsAsFactors=FALSE)
	names(grouping)<-c("group ID","peak IDs","charge")
	pattern[[3]]<-grouping
	########################################################################################
	# (4) Atom counts ######################################################################
	pattern[[4]]<-"no information"
	# (5) Count of pattern groups ##########################################################
	charge_count<-data.frame(charge_count,counted,stringsAsFactors=FALSE)
	names(charge_count)<-c("Charge level","Counts")
	pattern[[5]]<-charge_count
	# (6) Removal by rules #################################################################
	# (7) Number of peaks with pattern group overlapping ###################################
	pattern[[6]]<-numeric(0)
	pattern[[7]]<-"no information"
	# (8) Number of peaks per within-group interaction levels ##############################
	pattern[[8]]<-"no information"
	# (9) Counts of isotopes ###############################################################
	if(quick){
		pattern[[9]]<-"no information using quick"
	}else{
		isos<-rep(isotope_key,length(charge_key))
		chrgs<-rep(charge_key,length(isotope_key))
		chrgs<-chrgs[order(chrgs)]
		incr_count<-rep(0,length(chrgs))
		group_count<-rep(0,length(chrgs))
		element<-c()
		for(i in 1:length(isos)){
			element<-c(element,as.character(isotopes[isotopes[,2]==isos[i],1][1]))
		}
		# on mass increment counts
		for(i in 1:length(getit2)){
			if(getit2[i]!="none"){
				get1<-strsplit(getit2[i],"/")[[1]]
				get2<-strsplit(getit4[i],"/")[[1]]
				for(j in 1:length(get1)){
					incr_count[isos==get1[j] & chrgs==get2[j]]<-(incr_count[isos==get1[j] & chrgs==get2[j]]+1)
				}
			}
		}
		# on group counts
		for(i in 1:length(peakIDs)){
			peaks<-as.numeric(strsplit(peakIDs[i],",")[[1]])
			get3<-c()
			get4<-c()
			for(j in 1:length(peaks)){
				if(getit2[peaks[j]]!="none"){
					get3<-c(get3,strsplit(getit2[peaks[j]],"/")[[1]])
					get4<-c(get4,as.numeric(strsplit(getit4[peaks[j]],"/")[[1]]))
				}
			}
			get1<-unique(data.frame(get3,get4))
			for(j in 1:length(get1[,1])){
				group_count[isos==get1[j,1] & chrgs==get1[j,2]]<-(group_count[isos==get1[j,1] & chrgs==get1[j,2]]+1)
			}
		}
		counts<-data.frame(isos,chrgs,incr_count,group_count,element,stringsAsFactors=FALSE)
		names(counts)<-c("isotope","charge","peak counts","group counts","element")	
		counts<-counts[(counts[,3]!=0 | counts[,4]!=0),]
		pattern[[9]]<-counts
	}
	# (10) Elements ########################################################################
	if(quick){
		pattern[[10]]<-"no information"
	}else{
		pattern[[10]]<-as.character(unique(counts[,5]))
	}
	# (11) Charges #########################################################################
	pattern[[11]]<-use_charges
	# (12) Rule settings ###################################################################
	pattern[[12]]<-"no information"
	########################################################################################
	names(pattern)<-c("Patterns","Parameters","Peaks in pattern groups","Atom counts","Count of pattern groups",
	"Removals by rules","Number of peaks with pattern group overlapping",
	"Number of peaks per within-group interaction levels",
	"Counts of isotopes","Elements","Charges","Rule settings");
	cat(paste(" queries: ",retr_1," - done.\n",sep=""));
	return(pattern);
	########################################################################################

}






	
	
	
	
	
	
	
	
	
		
	 
 
 
