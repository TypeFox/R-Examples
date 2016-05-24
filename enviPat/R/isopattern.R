isopattern <-
function(
  isotopes,
  chemforms,
  threshold=.001,
  charge=FALSE,
  emass=0.00054858,
  plotit=FALSE,
  algo=1,
  rel_to=0,
  verbose=TRUE,
  return_iso_calc_amount=FALSE
  ){

    ############################################################################
    # (1) issue warnings #######################################################
    if(length(isotopes)!=5){stop("WARNING: invalid isotope list\n")}
    if(threshold>100 || threshold<0){stop("WARNING: invalid threshold; 0<=threshold<100.\n")}  
    if(plotit!="TRUE"&plotit!="FALSE"){stop("WARNING: plotit invalid. TRUE, FALSE.\n")}
    if(emass!=0.00054858){cat("NOTE: You are sure that is the mass of an electrone?")}
    if((length(charge)!=length(chemforms)) & length(charge)>1){stop("length of charge does not match number of chemforms!\n")}
    if(any(charge==0) & any(charge!=FALSE)){stop("WARNING: charge=0?")}
    if(any(is.numeric(charge)==FALSE) & any(charge!=FALSE)){stop("WARNING: charge either numeric or FALSE!")}
	if(length(charge)==1 & length(chemforms)>1){charge<-rep(charge,length(chemforms))}
    if(!any(algo==c(1,2))){stop("invalid algo argument!")}
    options(digits=10);
	if(!any(rel_to==c(0,1,2,3,4))){stop("invalid rel_to")}
	if(!is.logical(verbose)){stop("invalid verbose")}
    if(return_iso_calc_amount=="TRUE"){return_iso_calc_amount2=1}else{return_iso_calc_amount2=0}
	############################################################################
    # (2) run isotope pattern generator ######################################## 
    pattern<-list(0)
    for(i in 1:length(chemforms)){
      if(algo==1){ # = algo_1
        out <- .Call( "iso_pattern",
        s1 = as.character(chemforms[i]),   # chemical formula
        pl = as.integer(1E6),             # number of peaks to be reserved for
        t1 = as.double(threshold),        # relative intensity cutoff
        iso_list_elem = as.character(isotopes[,1]),  # isotope list: Element
        iso_list_iso = as.character(isotopes[,2]),  # isotope list: Isotope
        iso_list_mass = as.numeric(isotopes[,3]),  # isotope list: Isotope mass
        iso_list_abu = as.numeric(isotopes[,4]),  # isotope list: Isotope abundance
        rtm = as.integer(rel_to),     # 0:relative to highest, 1:relative to mono peak
        rica = as.integer(return_iso_calc_amount2),
        PACKAGE="enviPat"
		);
      }
      if(algo==2){ # = algo_4
        out <- .Call( "iso_pattern_4",
        s1 = as.character(chemforms[i]),   # chemical formula
        pl = as.integer(1E6),             # number of peaks to be reserved for
        t1 = as.double(threshold),        # relative intensity cutoff
        iso_list_elem = as.character(isotopes[,1]),  # isotope list: Element
        iso_list_iso = as.character(isotopes[,2]),  # isotope list: Isotope
        iso_list_mass = as.numeric(isotopes[,3]),  # isotope list: Isotope mass
        iso_list_abu = as.numeric(isotopes[,4]),  # isotope list: Isotope abundance
        rtm = as.integer(rel_to),     # 0:relative to highest, 1:relative to mono peak
        rica = as.integer(return_iso_calc_amount2),
        PACKAGE="enviPat"
		);
      }
      # parse output ###########################################################
      if(length(out[[1]])==0){
         pattern[[i]]<-"error";
      }else{
        if(return_iso_calc_amount2){
            pattern[[i]]<-out
        }else{
            out2<-out[order(out[,1],decreasing=FALSE),,drop=FALSE]
            colnames(out2)[1]<-"m/z"
            if(charge[i]!=FALSE){
                out2[,1]<-c(out2[,1]-(charge[i]*emass));  # electrone mass
                if(charge[i]!=1){
                    out2[,1]<-c(out2[,1]/abs(charge[i]));  # /charge=z
                }
            }
            pattern[[i]]<-out2
            if(plotit==TRUE){
                plot(out2[,1],out2[,2],type="h",
                xlab="m/z",ylab="Relative abundance",main=names(pattern)[i])
            }
        }
      }  
  }
  names(pattern)<-as.character(chemforms);
  if(verbose){cat(" done.");}
  ############################################################################ 
  # (3) output ###############################################################
  return(pattern)
  ############################################################################
}


