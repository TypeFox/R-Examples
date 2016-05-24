pattern.search <-
function(
	peaklist,
	iso,
	cutint=min(peaklist[,2]),
	rttol=c(-0.5,0.5),
	mztol=3,
	mzfrac=0.1,
    ppm=TRUE,
	inttol=0.5,
	rules=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
    deter=FALSE,
	entry=20
){

    ############################################################################
    # (0) check inputs #########################################################
    if(mzfrac>1 || mzfrac<=0){ stop("mzfrac must be >0 and <=1") };
    if(mztol<=0){warning("mztol should be >0!")};
    if(inttol>1 || inttol<0){ stop("inttol must be >0 and <=1") };
    if(length(rttol)!=2){stop("rttol must have a lower and an upper bound!")};
    if(rttol[1]>rttol[2]){stop("minimum > maximum for rttol!")};
    if(length(rules)<11){stop("wrong parameter setting: number of rules < 8!")}
	if(!is.data.frame(peaklist)){stop("peaklist must be a data.frame")}
	if(length(peaklist[1,])>3){stop("peaklist with > 3 columns not allowed")}
	if(!length(peaklist[,1])>1){stop("peaklist with one entry - doesn`t make sense ...")}
	if(!is.numeric(peaklist[,1]) || !is.numeric(peaklist[,2]) || !is.numeric(peaklist[,3]) ){stop("peaklist columns not numeric")}
    if(rules[4]==TRUE & any(iso$elements=="C")==FALSE & deter!=TRUE){stop("How is rule #7 supposed to work if carbon is not part of the iso argument? Include carbon or set rules[7] to FALSE.")}
	############################################################################
    cat("\n (1) Assemble lists ... ");
    # (1) define parameters / lists / matrices / ... ###########################
    # (1.1) sort peaklist
    getback<-order(peaklist[,3],peaklist[,1],decreasing=FALSE);
    samples<-peaklist[getback,];
    # retrieve data from isos
    isomat<-iso[[2]];
    isos<-iso[[1]];
    charges<-abs(iso[[3]]);
    manyisos<-iso[[4]];
    elements<-iso[[5]];
    if(deter==FALSE){
		for(i in 1:length(elements)){
			if(any(iso[[1]][,1]==elements[i])!=TRUE){
				stop(paste("Element ",elements[i]," not found in iso[[5]]!",sep=""))
			};
		};
    }else{
		rules=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE);
    };
    # (1.2) objects required for the screening step ...
    # (1.2.1) ... 1: find mass increments
    alls<-length(samples[,1]);
	ID<-getback; # ID<-seq(1:alls);
    getit1<-rep("none",alls);   # (1) which isotope?
    getit2<-rep("0",alls);      # (2) from which peak?
    getit4<-rep("0",alls);      # (3) to which peak?
    getit5<-rep("0",alls);      # (4) within [1] large or [2] small mass tolerance?
    getit6<-rep("0",alls);      # (5) with which charge?
    # (1.2.2) ... 2: validate increments
    mpoldnew<-max(isomat[,3]);
    countrem1<-c(0);
    countrem2<-c(0);
    countrem3<-c(0);
    countrem4<-c(0);
    countrem5<-c(0);
    countrem6<-c(0);
    countrem7<-c(0);
    countrem8<-c(0);
    countrem9<-c(0);
    countrem11<-c(0);    
    # (1.2.3) ...3: grouping & estimating atom numbers per element!
    groupcount<-c(1);
    group1<-rep("0",alls);      # which group? per charge level!
    group2<-rep("0",alls);      # which tree level?
    groupinfo<-rep();           # how many atoms?
    group3<-c();                # number of group ...
    group4<-c();                # ... and ID of peaks in that group!
    group5<-rep("0",alls);      # store charge level
    group6<-c();                # store charge level
    #
    cat("done.");
    ############################################################################

    ############################################################################
    cat("\n (2) Screen for mass increments ... ");
    # (2) screen for isotope dmass #############################################
    #dyn.load(paste(.libPaths(),"/nontarget/temp/massCpp.dll",sep=""));
    if(ppm==TRUE){ppm2=1}else{ppm2=2};
    getit1a<-rep(0,length(samples[,1])*entry);
    getit2a<-rep(0,length(samples[,1])*entry);
    getit4a<-rep(0,length(samples[,1])*entry);
    getit5a<-rep(0,length(samples[,1])*entry);
    getit6a<-rep(0,length(samples[,1])*entry);
    maxmass=max(isomat[,2]);
    len3<-length(samples[,1])
    result<-.C("mass",
		as.double(samples[,1]),
		as.double(samples[,3]),
		as.integer(len3),
		as.double(mztol*2),
		as.double(mzfrac*2),
		as.double(rttol[1]),
		as.double(rttol[2]),
		as.integer(manyisos), # 8
		as.double(isomat[,2]),
		as.integer(isomat[,4]),
		as.double(maxmass),
		as.integer(isomat[,7]), # 12
		as.integer(entry),
		as.integer(ppm2), # 14
		as.integer(getit1a),
		as.integer(getit2a),
		as.integer(getit4a),
		as.integer(getit5a),
		as.integer(getit6a), # 19
		PACKAGE="nontarget"
    )
    # generate outputs: ########################################################
    isomat[,4]<-result[10];
    # (1) which isotope? #######################################################
    for(i in 1:(alls)){
		if(any(result[15][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0)){
			getit1[i]<-paste(getit1[i],"/",paste0(result[15][[1]][((i-1)*entry+1):((i-1)*entry+entry)][result[15][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0],collapse="/"),sep="")
		}
	}
    # (2) from which peak? #####################################################
    for(i in 1:(alls)){
		if(any(result[16][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0)){
			getit2[i]<-paste(getit2[i],"/",paste0(result[16][[1]][((i-1)*entry+1):((i-1)*entry+entry)][result[16][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0],collapse="/"),sep="")
		}
	}
    # (3) to which peak? #######################################################
    for(i in 1:(alls)){
		if(any(result[17][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0)){
			getit4[i]<-paste(getit4[i],"/",paste0(result[17][[1]][((i-1)*entry+1):((i-1)*entry+entry)][result[17][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0],collapse="/"),sep="")
		}
	}
    # (4) tolerance: small or large? ###########################################
    for(i in 1:(alls)){
		for(j in 1:entry){
			if(result[18][[1]][(i-1)*entry+j]==1){
				getit5[i]<-paste(getit5[i],"small",sep="/")
			};
			if(result[18][[1]][(i-1)*entry+j]==2){
				getit5[i]<-paste(getit5[i],"large",sep="/")
			};
		}
	};
    # (5) charge level: ########################################################
    for(i in 1:(alls)){
		if(any(result[19][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0)){
			getit6[i]<-paste(getit6[i],"/",paste0(result[19][[1]][((i-1)*entry+1):((i-1)*entry+entry)][result[19][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0],collapse="/"),sep="")
		}
	}	
    if(result[13][[1]]!=entry){cat("WARNING: entry overflow -> links missing! Decrease mztol? Increase entry argument?")};
	rm(result);
    #dyn.unload(paste(.libPaths(),"/nontarget/temp/massCpp.dll",sep=""));
    #data.frame(samples[,1],samples[,3],getit4,getit2,getit1,getit5,getit6);
    cat("done.");
    ############################################################################

    ############################################################################
    if(deter==FALSE){
    cat("\n (3) Check plausibility ... ");
    # (3) remove invalid dmass-links based on rules (1) to (3) #################
    for(i in 1:length(getit4)){
    if(getit4[i]!="0"){ # anything to check?
      #
      ##########################################################################
      # (3.1) RULE1: intensity ratio check over ALL isotopes at charge level ###
      that1<-c(strsplit(getit4[i],"/")[[1]][-1]);
      that2<-rep(TRUE,length(that1));
      if(rules[1]=="TRUE"){
        that3<-c(strsplit(getit1[i],"/")[[1]][-1]);
        that15<-c(strsplit(getit5[i],"/")[[1]][-1]);
        that100<-c(as.numeric(strsplit(getit6[i],"/")[[1]][-1]));
        for(j in 1:length(that1)){
        if(isomat[as.numeric(that3[j]),5]==1){ # exclude double-distanced peaks
          # including intensity tolerances "inttol"
          possnumb<-c(((samples[as.numeric(that1[j]),2]*(1-inttol))/(samples[i,2]*(1+inttol)))*mpoldnew); # based on isotope with largest mpoldnew
          if(possnumb<1){possnumb<-c(1)};
          if((possnumb*3.0078)>(samples[i,1]*max(that100))){ # largest charge z, minimum mass of hydrogen
              that2[j]<-FALSE;
              countrem1<-c(countrem1+1);
          }; # if ...
        }; # if single-distanced ...
        }; # for j ...
        # remove entries in getit4[i] ("to"), getit2[->i] ("from") and isomat!
        that2[isomat[as.numeric(that3),5]!=1]<-TRUE; # if double-distanced
        for(m in 1:length(that2)){
          if(that2[m]==FALSE){
            getit2[as.numeric(that1[m])]<-sub(paste("/",i,sep=""),"",getit2[as.numeric(that1[m])])
          };
        };
        if(any(that2==FALSE)){
          isomat[as.numeric(that3[that2==FALSE]),4]<-isomat[as.numeric(that3[that2==FALSE]),4]-1;
        };
        that1<-that1[that2];
        that3<-that3[that2];
        that15<-that15[that2];
        that100<-that100[that2];
        if(any(that2)){
          that4<-c();
          that5<-c();
          that20<-c();
          that200<-c();
          for(n in 1:length(that1)){
            that4<-paste(that4,that1[n],sep="/");
            that5<-paste(that5,that3[n],sep="/");
            that20<-paste(that20,that15[n],sep="/");
            that200<-paste(that200,that100[n],sep="/");
          };
          getit4[i]<-paste("0",that4,sep="");
          getit1[i]<-paste("0",that5,sep="");
          getit5[i]<-paste("0",that20,sep="");
          getit6[i]<-paste("0",that200,sep="");
        }else{
          getit4[i]<-"0";
          getit1[i]<-"none"; # (1) which isotope?
          getit5[i]<-"0";
          getit6[i]<-"0";
        };
      } # RULE1
      #
      ##########################################################################
      # (3.2) RULE2: intensity ratio check with LARGE mass tolerance ###########
      if(any(that2)){ # anything left to check?
      if(rules[2]=="TRUE"){
          that1<-c(strsplit(getit4[i],"/")[[1]][-1]);
          that2<-rep(TRUE,length(that1));
          that3<-c(strsplit(getit1[i],"/")[[1]][-1]);
          that6<-c(that1[duplicated(that1)==FALSE]); 
          that7<-c(isomat[as.numeric(that3),5]==1);
          #that15<-c(strsplit(getit5[i],"/")[[1]][-1]);
          that100<-c(as.numeric(strsplit(getit6[i],"/")[[1]][-1]));
          for(j in 1:length(that6)){
          if(any(that7[that1==that6[j]])){ # only on single-distanced peaks!
            mpoldnew2<-min(1/isomat[as.numeric(that3[that1==that6[j]&that7==TRUE]),3]);
            ratioC<-max(isomat[as.numeric(that3[that1==that6[j]&that7==TRUE]),6]);
            possnumb<-min(((samples[as.numeric(that1[that1==that6[j]]),2]*(1-inttol))/(samples[i,2]*(1+inttol)))*mpoldnew2);
            if(possnumb<1){possnumb<-c(1)};
            possmass<-min(isos[as.logical(match(isos[,1],isos[as.logical(match(isos[,2],isomat[as.numeric(that3[that1==that6[j]&that7==TRUE]),1],nomatch=FALSE)),1],nomatch=FALSE)),3])
            if(ratioC!=0){
              if((possnumb*possmass+((possnumb/ratioC)*12))>(samples[i,1]*max(isomat[as.numeric(that3[that1==that6[j]&that7==TRUE]),7]))){ # include ratios to C! BEWARE if ration=0->Inf->alway>mass
                  that2[as.numeric(that1)==as.numeric(that6[j])]<-FALSE;
                  countrem2<-c(countrem2+1);
              }; # if ...
            }else{
              if((possnumb*possmass)>(samples[i,1]*max(isomat[as.numeric(that3[that1==that6[j]&that7==TRUE]),7]))){ # include ratio to C
                  that2[that1==that6[j]]<-FALSE;
                  countrem2<-c(countrem2+1);
              }; # if ...
            };
          }; #
          }; # for
          # remove entries in getit2[i] (i.e."to") only  / NOT "from!"
          # reset: keep non-single distanced peaks!
          that2[isomat[as.numeric(that3),5]!=1]<-TRUE;
          for(m in 1:length(that2)){
            if(that2[m]==FALSE){
              getit2[as.numeric(that1[m])]<-sub(paste("/",i,sep=""),"",getit2[as.numeric(that1[m])])
            };
          };
          if(any(that2==FALSE)){
            isomat[as.numeric(that3[that2==FALSE]),4]<-isomat[as.numeric(that3[that2==FALSE]),4]-1
          };
          that1<-that1[that2];
          that3<-that3[that2];
          that15<-that15[that2];
          that100<-that100[that2];
          if(any(that2)){
            that4<-c();
            that5<-c();
            that20<-c();
            that200<-c();
            for(n in 1:length(that1)){
              that4<-paste(that4,that1[n],sep="/");
              that5<-paste(that5,that3[n],sep="/");
              that20<-paste(that20,that15[n],sep="/");
              that200<-paste(that200,that100[n],sep="/");
            };
            getit4[i]<-paste("0",that4,sep="");
            getit1[i]<-paste("0",that5,sep="");
            getit5[i]<-paste("0",that20,sep="");
            getit6[i]<-paste("0",that200,sep="");
          }else{
            getit4[i]<-"0";
            getit1[i]<-"none"; # (1) which isotope?
            getit5[i]<-"0";
            getit6[i]<-"0";
          };
      } # if
      } # RULE2
      #
      ##########################################################################
      # (3.3) RULE3: intensity ratio check with SMALL mass tolerance ###########
      if(rules[3]=="TRUE"){
      if(any(that2)){# anything left to check?
          that1<-c(strsplit(getit4[i],"/")[[1]][-1]);
          that2<-rep(TRUE,length(that1));
          that3<-c(strsplit(getit1[i],"/")[[1]][-1]);
          that7<-c(isomat[as.numeric(that3),5]==1);
          that15<-c(strsplit(getit5[i],"/")[[1]][-1]);
          that100<-c(as.numeric(strsplit(getit6[i],"/")[[1]][-1]));
          for(j in 1:length(that1)){
          if( that15[j]=="small" && that7[j]==TRUE ){ # within small mass tolerance? single distanced?
            mpoldnew2<-c(1/isomat[as.numeric(that3[j]),3]);
            ratioC<-c(isomat[as.numeric(that3[j]),6]);
            possnumb<-c(((samples[as.numeric(that1[j]),2]*(1-inttol))/(samples[i,2]*(1+inttol)))*mpoldnew2);
            if(possnumb<1){possnumb<-c(1)};
            possmass<-min(isos[as.logical(match(isos[,1],isos[as.logical(match(isos[,2],isomat[as.numeric(that3[j]),1],nomatch=FALSE)),1],nomatch=FALSE)),3])
            if(ratioC!=0){
              if((possnumb*possmass+((possnumb/ratioC)*12))>(samples[i,1]*c(isomat[as.numeric(that3[j]),7]))){ # include ratios to C! BEWARE if ration=0->Inf->alway>mass
                  that2[j]<-FALSE;
              }; # if ...
            }else{
              if((possnumb*possmass)>(samples[i,1]*c(isomat[as.numeric(that3[j]),7]))){ # include ratio to C
                  that2[j]<-FALSE;
              }; # if ...
            };
          };    # if with small mass range
          };    # for all possible isotopes
          # remove entries in getit2[i] (i.e."to") only
          # reset: keep non-single distanced peaks!
          that2[isomat[as.numeric(that3),5]!=1]<-TRUE;
          for(m in 1:length(that2)){ 
            if(that2[m]==FALSE){
              #if(length(that1[that1==that1[m]])==1){ # if several hits: do not remove entry but set to "large"
                getit2[as.numeric(that1[m])]<-sub(paste("/",i,sep=""),"",getit2[as.numeric(that1[m])]);
                countrem3<-c(countrem3+1);
              #}else{
                #that2[m]<-TRUE;
                #that15[m]<-"large";
                #countrem3<-c(countrem3+1);
              #}
            }
          }
          if(any(that2==FALSE)){isomat[as.numeric(that3[that2==FALSE]),4]<-isomat[as.numeric(that3[that2==FALSE]),4]-1};
          that1<-that1[that2];
          that3<-that3[that2];
          that15<-that15[that2];
          that100<-that100[that2];
          if(any(that2)){
            that4<-c();
            that5<-c();
            that20<-c();
            that200<-c();
            for(n in 1:length(that1)){
              that4<-paste(that4,that1[n],sep="/");
              that5<-paste(that5,that3[n],sep="/");
              that20<-paste(that20,that15[n],sep="/");
              that200<-paste(that200,that100[n],sep="/");
            };
            getit4[i]<-paste("0",that4,sep="");
            getit1[i]<-paste("0",that5,sep="");
            getit5[i]<-paste("0",that20,sep="");
            getit6[i]<-paste("0",that200,sep="");
          }else{
            getit4[i]<-"0";
            getit1[i]<-"none"; # (1) which isotope?
            getit5[i]<-"0";
            getit6[i]<-"0";
          };
      }# if any in that2
      }# rule 3
      ##########################################################################
      # (3.4) RULE 4: minimum intensity for all used isotopes reached? #########
      if(rules[4]=="TRUE"){
      if(any(that2)){# anything left to check?
        that1<-c(strsplit(getit4[i],"/")[[1]][-1]);
        that2<-rep(TRUE,length(that1));
        that3<-c(strsplit(getit1[i],"/")[[1]][-1]);
        that15<-c(strsplit(getit5[i],"/")[[1]][-1]);
        that100<-c(as.numeric(strsplit(getit6[i],"/")[[1]][-1]));
        for(j in 1:length(that1)){
          if(isomat[as.numeric(that3[j]),5]==1){ # exclude double-distanced peaks
            if(((samples[as.numeric(that1[j]),2]+(inttol*samples[as.numeric(that1[j]),2]))/(samples[i,2]-(inttol*samples[i,2])))<min(isomat[,3])){
			  that2[j]<-FALSE;
              countrem4<-c(countrem4+1);        
            }
          }
        }
        # remove entries in getit4[i] ("to"), getit2[->i] ("from") and isomat!
        that2[isomat[as.numeric(that3),5]!=1]<-TRUE; # if double-distanced
        for(m in 1:length(that2)){
          if(that2[m]==FALSE){
            getit2[as.numeric(that1[m])]<-sub(paste("/",i,sep=""),"",getit2[as.numeric(that1[m])])
          };
        };
        if(any(that2==FALSE)){
          isomat[as.numeric(that3[that2==FALSE]),4]<-isomat[as.numeric(that3[that2==FALSE]),4]-1;
        };
        that1<-that1[that2];
        that3<-that3[that2];
        that15<-that15[that2];
        that100<-that100[that2];
        if(any(that2)){
          that4<-c();
          that5<-c();
          that20<-c();
          that200<-c();
          for(n in 1:length(that1)){
            that4<-paste(that4,that1[n],sep="/");
            that5<-paste(that5,that3[n],sep="/");
            that20<-paste(that20,that15[n],sep="/");
            that200<-paste(that200,that100[n],sep="/");
          };
          getit4[i]<-paste("0",that4,sep="");
          getit1[i]<-paste("0",that5,sep="");
          getit5[i]<-paste("0",that20,sep="");
          getit6[i]<-paste("0",that200,sep="");
        }else{
          getit4[i]<-"0";
          getit1[i]<-"none"; # (1) which isotope?
          getit5[i]<-"0";
          getit6[i]<-"0";
        };
      }# if any in that2
      } # rule 4
      ##########################################################################
      # (3.5) RULE 5: minimum intensity for isotopes within large reached? #####
      if(rules[5]=="TRUE"){
      if(any(that2)){ # anything left to check?
        that1<-c(strsplit(getit4[i],"/")[[1]][-1]);
        that2<-rep(TRUE,length(that1));
        that3<-c(strsplit(getit1[i],"/")[[1]][-1]);
        that15<-c(strsplit(getit5[i],"/")[[1]][-1]);
        that100<-c(as.numeric(strsplit(getit6[i],"/")[[1]][-1]));
        for(j in 1:length(that1)){
          if(isomat[as.numeric(that3[j]),5]==1){ # exclude double-distanced peaks
            if(
            ((samples[as.numeric(that1[j]),2]*(1+inttol))/(samples[i,2]*(1-inttol)))
            < min(isomat[as.numeric(that3[that1==that1[j]]),3])
            ){
              that2[j]<-FALSE;
              countrem5<-c(countrem5+1);        
            }
          }
        }
        # remove entries in getit4[i] ("to"), getit2[->i] ("from") and isomat!
        that2[isomat[as.numeric(that3),5]!=1]<-TRUE; # if double-distanced
        for(m in 1:length(that2)){
          if(that2[m]==FALSE){
            getit2[as.numeric(that1[m])]<-sub(paste("/",i,sep=""),"",getit2[as.numeric(that1[m])])
          };
        };
        if(any(that2==FALSE)){
          isomat[as.numeric(that3[that2==FALSE]),4]<-isomat[as.numeric(that3[that2==FALSE]),4]-1;
        };
        that1<-that1[that2];
        that3<-that3[that2];
        that15<-that15[that2];
        that100<-that100[that2];
        if(any(that2)){
          that4<-c();
          that5<-c();
          that20<-c();
          that200<-c();
          for(n in 1:length(that1)){
            that4<-paste(that4,that1[n],sep="/");
            that5<-paste(that5,that3[n],sep="/");
            that20<-paste(that20,that15[n],sep="/");
            that200<-paste(that200,that100[n],sep="/");
          };
          getit4[i]<-paste("0",that4,sep="");
          getit1[i]<-paste("0",that5,sep="");
          getit5[i]<-paste("0",that20,sep="");
          getit6[i]<-paste("0",that200,sep="");
        }else{
          getit4[i]<-"0";
          getit1[i]<-"none"; # (1) which isotope?
          getit5[i]<-"0";
          getit6[i]<-"0";
        };
      }# if any in that2
      } # rule 5
      ##########################################################################
      # (3.6) RULE 6: minimum intensity for isotopes within small reached? #####
      if(rules[6]=="TRUE"){
      if(any(that2)){# anything left to check?
          that1<-c(strsplit(getit4[i],"/")[[1]][-1]);
          that2<-rep(TRUE,length(that1));
          that3<-c(strsplit(getit1[i],"/")[[1]][-1]);
          that7<-c(isomat[as.numeric(that3),5]==1);
          that15<-c(strsplit(getit5[i],"/")[[1]][-1]);
          that100<-c(as.numeric(strsplit(getit6[i],"/")[[1]][-1]))
          for(j in 1:length(that1)){
            if( that15[j]=="small" && that7[j]==TRUE ){
              if(
                ((samples[as.numeric(that1[j]),2]*(1+inttol))/(samples[i,2]*(1-inttol)))
                < (isomat[as.numeric(that3[j]),3])
              ){
                that2[j]<-FALSE;
                countrem6<-c(countrem6+1);        
              }
            }
          }
          # remove entries in getit2[i] (i.e."to") only
          # reset: keep non-single distanced peaks!
          that2[isomat[as.numeric(that3),5]!=1]<-TRUE;
          for(m in 1:length(that2)){ 
            if(that2[m]==FALSE){ 
              #if(length(that1[that1==that1[m]])==1){ # if several hits: do not remove entry but set to "large"
                getit2[as.numeric(that1[m])]<-sub(paste("/",i,sep=""),"",getit2[as.numeric(that1[m])]);
                countrem6<-c(countrem6+1);
              #}else{
                #that2[m]<-TRUE;
                #that15[m]<-"large";
                #countrem6<-c(countrem6+1);
              #}
            }
          }
          if(any(that2==FALSE)){
            isomat[as.numeric(that3[that2==FALSE]),4]<-isomat[as.numeric(that3[that2==FALSE]),4]-1
          };
          that1<-that1[that2];
          that3<-that3[that2];
          that15<-that15[that2];
          that100<-that100[that2];
          if(any(that2)){
            that4<-c();
            that5<-c();
            that20<-c();
            that200<-c();
            for(n in 1:length(that1)){
              that4<-paste(that4,that1[n],sep="/");
              that5<-paste(that5,that3[n],sep="/");
              that20<-paste(that20,that15[n],sep="/");
              that200<-paste(that200,that100[n],sep="/");
            };
            getit4[i]<-paste("0",that4,sep="");
            getit1[i]<-paste("0",that5,sep="");
            getit5[i]<-paste("0",that20,sep="");
            getit6[i]<-paste("0",that200,sep="");
          }else{
            getit4[i]<-"0";
            getit1[i]<-"none"; # (1) which isotope?
            getit5[i]<-"0";
            getit6[i]<-"0";
          };
      } # if any in that2
      } # rule 6
      ##########################################################################
      # (3.7) RULE 7: carbon ratio check with SMALL mass tolerance #############
      if(rules[7]=="TRUE"){
      if(any(that2)){ # anything left to check?
          that1<-c(strsplit(getit4[i],"/")[[1]][-1]);
          that2<-rep(TRUE,length(that1));
          that3<-c(strsplit(getit1[i],"/")[[1]][-1]);
          that7<-c(isomat[as.numeric(that3),5]==1);
          that15<-c(strsplit(getit5[i],"/")[[1]][-1]);
          that100<-c(as.numeric(strsplit(getit6[i],"/")[[1]][-1]));
          for(j in 1:length(that1)){
              if(that15[j]=="small" && that7[j]==TRUE && isomat[as.numeric(that3[j]),1]!="13C" && isomat[as.numeric(that3[j]),6]!=0 ){
                mpoldnew2<-c(1/isomat[as.numeric(that3[j]),3]);
                ratioC<-c(isomat[as.numeric(that3[j]),6]);
                possnumb<-c(((samples[as.numeric(that1[j]),2]*(1-inttol))/(samples[i,2]*(1+inttol)))*mpoldnew2);
                if(possnumb<1){possnumb<-c(1)}; # must have at least one atom for a peak signal
                numbC<-c(possnumb/ratioC);
                if(numbC<1){numbC<-c(1)}; # must have at least one carbon atom!
                if(samples[i,2]*isomat[isomat[,1]=="13C",3][1]*(1/numbC)>cutint){
                  if(rules[10]==FALSE){ # beware of interfering peaks?
                    if(any(isomat[as.numeric(that3),1]=="13C")==FALSE){
                      that2[j]<-FALSE;
                      countrem7<-c(countrem7+1);
                    }
                  }else{
                    if(
                      any(
                        samples[,1]>=(samples[i,1]+(1.003355/isomat[as.numeric(that3[j]),7])-0.5) &
                        samples[,1]<=(samples[i,1]+(1.003355/isomat[as.numeric(that3[j]),7])+0.5) &
                        samples[,3]>=(samples[i,3]+rttol[1]) &
                        samples[,3]<=(samples[i,3]+rttol[2])
                       )==FALSE
                    ){
                      that2[j]<-FALSE;
                      countrem7<-c(countrem7+1);
                    }
                  }
                }
              }
          }
          # remove entries in getit2[i] (i.e."to") only !!!
          # reset: keep non-single distanced peaks!
          that2[isomat[as.numeric(that3),5]!=1]<-TRUE;
          for(m in 1:length(that2)){if(that2[m]==FALSE){getit2[as.numeric(that1[m])]<-sub(paste("/",i,sep=""),"",getit2[as.numeric(that1[m])])};}
          if(any(that2==FALSE)){isomat[as.numeric(that3[that2==FALSE]),4]<-isomat[as.numeric(that3[that2==FALSE]),4]-1};
          that1<-that1[that2];
          that3<-that3[that2];
          that15<-that15[that2];
          that100<-that100[that2];
          if(any(that2)){
            that4<-c();
            that5<-c();
            that20<-c();
            that200<-c();
            for(n in 1:length(that1)){
              that4<-paste(that4,that1[n],sep="/");
              that5<-paste(that5,that3[n],sep="/");
              that20<-paste(that20,that15[n],sep="/");
              that200<-paste(that200,that100[n],sep="/");
            };
            getit4[i]<-paste("0",that4,sep="");
            getit1[i]<-paste("0",that5,sep="");
            getit5[i]<-paste("0",that20,sep="");
            getit6[i]<-paste("0",that200,sep="");
          }else{
            getit4[i]<-"0";
            getit1[i]<-"none"; # (1) which isotope?
            getit5[i]<-"0";
            getit6[i]<-"0";
          };
      }# if any in that2
      }# rule 11
      ##########################################################################
    }; # if
    }; # for
    #it<-data.frame(ID,getit4,getit2,getit1,getit5,getit6);
    #names(it)<-c("ID","to","from","isotope","tol","charge")
    cat("done.");
    }else{
    cat("\n (3) Plausibility tests skipped. ");
    } # if deter == FALSE
    ############################################################################
    ############################################################################

    ############################################################################
    ############################################################################
    cat("\n (4) Group peaks within charge levels: ");
    # (4) group! ###############################################################
    along<-order(samples[,1],decreasing=FALSE)
    for(z in 1:length(charges)){
		group1b<-rep("0",alls);      # which group? per charge level! renew per charge level!
		group2b<-rep("0",alls);      # which interaction level?
		group5b<-rep("0",alls);      # ... and which charge?
		i<-c(1);
		while(i<alls){
			# correct entry, if peak points at itself! ###############################
			these1<-c(as.numeric(strsplit(getit4[along[i]],"/")[[1]][-1]));
			if(any(these1==along[i])){ # remove any self-reference
				these1<-these1[these1!=along[i]];
				if(length(these1)==0){
					getit1[along[i]]<-"none";
					getit4[along[i]]<-"0";
					getit5[along[i]]<-"0";
					getit6[along[i]]<-"0";
				}else{
					these1<-strsplit(getit1[along[i]],"/")[[1]][-1];
					these4<-as.numeric(strsplit(getit4[i],"/")[[1]][-1]);
					these5<-strsplit(getit5[along[i]],"/")[[1]][-1];
					these6<-strsplit(getit6[along[i]],"/")[[1]][-1];
					these1<-these1[these4!=along[i]];
					these5<-these5[these4!=along[i]];
					these6<-these6[these4!=along[i]];
					these4<-these4[these4!=along[i]];
					getit1[along[i]]<-"none";
					getit4[along[i]]<-"0";
					getit5[along[i]]<-"0";
					getit6[along[i]]<-"0";
# collapes by paste, not loop
					for(j in 1:length(these1)){getit1[along[i]]<-paste(getit1[along[i]],"/",these1[j],sep="");};
					for(j in 1:length(these4)){getit4[along[i]]<-paste(getit4[along[i]],"/",these4[j],sep="");};
					for(j in 1:length(these5)){getit5[along[i]]<-paste(getit5[along[i]],"/",these5[j],sep="");};
					for(j in 1:length(these6)){getit6[along[i]]<-paste(getit6[along[i]],"/",these6[j],sep="");};
				};
			};
			if( (getit4[along[i]]!="0") & (group1b[along[i]]==0) & (grepl(as.character(charges[z]),getit6[along[i]])) ){  # group1b: schon als M+X erfasst
				########################################################################
				these1<-c(along[i],as.numeric(strsplit(getit4[along[i]],"/")[[1]][-1]));
				# these1<-c(i,as.numeric(strsplit(getit4[i],"/")[[1]][-1]));
				these5<-as.numeric(strsplit(getit6[along[i]],"/")[[1]])[-1];
				these1<-these1[c(TRUE,these5==charges[z])];
				these1<-as.numeric(levels(as.factor(these1))); # remove double entries
				group2b[along[i]]<-c("1/0");
				group2b[these1[these1!=along[i]]]<-paste("2",group2b[these1[these1!=along[i]]],sep="/");
				allpeaks<-these1;
				if(length(these1)>1){
					newpeaks1<-these1[these1!=along[i]];
				}else{
					newpeaks1<-c()
				};
				level<-c(3);
				while(length(newpeaks1)>0){
					newlevel<-c();
					newpeaks2<-c();
					for(m in 1:length(newpeaks1)){
						these1<-c(as.numeric(strsplit(getit4[newpeaks1[m]],"/")[[1]][-1]));
						these5<-as.numeric(strsplit(getit6[newpeaks1[m]],"/")[[1]])[-1];
						these1<-c(these1[these5==charges[z]]);
						if(length(these1)>0){
							for(n in 1:length(these1)){
								if(any(these1[n]==allpeaks)!=TRUE){
									newpeaks2<-c(newpeaks2,these1[n]);
									newlevel<-c(newlevel,level);
								}; # if
							}; # for
						}; # if
					};
					# remove double entries in [newpeaks2, newlevel] #####################
					if(length(newpeaks2)>1){
						bad<-c();
						unt<-length(newpeaks2)
						for(l in 1:(unt-1)){
							if( any(  (newpeaks2[l]==newpeaks2[(l+1):unt]) & (newlevel[l]==newlevel[(l+1):unt])  ) ){bad<-c(bad,l)}
						}
						if(length(bad)>0){
							newpeaks2<-newpeaks2[-bad];
							newlevel<-newlevel[-bad];
						}
					}
					if(length(newpeaks2)>0){
						for(l in 1:length(newpeaks2)){
							group2b[newpeaks2[l]]<-paste(newlevel[l],group2b[newpeaks2[l]],sep="/");
						}
					}
					# clean for new round ################################################
					allpeaks<-c(allpeaks,newpeaks1,newpeaks2);
					allpeaks<-as.numeric(levels(as.factor(allpeaks)));
					newpeaks1<-newpeaks2;
					newpeaks1<-as.numeric(levels(as.factor(newpeaks1)));
					level<-c(level+1);
				}; # while
				these1<-allpeaks;
				these1<-these1[these1!=0]; # dispensable?
				########################################################################
				# RULE8: calculate feasible mass range -> apply as filter after grouping
				# shift to these 2, keep these 1 for rules[9]
				# these2: defines upper bound
				if(rules[8]=="TRUE"){
					topint<-samples[along[i],2];
					topmass<-samples[along[i],1];
					topcount<-ceiling(topmass/max(isomat[,6]));
					topput<-c(max(isomat[,2])/charges[z])
					mpoldnew<-min(1/isomat[,3]);
					toprep<-c(1)
					while(topint>0 && topcount>0){ #topint>cutint->what if these[1]<cutint? 
						topmass<-c(topmass+topput);
						topint<-c(topint*(1/mpoldnew)*(topcount/toprep));
						toprep<-c(toprep+1);
						topcount<-c(topcount-1)
					};
					getit<-c(samples[these1,1]<=topmass);
					if(any(getit==FALSE)){countrem8<-countrem8+1}
					getit[c(1,2)]<-TRUE;
					these2<-these1[getit];
				}else{
					these2<-these1;
				};
				########################################################################
				# RULE6 + RULE7: pattern plausibility ##################################
				plaus<-TRUE; # per se before Rule 6 is evaluated!
				if(rules[9]=="TRUE"){
					dat2<-samples[these1,];
					dat4<-seq(1:length(these1));
					these4<-these1;
					dat4<-dat4[order(dat2[,1],decreasing=FALSE)];
					these4<-these4[order(dat2[,1],decreasing=FALSE)];
					dat2<-dat2[order(dat2[,1],decreasing=FALSE),];
					monomass<-dat2[1,1];
					monointens<-dat2[1,2];
					this8<-as.numeric(strsplit(getit1[along[i]],"/")[[1]])[-1];                # this isotope in "isomat" ...
					this8b<-c(strsplit(getit5[along[i]],"/")[[1]])[-1];
					this8c<-as.numeric(strsplit(getit6[along[i]],"/")[[1]])[-1];
					this10<-as.numeric(strsplit(getit4[along[i]],"/")[[1]])[-1];               # ... for this daughter peak ID
					this11<-seq(1:length(this8));
					this8<-this8[this8b=="small" & this8c==charges[z]];
					this10<-this10[this8b=="small" & this8c==charges[z]];
					this11<-this11[this8b=="small" & this8c==charges[z]];
					this9<-c();                                                         # ... with this entry in peak group "dat2"
					for(y in 1:length(this10)){
						this9<-c(this9,dat4[these4==this10[y]])
					};
					if(length(this8)>0){
						if(rules[10]==FALSE){                        # problem: e.g. 37Cl shading 2*13C
							if(ppm==TRUE){
								mztol2<-c((mztol*dat2[1,1]/1e6));
							}else{
								mztol2<-mztol;
							};
						}else{
							mztol2<-0.5
						}
						plaus<-rep(TRUE,length(this8));
						for(l in 1:length(this8)){ # over all matches found in "isomat"
							# get number of atoms per dmass
							numb<-floor((dat2[this9[l],2]*(1-inttol))/(monointens*(1+inttol)*isomat[this8[l],3]));
							# more such peaks expected?
							if(numb>2){
								mass3<-c(dat2[1,1]+isomat[this8[l],2]); # start with M not M+1 mass!
								int3<-c(dat2[this9[l],2]);
								count<-c(2);
								while((int3[length(int3)]*(1-inttol))>=(cutint*(1+inttol)) & count<numb){
									mass3<-c(mass3,(mass3[length(mass3)]+isomat[this8[l],2]));
									int3<-c(int3,int3[length(int3)]*isomat[this8[l],3]*((numb-count+1))/((count)));
									count<-c(count+1);
								};
								if(length(mass3)>1){
									mass3<-mass3[-length(mass3)];
									int3<-int3[-length(int3)]
								};
								for(j in 1:length(mass3)){ # do not check intensities - only if masses exist
									if(any( dat2[,1]<=(mass3[j]+(mztol2)) & dat2[,1]>=(mass3[j]-(mztol2)) )){      # code can be improved: search at generation level + link!
										plaus[l]<-TRUE
									}else{
										plaus[l]<-FALSE
									};
								}; # with mztol*3!
							}; # if still plausible ...
						}; # for ...
					}; # if ...
				}; # on RULE 6
				# to have more rules inserted, the following is separated from RULE 6:
				if(all(plaus)){
					group1b[these2]<-paste(groupcount,group1b[these2],sep="/");
					group5b[these2]<-paste(charges[z],group5b[these2],sep="/");
					# estimate number of atoms! ##########################################
					these16<-c();
					if(deter==FALSE){
						for(b in 1:length(these2)){
							that1<-c(strsplit(getit4[these2[b]],"/")[[1]][-1]);
							that3<-c(strsplit(getit1[these2[b]],"/")[[1]][-1]);
							that15<-c(strsplit(getit5[these2[b]],"/")[[1]][-1]);
							that100<-c(as.numeric(strsplit(getit6[these2[b]],"/")[[1]][-1]));
							if(any(that15=="small")){
								for(d in 1:length(that15)){
									if(that15[d]=="small" && that100[d]==charges[z]){
										count<-floor((samples[as.numeric(that1[d]),2]*(1-inttol))/(samples[as.numeric(these2[b]),2]*(1+inttol)*isomat[as.numeric(that3[d]),3]));
										if(count<1){
											count<-c(1)
										}; # at least one atom for a peak
										these16<-paste(these16,"/",isos[isos[,2]==isomat[as.numeric(that3[d]),1],1],":",count,sep="");
									}
								}
							}
						} # for b
					} # deter
					groupinfo<-c(groupinfo,paste("minimum atom counts: no information",these16,sep=""));
					group3<-c(group3,groupcount);
					group6<-c(group6,charges[z]);
					that16<-c(these2[1]);
					for(f in 2:length(these2)){
						that16<-paste(that16,",",these2[f],sep="")
					}
					group4<-c(group4,that16);
					groupcount<-c(groupcount+1);
					i<-c(i+1);				
				}else{ # if not everything is plausible ...
					countrem9<-c(countrem9+1);
					# remove affected links:
					this11<-this11[plaus==FALSE];       # from rule 6!
					this8<-this8[plaus==FALSE];         # from rule 6!
					if(length(this11)==length(strsplit(getit1[along[i]],"/")[[1]][-1])){
						getit1[along[i]]<-"none";
						getit4[along[i]]<-"0";
						getit5[along[i]]<-"0";
						getit6[along[i]]<-"0";
						isomat[this8,4]<-isomat[this8,4]-1;
						group2b[these1]<-"0";
						i<-c(i+1);
					}else{
						getit1a<-strsplit(getit1[along[i]],"/")[[1]][-1][-this11];
						getit1[along[i]]<-"0";for(y in 1:length(getit1a)){getit1[along[i]]<-paste(getit1[along[i]],"/",getit1a[y],sep="")};
						getit4a<-strsplit(getit4[along[i]],"/")[[1]][-1][-this11];
						getit4[along[i]]<-"0";for(y in 1:length(getit4a)){getit4[along[i]]<-paste(getit4[along[i]],"/",getit4a[y],sep="")};
						getit5a<-strsplit(getit5[along[i]],"/")[[1]][-1][-this11];
						getit5[along[i]]<-"0";for(y in 1:length(getit5a)){getit5[along[i]]<-paste(getit5[along[i]],"/",getit5a[y],sep="")};
						getit6a<-strsplit(getit6[along[i]],"/")[[1]][-1][-this11];
						getit6[along[i]]<-"0";for(y in 1:length(getit6a)){getit6[along[i]]<-paste(getit6[along[i]],"/",getit6a[y],sep="")};
						isomat[this8,4]<-isomat[this8,4]-1;
						group2b[these1]<-"0";
						if(i>1){ i<-c(i-1);} # still peaks left for grouping? reset to previous peak. 
					};
				};
			}else{
				i<-c(i+1)
			}; # if conditions 1-3
			##########################################################################
		}; # while i
		############################################################################
		# merge results from different charge levels! ##############################
		for(x in 1:alls){ 
			if(group1b[x]!="0"){
				group1[x]<-paste(group1[x],group1b[x],sep="/");
				group1[x]<-sub("/0/","/",group1[x]);
				group2[x]<-paste(group2[x],group2b[x],sep="/");
				group2[x]<-sub("/0/","/",group2[x]);
				group5[x]<-paste(group5[x],group5b[x],sep="/");
				group5[x]<-sub("/0/","/",group5[x]);
			}
		};
		cat(paste(charges[z],"/",sep=""));
    } # for z = charge level ###################################################
    ############################################################################
    for(x in 1:alls){ 
		if(group1[x]!="0"){
			group1[x]<-substr(group1[x],1,nchar(group1[x])-2)
			group2[x]<-substr(group2[x],1,nchar(group2[x])-2)
			group5[x]<-substr(group5[x],1,nchar(group5[x])-2)
		}
	};
    #data.frame(group1,group2,group5);
    ############################################################################
    ############################################################################
    # Rule 11: remove nested groups = merge with larger / equal sized group #####
    if(rules[11]==TRUE && length(charges)>1){
      removals1<-c();
      removals2<-c();
      removals3<-c();
      removals4<-c();
      countrem11<-c(0);
      for(i in 1:alls){
        if(group1[i]!="0"){
          this1<-as.numeric(strsplit(group1[i],"/")[[1]][-1])
          if(length(this1)>1){
            for(n in 1:(length(this1)-1)){
              for(m in (n+1):length(this1)){
                  if(group6[this1[n]]!=group6[this1[m]]){ # not on same charge level (anyway impossible after grouping)!
                      this2<-as.numeric(strsplit(group4[this1[n]],",")[[1]]);
                      this3<-as.numeric(strsplit(group4[this1[m]],",")[[1]]);
                      if( any(is.na(match(this2,this3)))!=TRUE ){
                            removals1<-c(removals1,this1[n]); # keep
                            removals2<-c(removals2,this1[m]); # bed into
                            removals3<-c(removals3,i);
                            removals4<-c(removals4,group6[this1[n]]);
                        }
                      if( length(this2)!=length(this3) ){ # otherwise done above
                      if( any(is.na(match(this3,this2)))!=TRUE ){
                            removals1<-c(removals1,this1[m]);
                            removals2<-c(removals2,this1[n]);
                            removals3<-c(removals3,i);
                            removals4<-c(removals4,group6[this1[m]]);
                        }
                      }
                      }
                }
              }
            }
          } # if
        } # for
      # remove double entries ...
      bad<-c();
      unt<-length(removals1)
      for(i in 1:(length(removals1)-1)){
        if(any(removals1[i]==removals1[(i+1):unt]) & any(removals2[i]==removals2[(i+1):unt])){ bad<-c(bad,i)}
      }
      if(length(bad)>0){
        removals1<-removals1[-bad];   # group to be merged
        removals2<-removals2[-bad];   # to be kept
        removals3<-removals3[-bad];   # i
        removals4<-removals4[-bad];   # charge level to be merged
        countrem11<-length(removals4); #
        #data.frame(removals1,removals2,removals4,removals3);
        # correct group1 and group6: merge and delete!
        this1<-rep(TRUE,length(group4));
        for(i in 1:length(removals3)){
            group3[removals2[i]]<-paste(group3[removals2[i]],removals1[i],sep="/");
            group6[removals2[i]]<-paste(group6[removals2[i]],removals4[i],sep="/");
            this1[removals1[i]]<-FALSE
        };
        group3<-group3[this1];
        group4<-group4[this1];
        group6<-group6[this1];
      }else{
        countrem11<-c(0);
      }; # if bad>0
    }; # rule 8
    # data.frame(group3,group4,group6);
    # data.frame(ID,getit4,getit1,getit5,getit6);
    cat("done.");
    ############################################################################

    ############################################################################
    # only keep groups of size >=minpeaks ? ####################################
    #excl<-c()
    #if(minpeaks!=FALSE & length(group3)>0){
    #    for(i in 1:length(group4)){
    #      that<-as.numeric(strsplit(group4[i],",")[[1]])
    #      if(length(that)<minpeaks){excl<-c(excl,i)}
    #    }
    #group3<-group3[-excl];
    #group4<-group4[-excl];
    #group6<-group6[-excl];
    #}
    ############################################################################

    ############################################################################
    cat("\n (5) Create output... ");
    ############################################################################
    overlap<-rep(0,100);
    for(i in 1:alls){
        if(group1[i]!="0"){
          this11<-strsplit(group1[i],"/")[[1]];
          this11<-this11[-length(this11)];
          overlap[length(this11)]<-c(overlap[length(this11)]+1);
        }
    }
    overlap<-overlap[overlap!=0]
    if(length(overlap)>0){
      this11<-data.frame(seq(1:length(overlap)),overlap)
      names(this11)<-c("Number of groups in overlap","Peak counts")
    }else{
     this11<-"No overlaps detected"
    }
    ############################################################################
    deep<-rep(0,1000);
    for(i in 1:length(group2)){
		if(group2[i]!="0"){deep[as.numeric(strsplit(group2[i],"/")[[1]])]<-c((deep[as.numeric(strsplit(group2[i],"/")[[1]])])+1);
    }};
    deep<-deep[deep!=0];
    if(length(deep)>0){
        deep<-data.frame(seq(1:length(deep)),deep);
        names(deep)<-c("interaction level","peak counts");
    }else{
        deep<-"No groups formed"
    }
    ############################################################################
    hits<-data.frame(isomat[,c(1,7,4)],rep(0,length(isomat[,1])),rep("0",length(isomat[,1])),stringsAsFactors=FALSE);
    names(hits)<-c("isotope","charge","peak counts","group counts","element");
    # increment counts
    for(j in 1:length(getit1)){
		if(getit1[j]!="none"){
			this1<-as.numeric(strsplit(getit1[j],"/")[[1]][-1]);
			for(n in 1:length(this1)){
				hits[this1[n],3]<-c( hits[this1[n],3]+ 1)
			}
		}
    }
    # group counts
    if(length(group4)>0){ # anything found at all?
      for(j in 1:length(group4)){
        hit<-c();
        this<-as.numeric(strsplit(as.character(group4[j]),",")[[1]]);
        that<-as.numeric(strsplit(as.character(group6[j]),"/")[[1]]);
        for(n in 1:length(this)){
          this1<-as.numeric(strsplit(getit1[this[n]],"/")[[1]][-1]);
          this2<-as.numeric(strsplit(getit4[this[n]],"/")[[1]][-1]);
          this3<-as.numeric(strsplit(getit6[this[n]],"/")[[1]][-1]);
          if(length(this1)>0){
          for(m in 1:length(this1)){
            if(any(this==this2[m]) & any(that==this3[m])){
               hit<-c(hit,this1[m])
            }
          }
          }
        }
        hit<-as.numeric(levels(as.factor(hit)))
        hits[hit,4]<-c( hits[hit,4]+ 1)
      }
      hits<-hits[order((hits[,1]),(hits[,2]),decreasing=FALSE),];
      hits[,5]<-as.character(hits[,5]);
    }
    if(deter==FALSE){
      hits[,5]<-as.character(hits[,5]);
      for(i in 1:length(hits[,1])){
        hits[i,5]<-as.character(iso[[1]][,1][as.character(iso[[1]][,2])==as.character(hits[i,1])])[1]
      }
    }
    ############################################################################
    removals<-data.frame(seq(1:length(rules)),rep(0,length(rules)),stringsAsFactors=FALSE);
    names(removals)=c("Rule","Counts");
    removals[1,2]=countrem1;
    removals[2,2]=countrem2;
    removals[3,2]=countrem3;
    removals[4,2]=countrem4;
    removals[5,2]=countrem5;
    removals[6,2]=countrem6;
    removals[7,2]=countrem7;
    removals[8,2]=countrem8;
    removals[9,2]=countrem9;    
    removals[11,2]=countrem11;
    removals<-removals[-c(10),]; 
    if(rules[10]==TRUE){removals[7,1]<-"7/10";removals[9,1]<-"9/10"};
    ############################################################################
    # correct entries:
    for(i in 1:alls){
      if(getit2[i]!="0"){getit2[i]<-substr(getit2[i],3,nchar(getit2[i]))};
      if(getit4[i]!="0"){
		getit4[i]<-substr(getit4[i],3,nchar(getit4[i]))
		this30<-as.numeric(strsplit(getit4[i],"/")[[1]])
		this31<-as.character(ID[this30[1]])
		if(length(this30)>1){
			for(j in 2:length(this30)){
				this31<-paste(this31,"/",ID[this30[j]],sep="")
			}
		}
		getit4[i]<-this31;
	  };
      if(getit5[i]!="0"){getit5[i]<-sub("0/","",getit5[i])};
      if(getit6[i]!="0"){getit6[i]<-sub("0/","",getit6[i])};
      if(getit1[i]!="none"){
        this12<-as.numeric(strsplit(getit1[i],"/")[[1]][-1]);
        if(length(this12)==1){
          getit1[i]<-isomat[this12,1];
        }else{
          getit1[i]<-isomat[this12[1],1]
          for(j in 2:length(this12)){getit1[i]<-paste(getit1[i],isomat[this12[j],1],sep="/");}
        }
      };
    if(group1[i]!="0"){group1[i]<-substr(group1[i],3,(nchar(group1[i])))};
    if(group5[i]!="0"){group5[i]<-substr(group5[i],3,(nchar(group5[i])))};
    if(group2[i]!="0"){group2[i]<-substr(group2[i],3,(nchar(group2[i])))};
    }
    if(length(groupinfo)>0){
      for(i in 1:length(groupinfo)){
         if(groupinfo[i]!="minimum atom counts: no information"){groupinfo[i]<-sub("no information/","",groupinfo[i]);}
      }
    }else{
      groupinfo<-"No groups detected"
    }
    ############################################################################
    groupcount<-data.frame(charges,rep(0,length(charges)),stringsAsFactors=FALSE);
    names(groupcount)<-c("Charge level","Counts");
    that2<-seq(1:length(charges));
    if(length(group6)>0){
      for(i in 1:length(group6)){
            that1<-as.numeric(strsplit(as.character(group6[i]),"/")[[1]]);
            that3<-c();
            for(j in 1:length(that1)){that3<-c(that3,that2[charges==that1[j]])}
            groupcount[that3,2]<-c(groupcount[that3,2]+1);
      }
    }else{
      groupcount<-"No groups detected"
    }
    ############################################################################
    grouped_samples<-data.frame(samples,ID,group1,group2,getit4,getit1,getit5,getit6,stringsAsFactors=FALSE);
    grouped_samples<-grouped_samples[order(ID,decreasing=FALSE),]
	names(grouped_samples)<-c(names(samples),"peak ID","group ID","interaction level","to ID","isotope(s)","mass tolerance","charge level")
    #
    parameters<-data.frame(rttol[1],rttol[2],mztol,mzfrac,ppm,inttol,cutint,deter,stringsAsFactors=FALSE)
	#
    if(length(group4)>0){
		for(k in 1:length(group3)){
			group3[k]<-paste("/",group3[k],"/",sep="")
		}
		for(k in 1:length(group4)){
			this30<-as.numeric(strsplit(group4[k],",")[[1]])
			this31<-as.character(ID[this30[1]])
			if(length(this30)>1){
				for(j in 2:length(this30)){
					this31<-paste(this31,",",ID[this30[j]],sep="")
				}
			}
			group4[k]<-this31;
		}
		grouping<-data.frame(group3,group4,group6);
		names(grouping)<-c("group ID","peak IDs","charge level");
    }else{
		grouping<-"no groups assembled"
    }
    #
    pattern<-list(grouped_samples,parameters,grouping,groupinfo,groupcount,removals,this11,deep,hits,elements,iso[[3]],rules);
    #
    names(pattern)<-c(
		"Patterns",
		"Parameters",
		"Peaks in pattern groups",
		"Atom counts",
		"Count of pattern groups",
		"Removals by rules",
		"Number of peaks with pattern group overlapping",
		"Number of peaks per within-group interaction levels",
		"Counts of isotopes",
		"Elements",
		"Charges",
		"Rule settings"
	);
    ############################################################################
    cat("done.\n\n");
    # clean up! ##################################################################
    #rm(samples,getback,isos,isomat,manyisos,alls,ID,getit1,getit2,getit4,rtlow,rtup,mpoldnew,countrem1,countrem2,countrem3,
    #uptol,lowtol,findsam,howmany,that1,that2,that3,that6,that7,possnumb,possmass,minmass,groupcount,group1,group2,
    #countrem6,leng,this8,this9,these1,these2,level,groupcount,plaus,dat2,monomass,monointens,mas3,int3,count,this10,
    #this11,this12);
    ############################################################################
    return(pattern);
    ############################################################################

}
