combine <-
function(
	pattern,
	adduct,
	homol=FALSE,
	rules=c(FALSE,FALSE,FALSE),
	dont=FALSE
){

    ############################################################################
    # rules[1]: deal with interfering peaks also seperately ? ##################
    # rules[2]: remove single-peaked components? ###############################
    # rules[3]: remove components not part of a homologue series ? #############
    # warning(1): only one adduct found per main isotope group? ################
    # warning(2): adduct ranks ok? #############################################
    # warning(3): any interferring peaks? ######################################
    ############################################################################
    # (0.1) define parameters & check ##########################################
    cat("\n (1) Assemble lists...");
    comp1a<-c();    # = ID of pattern group 
    comp1b<-c();    # = IDs of peaks in main pattern group
    comp2a<-c();    # = ID of adduct group(s)
    comp2b<-c();    # = IDs of peaks in adduct group
    comp2c<-c();    # = main adduct, i.e. that of comp1b
    comp2d<-c();    # = other adducts, i.e. those of comp2b
    comp3<-c();     # = ID of homologue serie(s) 
    comp4<-c();     # = IDs of interfering peaks
    comp5<-c();     # = IDs of interfering pattern groups
    comp6<-c();     # = IDs of interfering adduct groups
    comp7<-c();     # = consistent? -> ok/(1,2,3,4,5,6)
    comp8<-c();     # = isotope relations
    comp9<-c();     # = charge levels
    if(length(pattern)>1 & length(adduct)>1){
		if(length(pattern[[1]][,1])!=length(adduct[[1]][,1])){stop(("Different data sets pattern<->adduct used for combining!"))}
		if(any(all(pattern[[1]][,1]==adduct[[1]][,1])!=TRUE)){stop("Different data sets pattern<->adduct used for combining!")};
    }
    if(length(homol)>1 & length(adduct)>1){
		if(length(homol[[1]][,1])!=length(adduct[[1]][,1])){stop(("Different data sets homol<->adduct used for combining!"))}
		if(any(all(homol[[1]][,1]==adduct[[1]][,1])!=TRUE)){stop("Different data sets pattern<->adduct used for combining!")}
    };
    if(length(homol)>1 & length(pattern)>1){
		if(length(homol[[1]][,1])!=length(pattern[[1]][,1])){stop(("Different data sets homol<->pattern used for combining!"))}    
		if(any(all(homol[[1]][,1]==pattern[[1]][,1])!=TRUE)){stop("Different data sets pattern<->adduct used for combining!")}
    };
    if(rules[3]==TRUE & length(homol)<1){stop("rule3 = TRUE not applicable if no correct homologue series are provided!")};
    if(dont[1]!=FALSE){
		if(length(dont)>4){stop("Invalid dont argument")};
		if(any(dont>4)){stop("Invalid dont argument")};
		if(any(dont<1)){stop("Invalid dont argument")};
    }
    if(length(pattern)>1 & length(adduct)>1){
		if(pattern[[11]][1]!="FALSE"){
			if((adduct[[2]][4]=="negative" & any(as.numeric(pattern[[11]])>0)) || (adduct[[2]][4]=="positive" & any(as.numeric(pattern[[11]])<0))){
				warning("Are charges of adduct vs. pattern groups consistent?")
			}
		}
	}
    ###########################################################################
    # (0.2) local function definitions ########################################
    rankit<-function(a,b,inttol){
		# on input vector a ###################################################
		a_low<-a*(1-(inttol*2))
		a_up<-a*(1+(inttol*2))
		orda<-order(a,decreasing=TRUE)
		ord_a<-rep(0,length(a))
		dor<-c(length(a));
		ord_a[orda[1]]<-c(length(a));
		for(j in 2:length(a)){
           if(a_up[orda[j]]<a_low[orda[(j-1)]]){dor<-c(dor-1);}
           ord_a[orda[j]]<-dor;
		}
		# on input vector b ###################################################
		b_low<-b*(1-(inttol*2))
		b_up<-b*(1+(inttol*2))
		ordb<-order(b,decreasing=TRUE)
		ord_b<-rep(0,length(b))
		dor<-c(length(b));
		ord_b[ordb[1]]<-c(length(b));
		for(j in 2:length(b)){
           if(b_up[ordb[j]]<b_low[ordb[(j-1)]]){dor<-c(dor-1);}
           ord_b[ordb[j]]<-dor;
		}
		#######################################################################
		if(sum(ord_a)>=sum(ord_b)){
			return(all(ord_a>=ord_b))
		}else{
			return(all(ord_a<=ord_b))
		}
    }
    cat("done.");
    ##########################################################################

    ##########################################################################
    # combine # 1 ############################################################
    cat("\n (2) Combine...");
    if(length(pattern[[1]])>1){
		no3<-rep(TRUE,length(pattern[[1]][,1]));
		intord<-order(pattern[[1]][,2],decreasing=TRUE);
    }else{
		no3<-rep(TRUE,length(adduct[[1]][,1]));
		intord<-order(adduct[[1]][,2],decreasing=TRUE);
    }
    if(length(homol[[1]])>1){ # save homologue IDs per peak ID ###############
		getit<-rep(0,length(no3));
		for(i in 1:length(homol[[3]][,1])){
			compit<-as.numeric(strsplit(as.character(homol[[3]][i,2]),",")[[1]]);
			for(m in 1:length(compit)){
				getit[compit[m]]<-paste(getit[compit[m]],"/",i,sep="");
			}
		}
    }
    ##########################################################################
    for(i in 1:length(intord)){ 
    if(no3[intord[i]]!=FALSE){ # if not yet included in some group ###########
        ######################################################################
        comp7<-c(comp7,"-");
        ######################################################################
        # on isotope pattern #################################################
        ######################################################################
        if(length(pattern[[1]])>1){ # if available ...
            if(pattern[[1]][intord[i],5]!="0"){
                allpeaks<-c(intord[i]);
                get1<-as.numeric(strsplit(as.character(pattern[[1]][allpeaks,5]),"/")[[1]]);
                get2<-c()
                get2a<-as.numeric(strsplit(as.character(pattern[[3]][grepl(paste("/",get1[1],"/",sep=""),
                as.character(pattern[[3]][,1]),fixed=TRUE),2]),",")[[1]]);
                get2<-c(get2,get2a);
                comp1a<-c(comp1a,as.character(pattern[[1]][allpeaks,5]));       # = ID of pattern group 
                compit<-as.character(get2[1]);
                for(m in 2:length(get2)){compit<-paste(compit,",",get2[m],sep="")};
                comp1b<-c(comp1b,compit);                                       # = IDs of peaks in pattern group
                allpeaks<-c(allpeaks,get2);
                allpeaks<-as.numeric(levels(as.factor(allpeaks)));
                no3[allpeaks]<-FALSE;
                get6<-c()
                for(k in 1:length(allpeaks)){
                    if(pattern[[1]][allpeaks[k],7]!="0"){
						get3<-as.numeric(strsplit(as.character(pattern[[1]][allpeaks[k],7]),"/")[[1]]);
						get4<-strsplit(as.character(pattern[[1]][allpeaks[k],8]),"/")[[1]];
						get5<-strsplit(as.character(pattern[[1]][allpeaks[k],9]),"/")[[1]];
						for(y in 1:length(get3)){
							if(any(allpeaks==get3[y])){
								get6<-c(get6,paste(get4[y],"(",get5[y],")",sep=""));
							}
						}
                    }
                }
                if(length(get6)>0){
					if(length(get6)>1){get6<-as.character(levels(as.factor(get6)))};
					get7<-get6[1]
					if(length(get6)>1){
						for(z in 2:length(get6)){
							get7<-paste(get7,"/",get6[z],sep="");                                      
						}
					}
                }
                comp8<-c(comp8,get7);
                comp9<-c(comp9,as.character(pattern[[3]][grepl(paste("/",get1[1],"/",sep=""),as.character(pattern[[3]][,1]),fixed=TRUE),3]));
                if(length(strsplit(as.character(pattern[[3]][grepl(paste("/",get1[1],"/",sep=""),as.character(pattern[[3]][,1]),fixed=TRUE),3]),"/")[[1]])>1){
					comp7[length(comp7)]<-paste(comp7[length(comp7)],",4",sep="")
                }      
            }else{
				comp1a<-c(comp1a,"-");   # = ID of pattern group 
                comp1b<-c(comp1b,as.character(intord[i]));   # = IDs of peaks in pattern group
                comp8<-c(comp8,"-"); # = isotope relations found
                comp9<-c(comp9,"-"); # = isotope relations found
                allpeaks<-c(intord[i]);
                no3[intord[i]]<-FALSE;
            }
        }else{
            comp1a<-c(comp1a,"-");   # = ID of pattern group 
            comp1b<-c(comp1b,as.character(intord[i]));   # = IDs of peaks in pattern group
            comp8<-c(comp8,"-"); # = isotope relations found
            comp9<-c(comp9,"-"); # = isotope relations found
            allpeaks<-c(intord[i]);
            no3[intord[i]]<-FALSE;
        }
        ######################################################################
        # on adduct pattern ##################################################
        ######################################################################
        if(length(adduct[[1]])>1){ # if available ...
            thesepeaks<-c();                                         
            addID<-c();
            addfrom<-c();
            addto<-c();
            for(m in 1:length(allpeaks)){
				get1<-as.numeric(strsplit(as.character(adduct[[1]][allpeaks[m],5]),"/")[[1]]);
				get3<-strsplit(as.character(adduct[[1]][allpeaks[m],7]),"//")[[1]];
				if(get1[1]!=0){
					for(n in 1:length(get1)){
						get2<-as.numeric(strsplit(as.character(adduct[[3]][get1[n],2]),",")[[1]]);
						for(k in 1:length(get2)){
							if(any(allpeaks==get2[k])==FALSE){
								thesepeaks<-c(thesepeaks,get2[k]);
							};
						};
						addID<-c(addID,get1[n]);
					};
					for(n in 1:length(get3)){
						get4<-strsplit(get3[n],"<->")[[1]];
						addfrom<-c(addfrom,get4[1]);
						addto<-c(addto,get4[2]);
					};
				};
            };
            if(length(thesepeaks)>0){
                thesepeaks<-as.numeric(levels(as.factor(thesepeaks)));
                addID<-as.numeric(levels(as.factor(addID)));
                # store adduct group IDs #######################################
                compit<-as.character(addID[1]);
                if(length(addID)>1){
					for(m in 2:length(addID)){
						compit<-paste(compit,"/",as.character(addID[m]),sep="");
					}
                }
                comp2a<-c(comp2a,compit);   # = ID of adduct group(s)
                # store adduct peak IDs ########################################
                compit<-as.character(thesepeaks[1]);
                if(length(thesepeaks)>1){
					for(m in 2:length(thesepeaks)){
						compit<-paste(compit,",",as.character(thesepeaks[m]),sep="");
					}
                }
                comp2b<-c(comp2b,compit);   # = IDs of peaks in adduct group
                # store adduct from ############################################
                addfrom<-as.character(levels(as.factor(addfrom)));               
                compit<-as.character(addfrom[1]);
                if(length(addfrom)>1){
					for(m in 2:length(addfrom)){
						compit<-paste(compit,",",as.character(addfrom[m]),sep="");
					}
                }
                comp2c<-c(comp2c,compit);   # = main adduct, i.e. that of comp1b
                if(length(addfrom)>1){comp7[length(comp7)]<-paste(comp7[length(comp7)],",1",sep="")}
                # store addut to ###############################################
                addto<-as.character(levels(as.factor(addto)));
                compit<-as.character(addto[1]);
                if(length(addto)>1){
					for(m in 2:length(addto)){
						compit<-paste(compit,",",as.character(addto[m]),sep="");
					}
                }
                comp2d<-c(comp2d,compit);   # = other adducts, i.e. those of comp2b
                ################################################################
                # check adduct plausbility #####################################
                # (1) same adducts? ... done ###################################
                if( length(addfrom)==1 & length(allpeaks)>1 ){ 
					# allpeaks>1 also ensures that pattern data set available! #
					# (0) derive raltion matrix ################################
					mat<-matrix(ncol=length(allpeaks),nrow=(length(addto)+1),0);
					colnames(mat)<-allpeaks;
					rownames(mat)<-c(addfrom,addto);
					for(z in 1:length(allpeaks)){
						mat[1,z]<-adduct[[1]][allpeaks[z],2];
						geta<-as.numeric(strsplit(as.character(adduct[[1]][allpeaks[z],6]),"/")[[1]]);
						getb<-strsplit(as.character(adduct[[1]][allpeaks[z],7]),"//")[[1]];
						if(geta[1]!=0){
							for(y in 1:length(geta)){
								mat[rownames(mat)==strsplit(getb[y],"<->")[[1]][2],z]<-adduct[[1]][geta[y],2];
							};
						};
					};
					cutint<-pattern[[2]][[7]];
					inttol<-pattern[[2]][[6]];
					mat[(mat*(1-2*inttol))<cutint]<-0;
					# (2) adduct intensity ranks ok? ###########################
					for(z in 2:length(mat[,1])){
						if(rankit(mat[1,],mat[z,],inttol)==FALSE){
							comp7[length(comp7)]<-paste(comp7[length(comp7)],",2",sep="")
						}
					}
					# (3) relative adduct intensities ok? ######################
					# ... skipped ##############################################
                }
                ################################################################
                allpeaks<-c(allpeaks,thesepeaks);
                allpeaks<-as.numeric(levels(as.factor(allpeaks)));
                no3[allpeaks]<-FALSE; 
            }else{
                comp2a<-c(comp2a,"-");   # = ID of adduct group(s)
                comp2b<-c(comp2b,"-");   # = IDs of peaks in adduct group
                comp2c<-c(comp2c,"-");   # = main adduct, i.e. that of comp1b
                comp2d<-c(comp2d,"-");   # = other adducts, i.e. those of comp2b         
            }
        }else{
            comp2a<-c(comp2a,"-");   # = ID of adduct group(s)
            comp2b<-c(comp2b,"-");   # = IDs of peaks in adduct group
            comp2c<-c(comp2c,"-");   # = main adduct, i.e. that of comp1b
            comp2d<-c(comp2d,"-");   # = other adducts, i.e. those of comp2b         
        }
        ######################################################################
        # on interfering peaks ###############################################
        allpeaks<-as.numeric(levels(as.factor(allpeaks)));
        oldpeaks<-allpeaks;
        if(length(allpeaks)>1){
            IDpeak<-c("");
            IDpat<-c("");
            IDadd<-c("");              
            newpeaks1<-c(allpeaks);
            newpeaks2<-c();
            while(length(newpeaks1)>0){
				for(j in 1:length(newpeaks1)){
					##################################################################
					# find other pattern groups ######################################
					if(length(pattern)>1){
						that<-strsplit(as.character(pattern[[1]][newpeaks1[j],5]),"/")[[1]];
						if(that[1]!="0"){
							for(n in 1:length(that)){
								if( grepl(paste("/",that[n],"/",sep=""),as.character(comp1a[i]))!=TRUE & grepl(paste("/",that[n],"/",sep=""),IDpat)!=TRUE ){
									this<-as.numeric(strsplit(as.character(pattern[[3]][grepl(paste("/",that[n],"/",sep=""),pattern[[3]][,1]),2]),",")[[1]])
									for(m in 1:length(this)){
										doit<-FALSE;
										if(any(oldpeaks==this[m])!=TRUE){
											IDpeak<-paste(IDpeak,",",as.character(this[m]),sep="");
											newpeaks2<-c(newpeaks2,this[m]);
											doit<-TRUE;
										}
										if(doit==TRUE){
											IDpat<-paste(IDpat,paste("/",as.character(that[n]),"/",sep=""),sep="");
										}
									}
								}
							}
						}
					}
					# find other adduct groups #######################################
					if(length(adduct)>1){
						that<-strsplit(as.character(adduct[[1]][newpeaks1[j],5]),"/")[[1]];
						if(that[1]!="0"){
							for(n in 1:length(that)){
								if( grepl(paste("/",that[n],"/",sep=""),as.character(comp2a[i]))!=TRUE & grepl(paste("/",that[n],"/",sep=""),IDadd)!=TRUE ){
									this<-as.numeric(strsplit(as.character(adduct[[3]][grepl(paste("/",that[n],"/",sep=""),adduct[[3]][,1]),2]),",")[[1]])
									for(m in 1:length(this)){
											doit<-FALSE;
											if(any(oldpeaks==this[m])!=TRUE){
												IDpeak<-paste(IDpeak,",",as.character(this[m]),sep="");
												newpeaks2<-c(newpeaks2,this[m]);
												doit<-TRUE;
											}
											if(doit==TRUE){
												IDadd<-paste(IDadd,paste("/",as.character(that[n]),"/",sep=""),sep="");
											}
									}
								}
							}
						}
					}
					##################################################################
				}		
				newpeaks2<-as.numeric(as.character(levels(as.factor(newpeaks2))));
				oldpeaks<-c(oldpeaks,newpeaks2)
				newpeaks1<-newpeaks2;
				newpeaks2<-c();
            } # while newpeaks1
            if(nchar(IDpeak)>0){
				comp7[length(comp7)]<-paste(comp7[length(comp7)],",3",sep="")
				comp4<-c(comp4,substr(IDpeak,2,nchar(IDpeak)));
            }else{
				comp4<-c(comp4,"-");
            }
			if(nchar(IDpat)>0){
				comp5<-c(comp5,IDpat);
            }else{
				comp5<-c(comp5,"-");
            }
            if(nchar(IDadd)>0){
				comp6<-c(comp6,IDadd);
            }else{
				comp6<-c(comp6,"-");
            }  
        }else{
			comp4<-c(comp4,"-");
            comp5<-c(comp5,"-");
            comp6<-c(comp6,"-");
        }
        ######################################################################
        # rules[1]: deal with interfering peaks also seperately ? ############
        if(rules[1]!=TRUE){
			allpeaks<-oldpeaks;
			no3[allpeaks]<-FALSE;
        }
        ######################################################################
        # on homologue series ################################################          
        if(length(homol[[1]])>1){ # if available ...
            compit<-c();
            for(m in 1:length(allpeaks)){
				if(getit[allpeaks[m]]!="0"){
					compit<-paste(compit,substr(getit[allpeaks[m]],2,nchar(getit[allpeaks[m]])),sep="");
				}
            }
            if(length(compit)>0){
				comp3<-c(comp3,compit);
            }else{
				comp3<-c(comp3,"-");    # = ID of homologue serie(s)             
            }
        }else{
			comp3<-c(comp3,"-");    # = ID of homologue serie(s) 
        }
        ######################################################################
    }
    }
    # data.frame(comp1a,comp1b,comp2a,comp2b,comp3)
    # data.frame(comp1a,comp1b,comp2a,comp2b,comp3,comp4,comp5,comp6,comp7)
    cat("done.");
    ############################################################################
        
    ############################################################################
    cat("\n (3) Apply rules 2&3...");
    # remove single-peaked components? #########################################
    if(rules[2]==TRUE){
		getit<-c()
		for(i in 1:length(comp1b)){
			if(length(strsplit(comp1b[i],",")[[1]])==1){
				if(comp2b[i]=="-" & comp3[i]=="-" ){
					getit<-c(getit,i);
				}
			}
		}
		comp1a<-comp1a[-getit]    # = ID of pattern group 
		comp1b<-comp1b[-getit]    # = IDs of peaks in pattern group
		comp2a<-comp2a[-getit]    # = ID of adduct group(s)
		comp2b<-comp2b[-getit]    # = IDs of peaks in adduct group
		comp2c<-comp2c[-getit]    # = main adduct, i.e. that of comp1b
		comp2d<-comp2d[-getit]    # = other adducts, i.e. those of comp2b
		comp3<-comp3[-getit]      # = ID of homologue serie(s) 
		comp4<-comp4[-getit]      # = IDs of interfering peaks
		comp5<-comp5[-getit]      # = IDs of interfering pattern groups
		comp6<-comp6[-getit]      # = IDs of interfering adduct groups
		comp7<-comp7[-getit]      # = consistent? -> ok/(1,2,3,4,5,6)
		comp8<-comp8[-getit]      # = isotope relations in main pattern group
		comp9<-comp9[-getit]      # = charge level
    }
    ############################################################################
    # remove components not part of a homologue series ? #######################
    if(rules[3]==TRUE & length(homol)>1){
		getit<-c()
		for(i in 1:length(comp3)){
			if(comp3[i]!="-"){
				getit<-c(getit,i);
			};
		};
		comp1a<-comp1a[getit]    # = ID of pattern group 
		comp1b<-comp1b[getit]    # = IDs of peaks in pattern group
		comp2a<-comp2a[getit]    # = ID of adduct group(s)
		comp2b<-comp2b[getit]    # = IDs of peaks in adduct group
		comp2c<-comp2c[getit]    # = main adduct, i.e. that of comp1b
		comp2d<-comp2d[getit]    # = other adducts, i.e. those of comp2b      
		comp3<-comp3[getit]      # = ID of homologue serie(s) 
		comp4<-comp4[getit]      # = IDs of interfering peaks
		comp5<-comp5[getit]      # = IDs of interfering pattern groups
		comp6<-comp6[getit]      # = IDs of interfering adduct groups
		comp7<-comp7[getit]      # = consistent? -> ok/(1,2,3,4,5,6)
		comp8<-comp8[getit]      # = isotope relations in main pattern group
		comp9<-comp9[getit]      # = charge level
    };
    # data.frame(comp1a,comp1b,comp2a,comp2b,comp3,comp4,comp5,comp6,comp7)
    ############################################################################
    if(rules[2]==FALSE & rules[3]==FALSE){cat("skipped.")}else{cat("done.")}
    ############################################################################

    ############################################################################
    # summarize / return #######################################################
    cat("\n (4) Generate output...");
    # main list ################################################################
    comp7[comp7!="-"]<-substr(comp7[comp7!="-"],3,nchar(comp7[comp7!="-"]))
    comps<-data.frame(seq(1,length(comp1a),1),comp1a,comp1b,comp2a,comp2b,comp3,comp4,comp5,comp6,comp7,comp2c,comp2d)
    names(comps)<-c("Component ID |","ID pattern group |","ID pattern peaks |","ID adduct group(s) |","ID adduct peaks |",
    "ID homologue series |","ID interfering peaks |","ID interfering pattern group(s) |",
    "ID interfering adduct group(s) |","Warnings |","pattern group adduct|","adduct group adduct(s) |");
    # exclude components with warnings in comp 7 ###############################
    if(length(dont)>0){
		if(dont[1]!=FALSE){
			excl<-rep(TRUE,length(comp7));
			for(i in 1:length(comp7)){
				if(comp7[i]!="-"){  
					this<-as.numeric(strsplit(as.character(comp7[i]),",")[[1]])
					for(j in 1:length(this)){
						if(any(dont==this[j])){
						excl[i]<-FALSE;
					}
				}
			}
		}
		comps<-comps[excl,]
		comp8<-comp8[excl]
		comp9<-comp9[excl]
		} # if done
    } # if done
    # sort by decreasing intensity of ALL peaks and ############################ 
    # calculate mean component size ############################################
    # index peaks in components ################################################
    int<-c();
    intID<-c();
    intmz<-c();
    intrt<-c();
    numb<-c();
    allin<-c();
    partin<-c();
    for(i in 1:length(comps[,1])){
		this<-c();
		if(comps[i,3]!="-"){
			this<-c(this,as.numeric(strsplit(as.character(comps[i,3]),",")[[1]]));
		}
		if(comps[i,5]!="-"){
			this<-c(this,as.numeric(strsplit(as.character(comps[i,5]),",")[[1]]));
		}
		if(comps[i,7]!="-"){
			this<-c(this,as.numeric(strsplit(as.character(comps[i,7]),",")[[1]]));
			partin<-c(partin,as.numeric(strsplit(as.character(comps[i,7]),",")[[1]]));
		}
		this<-as.numeric(as.character(levels(as.factor(this))));
		numb<-c(numb,length(this));
		if(length(pattern)>1){
			int<-c(int,max(pattern[[1]][this,2])[1]);
			intID<-c(intID,pattern[[1]][this,4][pattern[[1]][this,2]==max(pattern[[1]][this,2])[1]][1]);
			intmz<-c(intmz,pattern[[1]][this,1][pattern[[1]][this,2]==max(pattern[[1]][this,2])[1]][1]);
			intrt<-c(intrt,pattern[[1]][this,3][pattern[[1]][this,2]==max(pattern[[1]][this,2])[1]][1]);        
		}else{
			int<-c(int,max(adduct[[1]][this,2])[1]);
			intID<-c(intID,adduct[[1]][this,4][adduct[[1]][this,2]==max(adduct[[1]][this,2])[1]][1]);
			intmz<-c(intmz,adduct[[1]][this,1][adduct[[1]][this,2]==max(adduct[[1]][this,2])[1]][1]);
			intrt<-c(intrt,adduct[[1]][this,3][adduct[[1]][this,2]==max(adduct[[1]][this,2])[1]][1]);        
		}
		allin<-c(allin,this);
    }
    comps<-cbind(comps,intID,intmz,int,intrt,comp8,comp9)
    names(comps)[13:18]<-c("HI peak ID |","HI m/z |","Highest intensity (HI) |","HI RT |","Isotope(s) d(m/z) |","z")
    comps<-comps[order(int,decreasing=TRUE),]
    comps[,1]<-seq(1,length(comps[,1]),1);
    numb<-mean(numb);
    if(length(pattern)>1){
		no1<-rep(FALSE,length(pattern[[1]][,1]));
    }else{
		no1<-rep(FALSE,length(adduct[[1]][,1]));
    }
    no1[allin]<-TRUE;
    if(length(pattern)>1){
		no2<-rep(FALSE,length(pattern[[1]][,1]));
    }else{
		no2<-rep(FALSE,length(adduct[[1]][,1]));
    }
    no2[partin]<-TRUE;
    out<-c(round(length(no1[no1])/length(no1),digits=3),round(length(no1[no1]),digits=0),round(length(no2[no2]),digits=0),round(numb,digits=2));
    names(out)<-c("Fraction of peaks in components","Number of peaks in components","Number of interfering peaks","Mean number of peaks in components");
    ############################################################################
    # get parameters ###########################################################
    if(length(pattern)>1){
		param<-pattern[[2]];
		param<-param[-1];
		if(length(adduct)>1){
			param[1]<-max(adduct[[2]][1],pattern[[2]][1],pattern[[2]][2]);
		}else{
			param[1]<-max(pattern[[2]][1],pattern[[2]][2]);
		}
    }else{
		param<-adduct[[2]];
    }
    ############################################################################
    comp<-list(comps,FALSE,FALSE,FALSE,no1,out,param);
    if(length(pattern)>1){comp[[2]]<-pattern[[1]]}else{comp[[2]]<-"No isotope pattern groups available."};
    if(length(adduct)>1){comp[[3]]<-adduct[[1]]}else{comp[[3]]<-"No adduct groups available."};
    if(length(homol)>1){comp[[4]]<-homol[[1]]}else{comp[[4]]<-"No homologue series available."};
    names(comp)<-c("Components","pattern peak list","adduct peak list","homologue list","Peaks in components","Summary","Parameters");
    cat("done.\n\n");
    return(comp);
    ############################################################################
}
