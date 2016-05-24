homol.search <-
function(
	peaklist,
	isotopes,	
	elements=c("C","H","O"),
	use_C=FALSE,
	minmz=5,
	maxmz=120,
	minrt=2,
	maxrt=2,
	ppm=TRUE,
	mztol=3.5,
    rttol=0.5,
	minlength=5,
	mzfilter=FALSE,
	vec_size=3E6,
	mat_size=3,
	R2=.98,
	spar=.45,
	plotit=FALSE,
	deb=0
){
  
    ##########################################################################
    # (0) Issue warnings: check arguments ####################################
	if(!is.logical(use_C)){stop("use_C must be logial.")}
    if(ppm==TRUE & mztol>10){cat("Too big mztol?")};
    if(minlength<3){stop("invalid minlen argument. set minlen >= 3")};
    if(mztol<=0){warning("mztol should be >0!")};
	if(length(elements)<1){stop("specify elements")}
	if(!is.numeric(mzfilter) & mzfilter[1]!="FALSE"){stop("mzfilter must either be a numeric vector or set to FALSE")}
	if(elements[1]!="FALSE"){if(any(is.na(match(elements,isotopes[,1])))){ stop("unknown elements") }}
	if(length(peaklist[1,])>3){stop("peaklist with > 3 columns not allowed")}
	if(!is.numeric(peaklist[,1]) || !is.numeric(peaklist[,2]) || !is.numeric(peaklist[,3]) ){stop("peaklist columns not numeric")}
	if(R2!=FALSE){if((R2<=0)||(R2>1)){stop("R2 must be either FALSE or 0<R2<=1")}}
	if(R2!=FALSE){if((spar<=0)||(spar>1)){stop("R2 must be either FALSE or 0<R2<=1")}}	
	if(mzfilter[1]!="FALSE"){
		if(!is.numeric(mzfilter)){stop("mzfilter must either be FALSE or (a vector of) numeric")}
		mzfilter<-mzfilter[order(mzfilter)]
		if(any(mzfilter<0)){stop(" Negative mzfilter values found - abort.")}
		if(min(mzfilter)<minmz){stop(" Minimum mzfilter value smaller than minmz - abort.")}
		if(max(mzfilter)>maxmz){stop(" Maximum mzfilter value larger than maxmz - abort.")}	
	}
	if(!is.data.frame(peaklist) & !is.matrix(peaklist) ){stop("peaklist must be a numeric data.frame")}
	if(!is.matrix(peaklist)){
		peaklist<-data.matrix(peaklist)
	}
	if(length(mztol)==1){
		if(ppm==TRUE){
			delmz<-c(mztol*max(peaklist[,1])/1e6);
			peaklist4<-c(mztol*peaklist[,1]/1e6);
		}else{
			delmz<-mztol;
			peaklist4<-rep(mztol,length(peaklist[,1]))
		}
	}else{
		if(length(mztol)!=length(peaklist[,1])){
			peaklist4<-mztol
		}else{
			stop("mztol: must be either one value or of length peaklist[,1]")
		}
	}
	inter<-interactive()
	##########################################################################
    # (1) retrieve feasible mass differences & all combinations thereof ######
	# (1.1) upper & lower mass defect / mass bound ###########################
 	##########################################################################
	cat("\n(1-3) Build bounds, trees & nearest neighbour path ... ");	
	if(elements[1]==FALSE){
		elements<-unique(as.character(isotopes[,1]))
	}
	delmass<-c();
	c_ratio<-c();
	for(i in 1:length(elements)){
		isos<-isotopes[isotopes[,1]==elements[i],]
		delmass<-c(delmass,
			isos[isos[,3]==min(isos[,3]),3]
		);
		c_ratio<-c(c_ratio,unique(isos[,5]));
	}	
	c_ratio[c_ratio==0]<-Inf # well...
	if(use_C){
		delmass<-( c( delmass-round(delmass,digits=0) ) / (delmass+(1/c_ratio*12)) );		
	}else{
		delmass<-( c( delmass-round(delmass,digits=0) ) / delmass );
	}
	maxup<-c(max(delmass));
	maxdown<-c(abs(min(delmass)));
	##########################################################################
    # (2) Set internal data & parameters #####################################
	maxmove<-c(max(maxup,maxdown))
	shift<-c(maxmz*(maxdown+maxup))
	mass_def<-c(peaklist[,1]-round(peaklist[,1]))		# calculate mass defect
	peaklist2<-peaklist[,-2]
	uplim<-c(mass_def-(peaklist2[,1]*maxup))				# upward mass defect shift - (m/z)/mass defect lower intercept
	uplim_tol<-c(max(peaklist4)*maxup*2)					# upward mass defect shift - maximum tolerance
	downlim<-c(mass_def+(peaklist2[,1]*maxdown))			# downward mass defect shift - (m/z)/mass defect upper intercept
	downlim_tol<-c(max(peaklist4)*maxdown*2)				# downward mass defect shift - maximum tolerance
	peaklist2<-cbind(peaklist2,uplim,downlim)
	peaklist2<-as.matrix(peaklist2)
	max_delmz<-c(4*max(peaklist4)) 						# maximum m/z-distance gap - used for early aborting
	scaled<-c(abs(maxmz-minmz),abs(maxrt+minrt),shift)	# scaling used for nearest neighbour-search
	if(any(scaled==0)){stop("debug me on issue #5: scaled entry ==0 -> wrong parameters=")}
	peaklist3<-cbind(peaklist2[,c(1,2)],mass_def)
	peaklist3<-as.matrix(peaklist3)
	bounds<-matrix(ncol=2,nrow=4,0)						# store search bounds
	marked<-matrix(ncol=3,nrow=length(peaklist2[,1]),0) # first sweep = 1 -> do not set to 0!
	marked[,1]<--1;
	marked[,3]<-c(1:length(peaklist2[,1]))
	colnames(marked)<-c("sweep","done","ID")
	tupels<-matrix(nrow=vec_size,ncol=7,0)				# initially used to store triplets, then any n-tupels
	tupeldo<-c(1);											# where to write into tupels
	new_found<-rep(0,length(peaklist2[,1]))				# store newly detected m/z-differences
	new_found_ceiling<-rep(0,length(peaklist2[,1])) 	# if mass defect rounding hits ceiling
	ceiled<-FALSE
	new_found_floor<-rep(0,length(peaklist2[,1]))		# if mass defect rounding hits floor
	floored<-FALSE
	dist_ID<-c()										# point ID
	dist_dist<-c()										# store distance
	mz_last<-c(0)										# store last sweep`s m/z-value for correction
	a<-c(0)												# track over sweeps: how many points within bounds?
	use<-c(1)											# for nearest neighbour search
	along<-rep(0,length(peaklist2[,1]))					# for nearest neighbour search
	along[1]<-1											# for nearest neighbour search
	##########################################################################
	# (3) build kd-trees & run nearest-neighbour search ######################
	peakTree<-.Call("kdtree", 
			as.matrix(peaklist2),
			PACKAGE="nontarget"
	);
	if(!any(deb==4)){ # along nearest neighbour path ...
		peakTree3<-.Call("kdtree", 
				as.matrix(peaklist3),
				PACKAGE="nontarget"
		);			
		for(i in 2:length(peakTree3[,1])){	
			use2<-.Call("search_kdtree3", 
				peaklist3, 	# rows: c(m/z,RT,UB,LB)
				peakTree3,  # peaks - search tree
				use,
				scaled,	
				PACKAGE="nontarget"
			)[,1]
			if(use2<1){stop("debug me on issue #4 - use2=-1,<1")}
			.Call("node_delete", 
				use,
				peaklist3,
				peakTree3,
				PACKAGE="nontarget"
			);
			along[i]<-use2
			use<-use2
		}
		if(sum(peakTree3[,3]!=0)>1){stop("debug me on issue #5: peakTree3 was not fully emptied in NN-search!")}
		if(any(peakTree3[,3]>1)){stop("debug me on issue #6: peakTree3 corrupted!")}		
	}else{ # ... or, for comparison, via random sampling
		along<-sample(
			seq(1,length(peaklist[,1]),1),
			length(peaklist[,1]),
			replace=FALSE
		)
	}
	cat("done.");
	##########################################################################
	# (4) Sweep through nearest neighbour path ###############################
	cat("\n(4) Triplet extraction \n");	
	if(inter){pBar <- txtProgressBar( min = 0, max = length(peaklist2[,1]), style = 3 )}
	for(i in 1:length(peaklist2[,1])){

		if(inter){setTxtProgressBar(pBar,i,title = NULL, label = NULL)}
		use<-along[i]
		######################################################################
		# upper area sweep ###################################################
		bounds[1,1]<-(peaklist2[use,1]+minmz)
		bounds[1,2]<-(peaklist2[use,1]+maxmz) 
		bounds[2,1]<-(peaklist2[use,2]-minrt)
		bounds[2,2]<-(peaklist2[use,2]+maxrt)
		bounds[3,1]<-((peaklist2[use,3])-shift-uplim_tol)	# maxup LB
		bounds[3,2]<-(peaklist2[use,3]+uplim_tol) 			# maxup UB
		bounds[4,1]<-(peaklist2[use,4]-downlim_tol)			# maxdown LB
		bounds[4,2]<-(peaklist2[use,4]+shift+downlim_tol) 	# maxdown UB
		.Call("search_kdtree_homol", 
			peaklist2, 	# rows: c(m/z,RT,UB,LB)
			peakTree,  	# peaks - search tree
			bounds,
			marked,		# only write to new_found those of marked[,1]!=(i-1) to omit last sweep`s new_found
			i,
			new_found,	
			1, 			# clean new_found?
			PACKAGE="nontarget"
		)
		# -> mass defect rounding hits ceiling? ##############################
		if((mass_def[use]+(maxmove*maxmz))>=.5){
			bounds[3,1]<-((peaklist2[use,3])-shift-1-uplim_tol)	# maxup LB
			bounds[3,2]<-(peaklist2[use,3]-1+uplim_tol) 		# maxup UB
			bounds[4,1]<-(peaklist2[use,4]-1-downlim_tol)		# maxdown LB
			bounds[4,2]<-(peaklist2[use,4]+shift-1+downlim_tol) # maxdown UB
			.Call("search_kdtree_homol", 
				peaklist2, 	# rows: c(m/z,RT,UB,LB)
				peakTree,  	# peaks - search tree
				bounds,
				marked,		# only write to new_found those of marked[,1]!=(i-1) to omit last sweep`s new_found
				i,
				new_found_ceiling,	
				1, 			# clean new_found?
				PACKAGE="nontarget"
			)
			if(new_found_ceiling[1]!=0){ceiled<-TRUE}
		}
		# lower area sweep ###################################################
		bounds[1,1]<-(peaklist2[use,1]-maxmz)
		bounds[1,2]<-(peaklist2[use,1]-minmz) 
		bounds[2,1]<-(peaklist2[use,2]-maxrt)
		bounds[2,2]<-(peaklist2[use,2]+minrt)
		bounds[3,1]<-((peaklist2[use,3])-uplim_tol)			# maxup LB
		bounds[3,2]<-((peaklist2[use,3])+shift+uplim_tol) 	# maxup UB
		bounds[4,1]<-((peaklist2[use,4])-shift-downlim_tol)	# maxdown LB
		bounds[4,2]<-(peaklist2[use,4]+downlim_tol) 		# maxdown UB
		.Call("search_kdtree_homol", 
			peaklist2, 	# rows: c(m/z,RT,UB,LB)
			peakTree,  	# peaks - search tree
			bounds,
			marked,		# only write to new_found those of marked[,1]!=(i-1) to omit last sweep`s new_found
			i,
			new_found,	
			0, 			# clean new_found?
			PACKAGE="nontarget"
		)	
		# -> mass defect rounding hits floor? ################################
		if((mass_def[i]-(maxmove*maxmz))<=-.5){
			bounds[3,1]<-((peaklist2[use,3])+1-uplim_tol)			# maxup LB
			bounds[3,2]<-((peaklist2[use,3])+shift+1+uplim_tol) 	# maxup UB
			bounds[4,1]<-((peaklist2[use,4])-shift+1-downlim_tol)	# maxdown LB
			bounds[4,2]<-(peaklist2[use,4]+1+downlim_tol) 			# maxdown UB
			.Call("search_kdtree_homol", 
				peaklist2, 	# rows: c(m/z,RT,UB,LB)
				peakTree,  	# peaks - search tree
				bounds,
				marked,		# only write to new_found those of marked[,1]!=(i-1) to omit last sweep`s new_found
				i,
				new_found_floor,	
				1, 			# clean new_found?
				PACKAGE="nontarget"
			)
			if(new_found_floor[1]!=0){floored<-TRUE}
		}
		######################################################################
		# delete or update old points ########################################
		if(length(dist_ID)>0){
			stay<-(marked[dist_ID,1]==i)
			dist_ID<-dist_ID[stay]	
			dist_dist<-dist_dist[stay]		
			dist_dist<-(dist_dist+mz_last-peaklist2[use,1])
		}
		mz_last<-peaklist2[use,1] # for next=i+1 update
		######################################################################
		# add new points #####################################################
		if(new_found[1]!=0){
			many<-(sum(new_found!=0))
			a<-(a+many)
			dist_ID<-append(dist_ID,new_found[1:many])
			dist_dist<-append(dist_dist,
				peaklist2[new_found[1:many],1]-peaklist2[use,1]	
			)
		}
		if(ceiled){
			many<-sum(new_found_ceiling!=0)
			dist_ID<-append(dist_ID,new_found_ceiling[1:many])
			dist_dist<-append(dist_dist,
				peaklist2[new_found_ceiling[1:many],1]-peaklist2[use,1]	
			)	
			new_found_ceiling[]<-0;
			ceiled<-FALSE;
			a<-(a+many)
		}
		if(floored){
			many<-sum(new_found_floor!=0)
			dist_ID<-append(dist_ID,new_found_floor[1:many])
			dist_dist<-append(dist_dist,
				peaklist2[new_found_floor[1:many],1]-peaklist2[use,1]	
			)	
			new_found_floor[]<-0;
			floored<-FALSE;		
			a<-(a+many)
		}
		if(length(dist_ID)>1){
			##################################################################
			# resort to new centre point i ###################################
			ord<-order(abs(dist_dist)) # decreasing=FALSE!
			dist_ID<-dist_ID[ord]
			dist_dist<-dist_dist[ord]
			##################################################################
			# find triplets ##################################################
			.Call("homol_triplet", 
				peaklist3,
				dist_ID,
				dist_dist,
				tupels,
				tupeldo,
				peaklist4,
				use, # = current point
				max_delmz,
				rttol,
				PACKAGE="nontarget"
			);
			if(tupeldo>vec_size){stop("\n tupels out of bounds. increase vec_size")}
			##################################################################
		}	
		######################################################################
	}
	if(inter){close(pBar)}
	tupeldo<-(tupeldo-1)
	if(tupeldo==0){stop("no series detected")}
	tupels<-tupels[1:tupeldo,,drop=FALSE]
	tupels<-tupels[tupels[,1]!=0,,drop=FALSE] # backup - if tupeldo fails 
	tupels<-tupels[order(tupels[,4]),,drop=FALSE]
	if(plotit==1){
		######################################################################	
		# path in m/z vs. mass defect ########################################
		plot.new()
		plot.window(xlim=c(min(peaklist2[,1]),max(peaklist2[,1])),ylim=c(min(mass_def),max(mass_def)))
		points(
			peaklist2[,1],mass_def,
			pch=19,cex=.1,col="lightgrey"
		)
		title(xlab="m/z",ylab="Mass defect",main="Nearest neighbour path")
		box();axis(1);axis(2)
		for(j in 2:length(along)){
			lines(
				c(peaklist2[along[j],1],peaklist2[along[j-1],1]),c(mass_def[along[j]],mass_def[along[j-1]]),
				,col="black",lwd=.5)
		}
		Sys.sleep(5)
		######################################################################	
		plot.new()
		plot.window(xlim=c(min(peaklist2[,1]),max(peaklist2[,1])),ylim=c(min(peaklist2[,2]),max(peaklist2[,2])))
		points(
			peaklist2[,1],peaklist2[,2],
			pch=19,cex=.1,col="lightgrey"
		)
		title(xlab="m/z",ylab="RT",main="Nearest neighbour path")
		box();axis(1);axis(2)
		for(j in 2:length(along)){
			lines(
				c(peaklist2[along[j],1],peaklist2[along[j-1],1]),c(peaklist2[along[j],2],peaklist2[along[j-1],2]),
				,col="black",lwd=.5)
		}
		Sys.sleep(5)
		######################################################################	
	}
	cat("   done.");
	##########################################################################
	# (5) Apply mzfilter #####################################################
	if( mzfilter[1]!="FALSE" ){
		cat("\n(5) Apply mzfilter ... ");	
		keep<-rep(FALSE,tupeldo)
		for(k in 1:length(mzfilter)){
			keep[	
				tupels[,4]<=mzfilter[k] &
				tupels[,5]>=mzfilter[k]	
			]<-TRUE
		}
		if(!any(keep)){stop("No homologues detected after mzfilter application")}
		tupels<-tupels[keep,,drop=FALSE]
		cat("done.");
	}else{
		cat("\n(5) Skip mzfilter. ");	
	}
	##########################################################################
	# (6) Combine triplets to n-tupel - possibly check for smoothness ########
	cat("\n(6) Combine n-tupels, n(rejects):");	
	HS<-list();
	HS_length<-3;
	found<-0;
	tupels<-cbind(tupels[,],seq(1,length(tupels[,1]),1)) # dummy column
	while(any(tupels[,1]!=0)){
	    cat(paste(" ",HS_length," - ",length(tupels[,1]),sep=""));
		keeper<-rep(0,length(tupels[,1]))
		if(length(tupels[,1])==1){break}
		if(length(tupels[,1])==0){stop("\n debug me, issue #1")}
		if(any(tupels[,1]<1)){stop("\n debug me, issue #3")}
		merged_tupels<-.Call("combine_tuple", 
			tupels[,1:HS_length],
			tupels[,(HS_length+1):(HS_length+4)],
			keeper,
			mat_size
		)	
		merged_tupels<-merged_tupels[merged_tupels[,1]!=0,,drop=FALSE]	
		if(length(merged_tupels[,1])==0){
			HS[[HS_length]]<-tupels[keeper==0,,drop=FALSE]
			if(HS_length>=minlength){
				found<-c(found+length(HS[[HS_length]][,1]))
			}
			break;
		}
		lang<-length(merged_tupels[1,])
		keeper_2<-rep(1,length(merged_tupels[,1]))
		if(R2!=FALSE){
			reject<-0;
			for(i in 1:length(merged_tupels[,1])){
				x<-peaklist3[
					merged_tupels[i,1:(HS_length+1)],1
				]
				if(any(duplicated(x)) & (HS_length<4)){next} # bypass fitting issues
				y<-peaklist3[
					merged_tupels[i,1:(HS_length+1)],2
				]
				if(any(duplicated(y)) & (HS_length<4)){next} # bypass fitting issues
				b <- predict(stats::smooth.spline(x,y,cv=FALSE,spar=.45))
				SS_tot<-sum((y-mean(y))^2);
				SS_res<-sum((y-b$y)^2);
				R2_i<-(1-(SS_res/SS_tot));
				if(R2_i<R2){
					keeper[merged_tupels[i,lang]]<-(keeper[merged_tupels[i,lang]]-1)
					keeper[merged_tupels[i,(lang-1)]]<-(keeper[merged_tupels[i,(lang-1)]]-1)
					keeper_2[i]<-0
					reject<-(reject+1)
				}
				if(plotit==2 & (HS_length>=minlength)){
					plot(y,x,pch=19,main=R2_i,xlab="m/z",ylab="RT");
					points(b$y,b$x,type="l",col="red");	
					Sys.sleep(1.4)
				}
			}
			cat(paste("(",reject,")",sep=""));
		}else{
			cat(paste("(0)",sep=""));			
		}
		HS[[HS_length]]<-tupels[keeper==0,,drop=FALSE]
		cat(", ",sum(keeper!=0),sep="")
		if(HS_length>=minlength){
			found<-c(found+length(HS[[HS_length]][,1]))
		}
		HS_length<-(HS_length+1)
		merged_tupels<-merged_tupels[keeper_2==1,,drop=FALSE]
		merged_tupels<-merged_tupels[order(merged_tupels[,(HS_length+1)]),,drop=FALSE]
		if( length(merged_tupels[,1])==1 ){ # single-rowed?
			HS[[HS_length]]<-merged_tupels;
			if(HS_length>=minlength){
				found<-c(found+length(HS[[HS_length]][,1]))
			}			
			break;
		}
		tupels<-merged_tupels;	
	}
	if(any(deb==1)){return(HS)}
	if(length(HS)<minlength){
		stop("\n No homologue series detected. Check parameters (e.g., minlength) and units?");
	}
	cat(" - done.");
	##########################################################################	
	# (7) remove nested HS ###################################################
	if(!any(deb==2)){
		cat("\n(7) Remove nested HS:  \n");
		reject<-0;
		# membership of peaks in HS #########################################
		peakHS<-list(0)
		for(i in 1:length(peaklist[,1])){
			peakHS[[i]]<-matrix(ncol=2,nrow=0)
		}
		if(length(HS)>minlength){ # anything to do?		
			for(i in minlength:length(HS)){		
				if(length(HS[[i]])>0){
					HS[[i]][,length(HS[[i]][1,])]<-1 # use last column as marker for keeping a HS
					for(a in 1:length(HS[[i]][,1])){
						for(b in 1:i){
							peakHS[[ HS[[i]][a,b] ]]<-(
								rbind(peakHS[[ HS[[i]][a,b] ]],c(i,a))
							)
						}
					}
				}
			}
		}
		# check peak memberships of HS for subsets ##########################
		if(inter){pBar <- txtProgressBar( min = 0, max = length(peaklist2[,1]), style = 3 )}
		for(i in 1:length(peakHS)){
			if(inter){setTxtProgressBar(pBar,i,title = NULL, label = NULL)}
			if(length(peakHS[[i]])>2){
				ns<-unique(peakHS[[i]][,1])
				if(length(ns)>1){
					for(a in 1:(length(peakHS[[i]][,1])-1)){
						n<-peakHS[[i]][a,1]
						for(b in (a+1):length(peakHS[[i]][,1])){
							m<-peakHS[[i]][b,1]
							if(n<m){
								if( m>=((n*2)-1) ){ # series nesting requirement
									this<-HS[[n]][peakHS[[i]][a,2],1:n]
									that<-HS[[m]][peakHS[[i]][b,2],1:m]
									if(all(!is.na(match(this,that)))){
										k<-peakHS[[i]][a,2]
										HS[[n]][k,length(HS[[n]][1,])]<-0
										reject<-(reject+1)
										#stop()	
									}
									
								}
							}
						}
					}	
				}
			}
		}
		if(inter){close(pBar)}
		# remove unmarked HS ################################################
		found<-0
		if(length(HS)>minlength){ # anything to do?		
			for(i in minlength:length(HS)){		
				if(length(HS[[i]])>0){
					HS[[i]]<-(HS[[i]][HS[[i]][,length(HS[[i]][1,])]==1,,drop=FALSE])
					found<-(found+length(HS[[i]][,1]))
				}
			}
		}
		#####################################################################
		cat(paste(reject," cases - done.",sep=""));
	}else{
		cat("\n(7) Removal of nested HS skipped");	
	}
	##########################################################################	
	# (8) cluster similar HS #################################################
	if(any(deb==3)){
		cat("\n(8) Grouping of superjacent HS:  \n");
		# HS_IDs: table to translate HS list to a continous ID and back ######
		# use 2nd last column in HS to place a unique ID #####################
		# peakHS: membership of peaks in HS ##################################
		if(length(HS)<minlength){stop("length(HS)<minlength")}
		peakHS<-list(0)
		for(i in 1:length(peaklist[,1])){
			peakHS[[i]]<-numeric(0)
		}
		number<-1
		for(i in minlength:length(HS)){		
			if(length(HS[[i]])>0){
				HS[[i]][,i+5]<-seq(number,(number+length(HS[[i]][,1])-1),1)
				number<-(number+length(HS[[i]][,1]))
				for(a in 1:length(HS[[i]][,1])){
					for(b in 1:i){
						peakHS[[ HS[[i]][a,b] ]]<-(
							c(
								peakHS[[ HS[[i]][a,b] ]],
								HS[[i]][a,i+5]
							)
						)
					}
				}
			}	
		}
		HS_IDs<-matrix(nrow=(number-1),ncol=4,0)
		number<-1
		for(i in minlength:length(HS)){		
			if(length(HS[[i]])>0){
				HS_IDs[number:(number+length(HS[[i]][,1])-1),1]<-HS[[i]][,i+5];
				HS_IDs[number:(number+length(HS[[i]][,1])-1),2]<-i;
				HS_IDs[number:(number+length(HS[[i]][,1])-1),3]<-seq(1,length(HS[[i]][,1]),1)
				number<-(number+length(HS[[i]][,1]))
			}	
		}
		if(any(duplicated(HS_IDs))){stop("debug me on issue #6: duplicated IDs in similarity grouping of HS!")}
		# get matrix with unique HSpairs #################################### 				
		at<-1;
		leng<-5000;
		from<-rep(0,leng);
		to<-rep(0,leng);
		for(i in 1:length(peakHS)){
			if(length(peakHS[[i]])>1){
				these<-peakHS[[i]]
				these<-these[order(these)]
				for(a in 1:(length(these)-1)){
					for(b in (a+1):length(these)){
						from[at]<-these[a]
						to[at]<-these[b]
						at<-(at+1)
						if(at>leng){ # extend the vector
							from<-c(from,rep(0,leng))
							to<-c(to,rep(0,leng))
							leng<-(leng+leng)
						}
					}
				}
			}
		}
		from<-from[from!=0]
		to<-to[to!=0]
		HSpairs<-unique(cbind(from,to))
		# derive similarities ###############################################
		sim<-rep(0,length(HSpairs[,1]))
		if(inter){pBar <- txtProgressBar( min = 0, max = length(HSpairs[,1]), style = 3 )}
		for(i in 1:length(HSpairs[,1])){
			if(inter){setTxtProgressBar(pBar,i,title = NULL, label = NULL)}
			# get model for peak set a
			n_a<-HS_IDs[HSpairs[i,1],2]
			m_a<-HS_IDs[HSpairs[i,1],3]
			peaks_a<-HS[[n_a]][m_a,1:n_a]
			x_a<-peaklist[peaks_a,1]
			if(any(duplicated(x_a)) ){next} # bypass fitting issues
			y_a<-peaklist[peaks_a,3]
			if(any(duplicated(y_a))){next} # bypass fitting issues			
			model_a <- smooth.spline(x_a,y_a,cv=FALSE,spar=.45)
			# get model for peak set b				
			n_b<-HS_IDs[HSpairs[i,2],2]
			m_b<-HS_IDs[HSpairs[i,2],3]
			peaks_b<-HS[[n_b]][m_b,1:n_b]
			x_b<-peaklist[peaks_b,1]
			if(any(duplicated(x_b)) ){next} # bypass fitting issues
			y_b<-peaklist[peaks_b,3]
			if(any(duplicated(y_b))){next}  # bypass fitting issues			
			model_b <- smooth.spline(x_b,y_b,cv=FALSE,spar=.45)		
			# peaks a done by model_b			
			predict_a<-predict(model_b,x_a)			
			SS_tot_a<-sum((y_a-mean(y_a))^2);			
			SS_res_a<-((y_a-predict_a$y)^2);	
			SS_res_a[!is.na(match(peaks_a,peaks_b))]<-0 # replace full peak matches
			SS_res_a<-sum(SS_res_a);
			R2_a<-(1-(SS_res_a/SS_tot_a));
			# peaks b done by model_a			
			predict_b<-predict(model_a,x_b)			
			SS_tot_b<-sum((y_b-mean(y_b))^2);			
			SS_res_b<-((y_b-predict_b$y)^2);		
			SS_res_b[!is.na(match(peaks_b,peaks_a))]<-0	# replace full peak matches
			SS_res_b<-sum(SS_res_b);
			R2_b<-(1-(SS_res_b/SS_tot_b));
			# store maximum similarity
			sim[i]<-max(c(R2_a,R2_b))
			# plot it!
			if(plotit==3){	
				plot.new()
				plot.window(xlim=c(min(c(x_a,x_b)),max(c(x_a,x_b))),ylim=c(min(c(y_a,y_b)),max(c(y_a,y_b))))
				box();axis(1);axis(2);
				points(x_a,y_a,pch=19,cex=.5,col="red")
				points(x_a,y_a,type="l",col="red")
				points(x_b,y_b,pch=19,cex=.5,col="blue")
				points(x_b,y_b,type="l",col="blue")	
				title(main=paste(round(R2_a,2)," vs. ",round(R2_b,2)))
				Sys.sleep(.3)
			}
		}
		if(inter){close(pBar)}
		HSpairs<-HSpairs[sim>=.9,,drop=FALSE]
		if(length(HSpairs[,1])>0){ # any pairs formed? ...
			from<-HSpairs[,1]
			to<-HSpairs[,2]
			# add non-paired HS to have them assigned as a group, too
			gotthose<-unique(c(from,to))
			singled<-HS_IDs[-gotthose,1]
			use<-c(rep(1,length(from)),rep(2,length(to)),rep(1,length(singled)))
			from2<-c(from,to,singled)
			to2<-c(to,from,singled)
			use<-use[order(from2,decreasing=FALSE)]
			to2<-to2[order(from2,decreasing=FALSE)]
			from2<-from2[order(from2,decreasing=FALSE)]
			groups<-.Call("metagroup",
				as.integer(from2),
				as.integer(to2),
				PACKAGE="nontarget" 
			)
			to2<-to2[use==1]
			from2<-from2[use==1]
			groups<-groups[use==1]
			HS_IDs[to2,4]<-groups
			HS_IDs[from2,4]<-groups
			if(any(HS_IDs[,4]==0)){stop("debug me on issue #7: ungrouped HS formed.")}
		}else{ # otherwise, insert a consecutive grouping ID
			# insert a consecutive grouping ID
			HS_IDs[,4]<-seq(1,length(HS_IDs[,4]),1);
			groups<-length(HS_IDs[,4]);
		}
		cat(paste(max(groups)," HS cluster formed - done.",sep=""));
	}else{
		# insert a consecutive grouping ID
		number<-1
		for(i in minlength:length(HS)){		
			if(length(HS[[i]])>0){
				HS[[i]][,i+5]<-seq(number,(number+length(HS[[i]][,1])-1),1)
				number<-(number+length(HS[[i]][,1]))
			}	
		}
		HS_IDs<-matrix(nrow=(number-1),ncol=4,0)
		number<-1
		for(i in minlength:length(HS)){		
			if(length(HS[[i]])>0){
				HS_IDs[number:(number+length(HS[[i]][,1])-1),1]<-HS[[i]][,i+5];
				HS_IDs[number:(number+length(HS[[i]][,1])-1),2]<-i;
				HS_IDs[number:(number+length(HS[[i]][,1])-1),3]<-seq(1,length(HS[[i]][,1]),1)
				number<-(number+length(HS[[i]][,1]))
			}	
		}
		if(any(duplicated(HS_IDs))){stop("debug me on issue #6: duplicated IDs in similarity grouping of HS!")}
		HS_IDs[,4]<-seq(1,length(HS_IDs[,4]),1);
	}
	##########################################################################	
	# (9) Generate data output ###############################################
	cat(paste("\n(9) Parse output for ",max(HS_IDs[,4])," homologue series and ",found," cluster ... ",sep=""));
	# generate HS group lists ################################################
	HS_groups<-list(0)
	for(i in 1:length(HS_IDs[,4])){
		HS_groups[[HS_IDs[i,4]]]<-numeric(0)
	}
	for(i in 1:length(HS_IDs[,4])){
		HS_groups[[HS_IDs[i,4]]]<-c(HS_groups[[HS_IDs[i,4]]],HS_IDs[i,1])
	}
	if(plotit==4){
		for(i in 1:length(HS_groups)){
			those<-HS_IDs[HS_groups[[i]],,drop=FALSE]
			peaks<-c()
			delmz<-c()
			for(j in 1:length(those[,1])){
				peaks<-c(peaks,
					HS[[those[j,2]]][those[j,3],1:those[j,2]]
				)
				delmz<-c(delmz,
					mean(HS[[those[j,2]]][those[j,3],(those[j,2]+1):(those[j,2]+2)])
				)
			}
			min_mz<-min(peaklist[peaks,1])
			max_mz<-max(peaklist[peaks,1])
			min_RT<-min(peaklist[peaks,3])
			max_RT<-max(peaklist[peaks,3])
			plot.new();
			plot.window(xlim=c(min_mz,max_mz),ylim=c(min_RT,max_RT));			
			title(xlab="m/z",ylab="RT",main=paste(length(those[,1]),"HS"))
			box();axis(1);axis(2);
			this<-round(delmz,digits=2);
			that<-levels(as.factor(this));
			colo<-rainbow(length(that))
			for(j in 1:length(those[,1])){
				peaks<-HS[[those[j,2]]][those[j,3],1:those[j,2]]
				points(peaklist[peaks,1],peaklist[peaks,3],
					col=colo[that==this[j]],pch=19,cex=0.5)
				points(peaklist[peaks,1],peaklist[peaks,3],
					col=colo[that==this[j]],type="l")			
			}
			Sys.sleep(.8)
		}
	}
    # generate peaklist with links & #########################################
    # generate component list with relevant m/z & RT increments ##############
    group1<-c();  # HS ID
    group2<-c();  # peak IDs
    group3<-c();  # mzincr
    group4<-c();  # retincr
    group5<-c();  # retmin
    group6<-c();  # retmax
	group7<-c();  # retdel
	group8<-c();  # HS cluster ID
	len<-length(peaklist[,1])
    getit1<-rep("0",len);      	# (1) HS ID
    getit2<-rep("0",len);      	# (2) level
    getit3<-rep("0",len);   	# (3) to which peak
    getit4<-rep("none",len);   	# (4) mean dmass
    getit5<-rep("none",len);    # (5) mean RT
    getit6<-rep("0",len);      	# (6) HS cluster ID	
	peakID<-list(0);
	atgroup<-1;
	for(a in minlength:length(HS)){
		if(length(HS[[a]])>0){ 		
			for(b in 1:length(HS[[a]][,1])){
				peakID[[atgroup]]<-numeric(0)
				listit<-list(0)
				atthat<-""
				meanmz<-mean(HS[[a]][b,((a+1):(a+2))]);
				meanRT<-mean(HS[[a]][b,((a+3):(a+4))]);
				minRT<-min(peaklist[HS[[a]][b,1:a],3]);
				maxRT<-max(peaklist[HS[[a]][b,1:a],3]);
				difRT<-(maxRT-minRT);
				for(k in 1:a){
					getit1[HS[[a]][b,k]]<-paste(getit1[HS[[a]][b,k]],"/",as.character(atgroup),sep="")
					getit2[HS[[a]][b,k]]<-paste(getit2[HS[[a]][b,k]],"/",as.character(k),sep="")
					if(k<a){
						getit3[HS[[a]][b,k]]<-paste(getit3[HS[[a]][b,k]],"/",as.character(HS[[a]][b,(k+1)]),sep="")
					}
					getit4[HS[[a]][b,k]]<-paste(getit4[HS[[a]][b,k]],"/",as.character(round(meanmz,digits=4)),sep="")
					getit5[HS[[a]][b,k]]<-paste(getit5[HS[[a]][b,k]],"/",as.character(round(meanRT,digits=4)),sep="")	
					getit6[HS[[a]][b,k]]<-paste(getit6[HS[[a]][b,k]],"/",as.character(HS_IDs[HS[[a]][b,a+5],4]),sep="")
					listit[[k]]<-HS[[a]][b,k];
					atthat<-paste(atthat,",",as.character(HS[[a]][b,k]),sep="");
				}
				peakID[[atgroup]]<-listit;
				group1<-c(group1,atgroup); 
				group2<-c(group2,substr(atthat,2,nchar(atthat)));
				group3<-c(group3,meanmz);
				group4<-c(group4,meanRT);
				group5<-c(group5,minRT);				
				group6<-c(group6,maxRT);				
				group7<-c(group7,difRT);
				group8<-c(group8,HS[[a]][b,a+5]);
				atgroup<-(atgroup+1);
			}
		}
	}
	for(i in 1:length(getit1)){
		if(getit1[i]!="0"){
			getit1[i]<-substr(getit1[i],3,nchar(getit1[i]))
			getit2[i]<-substr(getit2[i],3,nchar(getit2[i]))
			getit4[i]<-substr(getit4[i],6,nchar(getit4[i]))
			getit5[i]<-substr(getit5[i],6,nchar(getit5[i]))
			getit6[i]<-substr(getit6[i],3,nchar(getit6[i]))			
		}
		if(getit3[i]!="0"){
			getit3[i]<-substr(getit3[i],3,nchar(getit3[i]))
		}
	}
	grouped_samples<-data.frame(peaklist,seq(1,len,1),getit1,getit2,getit3,getit4,getit5,getit6,stringsAsFactors=FALSE);
	names(grouped_samples)<-c("mz","intensity","RT","peak ID","HS IDs","series level","to ID","m/z increment","RT increment","HS cluster");
	grouping<-data.frame(group1,group2,group3,group4,group5,group6,group6-group5,group8,stringsAsFactors=FALSE);
	names(grouping)<-c("HS IDs","peak IDs","m/z increment","RT increment","min. RT in series","max. RT in series","max.-min. RT","HS cluster");
	# store parameters used ##################################################
	parameters<-list(0)
	parameters[[1]]<-elements;
	parameters[[2]]<-use_C
	parameters[[3]]<-minmz
	parameters[[4]]<-maxmz
	parameters[[5]]<-minrt
	parameters[[6]]<-maxrt
	parameters[[7]]<-ppm	
	parameters[[8]]<-mztol
	parameters[[9]]<-rttol
	parameters[[10]]<-minlength
	parameters[[11]]<-vec_size
	parameters[[12]]<-spar
	parameters[[13]]<-R2
	parameters[[14]]<-plotit	
	names(parameters)<-c("elements","use_C","minmz","maxmz","minrt","maxrt","ppm","mztol","rttol","minlength","vec_size","spar","R2","plotit")
	cat("done.\n");
	##########################################################################		
	# store HS memberships per peak ##########################################
	homol<-list(grouped_samples,parameters,grouping,mzfilter,peakID,HS_groups);
	names(homol)<-c("Peaks in homologue series","Parameters","Homologue Series","m/z-Filter used","Peaks per level","HS cluster")
	##########################################################################
	return(homol); ###########################################################
	##########################################################################
	
}















