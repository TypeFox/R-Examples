#### This will be the elec.strat package 
#### Created by Mike Higgins and Hua Yang
#### Obtains p-values for audits using
#### stratified random samples.
####  
#### Has options to give improved
#### sample sizes within each strata.
####
#### Package ought to include CA house races!!!!!
#### and Minnesota 2006 senate race!!!!!
####
####Uses elec package.

#makeStratObj requires as input a  
#valid elec.data object Z with a 
#1) Z$V$CID column OR
#2) with a strat.col name given OR
#3) A CID vector of length(nrow(votes))
#   specifying the strata ID of the corresponding
#   batches.

#When both strat.col and CID are included
#priority is given to strat.col, then to CID

#makeStratObj creates a strat.elec.data object.

#The CID column will need for each batch 
#a label identifing to which stratum the batch belongs.

#If audit information is availible, 
#and not included in Z$audit,
#feed it in using auditTable argument.
#auditTable will require a data.frame of dimensions
#(no. of strata)x(2).
#Each row in the auditInfo data.frame is comprised of:
#1) A unique CID label
#2) The number of batches audited within the
#	stratum corresponding to that CID label

#The audit column is the 
#second column of the auditTable matrix

#If no audit information is available, 
#the audit column defaults to a column of zeroes.

#CIDnum assigns to each column 
#a number ranging from  1 to the no. of strata

#n denotes the number of batches within a strata.

#The Z$strat data.frame includes the columns (in order)
#CID
#CIDnum
#n
#audit

#Will also create 
#Z$CID.col: "CID"
#Z$CIDnum: The CID column with CIDnum instead of the CID label
#Z$V$e.max: maximumMarginBound(Z)

#If Z$audit[Z$PID.col] != NULL, create
#Z$audit$e.max: maximumMarginBound(Z) for PIDs in Z$audit

makeStratObj <- function(Z, strat.col = NULL, CID = NULL, auditTable = NULL){
	
	#Check if Z is elec.data object
	if(!is.elec.data(Z)){
		stop("makeStratObj requires an elec.data object.")
	}
	
	#Check if Z$V$CID column is null
	if(is.null(Z$V$CID) && is.null(strat.col) && is.null(CID)){
		stop("makeStratObj needs argument that identifies name of strata column.")
	}
	
	#If one of these is not null
	else{ 
		
		#If strat.col not null
		if(!is.null(strat.col)){
			
			#Check if votes column has strat.col
			if ( !(strat.col %in% names(Z$V) ) ) {
				stop("Name of column'", 
					 strat.col,"' not in vote data object Z$V." )
			}
			
			#Create CID column
			CID = Z$V[[strat.col]]
			Z$V$CID <- CID
			
			#If information about audit included in Z$audit, create auditTable from that.
			if(!is.null(Z$audit)){
				
				#Check if PID.col is in the audit data frame
				if(!(strat.col %in% names(Z$audit))){
					
					warning("Name of column'", 
						 strat.col,"' not in audit data object Z$audit. Will not use Z$audit information to create Z$strat" )
				}
				else{
					
					#Create auditTable
					CIDcol <- unique(Z$audit[[strat.col]])
					AUDcol <- sapply(CIDcol,function(stratum){
									 length(which(Z$audit[,strat.col] == stratum))})
					
					auditTable <- data.frame(CIDcol,AUDcol)
				}
			}
		}
		
		#If CID is not null
		else if(!is.null(CID)){
			
			#Check if CID is of right the length
			if(length(CID) != nrow(Z$V)){
			   stop(paste("Length of CID vector not equal",
					  "to no. of rows in vote data object."))
			}
			
			Z$V$CID <- CID
		}	
	}	
	
	#Check valid audit table	
	if(!is.null(auditTable)){
		
		#Check if dimensions of audit table are correct
		if(is.null(nrow(auditTable)) || 
		   is.null(ncol(auditTable)) ||
		   is.na(nrow(auditTable)) ||
		   is.na(ncol(auditTable)) ||
		   nrow(auditTable) != length(unique(Z$V$CID)) ||
		   ncol(auditTable) != 2){
			stop(paste("Audit table is not of correct dimensions.",  
				 "Audit table needs to be a data.frame containing two columns of length = no. of stratum."))
		}
		
		#Check if audit column is valid
		if(!is.numeric(auditTable[,2])){
			stop("Audit column not numeric.")
		}
		
		#Check if CID column matches CIDs in Z
		#Check unmatched auditCIDs or unmatched Z$V$CID
		if(sum(is.na(match(auditTable[,1],unique(Z$V$CID)))) != 0 ||
		   sum(is.na(match(unique(Z$V$CID),auditTable[,1]))) != 0){
			stop(paste("Stratum labels in audit table do not", 
					   "match stratum labels in Z$V"))
		}
	}
	
	#Get strata.  Will be CID.
	strata <- unique(Z$V$CID)
	
	#Get Num_Pct.  Will be n.
	Num_Pct <-  sapply(strata,function(stratum){
					   length(which(Z$V$CID == stratum))})
	
	#Get audit column.  Will be audit.
	if(!is.null(auditTable)){
		auditCol <-  sapply(strata,function(stratum){
							auditTable[which(auditTable[,1] == stratum),2]})
	}
	else{
		auditCol = rep(0,length(strata))
	}
	
	#Get CIDn.  Will become CIDnum.
	CIDn=rank(-Num_Pct, ties.method = "first")
	
	#For Z$V$CIDnum
	ID <- rep.int(0,NROW(Z$V))
	for (i in 1:length(strata)){
		ID[which(Z$V$CID == strata[i])] = CIDn[i]
	}
	Z$V$CIDnum <- ID
	
	Z$strat <- data.frame(CID = strata, 
						  CIDnum=CIDn,n=Num_Pct,audit = auditCol)
	
	Z$V <- Z$V[order(Z$V$CIDnum), ]
	Z$strat <- Z$strat[order(Z$strat$CIDnum), ]
	Z$CID.col <- "CID"
	
	#Make Z$V$e.max column
	Z$V$e.max <- maximumMarginBound(Z)
	
	#If Z$audit[Z$PID.col] != NULL, create
	#Z$audit$e.max: maximumMarginBound(Z) for PIDs in Z$audit
	if( !is.null(Z$audit[[Z$PID.col]]) ){
		
		#aud.ebs holds error bounds for audit.
		aud.ebs <- length(Z$audit[[Z$PID.col]])
		
		#For each PID in audit, find PID in Z$V, and get the error bound
		i = 1
		for(PIDs in Z$audit[[Z$PID.col]]){
			
			#Is the PID in Z$V?
			if(!(PIDs %in% (Z$V[[Z$PID.col]]) ) ){
				stop("Precinct ID in audit not found in votes data.frame")
			}
			
			#Is the PID in Z$V multiple times?
			if(length(which(Z$V[[Z$PID.col]] == PIDs)) > 1){
				stop("Precinct ID not unique in votes data.frame")	
			}
			
			#Get the row in votes data.frame corresponding to PIDs, and get the e.max value
			#Store in aud.ebs
			myPID <- Z$V[which(Z$V[[Z$PID.col]] == PIDs),]
			aud.ebs[i] <- myPID$e.max
			
			i = i + 1
		}
		
		#Make Z$audit$e.max column
		Z$audit$e.max <- aud.ebs
	
	}
	
	class(Z) = c(class(Z),"strat.elec.data")
	
	Z
}

# strat.elec.data creates a strat.elec.data object:
# an elec.data object with additional entries to allow for
# simple incorporation with the rest of the elec.strat package

#Calls elec.data function from elec package,
#see comments there for implementation

#Two ways to specify stratification:
#1) Give to strat.col argument the name 
#   of the column that contains strata information
#2) Give to CID argument a vector of length(nrow(votes))
#   specifying the strata ID of the corresponding
#   batches.
#Priority given to 1) then to 2), in case both arguments
#are filled.

#If audit is not null, and strat.col given,
#strat.elec.data will try to find strat.col in the audit dataframe
#to create an audit table.

#If strat.col not included, CID included, 
#and audit information is included,
#strat.elec.data will first try to find 
#the PID.col in the audit column,
#and use this information to generate the auditTable.
#If this column is not in the audit data frame,
#then auditTable will need to be given,
#otherwise Z$strat$audit defaults to a zero vector.

#auditTable will require a data.frame of dimensions
#(no. of strata)x(2).
#Each row in the auditInfo data.frame is comprised of:
#1) A unique CID label
#2) The number of batches audited within the
#	stratum corresponding to that CID label


#Creates Z$V$e.max column, filled by maximumMarginBound(Z),

#If no strata are specified, only an elec.data object is created.

strat.elec.data <- function(V, C.names=names(V)[2:length(V)], f = 1, 
audit=NULL, pool=TRUE, tot.votes.col="tot.votes", PID.col="PID",
strat.col = NULL, CID = NULL, auditTable = NULL){
	
	#Create the elec.data object
	Z = elec.data(V, C.names, f, audit, pool, tot.votes.col, PID.col)
	
	if(is.null(strat.col) && is.null(CID)){	
		warning(paste("No information identifying strata column included.",  
					  "Only creating elec.data object."))	
	}
	
	#Else, check if strat.col is null
	else{ 
		if(!is.null(strat.col)){
			
			#Check if votes column has strat.col
			if ( !(strat.col %in% names(V) ) ) {
				stop("Name of column'", 
					 strat.col,"' not in vote data object." )
			}
			CID = V[,strat.col]
			Z$V$CID <- CID
			if(!is.null(audit)){
				
				#Check if audit column has strat.col
				if ( !(strat.col %in% names(audit) ) ) {
					stop("Name of column'", 
						 strat.col,"' not in audit data object." )
				}
				
				#Create auditTable
				CIDcol <- unique(audit[,strat.col])
				AUDcol <- sapply(CIDcol,function(stratum){
						  length(which(audit[,strat.col] == stratum))})
				auditTable <- data.frame(CIDcol,AUDcol)
			}
		}
		
		#Else, must have non-null CID entry.
		else{
			
			#Check if CID is of right length
			if(length(CID) != nrow(V)){
				stop(paste("Length of CID vector not equal",
						   "to no. of rows in vote data object."))
			}			
			if(!is.null(audit)){
				
				#Check if PID.col is in the audit data frame
				if(Z$PID.col %in% names(audit)){
					
					#Initialize audit table.
					CIDcol <- unique(CID)
					AUDcol <- rep(0,length(CIDcol))
					
					for(PIDs in Z$audit[[Z$PID.col]]){
						
						#Is the PID in Z$V?
						if(!(PIDs %in% (Z$V[[Z$PID.col]]) ) ){
							stop("Precinct ID in audit not found in votes data.frame")
						}
						
						#Is the PID in Z$V multiple times?
						if(length(which(Z$V[[Z$PID.col]] == PIDs)) > 1){
							stop("Precinct ID not unique in votes data.frame")	
						}
						
						#Get the CID corresponding to PID
						myCID <- CID[which(Z$V[[Z$PID.col]] == PIDs)]
						
						#Get the index of that CID in CIDcol
						myCIDind <- which(CIDcol == myCID)
						
						#Add one to AUDcol
						AUDcol[myCIDind] <- AUDcol[myCIDind] + 1
						
					}
					
					#Create auditTable
					AUDcol <- as.integer(AUDcol)
					auditTable = data.frame(CIDcol,AUDcol)
				}
				
			}
			Z$V$CID <- CID
		}
		
		#Make the strat object
		Z <- makeStratObj(Z, auditTable = auditTable)
	}
	
	Z
}

#Checks if object is a strat.elec.data object
is.strat.elec.data = function(Z) {
	inherits(Z, "strat.elec.data")
}

#Takes as input a vector of stratum sizes.
#Gives as output an order of counties;
#Sequentially sample from these counties to preserve a
#sample proportional to strata size. 
#Ties broken by sample size

#Assumes no stratum of 0 precincts.

#Can obtain stratum.sizes from Z$audit$n

getPropStratOrder <- function(stratum.sizes){
	numStrat = length(stratum.sizes)
	totPrecinct = sum(stratum.sizes)
	
	#propOrder is a 3 column matrix:
	#1) identifies stratum.
	#2) gives j*totPrecinct/stratum.sizes[i],
	#   which is sorted to give propStratOrder.
	#3) gives stratum sizes, used to break ties
	propOrder = NULL
	for( i in 1:numStrat){
		for( j in 0:(stratum.sizes[i]-1)){
			propOrder = 
			rbind(propOrder,c(i,j*totPrecinct/stratum.sizes[i],stratum.sizes[i]))
		}
	}
	
	#Sort to get propStratOrder
	sort.index = order(propOrder[,2],-propOrder[,3])
	propStratOrder = propOrder[sort.index,1]
	propStratOrder
}

#Removes precincts with MaximumMarginBound = 0
#Adjusts Z$strat appropriately
#Requires a strat.elec.data object

takeOutZeroMMB = function(Z){
	
	#Check if strat.elec.data object.
	if(!is.strat.elec.data(Z)){
		stop("takeOutZeroMMB requires a strat.elec.data object.")
	}
	
	#Only keep batches with maximumMarginBound > 0
	Z$V = Z$V[which(maximumMarginBound(Z) > 0),]
	
	#Adjust Z$strat$n vector
	for(i in 1:nrow(Z$strat)){
		Z$strat$n[i] = length(which(Z$V[Z$CID.col][,1] == Z$strat[Z$CID.col][i,1]))
	}
	Z	
}

#Obtains the working margin and maximum error bounds

#Takes an elec.data object and the maximum vote swing
#(either number of votes or as taint).
#asTaint is used to indicate if using taint.
#asNumber is used to indicate if t the maximum overstatment of the margin in votes.
#M is the margin of victory, as a real number.
#M defaults as 1, as according to the MRO discrepancy measure.

#For correctness, asNumber should only be used
#two candidates and one winner,
#or when running simulations.
#Otherwise, p-values may be conservative.

#Output is a list of two objects
#1) M is the working margin after subtracting sum(tminU)
#2) u is the vector of maximum error bounds

getEbsMargin = function(Z, t, asTaint = FALSE, asNumber = FALSE, M = NULL){
	
	#Check if elec.data object.
	if(!is.elec.data(Z)){
		stop("getEbsMargin requires an elec.data object.")
	}
	
	#If M is NULL, defaults to 1.
	if(is.null(M)){
		M=1
	}
	
	#Check if M is a number.
	if(!is.numeric(M)){
		stop("M must be a number")
	}
	
	#If Z has Z$V$e.max, use that.  
	#If not, use maximumMarginBound
	if(is.null(Z$V$e.max)){
		Z$V$e.max <- maximumMarginBound(Z)
	}
	eMax = Z$V$e.max
	
	#If as number, divide by margin
	if(asNumber){
		t = t/Z$margin
	}
	
	#If using taint...
	if(asTaint){
		u = eMax - t*eMax
		M = max(M-sum(t*eMax),0)
	}
	
	#Otherwise, if not using taint...
	else{
		tMinU <- pmin(eMax, t)
    	M <- max(M - sum(tMinU), 0)
      	u <- eMax - tMinU 
	}
	
	#Output M and u
	M.u <- list(M=M,u=u)
	M.u
}

#Obtain the vector q necessary for branch and bound algorithm
#Requires a strat.elec.data object

getQ = function(Z){
	
	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("getQ requries a strat.elec.data object.")
	}
	p <- 1
	#Creates one vector p of length N
	p <- do.call('c', lapply(1:nrow(Z$strat), function(c){
			
			#Precincts with largest e.max values are ranked first				 
			k<- rank(-Z$V$e.max[Z$V[[Z$CID.col]]== Z$strat[[Z$CID.col]][c]],
					 ties.method = "first")
							 
			#Compute the change in probability				 
			pmax((Z$strat$n[c] - Z$strat$audit[c] - k + 1) , 0) / (Z$strat$n[c] - k + 1)}))
	
	#Make transformation for branch and bound algorithm
  	q <- -log(p)
  	q
	
}

#Cleans strat.elec.data file to make
#algorithms more efficient.
#Z is a strat.elec.data object.
#M is the working margin
#u is a vector of maximum error bounds
#use takeOutZeroMMB if taking
#batches with a maximumMarginBound of zero.

#Used as part of LKP bound, 
#equal-value relaxation bound,
#and branch and bound.

stratClean = function(Z, M, u, takeOutZeroMMB = FALSE){
	
	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("stratClean requires a strat.data.elec object.")
	}
	
	if(takeOutZeroMMB){
		takeOutZeroMMB(Z)
	}
	
	#Find strata with zero batches audited.
	#Optimal solutions for algorithms will give
	#error to all batches within these strata.
	zeroStrat = Z$strat$CID[which(Z$strat$audit == 0)]
	
	#Eliminating and renaming strata
	if(length(zeroStrat) > 0){
		
		#Change M, and eliminate strata from Z$V
		for (i in zeroStrat){
			M = M - sum(u[which(Z$V[Z$CID.col]==i)])
			u = u[-which(Z$V[Z$CID.col]==i)]
			Z$V = Z$V[-which(Z$V[Z$CID.col]==i),]
		}
		
		#Eliminate stratum from Z$V$strat
		Z$strat = Z$strat[-which(Z$strat$audit == 0),]
		
		#Reorder CIDnum if nrow(Z$strat > 0)
		if(nrow(Z$strat) > 0){
			for(i in 1:nrow(Z$strat)){
				Z$V$CIDnum[which(Z$V[Z$CID.col]==Z$strat[Z$CID.col][i,])] <- i
				Z$strat$CIDnum[i] <- i
			}
		}
	}
	
	M = max(M,0)
	Z.M.u <- list(Z=Z,M=M,u=u)
	Z.M.u
}

#propSizes finds a sample proportional to stratum sizes.
#Requires a strat.elec.data object and the total number of batches
#to be sampled.

propSizes <- function(Z, n){
	
	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("next.r requires a strat.data.elec object.")
	}
	
	#myAudit holds sample
	myAudit <- rep(0,nrow(Z$strat))
	
	#myOrd holds strat order
	myOrd <- getPropStratOrder(Z$strat$n)
	
	#Build sample
	for(i in 1:n){
		myAudit[myOrd[i]] <-myAudit[myOrd[i]] + 1 
	}
	
	myAudit
}

#first.r computes a vector of sample sizes
#using the first.r algorithm:
#It computes the q values given an initial choice of sample sizes,
#finds the batch that corresponds to the largest q value,
#finds the stratum that corresponds to that batch,
#and adds one to the component of the sample vector
#that corresponds to that stratum.
#This process continues for n total iterations.

#Requires a strat.elec.data object and the total number of batches
#to be sampled.
#Other arguments to be used in getEbsMargin can be included.

#initSamp is an initial choice of sample sizes
#that first.r uses to build a sample from
#and is used in the get.first.r.samp function

#Avoids calling getQ to improve speed.

first.r <- function(Z, n, t = 0, asTaint = FALSE, asNumber = FALSE, M = NULL, initSamp = NULL){
	
	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("first.r.sample requires a strat.data.elec object.")
	}
	
	#Get u
	M.u <- getEbsMargin(Z,t, asTaint, asNumber, M)
	u <- M.u$u
	
	#Store number of strata
	noStrat <- nrow(Z$strat)
	
	#myAudit holds the sample vector.
	#If initSamp = NULL, myAudit is a 0 vector.
	#Otherwise, myAudit <- initSamp.  
	#If length(initSamp) != noStrat, there is an error, so stop.
	if(is.null(initSamp)){
	   myAudit <- rep(0, noStrat)
	}
	else{
		if(length(initSamp) != noStrat){
			stop("Length of vector initSamp needs to be the same as the number of stratum.")
		}
		myAudit <- initSamp
	}
	
	#Get largest values of u in each stratum
	#Get number of precincts in each stratum
	topn <- rep(0,noStrat)
	topU <- rep(0,noStrat)
	for(i in 1:noStrat){
		topU[i] <- max(u[which(Z$V$CIDnum == i)])
		topn[i] <- length(which(Z$V$CIDnum == i))
	}
	
	
	for(i in 1:n){
		#Get updated q values.
		#Avoid call to getQ by only considering only first value of q
		#within a stratum.
		topQ <- -log(pmax((topn - myAudit) , 0) / (topn))
		
		#Order q/u, 
		#find CIDnum corresponding to smallest value of q/u, 
		#add 1 to the corresponding stratum in sample vector
		#Ties are broken by stratum size
		myOrder <- order(topQ/topU, -topn)
		myAudit[myOrder[1]] <- myAudit[myOrder[1]] + 1
	}
	
	#Sync myAudit with Z$strat$CIDnum
	#(in case Z$strat gets permuted somehow)
	myAudit <- myAudit[rank(Z$strat$CIDnum, ties.method = "first")]
	
	myAudit
	
}

#next.r computes a vector of sample sizes
#using the next.r algorithm:
#It computes the q values given an initial choice of sample sizes,
#finds the batch that correspond to the largest q value,
#finds the stratum that corresponds to that batch,
#and adds one to the component of the sample vector
#that corresponds to that stratum.
#Then it removes that batch from consideration.
#This process continues for n total iterations.

#Requires a strat.elec.data object and the total number of batches
#to be sampled.
#Other arguments to be used in getEbsMargin can be included.

#initSamp is an initial choice of sample sizes
#that next.r uses to build a sample from
#and is used in the get.next.r.samp function

#Avoids calling getQ to improve speed.
next.r <- function(Z, n, t = 0, asTaint = FALSE, asNumber = FALSE, M = NULL, initSamp = NULL){

	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("next.r requires a strat.data.elec object.")
	}
	
	#Get u
	M.u <- getEbsMargin(Z,t, asTaint, asNumber, M)
	u <- M.u$u
	
	#Store number of strata
	noStrat <- nrow(Z$strat)
	
	#myAudit holds the sample vector.
	#If initSamp = NULL, myAudit is a 0 vector.
	#Otherwise, myAudit <- initSamp.  
	#If length(initSamp) != noStrat, there is an error, so stop.
	if(is.null(initSamp)){
		myAudit <- rep(0, noStrat)
	}
	else{
		if(length(initSamp) != noStrat){
			stop("Length of vector initSamp needs to be the same as the number of stratum.")
		}
		myAudit <- initSamp
	}
	
	#Get number of precincts in each stratum
	topn <- rep(0,noStrat)
	for(i in 1:noStrat){
		topn[i] <- length(which(Z$V$CIDnum == i))
	}
	
	
	for(i in 1:n){
		#Get updated q values.
		#Avoid call to getQ by only considering only the myAudit[c]th largest
		#value of q within a stratum.
		#Use max(q, 1/(myAudit + 1))
		#to avoid infinities
		
		topQ <- -log(pmax(pmax((topn - 2*myAudit + 1) , 0) / (topn - myAudit + 1), 1/(myAudit + 1)))
		
		#Get myAudit[c]th largest u values.
		#u value defaults to 0 if myAudit = 0 (q = 0 in this case)
		topU <- rep(1,noStrat)
		for(i in 1:noStrat){
			if(myAudit[i] != 0){
				topU[i] <- sort(u[which(Z$V$CIDnum == i)], decreasing = TRUE)[myAudit[i]]
			}
		}
		
		#Order q/u, 
		#find CIDnum corresponding to smallest value of q/u, 
		#add 1 to the corresponding stratum in sample vector
		#Ties are broken by stratum size
		myOrder <- order(topQ/topU, -topn)
		myAudit[myOrder[1]] <- myAudit[myOrder[1]] + 1
	}
	
	#Sync myAudit with Z$strat$CIDnum
	#(in case Z$strat gets permuted somehow)
	myAudit <- myAudit[rank(Z$strat$CIDnum, ties.method = "first")]
	
	myAudit
	
}
							  
#LKPBound computes the LKP bound;
#Gives an upper-bound of the p-value.
#If LKP.lower.bound = TRUE, it will also give
#a lower-bound of p-value.

#Requires a strat.elec.data object as input.
#t, asTaint, and M used in call of getEbsMargin.

#If t = NULL, Z$audit and compute.stark.t used to get value of t.
#Defaults to maximum observed difference, uses maximum observed taint if asTaint = TRUE.
#bound.col,  calc.e_p, w_p used in call of compute.stark.t

#Use takeOutZeroMMB = TRUE if not going to include
#batches with a maximumMarginBound of zero.

#Uses stratClean: any stratum with no audited batches
#have maximum error assigned to all batches within stratum 

#Other arugments

LKPBound <- function(Z, t = NULL, asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB = TRUE, LKP.lower.bound = FALSE, bound.col = "e.max",  calc.e_p=calc.pairwise.e_p, w_p = weight.function("no.weight")){
		
	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("LKPBound requires a strat.elec.data object.")
	}
	
	#Check if valid audit column...
	#Are the audit values integers?
	if(!isTRUE(all.equal(floor(Z$strat$audit), Z$strat$audit))){
		stop("Audit column must have integer values")
	}
	
	#Are there no negative audit values?
	if(min(Z$strat$audit) < 0){
		stop("No negative values of audit column allowed")
	}
	
	#Take out zeroMMBs...
	if(takeOutZeroMMB){
		takeOutZeroMMB(Z)
	}
	
	#Find t if none is given
	if(is.null(t)){
		
		#compute.stark.t will produce a t that is not asNumber
		asNumber = FALSE
		
		#If audit not included in Z,
		if(is.null(Z$audit)){
			stop("t is null and Z$audit not given.")	
		}
		
		#If asTaint = FALSE, compute t as using default parameters 
		if(!asTaint){
			t <- compute.stark.t(Z, bound.col, calc.e_p, w_p)
		}
		
		#If asTaint = TRUE, compute t using w_p = taint!
		else{
			t <- compute.stark.t(Z, bound.col, calc.e_p, w_p = weight.function("taint"))
		}
	}
	
	#Get Ebs, Margin, clean Z
	M.u <- getEbsMargin(Z,t,asTaint,asNumber,M)
	M <- M.u$M
	u <- M.u$u
	
	Z.M.u <- stratClean(Z,M,u)
	Z <- Z.M.u$Z
	M <- Z.M.u$M
	u <- Z.M.u$u
	
	#Trivial solution unless M > 0.
	if (M > 0) {
		
		#Trivial solution unless nrow(Z$strat) > 0
		if(nrow(Z$strat) > 0){
		    
			#Get q
			q <- getQ(Z)
			
			#Put in order of q/u
			ord <- order(q/u) 
			
			q <- q[ord]
			u <- u[ord]
			
			#urun is the running sum of the u
			#qrun is the running sum of the q
			#Find value of B, 
			#sum_1^B u > M
		    urun <- 0
			qrun <- 0
		    B <- 0
    
		    while(urun + u[B+1] <= M){
				B <- B+1
				urun <- urun + u[B]
				qrun <- qrun + q[B]
		    }
			
			#Compute upper and lower p-values using LKP bound
			D.upper <- qrun + q[B+1] * (M - urun) / u[B+1]  
			if(urun == M){
				D.lower <- D.upper
			}
			else{
				D.lower <- qrun + q[B+1]
			}
			bounds <- exp(-c(D.upper, D.lower)) 
		}
		#Trivial soloution.
		else{
			bounds <- c(0,0)
		}
	} 
	
	#Trivial solution.
	else {
		bounds <- c(1,1)
	}
  
	names(bounds) <- c("LKP.upper.bound", "LKP.lower.bound")
	if(!LKP.lower.bound){
		bounds <- bounds[1]
	}
	bounds
}
			
#breakValue finds the smallest number of batches
#required to make sum(u[Batches]) >= M
#Used in eqValBound and withReplaceBound

#Requires as input a vector of upper-bounds u
#and a margin M.

#Outputs a scalar B

breakValue <- function(u,M){
	
	#Find break value B
	#sortu is the vector sorted in decreasing order
	#runu is the running sum of u
	B <- 0
	sortu <- sort(u, decreasing = TRUE)
	runu <- 0
	while(runu < M){
		B <- B + 1
		runu <- runu + sortu[B]
	}
	B
}


#eqValBound computes the equal value bound;
#relaxes condition that the sum of the values is
#greater than M to the condition that
#the number of samples is less than or equal to the break value.

#Gives an upper-bound of the p-value.

#Requires a strat.elec.data object as input.
#t, asTaint, asNumber, and M used in call of getEbsMargin.

#If t = NULL, Z$audit and compute.stark.t used to get value of t.
#Defaults to maximum observed difference, uses maximum observed taint if asTaint = TRUE.
#bound.col,  calc.e_p, w_p used in call of compute.stark.t

#Use takeOutZeroMMB = TRUE if not going to include
#batches with a maximumMarginBound of zero.

#Uses stratClean: any stratum with no audited batches
#have maximum error assigned to all batches within stratum

eqValBound <- function(Z, t = NULL, asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB=TRUE, bound.col = "e.max",  calc.e_p=calc.pairwise.e_p, w_p = weight.function("no.weight")){
	
	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("eqValBound requires a strat.elec.data object.")
	}
	
	#Check if valid audit column...
	#Are the audit values integers?
	if(!isTRUE(all.equal(floor(Z$strat$audit), Z$strat$audit))){
		stop("Audit column must have integer values")
	}
	
	#Are there no negative audit values?
	if(min(Z$strat$audit) < 0){
		stop("No negative values of audit column allowed")
	}
	
	if(takeOutZeroMMB){
		takeOutZeroMMB(Z)
	}
	
	#Find t if none is given
	if(is.null(t)){
		
	#compute.stark.t will produce a t that is not asNumber
		asNumber = FALSE
		
		#If audit not included in Z,
		if(is.null(Z$audit)){
			stop("t is null and Z$audit not given.")	
		}
		
		#If asTaint = FALSE, compute t as using default parameters 
		if(!asTaint){
			t <- compute.stark.t(Z, bound.col, calc.e_p, w_p)
		}
		
		#If asTaint = TRUE, compute t using w_p = taint!
		else{
			t <- compute.stark.t(Z, bound.col, calc.e_p, w_p = weight.function("taint"))
		}
	}	
	
	#Get Ebs, Margin, clean Z
	M.u <- getEbsMargin(Z,t,asTaint,asNumber,M)
	M <- M.u$M
	u <- M.u$u
	
	Z.M.u <- stratClean(Z,M,u)
	Z <- Z.M.u$Z
	M <- Z.M.u$M
	u <- Z.M.u$u
	
	#Holds p.value
	eqVal.u <- NULL
	
	#Trivial solution if M <= 0.
	if(M <= 0){
		eqVal.u <- 1
	}
	
	else{
		
		#Trivial solution if sum(u) < M
		if(sum(u) < M){
			eqVal.u	<- 0
		}
		
		#Otherwise...
		else{
			
			#Get break value
			B <- breakValue(u,M)
			
			#Get q values
			q <- getQ(Z)
			
			#sum up first B values of q to get
			#negLog of p-value
			negLogEqVal.u <- sum(sort(q)[1:B])
			
			#Make transformation to get p-value
			eqVal.u <- exp(-negLogEqVal.u)
		}
	}
	
	names(eqVal.u) <- "eqVal.u"
	eqVal.u
}

#withReplaceBound computes the bound as in the
#Conservative Election Audits paper by Stark.
#Works by showing that stratified random sampling
#is no worse than unstratified random sampling with replacement.
#Gives an upper-bound of the p-value.

#Requires a strat.elec.data object as input.
#t, asTaint, asNumber, and M used in call of getEbsMargin.

#Use takeOutZeroMMB = TRUE if not going to include
#batches with a maximumMarginBound of zero.

#Uses stratClean: any stratum with no audited batches
#have maximum error assigned to all batches within stratum

withReplaceBound <- function(Z, t = NULL, asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB = TRUE, bound.col = "e.max",  calc.e_p=calc.pairwise.e_p, w_p = weight.function("no.weight")){
								 
	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("withReplaceBound requires a strat.elec.data object.")
	}
						 
	#Check if valid audit column...
	#Are the audit values integers?
	if(!isTRUE(all.equal(floor(Z$strat$audit), Z$strat$audit))){
		stop("Audit column must have integer values")
	}
	
	#Are there no negative audit values?
	if(min(Z$strat$audit) < 0){
		stop("No negative values of audit column allowed")
	}
	
	if(takeOutZeroMMB){
		takeOutZeroMMB(Z)
	}
	
	#Find t if none is given
	if(is.null(t)){
		
		#compute.stark.t will produce a t that is not asNumber
		asNumber = FALSE
		
		#If audit not included in Z,
		if(is.null(Z$audit)){
			stop("t is null and Z$audit not given.")	
		}
		
		#If asTaint = FALSE, compute t as using default parameters 
		if(!asTaint){
			t <- compute.stark.t(Z, bound.col, calc.e_p, w_p)
		}
		
		#If asTaint = TRUE, compute t using w_p = taint!
		else{
			t <- compute.stark.t(Z, bound.col, calc.e_p, w_p = weight.function("taint"))
		}
	}
	
	
	#Get Ebs, Margin, clean Z
	M.u <- getEbsMargin(Z,t,asTaint,asNumber,M)
	M <- M.u$M
	u <- M.u$u
	
	Z.M.u <- stratClean(Z,M,u)
	Z <- Z.M.u$Z
	M <- Z.M.u$M
	u <- Z.M.u$u
	
	#Holds p.value
	wRep.u <- NULL
	
	#Trivial solution if M <= 0.
	if(M <= 0){
		wRep.u <- 1
	}
	
	else{
		
		#Trivial solution if sum(u) < M
		if(sum(u) < M){
			wRep.u <- 0
		}
		
		#Otherwise...
		else{
			
			#Get break value
			B <- breakValue(u,M)
			
			#Get no. of batches,
			#no. of batches sampled
			noBatch <- sum(Z$strat$n)
			
			#Get smallq as in Stark
			smallq <- noBatch - B
			
			#Get smallest sampling fraction
			sampFrac <- min(Z$strat$audit/Z$strat$n)
			
			#Get new no. of samps
			noSamp <- floor(sampFrac*noBatch)
			
			#Get p-value
			wRep.u <- (smallq/noBatch)^(noSamp)
			
		}	
	}
	
	names(wRep.u) <- "wRep.u"
	wRep.u					 
}	
						 

#Calls the Branch and Bound .C function.
#Pre-processes data; 
#Branch and bound may not be needed.
#Requires a strat.elec.data object.

#t, asTaint, asNumber, and M used in call of getEbsMargin.

#Use takeOutZeroMMB = TRUE if not going to include
#batches with a maximumMarginBound of zero.

#Outputs is a p-value, and if give.strategy = TRUE, the solution 
#to the 0-1 KP.

#Solution is given by the number of batches assigned error
#to each stratum, where stratum are identified by CIDnum

#No call to stratClean needed!

BaB <- function(Z, t =  NULL, asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB=TRUE, 
give.strategy = FALSE, bound.col = "e.max",  calc.e_p=calc.pairwise.e_p, w_p = weight.function("no.weight")){

	#Check if strat.elec.data object
	if(!is.strat.elec.data(Z)){
		stop("BaB requires a strat.elec.data object.")
	}
	
	#Check if valid audit column...
	#Are the audit values integers?
	if(!isTRUE(all.equal(floor(Z$strat$audit), Z$strat$audit))){
		stop("Audit column must have integer values")
	}
	
	#Are there no negative audit values?
	if(min(Z$strat$audit) < 0){
		stop("No negative values of audit column allowed")
	}
	
	if(takeOutZeroMMB){
		Z = takeOutZeroMMB(Z)	
	}
	
	#Find t if none is given
	if(is.null(t)){
		
		#compute.stark.t will produce a t that is not asNumber
		asNumber = FALSE
		
		#If audit not included in Z,
		if(is.null(Z$audit)){
			stop("t is null and Z$audit not given.")	
		}
		
		#If asTaint = FALSE, compute t as using default parameters 
		if(!asTaint){
			t <- compute.stark.t(Z, bound.col, calc.e_p, w_p)
		}
		
		#If asTaint = TRUE, compute t using w_p = taint!
		else{
			t <- compute.stark.t(Z, bound.col, calc.e_p, w_p = weight.function("taint"))
		}
	}
	
	#Get Ebs, Margin, clean Z
	M.u <- getEbsMargin(Z,t,asTaint,asNumber,M)
	M <- M.u$M
	u <- M.u$u
	
	#Initialize p.value and strategy
	p.value = 0
	strategy = rep(0,nrow(Z$strat))
	
	#If audit table is identically zero:
	if(sum(Z$strat$audit) == 0){
		warning("Audit table is identically zero.")
		
		#If M- sum(u) is greater than zero, 
		#tainted everything and still not greater than margin.
		#p.value = 0
		#Else, p.value = 1	
		if(sum(u) - M > 0){
			p.value = 0
		}
		else{
			p.value = 1
		}
		
		#Get p.value into a list
		returnlist <- list(p.value = p.value)
		
		if(give.strategy){
			
			#strategy is strat$n
			strategy = Z$strat$n
			returnlist$strategy <- strategy
		}
	}
	
	else{
		#Get Q
		q = getQ(Z)
		q[is.infinite(q)] <- 1000000
	
		#Get CIDnum
		CIDnum = Z$V$CIDnum
	
		#Run BaB
		result <- runBaB(u,q,M,CIDnum)
	
		#Get p.value
		p.value <- exp(-result[[1]])
	
		#Get p.value into a list
		returnlist <- list(p.value = p.value)
	
		if(give.strategy){	
	
			#Get strategy
			prestrat <- result[[2]]*result[[5]]
			for(i in 1:nrow(Z$strat)){
				strategy[i] <- sum(prestrat == i)	
			}
			names(strategy) <- 1:nrow(Z$strat)
			
			#Get strategy into a list
			returnlist$strategy <- strategy
		}
	}	
		
	#Output returnlist
	returnlist
}	

#Runs the BaB function in C.
runBaB <- function(u,q,M,CIDnum){
	#Just to check that there are no infinities...
	q[is.infinite(q)] <- 1000000
	
	#Find r, order things in terms of r
	r <- q/u
	index <- order(r)
	u <- u[index]
	q <- q[index]
	CIDnum <- CIDnum[index]
	
	#Run C script
	result <-.C("BaB", as.double(sum(q)), as.integer(rep(0,length(u))), 
						as.double(u), as.double(q), as.integer(CIDnum), as.integer(length(unique(CIDnum))), 
						as.integer(length(u)), as.double(M))
	result
}

#get.first.r.samp finds a sample proportional to stratum sizes,
#given a fixed value of t and an alpha, 
#to obtain a p-value less than alpha
#
#bal = TRUE will give the expected number of ballots
#samp = TRUE will give the optimal sample
#initn gives the first number of samples to check.

get.prop.samp <- function(Z, alpha, t, bal=TRUE, numSamp = TRUE, initn = 1, asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB=TRUE){
	
	#Initialize
	#n <- initn - 1 to make while loop start out
	#at n = initn
	
	#go = TRUE for while loop
	#myStrat hold the sample
	go = TRUE
	n <- initn - 1
	myStrat <- 0
	
	#If taking out zero MMB...
	if(takeOutZeroMMB){
		Z = takeOutZeroMMB(Z)	
	}
	
	#Get order for proportional strat
	ordProp <- getPropStratOrder(Z$strat$n)
	
	#Find prop strat for n samples. 
	#Check if the strat gives
	#p-value less than alpha
	while(go){
		
		#Increment n
		n <- n + 1
		
		#Copy Z
		b <- Z
		
		#If n < length(b$strat$audit),
		#then table function wont work
		if(n < nrow(b$strat)){
			aud <- rep(0,nrow(b$strat))
			aud[ordProp[1:n]] <- 1
			names(aud) <- 1:nrow(Z$strat)
			b$strat$audit <- as.table(aud)
		}
		
		else{
			b$strat$audit <- table(ordProp[1:n])
		}
		
		#Run branch and bound 
		pval <- BaB(b,t,asTaint,asNumber,M,takeOutZeroMMB)
		
		#Check if less than alpha
		if(pval < alpha){
			
			#Stop loop.  Record sample.
			go = FALSE
			myStrat <- b$strat$audit
		}
		
	}
	
	#We will return result
	result <- list(samp = myStrat)
	
	#If bal = TRUE
	if(bal){		
		#For each strata, find average number of ballots in
		#precincts belonging to that county
		
		#Initialize
		avgincounty <- rep(0,nrow(Z$strat))
		
		#Find average number of ballots per precinct
		#for each stratum
		for(i in 1:length(Z$strat$CIDnum)){
			avgincounty[i] <- mean(Z$V$tot.votes[which(
								   Z$V$CIDnum == Z$strat$CIDnum[i])])
		}
		
		#Find the expected number of ballots 
		#and add it to results
		exBal <- sum(myStrat*avgincounty)
		result$exBal  <- exBal
	}
	
	if(numSamp){
		
		#Add the sample to the list
		result$numSamp <- sum(myStrat)	
	}
	
	result
}

#get.first.r.samp finds a sample using the first. r algorithm,
#given a fixed value of t and an alpha, 
#to obtain a p-value less than alpha
#
#bal = TRUE will give the expected number of ballots
#numSamp = TRUE will give the number of samples
#initn gives the first number of samples to check.

get.first.r.samp <- function(Z, alpha, t, bal=TRUE, numSamp = TRUE, initn = 1, asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB=TRUE){
	
	#Initialize
	go = TRUE
	myStrat <- 0
	
	#If taking out zero MMB...
	if(takeOutZeroMMB){
		Z = takeOutZeroMMB(Z)	
	}
	
	#Copy Z
	b <- Z
	
	#initSamp holds initial sample.
	initSamp <- first.r(b, initn, t, asTaint, asNumber, M)
	
	#Run branch and bound
	b$strat$audit <- initSamp
	pval <- BaB(b,t,asTaint,asNumber,M,takeOutZeroMMB)
	
	#Check if less than alpha
	if(pval < alpha){
		
		#Stop loop.  Record sample.
		go = FALSE
		myStrat <- b$strat$audit
	}
	
	#If initn does not give a sample less than alpha,
	#keep incrementing n until obtaining a sample that is less than alpha
	while(go){
		
		#Copy Z
		b <- Z
		
		#Find first.r sample
		initSamp <- first.r(b, 1, t, asTaint, asNumber, M, initSamp = initSamp)
		b$strat$audit <- initSamp
		
		#Run branch and bound
		pval <- BaB(b,t,asTaint,asNumber,M,takeOutZeroMMB)
		
		#Check if less than alpha
		if(pval < alpha){
			
			#Stop loop.  Record sample.
			go = FALSE
			myStrat <- b$strat$audit
		}
	}
	
	#We will return result
	names(myStrat) <- 1:length(myStrat)
	result <- list(samp = myStrat)
	
	#If bal = TRUE
	if(bal){
			
		#For each strata, find average number of ballots in
		#precincts belonging to that county
		
		#Initialize
		avgincounty <- rep(0,nrow(Z$strat))
		
		#Find average number of ballots per precinct
		#for each stratum
		for(i in 1:length(Z$strat$CIDnum)){
			avgincounty[i] <- mean(Z$V$tot.votes[which(
								   Z$V$CIDnum == Z$strat$CIDnum[i])])
		}
		
		#Find the expected number of ballots 
		#and add it to results
		exBal <- sum(myStrat*avgincounty)
		result$exBal  <- exBal
	}
	
	if(numSamp){
		
		#Add the sample to the list
		result$numSamp <- sum(myStrat)	
	}
	
	result
}

#get.next.r.samp finds a sample using the next.r algorithm,
#given a fixed value of t and an alpha, 
#to obtain a p-value less than alpha
#
#bal = TRUE will give the expected number of ballots
#numSamp = TRUE will give the number of samples
#initn gives the first number of samples to check.

get.next.r.samp <- function(Z, alpha, t, bal=TRUE, numSamp = TRUE, initn = 1, asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB=TRUE){
	
	#Initialize
	go = TRUE
	myStrat <- 0
	
	#If taking out zero MMB...
	if(takeOutZeroMMB){
		Z = takeOutZeroMMB(Z)	
	}
	
	#Copy Z
	b <- Z
	
	#initSamp holds initial sample.
	initSamp <- next.r(b, initn, t, asTaint, asNumber, M)
	
	#Run branch and bound
	b$strat$audit <- initSamp
	pval <- BaB(b,t,asTaint,asNumber,M,takeOutZeroMMB)
	
	#Check if less than alpha
	if(pval < alpha){
		
		#Stop loop.  Record sample.
		go = FALSE
		myStrat <- b$strat$audit
	}
	
	#If initn does not give a sample less than alpha,
	#keep incrementing n until obtaining a sample that is less than alpha
	while(go){
		
		#Copy Z
		b <- Z
		
		#Find first.r strat
		initSamp <- next.r(b, 1, t, asTaint, asNumber, M, initSamp = initSamp)
		b$strat$audit <- initSamp
		
		#Run branch and bound
		pval <- BaB(b,t,asTaint,asNumber,M,takeOutZeroMMB)
		
		#Check if less than alpha
		if(pval < alpha){
			
			#Stop loop.  Record sample.
			go = FALSE
			myStrat <- b$strat$audit
		}
		
	}
	
	#Return sample
	names(myStrat) <- 1:length(myStrat)
	result <- list(samp = myStrat)
	
	#If bal = TRUE
	if(bal){
		
		#For each strata, find average number of ballots in
		#precincts belonging to that county
		
		#Initialize
		avgincounty <- rep(0,nrow(Z$strat))
					
		#Find average number of ballots per precinct
		#for each stratum
		for(i in 1:length(Z$strat$CIDnum)){
			avgincounty[i] <- mean(Z$V$tot.votes[which(
													   Z$V$CIDnum == Z$strat$CIDnum[i])])
		}
		
		#Find the expected number of ballots 
		#and add it to results
		exBal <- sum(myStrat*avgincounty)
		result$exBal  <- exBal
	}
	
	if(numSamp){
		
		#Add the number of samples to the list
		result$numSamp <- sum(myStrat)	
	}
	
	result
}

#nextSample is used in the optStrat function.
#Finds the next sample in the brute force solution
#of optStrat
nextSample = function(currentSample){
	
	#Initialize
	enn = sum(currentSample)
	len = length(currentSample)
	nl = enn+len
	j = 1
	
	#newConfig is manipulated to find
	#the next sample sizes
	currentSample = currentSample + 1
	newConfig = c(cumsum(currentSample),1)
	
	#Done = true: We have exhausted all samples
	#a = true: We have found the next sample
	done = FALSE
	a = FALSE
	
	#Gives the next sample
	while(!a){
		if(newConfig[j] + 1 == newConfig[j+1]){
			newConfig[j] = j
			j = j+1
			if(j == len){
				a = TRUE
				done = TRUE
			}	
		}
		else{
			newConfig[j] = newConfig[j]+1
			a = TRUE	
		}
	}
	
	newConfig[len] = nl
	newConfig[len+1] = 0
	newConfig = sort(newConfig)
	
	#If there are no more samples,
	#output stop.
	#If there is another sample
	#output the new sample.
	if(done){
		newSample = "stop"
	}
	
	else{
		newSample = diff(newConfig) - 1
	}
	
	newSample	
}

#optStrat finds an optimal sample,
#given a fixed value of t and an alpha, 
#so that the p-value is less than alpha.
#Optimal is defined as requiring the fewest number of audited batches.
#
#optBal = TRUE will find the sample, given the minimum number of audited batches
#that gives the minimum expected number of audited ballots
#(which requires substantial computing time)
#
#numSamp = TRUE will give the number of samples
#bal = TRUE will give the expected number of audited ballots, given the sample obtained by optStrat

optStrat = function(Z,alpha, t,  bal=TRUE, optBal=FALSE, numSamp = TRUE, asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB=TRUE){
	
	#If taking out zero MMB...
	if(takeOutZeroMMB){
		Z = takeOutZeroMMB(Z)	
	}
		
	#For each strata, find average number of ballots in
	#precincts belonging to that county
	avgincounty <- rep(0,nrow(Z$strat))	
	
	for(i in 1:length(Z$strat$CIDnum)){
		avgincounty[i] <- mean(Z$V$tot.votes[which(
												   Z$V$CIDnum == Z$strat$CIDnum[i])])
	}
	
	#Initialize: use first.r to get an upper-bound on the number of samples needed
	initn.strat = get.first.r.samp(Z ,alpha, t, FALSE, TRUE, 1, asTaint, asNumber, M, takeOutZeroMMB)
	initn <- initn.strat$numSamp
	initStrat <- initn.strat$samp
	
	#Start checking for solutions with initn - 1 batches audited
	initn = initn - 1
		
	#toStop is true if we stop the while loop
	toStop = FALSE
	
	while(!toStop){
		
		#Initialize cycling through all possible ways of allocating sample
		#Move on to next sample...
		nextSamp = TRUE
		
		#InitSamp is first sample to check
		initSamp = rep(0,nrow(Z$strat))	
		initSamp[nrow(Z$strat)] = initn
		
		#enn counts iterations.
		enn = 0
		
		#optSamp stores optimal sample,
		#minp stores minimum pvalue
		optSamp <- initSamp
		minp <- 1
		
		#Run through all samples, reduce initn by one if a sample found that
		#has pvalue less than alpha
		while(nextSamp){
			
			#fix Z
			fixZ <- Z
			
			#Get pvalue
			fixZ$strat$audit <- initSamp
			pvalue  = BaB(fixZ, t, asTaint, asNumber, M = NULL, takeOutZeroMMB, 
						   give.strategy = FALSE)$p.value
			
			#If found a better pvalue, 
			#update minp, optSamp
			if(pvalue < minp){
				minp <- pvalue
				optSamp <-initSamp	
			}
			
			#Get next sample
			initSamp = nextSample(initSamp)
			
			#if pvalue < alpha 
			#Update initStrat, quit loop.
			if(pvalue <= alpha){
				nextSamp = FALSE
				initStrat = optSamp
			}
			
			#If we have gone through all possible samples and did not
			#Find a sample with pvalue less than alpha,
			#Exit both loops.
			if((length(initSamp) != nrow(Z$strat)) & nextSamp){
				nextSamp = FALSE	
				toStop = TRUE
			}
			
			#For fun.
			enn = enn+ 1
		}
		
		#If we exited first loop but did not exit second, reduce initn
		if(!toStop){
			initn = initn - 1	
		}
		
	}
	
	#Add one to initn to get total number of samples
	initn = initn + 1
	
	#Return n.audit
	names(initStrat) <- 1:length(initStrat)
	n.audit <- list(samp = initStrat)
	
	#If numSamp, add number of samples
	if(numSamp){
		n.audit$numSamp <- initn
	}
	
	#For each strata, find average number of ballots in
	#precincts belonging to that county
	
	#Initialize
	avgincounty <- rep(0,nrow(Z$strat))
	
	#Find average number of ballots per precinct
	#for each stratum
	for(i in 1:length(Z$strat$CIDnum)){
		avgincounty[i] <- mean(Z$V$tot.votes[which(
												   Z$V$CIDnum == Z$strat$CIDnum[i])])
	}
	
	#If bal = TRUE
	if(bal){
		
		#Find the expected number of ballots 
		#and add it to results
		exBal <- sum(initStrat*avgincounty)
		n.audit$exBal  <- exBal
	}
	
	#If finding the sample that minimizes audited ballots...
	if(optBal){
		#For each strata, find average number of ballots in
		#precincts belonging to that county
		
		#Initialize
		optSamp <- initSamp
		minBal <- sum(avgincounty*Z$strat$n)
		initSamp = rep(0,nrow(Z$strat))	
		initSamp[nrow(Z$strat)] = initn
		enn = 0
		
		#Go until reach "stop" sample.
		while(length(initSamp) == nrow(Z$strat)){
			
			#fix Z
			fixZ <- Z
			
			#Get pvalue
			fixZ$strat$audit <- initSamp
			pvalue  = BaB(fixZ, t, asTaint, asNumber, M = NULL, takeOutZeroMMB, 
						  give.strategy = FALSE)$p.value
			
			#if pvalue <= alpha 
			#If found a new minBal, update minBal, optSamp
			if(pvalue <= alpha){
				if(sum(avgincounty*initSamp) < minBal){
					minBal <- sum(avgincounty*initSamp)
					optSamp <- initSamp
				}
			}
				
			#Get next sample
			initSamp = nextSample(initSamp)
			
			#For fun.
			enn = enn+ 1	
		}
		
		#Update n.audit
		names(optSamp) <- 1:length(optSamp)
		n.audit$samp <- optSamp
		if(bal){
			n.audit$exBal <- minBal
		}
		
	}
	n.audit
}

	

