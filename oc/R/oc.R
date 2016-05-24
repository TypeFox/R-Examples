#.First.lib <- function(lib,pkg) {
#   library.dynam("oc",pkg,lib)
#}


###################################################################################
#checkMember():  Takes vote matrix and returns total number of retained            
#                votes for each member.  Used as wnominate() helper.               
###################################################################################
                                                                                  
checkMember <- function(votes) {
     apply(votes!=9,1,sum)
}

##################################################################################
#checkVote():   Takes vote matrix and returns total number of retained votes      
#               for each member.  Used as wnominate() helper.                     
##################################################################################
                                                                                 
checkVote <- function(votes) {
    minorityPercent <- pmin(apply(votes==1,2,sum),apply(votes==6,2,sum))
    minorityPercent <- minorityPercent/apply(votes!=9,2,sum)
    minorityPercent[is.nan(minorityPercent)]<-0    
    return(minorityPercent)
}

oc <- function(rcObject, dims=2, minvotes=20, lop=0.025, polarity, verbose=FALSE) {

    cat("\nPreparing to run Optimal Classification...\n\n\tChecking data...\n\n")
    start <- proc.time()
    
    # The following modifies an arugment, but I think we have to live with this 
    # one... If we don't, then we have to copy rcObjects$votes which is not a good 
    # thing memorywise.    
    if(!is.null(rcObject$codes)) {
        rcObject$votes[rcObject$votes %in% rcObject$codes$yea] <- 1
        rcObject$votes[rcObject$votes %in% rcObject$codes$nay] <- 6
        rcObject$votes[rcObject$votes %in% c(rcObject$codes$missing,
                        rcObject$codes$notInLegis)] <- 9
    }
    
    #Core Error Checking
    if(!class(rcObject)=="rollcall")
        stop("Input is not of class 'rollcall'.")

    if(dims>10) stop("'dims' cannot exceed 10\n") 
    if(dims<1) stop("'dims' cannot be negative\n") 

    polarityToUse <- polarity
    if(is.character(polarityToUse)) {
        if(length(setdiff(polarityToUse,rownames(rcObject$votes)))!=0) {
            cat("\tThe following legislators do not exist:\n\t\t",
            setdiff(polarityToUse,rownames(rcObject$votes)),"\n")
            stop("'polarity' is incorrectly specified\n")      
        }
        duplicates<-rownames(rcObject$votes)[duplicated(rownames(rcObject$votes))]
        if(length(intersect(duplicates,polarityToUse))!=0){
            cat("\tThe following legislators have duplicate names:\n\t\t",
                intersect(duplicates,polarityToUse),"\n")
            stop("'polarity' is incorrectly specified\n")      
        }
        polarityToUse<-match(polarityToUse,rownames(rcObject$votes))
    }
    
    if(is.list(polarityToUse)) {
        column<-which(colnames(rcObject$legis.data)==polarityToUse[[1]])
        if(length(column)==0) stop("Variable '", polarityToUse[[1]] ,"' does not exist in 'legis.data'")
        if(length(setdiff(polarityToUse[[2]],rcObject$legis.data[,column])==0)) {
            cat("\tThe following items from your variable do not exist:\n\t\t",
            setdiff(polarityToUse[[2]],rcObject$legis.data[,column]),"\n")
            stop("'polarity' is incorrectly specified\n") 
        }
        duplicates<-rcObject$legis.data[duplicated(rcObject$legis.data[,column]),column]
    if(length(intersect(duplicates,polarityToUse[[2]]))!=0) {
        cat("\tThe following legislators have duplicate items from your list:\n\t\t",
        intersect(duplicates,polarityToUse[[2]]),"\n")
        stop("'polarity' is incorrectly specified\n")      
        }
        polarityToUse <- match(polarityToUse[[2]], rcObject$legis.data[,column])
    }
    if(dims!=length(polarityToUse)) stop("'polarity' must be a vector of length 'dims'\n")
    if(max(polarityToUse)>rcObject$n | min(polarityToUse)<1) {
        cat("\tThe following legislators do not exist:\n\t\t",
            polarityToUse[which(polarityToUse>rcObject$n | polarityToUse<1)],"\n")
        stop("'polarity' is incorrectly specified\n") 
    }
    if((lop<0) | (lop>1)) stop("'lop' does not fall between 0 to 1.\n")
    entry<-as.numeric(sort(unique(as.integer(rcObject$votes)), na.last=TRUE))
        #as.numeric is required because there is an R bug at work here
    if(!(identical(entry,c(1,6,9)) | identical(entry,c(1,6)) )) {
        cat("\tData contains the following values: ", entry, "\n")
        stop("Data contains values other than 1 or 6 or 9.\n")
    }
    
    #Check legislators and members for minimum requirements until there is no further change
    memberVotes <- checkMember(rcObject$votes)
    minorityVoteShare <- checkVote(rcObject$votes)
    tempvotes<-rcObject$votes
    tempvotes[memberVotes<minvotes,]<-9
    tempvotes[,minorityVoteShare<=lop]<-9
    memberVotes2 <- checkMember(tempvotes)
    minorityVoteShare2 <- checkVote(tempvotes)

    while(!identical(memberVotes,memberVotes2)|!identical(minorityVoteShare,minorityVoteShare2)){
        memberVotes <- memberVotes2
        minorityVoteShare <- minorityVoteShare2
        tempvotes[memberVotes2<minvotes,] <- 9
        tempvotes[,minorityVoteShare2<=lop] <- 9
        memberVotes2 <- checkMember(tempvotes)
        minorityVoteShare2 <- checkVote(tempvotes)
    }

    members2dump <- memberVotes < minvotes 
    votes2dump <- minorityVoteShare <= lop
    N <- sum(!members2dump)
    M <- sum(!votes2dump)

    #Check that legislators used to specify polarities are not deleted, then 
    #adjust them for deleted legislators
    if(any(members2dump[polarityToUse])) {
        cat("\t\tThe following legislators fail minimum vote requirements:\n\t\t",
        rownames(rcObject$votes)[intersect(which(members2dump),polarityToUse)], "\n")
        stop("\t'polarity' is incorrectly specified\n")
    }

    for(i in 1:dims) polarityToUse[i] <- polarityToUse[i]-sum(members2dump[1:polarityToUse[i]])
    if(length(polarityToUse)==1) polarityToUse <- c(polarityToUse,1)

        #Print dumped legislators and dumped votes, then calculate percentages remaining.
        if(all(!members2dump)) {
            cat("\t\tAll members meet minimum vote requirements.\n\n")
        }   else {  
        tempvotes <- rcObject$votes[,!votes2dump]
        if (verbose) {
            cat("\t\tMembers dropped:\n")
            for(i in 1:sum(members2dump))
                cat("\t\t", rownames(rcObject$votes)[which(members2dump)[i]],
                "(voted on only", checkMember(tempvotes)[which(members2dump)[i]],
                "of", dim(tempvotes)[2], "retained votes).\n")
        }
            cat("\t\t...", sum(members2dump) ,"of", dim(tempvotes)[1], "total members dropped.\n\n")
        }

        if(all(!votes2dump)) {
            cat("\t\tAll votes meet minimum lopsidedness requirement.\n\n")
        }  else {
        tempvotes <- rcObject$votes[!members2dump,]
        cat("\t\tVotes dropped:\n")
        if (verbose) {
        for(i in 1:sum(votes2dump))
            cat(sprintf("\t\tNumber %.0f (minority size = %2.1f", which(votes2dump)[i],
            round(100*checkVote(tempvotes)[which(votes2dump)[i]],1)), "%)\n",sep="")
        }
        cat("\t\t...", sum(votes2dump) ,"of", dim(tempvotes)[2], "total votes dropped.\n\n")
        }

    cat("\tRunning Optimal Classification...\n\n")
    flush.console()


    res <- .Fortran("oc",
            #inputs 
            as.integer(rcObject$votes[!members2dump,!votes2dump]),      #columnwise vector of matrix
            as.integer(N),          #rows of votes
            as.integer(M),              #cols of votes
            as.integer(dims),               #number of dimensions
            as.integer(polarityToUse),      #vectors of members to constrain 
                                            #to be positive on each dimension.
            #outputs
            classify=double((N+M)*4),       #Vector of classifications
            fits=double(2),
            volume=double(N),
            idealpoints=double(dims*N),         #colwise vector of ideal point estimates matrix
            midpoints=double(M),                #colwise vector of votes by dims matrix of midpoints 
            normalvectors=double(M*dims),
            eigenvalues=double(N),          #vector of eigenvalues
            exitstatus=integer(1) )                       #Exit status flag (1=OK)

    if(res$exitstatus!=1) stop("\n\n====== oc.f did not execute properly ======\n\n")

#    return(res)
#    can delete the rest
    
    legClassify <- matrix(res$classify[1:(4*N)], nrow=N, ncol=4, byrow=TRUE)
    rcClassify <- matrix(res$classify[(4*N+1):length(res$classify)],
                    nrow=M, ncol=4, byrow=TRUE)
    rcPRE <- apply(cbind(rcClassify[,1]+rcClassify[,3],rcClassify[,2]+rcClassify[,4]),1, min)
    rcPRE <- (rcPRE-rcClassify[,2]-rcClassify[,3])/rcPRE 
    #Legislator PREs omitted because Keith says they are not meaningful

    tempLegis<-cbind(rank(matrix(res$idealpoints[1:N], nrow=N, ncol=1)),
                    legClassify,
                    matrix(res$volume, nrow=N, ncol=1),    
                    matrix(res$idealpoints, nrow=N, ncol=dims))

    tempRC<-cbind(rcClassify,
                  rcPRE,
                  matrix(res$normalvectors,nrow=M,ncol=dims),
                  matrix(res$midpoints,nrow=M,ncol=1))

    legislators<-matrix(NA,rcObject$n,6+1*dims)
    rollcalls<-matrix(NA,rcObject$m,6+1*dims)
    legislators[!members2dump,]<-tempLegis
    rollcalls[!votes2dump,]<-tempRC

    colnames(legislators)<-c("rank","correctYea","wrongYea","wrongNay","correctNay",
            "volume", paste("coord",1:dims,"D",sep=""))
    rownames(legislators)<-rownames(rcObject$votes)
    colnames(rollcalls)<-c("correctYea","wrongYea","wrongNay","correctNay","PRE",paste("normVector",1:dims,"D",sep=""),"midpoints")
    if(!is.null(rcObject$legis.data)) legislators <- cbind(as.data.frame(rcObject$legis.data),legislators)

    OCobject<-list(legislators=legislators, rollcalls=rollcalls, dimensions=dims, eigenvalues=res$eigenvalues,fits=res$fits)
    class(OCobject) <- c("OCobject")
    cat("\n\nOptimal Classification completed successfully.")    
    cat("\nOptimal Classification took", (proc.time()-start)[3],"seconds to execute.\n\n")     
    OCobject

}   #End OC


