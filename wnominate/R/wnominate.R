#.First.lib <- function(lib,pkg) {
#   library.dynam("wnominate",pkg,lib)
#}

# function: qnprob
# Note: Beta and dimweight are unused and only passed
#       to match nomprob function
qnprob <- function(yea,nay,ideal,Beta,dimweight,normal=1) {
    if (normal==1) {
        plink = pnorm
    } 
    else {
        plink = plogis
    }
    cons <- matrix(1,nrow(ideal),1) %*% t( apply(nay*nay,1,sum) - apply(yea*yea,1,sum) )
    discrim <- 2*(yea - nay)
    ystar <-  cons + ideal %*% t(discrim)    
    plink(ystar)
}

generateTestData <- function(legislators=20, rcVotes=100,
            yea=matrix(runif(rcVotes,min=-0.2,max=0.7),nrow=rcVotes),
            nay=matrix(runif(rcVotes,min=-0.7,max=0.2),nrow=rcVotes),
            ideal=matrix(rnorm(legislators),nrow=legislators),
            Beta=15, dimweight=0.5,normal=1, seed = NULL, utility='nominate') {
    if (!is.null(seed)) {
        set.seed(seed) 
    }                           
    if (utility == 'nominate') {
        probfunc <- nomprob        
    }
    else{
       if (utility == 'QN') {
        probfunc <- qnprob
       }    
       else {
          stop("'utility' argument must be 'nominate' or 'QN'")
       }
    }
    
    fakeData<-probfunc(yea,nay,ideal,Beta,dimweight,normal)
    fakeData <- (fakeData > matrix(runif(legislators*rcVotes),legislators,rcVotes))*1  
    fakeData[fakeData==0]<-6
    lopsided<-rep(FALSE,rcVotes)
    legis.names<-paste("Legislator", 1:legislators, sep="")
    icpsrState <- sample(c(10,56,63,94), size=legislators,replace=TRUE)
    state<-icpsrState
    state[state==10]<-"Alaska"
    state[state==56]<-"California"
    state[state==63]<-"Ohio"
    state[state==94]<-"Texas"
    cd<-sample(1:435, size=legislators, replace=TRUE)
    icpsrLegis <- sample(1:100, size=legislators, replace=TRUE)
    partyCode <- sample(c(200,100), size=legislators,replace=TRUE)
    partyName <- partyCode
    partyName[partyName==200]<-"Republican"
    partyName[partyName==100]<-"Democrat"
    codes=list(yea=1:3,nay=4:6,notInLegis=0,missing=7:9) 
   
    legdata<-data.frame(state=state, icpsrState=icpsrState, cd=cd,
                        icpsrLegis=icpsrLegis, party=partyName, party=partyCode)
    rownames(legdata) <- legis.names
    rownames(fakeData) <- legis.names
    colnames(fakeData) <- paste("V",1:dim(fakeData)[2],sep="")
    rcObject<-list(votes=fakeData, codes=codes, n=legislators, m=rcVotes,
        lopsided=lopsided, legis.data=legdata, vote.data=NULL, desc=NULL,source=NULL)
    
    class(rcObject) <- c("rollcall")
    rcObject
}


nomprob <- function(yea, nay, ideal, Beta, dimweight, normal=1) {
    res <- .C("nomprob",
                 as.double(t(yea)),
                 as.double(t(nay)),
                 as.double(t(ideal)),
                 as.double(Beta),
                 as.double(dimweight),
                 as.integer(dim(yea)[1]),
                 as.integer(dim(ideal)[1]),
                 as.integer(dim(ideal)[2]),
                 yeaProb = double( dim(ideal)[1]*dim(yea)[1]),
         as.integer(normal))
    matrix(res$yeaProb,nrow=dim(ideal)[1],byrow=FALSE)
}



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


wnominate <- function(rcObject, ubeta=15, uweights=0.5, dims=2, minvotes=20,
                        lop=0.025, trials=3, polarity, verbose=FALSE) {

    cat("\nPreparing to run W-NOMINATE...\n\n\tChecking data...\n\n")
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
    if(trials<1) stop("'trials' cannot be negative\n") 

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
        }   
        else {  
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

    cat("\tRunning W-NOMINATE...\n\n")
    flush.console()

    fakedims<-max(dims,2)
   
    res <- .Fortran("wnom",
          #inputs 
              as.integer(rcObject$votes[!members2dump,!votes2dump]),
              as.integer(N), 
              as.integer(M),  
              as.integer(fakedims),       
              as.integer(trials),     
              as.integer(polarityToUse),   
          as.integer(dims),#

              #input/output 
              ubeta=as.single(ubeta),       
              uweights=as.single(rep(uweights,dims)),
              #outputs
              classify=single((N+M)*4),
              fits=single(3*dims),
              gmp=single(N+M), 
              idealpoints=single(dims*N), 
              covariances=single(((dims*(dims+1))/2)*N), 
              midpoints=single(dims*M),      
              spreads=single(dims*M),       
              eigenvalues=single(N), 
              exitstatus=integer(1),

        #wnom() call continued, with dynamic memory allocation in R
                integer(N*M),             #LDATA2 
                single(N+M),              #GMPGMP(NUMMEMBERS+NUMVOTES)
                single(3*fakedims),           #XFITS(3*DIMS)
                single(N*N),              #ZMAT2(NUMMEMBERS,NUMMEMBERS)
                single(N),                #WVEC2(NUMMEMBERS)
                single(N*N),              #DSTAR(NUMMEMBERS,NUMMEMBERS)
                single(fakedims*N),           #XDATA(NUMMEMBERS,DIMS) 
                single(N),                #XXX(NUMMEMBERS)
                single(fakedims*N),           #XDATA3(NUMMEMBERS,DIMS) 
                single(fakedims*M),           #ZMID(NUMVOTES,DIMS) 
                single(fakedims*M),           #DYN(NUMVOTES,DIMS)
                single(4*N),              #XSAVE(NUMMEMBERS,2,2)
                single(4*M),              #ZSAVE(NUMVOTES,2,2)
                single(2*M),              #CSAVE(NUMVOTES,2)
                integer(M),               #KAV(NUMVOTES)
                integer(M),               #KAY(NUMVOTES)
                integer(M),               #KAN(NUMVOTES)
                single(N),                #XD(NUMMEMBERS)
                integer(N),               #ISENS(NUMMEMBERS)
                integer(M),               #LERIC(NUMVOTES) 
                single(M*N*2),            #PSI(NUMMEMBERS,NUMVOTES,2)
                single(N*2*fakedims),         #XMEANX(2*NUMMEMBERS,DIMS) 
                single(N*3*fakedims),         #STDDEV(3*NUMMEMBERS,DIMS) 
                single(N*2*fakedims*fakedims),    #COVX(2*NUMMEMBERS,DIMS,DIMS)
                single(N*2*fakedims*fakedims),    #COVX2(2*NUMMEMBERS,DIMS,DIMS) 
                single((M+N+111)*4),      #KPJP(NUMMEMBERS+NUMVOTES+111,4) 
                single(M+N+111),          #YBIGL(NUMMEMBERS+NUMVOTES+111) 
                single(M+N+111),          #YYBIGL(NUMMEMBERS+NUMVOTES+111) 
                single(M*2*fakedims),         #STDDVZ(NUMVOTES,2,DIMS) 
                single(N*fakedims),           #STDDVX(NUMMEMBERS,DIMS)
                single(M+N+111),          #LMO(NUMMEMBERS+NUMVOTES+111) 
                single(M*N),              #POOLE(NUMMEMBERS,NUMVOTES) 
                single(N*fakedims),           #TRUEX(NUMMEMBERS,DIMS)
                single(N*fakedims),           #TRUEX2(NUMMEMBERS,DIMS) 
                single(M*fakedims),           #TRUEZMID(NUMVOTES,DIMS) 
                single(M*fakedims),           #TRUEDYN(NUMVOTES,DIMS)
                single(M*N),              #PROBMAT(NUMMEMBERS,NUMVOTES) 
                single(N*fakedims),           #XTARGET(NUMMEMBERS,DIMS) 
                single(N*fakedims),           #XMAT0(NUMMEMBERS,DIMS)
                single(N*fakedims),           #XSAVE2(NUMMEMBERS,DIMS) 
                single(N*fakedims),           #XSAVE3(NUMMEMBERS,DIMS) 
                single(N*fakedims))           #XMAT(NUMMEMBERS,DIMS)
 
   
    if(res$exitstatus!=1) stop("\n\n====== wnominate9707.f did not execute properly ======\n\n")

    #Roll wnom() output into a nomObject
    #classify is a rowwise matrix that begins with legislators first
    legClassify <- matrix(res$classify[1:(4*res[[2]])], nrow=res[[2]], ncol=4, byrow=TRUE)
    rcClassify <- matrix(res$classify[(4*res[[2]]+1):length(res$classify)],
                    nrow=res[[3]], ncol=4, byrow=TRUE)

    ## These equations preceeded wnomiante 0.96 and are incorrect, according to May 29, 2011 email 
    #legPRE <- apply(cbind(legClassify[,1]+legClassify[,3],legClassify[,2]+legClassify[,4]),1, min)
    #legPRE <- (legPRE-legClassify[,2]-legClassify[,3])/legPRE

    rcPRE <- apply(cbind(rcClassify[,1]+rcClassify[,3],rcClassify[,2]+rcClassify[,4]),1, min)
    rcPRE <- (rcPRE-rcClassify[,2]-rcClassify[,3])/rcPRE

    ## Legislator Correct Classification replaces PRE from May 29, 2011
    ## CC = (correctYea+correctNay)/(correctYea+wrongYea+wrongNay+correctNay)
    legCC <- (legClassify[,1] + legClassify[,4])/apply(legClassify,1,sum)
    
    tempRC<-cbind(rcClassify,
                    res$gmp[(res[[2]]+1):length(res$gmp)],
                    rcPRE,
                    matrix(res$spread,nrow=res[[3]],ncol=dims),
                    matrix(res$midpoint,nrow=res[[3]],ncol=dims))

    tempLegis<-cbind(legClassify,
                    res$gmp[1:res[[2]]],
                    legCC,
                    matrix(res$idealpoints, nrow=res[[2]], ncol=dims),
                    matrix(res$covariance,nrow=res[[2]],ncol=dims+choose(dims,2)))

    legislators<-matrix(NA,rcObject$n,6+2*dims+choose(dims,2))
    rollcalls<-matrix(NA,rcObject$m,6+2*dims)
    legislators[!members2dump,]<-tempLegis
    rollcalls[!votes2dump,]<-tempRC

    if(choose(dims,2)==0){
        colnames(legislators)<-c("correctYea", "wrongYea", "wrongNay", "correctNay",
            "GMP", "CC", paste("coord",1:dims,"D",sep=""),
            paste("se",1:dims,"D",sep=""))
    } else {
        colnames(legislators)<-c("correctYea", "wrongYea", "wrongNay", "correctNay",
            "GMP", "CC", paste("coord",1:dims,"D",sep=""), paste("se",1:dims,
            "D",sep=""), paste("corr.",1:choose(dims,2),sep=""))
    }

    rownames(legislators)<-rownames(rcObject$votes)

    colnames(rollcalls)<-c("correctYea", "wrongYea", "wrongNay", "correctNay", "GMP", "PRE",
        paste("spread",1:dims,"D",sep=""), paste("midpoint",1:dims,"D",sep=""))

    if(!is.null(rcObject$legis.data)) legislators <- cbind(as.data.frame(rcObject$legis.data),legislators)

    fits <- as.numeric(res$fits)
    names(fits) <- c(paste("correctclass",1:dims,"D",sep=""),
paste("apre",1:dims,"D",sep=""),
paste("gmp",1:dims,"D",sep=""))

    nomObject<-list(legislators=as.data.frame(legislators), rollcalls=as.data.frame(rollcalls), dimensions=dims, eigenvalues=res$eigenvalues,beta=res$ubeta,weights=res$uweights,fits=fits)
    class(nomObject) <- c("nomObject")
    cat("W-NOMINATE estimation completed successfully.")    
    cat("\nW-NOMINATE took", (proc.time()-start)[3],"seconds to execute.\n\n")     
    nomObject
    
}   #end wnominate()
