#
# This is the CGHcall function from CGHcall package (GPL License).
# 
# Calls aberrations for array CGH data using a six state mixture model.
#
# @title Calling aberrations for array CGH tumor profiles.
# @param inputSegmented a list. See Details for more informations.
# @param prior Options are all, not all, or auto. See details for more information.
# @param nclass  The number of levels to be used for calling. Either 3 (loss, normal, gain), 4 (including amplifications), 5 (including double deletions).
# @param organism  Either human or other. This is only used for chromosome arm information when prior is set to all or auto (and samplesize > 20).
# @param cellularity  A vector of cellularities ranging from 0 to 1 to define the contamination of your sample with healthy cells (1 = no contamination). 
# See details for more information.
# @param robustsig  Options are yes or no. yes enforces a lower bound on the standard deviation of the normal segments.
# @param nsegfit Maximum number of segments used for fitting the mixture model. Posterior probabilities are computed for all segments0.
# @param maxnumseg  Maximum number of segments per profile used for fitting the model
# @param minlsforfit	Minimum length of the segment (in Mb) to be used for fitting the model
# @param build	Build of Humane Genome. Either GRCh37, GRCh36, GRCh35 or GRCh34.
# @param ncpus	Number of cpus used for parallel calling. Has a large effect on computing time. ncpus larger than 1 requires package snowfall.
# @param verbose if TRUE print some informations.
# 
# @details 
# inputSegmented is a list containing:
#  \describe{
#   \item{copynumber}{A matrix. Each column contains a signal of copynumber for a profile. Each row corresponds to a genomic position of a probe.}
#   \item{segmented}{A matrix of the same size as copynumber. It contains the segmented signals.}
#   \item{chromosome}{A vector of length nrow(copynumber) containing the number of the chromosome for each position.}
#   \item{startPos}{A vector of length nrow(copynumber) containing the starting genomic position of each probe.}
#   \item{featureNames}{A vector of length nrow(copynumber) containing the names of each probes.}
#   \item{sampleNames}{A vector of length ncol(copynumber) containing the names of each profiles.}
# }
# 
# 
# Please read the article and the supplementary information for detailed information on the algorithm. 
# The parameter prior states how the data is used to determine the prior probabilities.
#  When set to all, the probabilities are determined using the entire genome of each sample.
#  When set to not all probabilites are determined per chromosome for each sample when organism is set to other or per chromosome arm when organism is human.
#  The chromosome arm information is taken from the March 2006 version of the UCSC database.
#  When prior is set to auto, the way probabilities are determined depends on the sample size.
#  The entire genome is used when the sample size is smaller than 20, otherwise chromosome (arm) information is used. 
#  Please note that CGHcall uses information from all input data to determine the aberration probabilities.
#  When for example triploid or tetraploid tumors are observed, we advise to run CGHcall separately on those (groups of) samples. 
#  Note that robustsig = yes enforces the sd corresponding to the normal segments to be at least half times the pooled gain/loss sd. 
#  Use of nsegfit significantly lower computing time with respect to previous CGHcall versions without much accuracy loss. 
#  Moreover, maxnumseg decreases the impact on the results of profiles with inferior segmentation results. 
#  Finally, minlsforfit decreases the impact of very small aberations (potentially CNVs rather than CNAs) on the fit of the model. 
#  Note that always a result for all segments is produced. IN MOST CASES, CGHcall SHOULD BE FOLLOWED BY FUNCTION ExpandCGHcall.
#  
#  
# @author Sjoerd Vosse, Mark van de Wiel, Ilari Scheinin
# @references 
# Mark A. van de Wiel, Kyung In Kim, Sjoerd J. Vosse, Wessel N. van Wieringen, Saskia M. Wilting and Bauke Ylstra. CGHcall: calling aberrations for array CGH tumor profiles. Bioinformatics, 23, 892-894.
#
# @export
CGHcall <- function(inputSegmented, prior="auto", nclass=5, organism="human",cellularity=1, robustsig="yes",nsegfit=3000,maxnumseg=100,minlsforfit=0.5, build="GRCh37",ncpus=1,verbose=TRUE) 
{
   gc() #memory usage in Mb
    
 
    ## Author: Mark van de Wiel
    ## Maintainer: Mark van de Wiel
    ## Email: mark.vdwiel@vumc.nl
 
    
    ## Changes since previous version
    ## - input is of class cgh
   
    timeStarted <- proc.time()

    bw <- round(20*minlsforfit*nrow(inputSegmented$copynumber)/44000) 
    ncolscl     <- ncol(inputSegmented$copynumber)
    
    if (length(cellularity) < ncolscl) cellularity <- rep(cellularity, ncolscl) else cellularity <- cellularity[1:ncolscl]
 
    
    whna <- !is.na(inputSegmented$chromosome) & !is.na(inputSegmented$featureNames)
    datmat <- cbind(inputSegmented$copynumber, inputSegmented$segmented)[whna,]
    #datmat <- cbind(c(-1.2,2.1,-0.3,-0.4),c(0,0,0,0));chr<-c(1,1,1,1);nc<-1
    posit       <- inputSegmented$startPos[whna]
    chr         <- inputSegmented$chromosome[whna]
    naam        <- inputSegmented$featureNames[whna]
    rm(whna,inputSegmented);
    gc()
    nc    <- ncol(datmat)/2 
    nclone <-nrow(datmat)
    datareg  <- .MakeData(datmat[,-(1:nc),drop=FALSE],chr)  
   
    #rm(combined); 
    
    gc() #added 24/11/2009
    dataprob    <- data.frame(naam, chr, posit)
    
 
    
    ## Determine method for prior probabilities
    if (prior == 'auto') {
        if (nc >  20) prior <- "not all";
        if (nc <= 20) prior <- "all";
    }
    
    ### Convert to chromosome arm if neccessary
    if (prior != "all" && organism == "human") {
       # temp            <- datareg[[1]];
        dataprob    <- .convertChromosomeToArm(dataprob,build); 
       # datareg[[1]]    <- .convertChromosomeToArm(temp);
        chr <- dataprob[,2] #BUG;repaired 22/06/09
    }

    if(verbose)
      cat("\n EM algorithm started ... \n");
    
    dat         <- c()
    datsm       <- c()
    datall         <- c()
    datsmall       <- c()
    regions     <- c()
    regionsall  <- c()
    regionsdat  <- c()
    regionsdatall  <- c()
    nclones     <- length(chr)
     profile     <- c()
     profileall <- c()
    countreg    <- 0
    newregions <- list()
#    nclonesj <- 0
    thr <- 1.5
   
    #here is the selection of regions used for fitting. First select regions long enough and with smoothed signal <= thr. Then, maximize 
    #regions per prof to maxseg
    #nc<-2
    gc() 
    datall <- as.vector(datmat[,1:nc]) #changed 16/7/10
    datsmall      <- as.vector(datmat[,-(1:nc)])
    rm(datmat);gc();
  
    for (j in (1:nc)) {
        regions1    <- datareg[[j]]
        nreg1       <- nrow(regions1)
        profileall     <- c(profileall, rep(j,nreg1))
        regionsall     <- rbind(regionsall, regions1)
        toaddall       <- (j-1)*c(nclones, nclones)
        regionsdatall  <- rbind(regionsdatall, regions1[,-3]+toaddall)

        segl        <- regions1[,2]-regions1[,1]+1
        whseg       <- which(segl>=bw)
        regions2    <- regions1[whseg,,drop=FALSE] #selects regions long enough to participate in fitting
        regions2    <- regions2[which(abs(regions2[,3]) <= thr),,drop=FALSE] #selects regions with smoothed signal small enough
        nreg2       <- nrow(regions2)
        countreg <- countreg + min(maxnumseg,nreg2)
        newregions[[j]] <- regions2
        }
        fracforfit <- min(1,nsegfit/countreg)
        
        
        for (j in (1:nc)) {
        regions2 <- newregions[[j]]
        nreg2       <- nrow(regions2)
        ordsig      <- order(regions2[,3])
        regions2b   <- data.frame(regions2,genord = 1:nreg2)
        regions2c   <- regions2b[ordsig,,drop=FALSE]
        if(nreg2 >= maxnumseg | fracforfit<1){
            takeseg     <- floor(seq(1,nreg2,length.out=round(min(nreg2,maxnumseg)*fracforfit)))
            regions2c <- regions2c[takeseg,,drop=FALSE]
            }
        regions2d <- regions2c[order(regions2c[,4]),,drop=FALSE]
        
        nreg2       <- nrow(regions2d)
        profile    <- c(profile, rep(j,nreg2))
        regions     <- rbind(regions, regions2d)
        #takerow <- c();for(k in 1:nrow(regions2d)) {takerow <- c(takerow,regions2d[k,1]:regions2d[k,2])}
#        dat <- c(dat,datareg[[1]][takerow,(3+j)])
#        datsm <- c(datsm,datareg[[1]][takerow,(3+nc+j)])
        toadd       <- (j-1)*c(nclones, nclones)
        regionsdat  <- rbind(regionsdat, regions2d[,-(3:4)]+toadd) #delete smrat value and index
        #        nclonesj <- sum((regions2d[,2]+1)-regions2d[,1]) 
    }
    gc()

    takechr     <- function(reg, chrom) {
        chrom[reg[1]]
    }
    
    chrarmreg   <- as.vector(apply(as.matrix(regions), 1, takechr, chrom=chr))
    
    nreg        <- nrow(regionsdat)
    nregall     <- nrow(regionsdatall)
    if(verbose)
    {
      cat(paste("Total number of segments present in the data:",nregall,"\n"))
      cat(paste("Number of segments used for fitting the model:",nreg,"\n"))      
    }

    allnc       <- sapply(1:nreg, .countcl, regionsdat=regionsdat)
    allsum      <- sapply(1:nreg, .sumreg, dat=datall, regionsdat=regionsdat)
    allsumsq    <- sapply(1:nreg, .sumsqreg, dat=datall, regionsdat=regionsdat)
    varregall   <- sapply(1:nreg, .varregtimescount, counts=allnc, dat=datall, regionsdat=regionsdat)
    lev         <- 1:nc
    varprofnc   <- cbind(varregall, profile, allnc)
    varprof     <- sapply(lev, .varproffun, vcnmat=varprofnc, profile=profile)
    
    gc()

    selk <- function(k, varprof,prof) {  #changed 22/06/09
        profk <- prof[k]
        return(max(0.001,varprof[profk]))
    }
    varprofall <- sapply(1:nreg, selk, varprof=varprof,prof=profile)
    
    selkcell <- function(k, cellprof,prof) {  #changed 22/06/09
        profk <- prof[k]
        return(max(0.001,cellprof[profk]))
    }

    allmean     <- allsum/allnc
    allmeanl    <- allmean[allmean < -0.15 & allmean > -0.7]
    allmeang    <- allmean[allmean>0.15 & allmean<0.7]
    allmean0    <- allmean[abs(allmean) <= 0.15]
    sdlst       <- if(length(allmeanl) <= 1) 0.01 else mad(allmeanl)
    meanlst     <- if(length(allmeanl) <= 0) -0.3 else mean(allmeanl)
    sdgst       <- if(length(allmeang) <= 1) 0.01 else mad(allmeang)
    meangst     <- if(length(allmeang) <= 0) 0.3 else mean(allmeang)
    sd0st       <- if(length(allmean0) <= 1) 0.01 else mad(allmean0)
    #sd0st       <- if(length(allmean0) <= 1) 0.01 else sd(allmean0)

    profchrom   <- chrarmreg
    if(robustsig=="no") {
        bstart      <- c(-log(0.2), -log(-(meanlst+0.1)), sqrt(-log(0.5)), -log(meangst-0.1), -log(0.2), 0.0001, sdlst, sd0st, sdgst, 0.0001, 0.0001)
    } else {
        bstart      <- c(-log(0.2), -log(-(meanlst+0.1)), sqrt(-log(0.5)), -log(meangst-0.1), -log(0.2), 0.0001, sdlst, sd0st+max(0,sqrt(sdgst^2/8+sdlst^2/8)-sd0st), sdgst, 0.0001, 0.0001) #robust option 
    }
    
    allcell <- sapply(1:nreg, selkcell, cellprof=cellularity,prof=profile)

    mus         <- c(-0.10-exp(-bstart[2])-0.3 - exp(-bstart[1]), -0.10-exp(-bstart[2]), -0.05+0.1*exp(-(bstart[3])^2), 0.1 + exp(-bstart[4]), 2*(0.10 + exp(-bstart[4]))+0.3+exp(-bstart[5]))
    priorp      <- matrix(rep(c(0.01, 0.09, 0.8, 0.08, 0.01, 0.01), nreg), ncol=6, byrow=TRUE)
    alpha0      <- .alpha0all(nreg, profchrom, priorp, bstart, varprofall, allsum, allsumsq, allnc,robustsig,allcell,prior)
    
    maxiter     <- 10
    stop        <- 0
    iter        <- 1
    thrpara     <- 0.015 #changed 15/6/09
    
    gc()
    if(ncpus>1){
            trysnow <- try(library(snowfall))  
            if(class(trysnow)=="try-error") {
            print("You should install 'snowfall' if you want to use parallel computing...")
            print("Escaping to slower computation on 1 cpu")
            ncpus <- 1
            }
        }
        
    while (stop == 0 & iter <= maxiter) {
        gc()#print(gc())
        #cat("Calling iteration", iter, ":\n")
        posterior0  <- sapply(1:nreg, .posteriorp, priorp=alpha0, pm=bstart, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,allcell=allcell,robustsig=robustsig)
        likprev     <- .totallik(bstart, nreg=nreg, posteriorprev=posterior0, alphaprev=alpha0, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,allcell=allcell,robustsig=robustsig)
        

        if(ncpus>1){
          snowfall::sfInit(cpus=ncpus,parallel=T)
          snowfall::sfLibrary("MPAgenomics")
          snowfall::sfExport("bstart","nreg", "posterior0", "alpha0","varprofall", "allsum", "allsumsq", "allnc","allcell","robustsig","ncpus")
        }
        
        pmt<-proc.time()
        optres      <- optim(bstart,.totallik, nreg=nreg, posteriorprev=posterior0, alphaprev=alpha0, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,allcell=allcell,robustsig=robustsig,ncpus=ncpus)
        toptim<-proc.time()-pmt
               
        if(ncpus>1){
          snowfall::sfRemoveAll()
          snowfall::sfStop()            
        }
        
        bprev       <- bstart
        bstart      <- optres$par       
        alpha0      <- .alpha0all(nreg, profchrom, alpha0, bstart, varprofall, allsum, allsumsq, allnc, robustsig,allcell, prior)
        rl          <- .reallik4(nreg, alpha0, bstart, varprofall, allsum, allsumsq, allnc, robustsig,allcell)
       
        llmin       <- optres$value
        musprev     <- c(-0.10-exp(-bprev[2])-0.3 - exp(-bprev[1]), -0.10-exp(-bprev[2]), -0.05+0.1*exp(-(bprev[3])^2), 0.1 + exp(-bprev[4]),(log2(2)/log2(1.5))*(0.1 + exp(-bprev[4])), (log2(2)/log2(1.5))*(0.10 + exp(-bprev[4]))+0.3+exp(-bprev[5]))
        mus         <- c(-0.10-exp(-bstart[2])-0.3 - exp(-bstart[1]), -0.10-exp(-bstart[2]), -0.05+0.1*exp(-(bstart[3])^2), 0.1 + exp(-bstart[4]),(log2(2)/log2(1.5))*(0.1 + exp(-bstart[4])),(log2(2)/log2(1.5))*(0.10 + exp(-bstart[4]))+0.3+exp(-bstart[5]))
        param       <- c(mus, bstart[6], bstart[7], bstart[8], bstart[9], bstart[10],bstart[11])
        paramprev   <- c(musprev, bprev[6], bprev[7], bprev[8], bprev[9], bprev[10],bprev[11])
        
        #print("optim results")
        #print(paste("time:",round(toptim[3])))
        #print(paste("minimum:",optres$value))
        if(robustsig=="no") {
            printmat <- matrix(c(j,rl,mus,sqrt(bstart[7]^2+(bstart[6])^2 + 0.0001),bstart[7],bstart[8],bstart[9],sqrt((bstart[10])^2+(bstart[9])^2 + 0.0001),sqrt((bstart[11])^2 + (bstart[10])^2+ (bstart[9])^2+ 0.0001)),nrow=1)
            colnames(printmat) <- c("j","rl","mudl","musl","mun","mug","mudg","mua","sddl","sdsl","sdn","sdg","sddg","sda")
            #print(printmat)
        } else {
            printmat<-matrix(c(j,rl,mus,sqrt(bstart[7]^2+(bstart[6])^2 + 0.0001),bstart[7],sqrt(1/8*((bstart[7])^2 + 0.0001)+1/8*((bstart[9])^2 + 0.0001)+((bstart[8]^2)+0.0001)),bstart[9],sqrt((bstart[10])^2+(bstart[9])^2 + 0.0001),sqrt((bstart[11])^2 + (bstart[10])^2 + (bstart[9])^2 + 0.0001)),nrow=1)  #robust option
            colnames(printmat) <- c("j","rl","mudl","musl","mun","mug","mudg","mua","sddl","sdsl","sdn","sdg","sddg","sda")
            #print(printmat)
        }
        
        
        if (max(abs(param[c(2:4,8)]-paramprev[c(2:4,8)])) <= thrpara) stop <- 1 #changed 15/6/09; also sdn should converge
        iter <- iter+1
        }
        mudl<-printmat[3];mul <- printmat[4];mug<-printmat[6];mudg<- printmat[7];mua<-printmat[8];sdl<-printmat[10];sdg<-printmat[12] #added 17/12
        if(verbose)
          cat("EM algorithm done ...\n")
        best            <- bstart
    
    #now start computing posteriors for ALL data
    if(verbose)
      cat("Computing posterior probabilities for all segments ...\n")
    #nregall        <- nrow(regionsdatall)
    allncall       <- sapply(1:nregall, .countcl, regionsdat=regionsdatall)
    allsumall      <- sapply(1:nregall, .sumreg, dat=datall, regionsdat=regionsdatall)
    allmeanall     <- allsumall/allncall  
    allsumsqall    <- sapply(1:nregall, .sumsqreg, dat=datall, regionsdat=regionsdatall)
    varregallall   <- sapply(1:nregall, .varregtimescount, counts=allncall, dat=datall, regionsdat=regionsdatall)
    lev             <- 1:nc
    varprofncall   <- cbind(varregallall, profileall, allncall)
    varprof_all     <- sapply(lev,.varproffun, vcnmat=varprofncall, profile=profileall)
    varprof_allall <- sapply(1:nregall, selk, varprof=varprof_all,prof=profileall)
    gc()
    
    allcellall     <- sapply(1:nregall, selkcell, cellprof=cellularity,prof=profileall)
    profchromall <- as.vector(apply(as.matrix(regionsall), 1, takechr, chrom=chr))
       
    if(prior!="all"){
    chrarm_alpha0 <- cbind(chrarmreg,alpha0)
    chrarm_alpha0uni <- unique(chrarm_alpha0) 
    chrarmregall   <- as.vector(apply(as.matrix(regionsall), 1, takechr, chrom=chr))
    chrarmreguni <- unique(chrarmreg)
    chrarmreguniall <- unique(chrarmregall)
    differ <- setdiff(chrarmreguniall,chrarmreguni)
    if(length(differ)==0){
    alpha0_all      <- t(sapply(chrarmregall,function(x){chrarm_alpha0uni[chrarm_alpha0uni[,1]==x,-1]})) 
    } else {
    mn <- as.vector(apply(chrarm_alpha0uni[,-1],2,mean))
    alpha0_all      <- t(sapply(chrarmregall,function(x){if(!is.element(x,differ)) chrarm_alpha0uni[chrarm_alpha0uni[,1]==x,-1] else mn}))
    } 
    alpha0_all     <- .alpha0all(nregall, profchromall, alpha0_all, bstart, varprof_allall, allsumall, allsumsqall, allncall, robustsig,allcellall, prior)
    } else alpha0_all <- matrix(rep(alpha0[1,],nregall),byrow=T,nrow=nregall)
    
    posteriorfin    <- t(sapply(1:nregall, .posteriorp, priorp=alpha0_all, pm=best, varprofall=varprof_allall, allsum=allsumall, allsumsq=allsumsqall, allnc=allncall, allcell=allcellall, robustsig=robustsig))
    
    overrule <- function(id,mnall){
    #id <- 17;mnall<-allmeanall
    toret <- posteriorfin[id,]   
    mn <- mnall[id]
    row1 <- toret
    if((mn > mug+sdg) & (which.max(row1)<=3)) {
        sumr <- sum(c(row1[4],row1[5],row1[6])); #added 17/12/10; if max prob belongs to non-gain state, overrule; redistribute probability
        if(sumr > 10^(-10)) toret<- c(0,0,0,row1[4]/sumr,row1[5]/sumr,row1[6]/sumr) else {
            toret<- c(0,0,0,0,0,0)
            wm<-which.min(c(abs(mn-mug),abs(mn-mudg),abs(mn-mua)))+3
            toret[wm]<-1 #assigns probability 1 to status with mode closest to the mean segment value; only used if probability calculation fails.
        }
    }
    if((mn < mul-sdl) & (which.max(row1)>=3)) {
        sumr <- sum(c(row1[1],row1[2])); #added 17/12/10; if max prob belongs to non-loss state, overrule
         if(sumr > 10^(-10)) toret<- c(row1[1]/sumr,row1[2]/sumr,0,0,0,0) else {
            toret<- c(0,0,0,0,0,0)
            wm<-which.min(c(abs(mn-mudl),abs(mn-mul)))
            toret[wm]<-1 #assigns probability 1 to status with mode closest to the mean segment value; only used if probability calculation fails.
        } 
    }
    
     if(mn < 0) {
        sumr <- sum(c(toret[1],toret[2],toret[3])); #added 20/04/12; 
        toret<- c(toret[1]/sumr,toret[2]/sumr,toret[3]/sumr,0,0,0) 
    }
    
    if(mn >= 0) {
        sumr <- sum(c(toret[3],toret[4],toret[5],toret[6])); #added 20/04/12; 
        toret<- c(0,0,toret[3],toret[4],toret[5],toret[6])/sumr 
    }
    return(toret)
    }
    
    posteriorfin <- t(sapply(1:nregall,overrule,mnall=allmeanall)) #OVERRULE FUNCTION ADD 17/12/10 TO DEAL WITH EXTREMES

    
    #amporgain <- function(row) {
#        if (row[6] >= 0.5) return(c(row[1]+row[2], row[3], row[4]+row[5], row[6]))
#        else return(c(row[1]+row[2], row[3], row[4]+row[5], row[6]))
#    }
   
    if (nclass == 3) posteriorfin1 <- cbind(posteriorfin[,1]+posteriorfin[,2], posteriorfin[,3], posteriorfin[,4]+posteriorfin[,5]+posteriorfin[,6])
    if (nclass == 4) posteriorfin1 <- cbind(posteriorfin[,1]+posteriorfin[,2], posteriorfin[,3], posteriorfin[,4]+posteriorfin[,5],posteriorfin[,6])
    if (nclass == 5) posteriorfin1 <- cbind(posteriorfin[,1],posteriorfin[,2], posteriorfin[,3], posteriorfin[,4]+posteriorfin[,5],posteriorfin[,6])
   

    posteriorfin2   <- cbind(profileall, posteriorfin1)
    rm(posteriorfin,posteriorfin1);gc() #added24/11/2009
    regionsprof     <- cbind(profileall, regionsall)
    listcall <- list(posteriorfin2,nclone,nc,nclass,regionsprof,params = printmat[,-c(1,2)],cellularity = cellularity)
    #save(listcall,file="listcall.Rdata")
    gc()
    timeFinished <- round((proc.time() - timeStarted)[1] / 60)
    if(verbose)
      cat("Total time:", timeFinished, "minutes\n")
   
    return(listcall)
}


##################################################################################################   

#'
#' Launch the process of segmentation labeling. This function uses functions from CGHcall package developped by Sjoerd Vosse, Mark van de Wiel and Ilari Scheinin. See the CGHcall package for more details.
#'
#' @title Calling aberrations in segmented copy-number signal.
#' 
#' @param segmentData A list (see details).
#' @param nclass The number of levels to be used for calling. Either 3 (loss, normal, gain), 4 (including amplifications), 5 (including double deletions).
#' @param cellularity Proportion of tumor cells in the sample ranging from 0 to 1 (default=1). Reflects the contamination of the sample with healthy cells (1 = no contamination).
#' @param verbose If TRUE, print some details.
#' @param ... other options of CGHcall functions
#'
#' @return A list with the same element as segmentData list and
#' \describe{
#'   \item{calls}{A matrix, of the same size as segmentData$copynumber matrix, containing the label of each point.
#'   -2=double loss, -1=loss, 0=normal, 1=gain, 2=amplification.}
#'   \item{segment}{A data.frame that summarizes the different segments found.}
#'   \item{probdloss}{(if CGHcall was run with nclass=5) A matrix of the same size as segmentData$copynumber matrix. It contains the probability for each segmented copynumber to be a double loss.}
#'   \item{probloss}{A matrix of the same size as segmentData$copynumber matrix. It contains the probability for each segment to be a loss.}
#'   \item{probdnorm}{A matrix of the same size as segmentData$copynumber matrix. It contains the probability for each segment to be normal.}
#'   \item{probdgain}{A matrix of the same size as segmentData$copynumber matrix. It contains the probability for each segment to be a gain.}
#'   \item{probdamp}{(if CGHcall was run with nclass=4 or 5) A matrix of the same size as segmentData$copynumber matrix. It contains the probability for each segment to be an amplification.}
#' }
#'
#'
#' @details
#' segmentData is a list containing:
#'  \describe{
#'   \item{copynumber}{A matrix. Each column contains a signal of copynumber for a profile. Each row corresponds to a genomic position of a probe.}
#'   \item{segmented}{A matrix of the same size as copynumber. It contains the segmented signals.}
#'   \item{chromosome}{A vector of length nrow(copynumber) containing the studied chromosome (number) for each position.}
#'   \item{startPos}{A vector of length nrow(copynumber) containing the starting genomic position of each probe.}
#'   \item{featureNames}{A vector of length nrow(copynumber) containing the names of each probe.}
#'   \item{sampleNames}{A vector of length ncol(copynumber) containing the names of each profile.}
#' }
#' 
#'
#' @author Quentin Grimonprez
#' @export
#'  
callingProcess=function(segmentData,nclass=5,cellularity=1,verbose=TRUE,...)
{
  ####creating a data.frame that sum up the different segment
  segmentMean=apply(segmentData$segmented,2,FUN=function(x){list(unique(x))})#different mean in the signal
  #starting index of the segment
  indexStart=lapply(1:ncol(segmentData$segmented),FUN=function(x){unlist(lapply(segmentMean[[x]],FUN=function(y){match(y,segmentData$segmented[,x])}))})
  #ending index of the segment
  indexEnd=lapply(indexStart,FUN=function(x){c(x[-1]-1,nrow(segmentData$segmented))})  
  #number of the chromosome for each segment
  chrom=unlist(lapply(1:ncol(segmentData$segmented),FUN=function(x){segmentData$chromosome[indexStart[[x]]]}))
  #starting position of the segment
  chromStart=unlist(lapply(1:ncol(segmentData$segmented),FUN=function(x){segmentData$startPos[indexStart[[x]]]}))
  #ending position of the segment
  chromEnd=unlist(lapply(1:ncol(segmentData$segmented),FUN=function(x){segmentData$startPos[indexEnd[[x]]]}))
  #means of the segment
  means=unlist(lapply(1:ncol(segmentData$segmented),FUN=function(x){segmentData$segmented[indexStart[[x]],x]}))
  
  
  ##calling process
  cat("Calling process...")
  #normalization for better set the zero level
  postseg=postsegnormalize(segmentData)
  #launch calling
  calls=CGHcall(postseg, nclass=nclass,cellularity=cellularity,verbose=verbose,...) 
  #extract results
  res=ExpandCGHcall(calls,postseg,verbose=verbose)
  cat(" OK\n") 
 
  #convert the calls
  res$calls[res$calls==-2]="double loss"
  res$calls[res$calls==-1]="loss"
  res$calls[res$calls==0]="normal"
  res$calls[res$calls==1]="gain"
  res$calls[res$calls==2]="amplification"
  
  ##add the data.frame to the result
  #calls of the segment
  call=unlist(lapply(1:ncol(segmentData$segmented),FUN=function(x){res$calls[indexStart[[x]],x]}))

  res$segment= data.frame(chrom=paste0("chr",chrom),chromStart=chromStart,chromEnd=chromEnd,probes=unlist(indexEnd)-unlist(indexStart)+1,means=means,calls=call)
  rownames(res$segment)=NULL
  
  return(res)
}



#' @title Create the list of parameters for \link{segmentation} function
#' 
#' @description create the list of parameters for \link{segmentation} function
#' 
#' @param copynumber A vector containing the copy-number signal for one patient and one chromosome.
#' @param chromosome Chromosome associated with the copy-number signal.
#' @param position Position of the signal.
#' @param featureNames Names of the probes (not necessary).
#' @param sampleNames Name of the sample (not necessary).
#'
#' @return a list in the right format for \link{segmentation} function
#'  
#' @author Quentin Grimonprez
#'
#' @export
segmentationObject=function(copynumber,chromosome,position,featureNames,sampleNames)
{
  #copynumber : vector
  if(missing(copynumber))
    stop("copynumber is missing.")
  if(!is.numeric(copynumber) || !is.vector(copynumber))
    stop("copynumber must be a vector of real.")
  #position : vector
  if(missing(position))
    position=1:length(copynumber)
  if(!is.numeric(position) || !is.vector(position))
    stop("position must be a vector of real.")
  if(length(position)!=length(copynumber))
    stop("copynumber and position must have the same length.")
  
  #chromosome : integer
  if(missing(chromosome))
    stop("chromosome is missing")
  if(!is.double(chromosome))
    stop("chromosome must be an integer.")
  if(!is.wholenumber(chromosome))
    stop("chromosome must be an integer.")

  #featureNames : vector of char
  if(missing(featureNames))
  {
    if(!is.null(names(copynumber)))
      featureNames=names(copynumber)
    else
      featureNames=as.character(1:nrow(copynumber))
  }
  if(!is.vector(featureNames))
    stop("featureNames must be a vector of character.")
  if(!is.character(featureNames) || (length(featureNames)!=length(position)) )
    stop("featureNames must be a vector of character.")
  
  #sampleNames : char
  if(missing(sampleNames))
      sampleNames="Sample.1"
  
  if(!is.character(sampleNames))
    stop("sampleNames must be a character.")
  if(length(sampleNames)!=1)
    stop("sampleNames must be a character.")
  
  return(list(copynumber=copynumber,chromosome=chromosome,position=position,featureNames=featureNames,sampleNames=sampleNames))
}


#' @title Create the list of parameters for \link{callingProcess} function
#' 
#' @description create the list of parameters for \link{callingProcess} function
#' 
#' @param copynumber A matrix containing the copy-number signal. Each column is a different patient.
#' @param segmented A matrix containing the segmented copy-number signal. Matrix of the same size as copynumber.
#' @param chromosome Chromosome associated with the copy-number signal.
#' @param position Position of the signal.
#' @param featureNames Names of the probes (not necessary).
#' @param sampleNames Name of the sample (not necessary).
#'
#' @return a list in the right format for \link{callingProcess} function
#' 
#' @author Quentin Grimonprez  
#' 
#'@export
callingObject=function(copynumber,segmented,chromosome,position,featureNames,sampleNames)
{
  #copynumber : matrix
  if(is.vector(copynumber))
    copynumber=as.matrix(copynumber)
  if(!is.numeric(copynumber) || !is.matrix(copynumber))
    stop("copynumber must be a matrix of real.")
  #segmented : matrix
  if(is.vector(segmented))
    segmented=as.matrix(segmented)
  if(!is.numeric(segmented) || !is.matrix(segmented))
    stop("segmented must be a matrix of real.")
  if( (nrow(segmented)!=nrow(copynumber)) || (ncol(copynumber)!=ncol(copynumber)) )
    stop("segmented and copynumber have not the same dimension.")
    
  #position : vector
  if(!is.numeric(position) || !is.vector(position))
    stop("position must be a vector of real.")
  if(length(position)!=nrow(copynumber))
    stop("The length of position and the number of rows of copynumber don't match.")
  
  #chromosome : vector of integer
  if(!is.vector(chromosome))
    stop("chromosome must be a vector of integer.")
  if(length(chromosome)!=length(position))
    stop("chromosome and position must have the same length.")
  if(sum(is.wholenumber(chromosome))!=length(chromosome))
    stop("chromosome must be a vector of integer.")
  
  #featureNames : vector of char
  if(missing(featureNames))
  {
    if(!is.null(rownames(copynumber)))
      featureNames=rownames(copynumber)
    else
      featureNames=as.character(1:nrow(copynumber))
  }
  if(!is.vector(featureNames))
    stop("featureNames must be a vector of character.")
  if(!is.character(featureNames) || (length(featureNames)!=length(position)) )
    stop("featureNames must be a vector of character.")
  
  #sampleNames : char
  if(missing(sampleNames))
  {
    if(!is.null(colnames(copynumber)))
      sampleNames=colnames(copynumber)
    else
      sampleNames=paste0("Sample.",1:ncol(copynumber))
  } 
  if(!is.character(sampleNames))
    stop("sampleNames must be a character.")
  if(length(sampleNames)!=ncol(copynumber))
    stop("sampleNames must be a character.")
  
  return(list(copynumber=copynumber,segmented=segmented,chromosome=chromosome,startPos=position,featureNames=featureNames,sampleNames=sampleNames))
}

#' convert CNA object (output of the function segment from DNAcopy package) into a list for the argument segmentData of the function \link{callingProcess}.
#'
#' @title Convert CNAobject
#' 
#' @param CNAobject Output object of segment function from DNAcopy package
#' 
#' @return a list at the required format of \link{callingProcess}.
#' 
#' @seealso callingProcess
#' 
#' @author Quentin Grimonprez
#' 
#' @export
CNAobjectToCGHcallObject=function(CNAobject)
{
  #number of samples and names
  sampleNames=names(CNAobject$data)[-c(1:2)]
  nbSample=length(sampleNames)
  
  #creation of the segmented matrix
  segmented=matrix(nrow=length(CNAobject$data$chrom),ncol=nbSample)
  
  #fill the segmented matrix
  apply(CNAobject$output,1,FUN=function(segment)#loop on all the segment
    {
      #index of the sample : index of the column of segmented to fill
      ind=which(sampleNames==segment[1])
      
      #index of the chromosome to look for the position on the right chromosome
      indChrom=which(CNAobject$data$chrom==as.numeric(segment[2]))

      #start index of the segment
      indSegmentStart=which(CNAobject$data$maploc[indChrom]==as.numeric(segment[3]))[1]
      #end index of the segment
      indSegmentEnd=which(CNAobject$data$maploc[indChrom]==as.numeric(segment[4]))
      indSegmentEnd=indSegmentEnd[length(indSegmentEnd)]#if position are several times

      #fill the segmented matrix
      segmented[indChrom[indSegmentStart:indSegmentEnd],ind]<<-rep(as.numeric(segment[6]),indSegmentEnd-indSegmentStart+1)
      invisible(NULL)
    })
  
  
  list(copynumber=as.matrix(CNAobject$data[-c(1:2)]),segmented=segmented,chromosome=CNAobject$data$chrom,startPos=CNAobject$data$maploc,sampleNames=sampleNames,featureNames=paste0(1:length(CNAobject$data$maploc)))
}