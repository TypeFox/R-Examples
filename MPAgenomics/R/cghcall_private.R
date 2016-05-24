# All the functions of this file are from CGHcall and CGHbase packages created by Sjoerd Vosse 
# and Mark van de Wiel under License GPL

##private functions from package CGHcall
# .MakeData
# .countcl
# .sumreg
# .sumsqreg
# .varregtimescount
# .varproffun
# .alphafun
# .alpha0all
# .totallik
# .posteriorp
# .xgivenknotrunc
# .minusEloglikreg
# .reallik4
# .reallikk4
# .assignNames

#private functions from CGHbase
# .getCentromere
# .convertChromosomeToArm
# 

# @author Sjoerd Vosse & Mark van de Wiel


####CGHcall
.MakeData <- function(smratall,chrnum) 
{  #updated 16/07/10

    nc <- ncol(smratall)
    allregions  <- vector("list", nc)
    ls      <- length(chrnum)
    chrnumfw <- c(100,chrnum[-ls])
    for (j in 1:(nc)) {
        smrat   <- smratall[,j] #changed 19/06/2009, much faster
        smratshfw <- c(100,smrat[-ls])
        smratshbw <- c(smrat[-1],100)
        mult <- (smrat-smratshfw)*(smratshbw-smratshfw) + (chrnumfw-chrnum)
        wm <- which(mult!=0)
        smwh <- smrat[wm]
        wmend <- c(wm[-1]-1,ls)
        regions <- cbind(wm,wmend,smwh)
        allregions[[j]] <- regions

    rm(smrat,smratshfw,smratshbw,mult,smwh,wm);gc()
    }
    return(allregions)
}


.countcl <- function(k, regionsdat) {
    regionsdat[k,2]-regionsdat[k,1]+1
}

.sumreg <- function(k, dat, regionsdat) {
    sum(dat[regionsdat[k,1]:regionsdat[k,2]])
}

.sumsqreg <- function(k, dat, regionsdat) {
    sum((dat[regionsdat[k,1]:regionsdat[k,2]])^2)
}

.varregtimescount <- function(k, counts, dat, regionsdat) {
    var(dat[regionsdat[k,1]:regionsdat[k,2]])*counts[k]
}

.varproffun <- function(prof, vcnmat, profile) {
    vcnprof <- vcnmat[profile == prof & !is.na(vcnmat[,1]),]
    if(!is.null(dim(vcnprof))) {
        return(sum(vcnprof[,1])/sum(vcnprof[,3]))
    } else {
        return(vcnprof[1]/vcnprof[3])
    }
}

.alphafun <- function(k, profile, priorp, pm, varprofall, allsum, allsumsq, allnc=allnc, robustsig, allcell, prior) {
    prof        <- profile[k]
    regionsk    <- which(profile==prof)
    nregk       <- length(regionsk)
    if (prior=="all") {
        nregk       <- length(profile)
        regionsk    <- 1:nregk
    }
    totpost     <- rep(1,nregk)%*%t(sapply(regionsk, .posteriorp, priorp=priorp, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig,allcell=allcell))
    return(totpost/nregk)
}

.alpha0all <- function(nreg, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell, prior) {
    prevprof    <- 0
    alphaall    <- c()
    alpha1      <- .alphafun(1, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell, prior)
    for (i in (1:nreg)) {
        if (prior != "all") {
            curprof <- profile[i]
            if (curprof==prevprof) {
                newalpha <- oldalpha
            }
            if (curprof != prevprof) {
                newalpha <- .alphafun(i, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell, prior)
            }
            oldalpha <- newalpha
            prevprof <- curprof
        }
        if (prior == "all") {
            newalpha <- alpha1
        }
        alphaall <- rbind(alphaall,newalpha)
    }
    return(alphaall)
}

.totallik <- function(pm, nreg, posteriorprev, alphaprev, varprofall, allsum, allsumsq, allnc, robustsig,allcell,ncpus=1) {
    if(ncpus==1) {
    return(sum(sapply(1:nreg, .minusEloglikreg, pm=pm, posteriorprev=posteriorprev, alphaprev=alphaprev, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc,robustsig=robustsig,allcell=allcell)))
    } else {
    snowfall::sfExport("pm")
    return(sum(snowfall::sfSapply(1:nreg, .minusEloglikreg, pm=pm, posteriorprev=posteriorprev, alphaprev=alphaprev, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc,robustsig=robustsig,allcell=allcell)))
    } 
}

.totallik2 <- function(pm, nreg, posteriorprev, alphaprev, varprofall, allsum, allsumsq, allnc, robustsig,allcell,ncpus=1) {
  if(ncpus==1) {
    return(sum(sapply(1:nreg, .minusEloglikreg, pm=pm, posteriorprev=posteriorprev, alphaprev=alphaprev, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc,robustsig=robustsig,allcell=allcell)))
  } else {
    snowfall::sfExport("pm")
    return(sum(snowfall::sfSapply(1:nreg, .minusEloglikreg, pm=pm, posteriorprev=posteriorprev, alphaprev=alphaprev, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc,robustsig=robustsig,allcell=allcell)))
  } 
}


.posteriorp <- function(k, priorp, pm, varprofall, allsum, allsumsq, allnc,robustsig,allcell) {
   all3        <- sapply(c(1,2,3,4,5,6), .xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig,allcell=allcell)
    maxim       <- max(all3)
    all         <- exp(all3-maxim)
    allprior    <- all*priorp[k,]
    tot         <- all%*%priorp[k,]
    #TODO vérifier pertinence correction bug avec segment extreme d'un seul point
    if (tot==0)
      return (allprior)
     
    return(allprior/tot)
}


.xgivenknotrunc <- function(k, class, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell,betas=c(1,0)) {

    sumx    <- allsum[k]
    sumxsq  <- allsumsq[k]
    ncl     <- allnc[k]
    v1      <- 2*varprofall[k]
    sd1     <- sqrt(v1)
    cell <- allcell[k]

    if (class==6) {
        mu <- (log2(2)/log2(1.5))*(0.10 + exp(-pm[4]))+0.3+exp(-pm[5])
    } else if (class==5) {
        mu <- (log2(2)/log2(1.5))*(0.10 + exp(-pm[4]))  
    } else if (class==4) {
        mu <- 0.10 + exp(-pm[class])
    } else if (class==3) {
        mu <- -0.05+0.1*exp(-(pm[class])^2) 
    } else if (class==1) {
        mu <- -0.10-exp(-pm[2])-0.3 - exp(-pm[1])
    } else {
        mu <- -0.10-exp(-pm[class])
    }
    
    munorm <- -0.05+0.1*exp(-(pm[class])^2) 
    #NEW
    mu <- ((1-cell)*munorm + cell*mu)*(betas[1]+betas[2]*sd1)
    
    v2 <- if(class == 4 | class == 2) 2*((pm[5+class])^2 + 0.0001) 
    else {if(class==3) {if(robustsig=="yes") 1/4*((pm[7])^2 + 0.0001)+1/4*((pm[9])^2 + 0.0001)+2*((pm[8]^2)+0.0001)
    else {2*((pm[8])^2 + 0.0001)}}
    else {if(class==1) 2*(pm[7]^2+(pm[6])^2 + 0.0001) 
    else {if(class==5) 2*((pm[10])^2+(pm[9])^2 + 0.0001)
    else 2*((pm[11])^2 + (pm[10])^2 + (pm[9])^2 + 0.0001)}}} # for stability
    
    v2norm <- if(robustsig=="yes") {1/4*((pm[7])^2 + 0.0001)+1/4*((pm[9])^2 + 0.0001)+2*((pm[8]^2)+0.0001)}
    else {2*((pm[8])^2 + 0.0001)}
    
    #NEW
    v2 <- if(class==3) v2norm else {(1-cell)^2*v2norm + cell^2*v2}
    #v2 <- (1-cell)*v2norm + cell^2*v2
    
    A       <- (sumx*v2+mu*v1)/(ncl*v2+v1)
    B       <- v1*v2/(ncl*v2+v1)
    C       <- sumxsq/v1 + mu^2/v2 - A^2/B
    logD    <- log(B)/2 -log(v2)/2
    logres  <- - C + logD
    res     <- exp(logres)
    return(logres)
}

.minusEloglikreg <- function(k, posteriorprev, alphaprev, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell) {
    all3 <- sapply(c(1,2,3,4,5,6), .xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc, robustsig=robustsig, allcell=allcell)
    # patch for windows bug "unable to evaluate fonction with inital parameter"
    nonOind=which(alphaprev[k,]!=0)
    return(-all3[nonOind]%*%posteriorprev[,k][nonOind] - (log(alphaprev[k,][nonOind])%*%posteriorprev[,k][nonOind]))
}

.reallik4 <- function(nreg, alpha, pm, varprofall, allsum, allsumsq, allnc,robustsig,allcell,betas=c(1,0)) {
    return(-sum(sapply(1:nreg, .reallikk4, pm=pm, alpha=alpha, varprofall =varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig,allcell=allcell,betas=betas)))
}
 
.reallikk4 <- function(k, alpha, pm, varprofall, allsum, allsumsq, allnc,robustsig,allcell,betas=c(1,0)) {
#betas=c(1,0);allcell=allcell_pr*0.2;k=1;pm=bstart; varprofall =varprof_allall_pr;allsum= allsumall_pr;allsumsq= allsumsqall_pr;allnc= allncall_pr;robustsig= robustsig;
    all3    <-sapply(c(1,2,3,4,5,6), .xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, 
    allnc=allnc,robustsig=robustsig,allcell=allcell,betas=betas)
    maxim   <- max(all3)
    all     <- exp(all3-maxim)
    tot     <- maxim + log(all%*%alpha[k,])
    return(tot)
}

## pour fonction expandCGHcall
# .assignNames
.assignNames <- function(matrix, object) {
    colnames(matrix) <- object$sampleNames
    rownames(matrix) <- object$featureNames
    matrix
}

# .callFromSeg
#.callFromSeg <- function(seg, assayData) {
#    result <- new("cghCall", assayData=assayData,
#                            phenoData=phenoData(seg), 
#                            featureData=featureData(seg), 
#                            annotation=annotation(seg), 
#                            experimentData=experimentData(seg)
#                            )
#    result
#}



######### CGHbase  pour fonction CGHcall
#private functiosn from CGHbase
# .getCentromere
# .convertChromosomeToArm


.getCentromere <- function(build) {
    build <- as.integer(gsub('[^0-9]', '', build))
    centromere <- matrix(NA, 23, 2);
    if (build == 34 || build == 16) {
        centromere[1,] <- c(120093537, 122333537)
        centromere[2,] <- c(91810927, 94810927)
        centromere[3,] <- c(90425755, 93325755)
        centromere[4,] <- c(49575659, 52575659)
        centromere[5,] <- c(46451142, 49451142)
        centromere[6,] <- c(58877002, 61877002)
        centromere[7,] <- c(57832528, 60832528)
        centromere[8,] <- c(43856263, 46856263)
        centromere[9,] <- c(44443361, 47443361)
        centromere[10,] <- c(39208941, 41588941)
        centromere[11,] <- c(51602318, 54602318)
        centromere[12,] <- c(34747961, 36142961)
        centromere[13,] <- c(14900000, 16768000)
        centromere[14,] <- c(15070000, 18070000)
        centromere[15,] <- c(15260000, 18260000)
        centromere[16,] <- c(36366889, 38166889)
        centromere[17,] <- c(22408570, 25408570)
        centromere[18,] <- c(15398887, 16762885)
        centromere[19,] <- c(26923622, 29923622)
        centromere[20,] <- c(26314569, 29314569)
        centromere[21,] <- c(10260000, 13260000)
        centromere[22,] <- c(11330000, 14330000)
        centromere[23,] <- c(57548803, 60548803)
        # centromere[24,] <- c(9757849, 12757849)
    } else if (build == 35 || build == 17) {
        centromere[1,] <- c(121147476, 123387476)
        centromere[2,] <- c(91748045, 94748045)
        centromere[3,] <- c(90587544, 93487544)
        centromere[4,] <- c(49501045, 52501045)
        centromere[5,] <- c(46441398, 49441398)
        centromere[6,] <- c(58938125, 61938125)
        centromere[7,] <- c(57864988, 60864988)
        centromere[8,] <- c(43958052, 46958052)
        centromere[9,] <- c(46035928, 49035928)
        centromere[10,] <- c(39244941, 41624941)
        centromere[11,] <- c(51450781, 54450781)
        centromere[12,] <- c(34747961, 36142961)
        centromere[13,] <- c(16000000, 17868000)
        centromere[14,] <- c(15070000, 18070000)
        centromere[15,] <- c(15260000, 18260000)
        centromere[16,] <- c(35143302, 36943302)
        centromere[17,] <- c(22187133, 22287133)
        centromere[18,] <- c(15400898, 16764896)
        centromere[19,] <- c(26923622, 29923622)
        centromere[20,] <- c(26267569, 28033230)
        centromere[21,] <- c(10260000, 13260000)
        centromere[22,] <- c(11330000, 14330000)
        centromere[23,] <- c(58465033, 61465033)
        # centromere[24,] <- c(11237315, 12237315)
    } else if (build == 36 || build == 18) {
        centromere[1,] <- c(121236957, 123476957)
        centromere[2,] <- c(91689898, 94689898)
        centromere[3,] <- c(90587544, 93487544)
        centromere[4,] <- c(49354874, 52354874)
        centromere[5,] <- c(46441398, 49441398)
        centromere[6,] <- c(58938125, 61938125)
        centromere[7,] <- c(58058273, 61058273)
        centromere[8,] <- c(43958052, 46958052)
        centromere[9,] <- c(47107499, 50107499)
        centromere[10,] <- c(39244941, 41624941)
        centromere[11,] <- c(51450781, 54450781)
        centromere[12,] <- c(34747961, 36142961)
        centromere[13,] <- c(16000000, 17868000)
        centromere[14,] <- c(15070000, 18070000)
        centromere[15,] <- c(15260000, 18260000)
        centromere[16,] <- c(35143302, 36943302)
        centromere[17,] <- c(22187133, 22287133)
        centromere[18,] <- c(15400898, 16764896)
        centromere[19,] <- c(26923622, 29923622)
        centromere[20,] <- c(26267569, 28033230)
        centromere[21,] <- c(10260000, 13260000)
        centromere[22,] <- c(11330000, 14330000)
        centromere[23,] <- c(58598737, 61598737)
    } else { # 37 / 19
        centromere[1,] <- c(121535434, 124535434)
        centromere[2,] <- c(92326171, 95326171)
        centromere[3,] <- c(90504854, 93504854)
        centromere[4,] <- c(49660117, 52660117)
        centromere[5,] <- c(46405641, 49405641)
        centromere[6,] <- c(58830166, 61830166)
        centromere[7,] <- c(58054331, 61054331)
        centromere[8,] <- c(43838887, 46838887)
        centromere[9,] <- c(47367679, 50367679)
        centromere[10,] <- c(39254935, 42254935)
        centromere[11,] <- c(51644205, 54644205)
        centromere[12,] <- c(34856694, 37856694)
        centromere[13,] <- c(16000000, 19000000)
        centromere[14,] <- c(16000000, 19000000)
        centromere[15,] <- c(17000000, 20000000)
        centromere[16,] <- c(35335801, 38335801)
        centromere[17,] <- c(22263006, 25263006)
        centromere[18,] <- c(15460898, 18460898)
        centromere[19,] <- c(24681782, 27681782)
        centromere[20,] <- c(26369569, 29369569)
        centromere[21,] <- c(11288129, 14288129)
        centromere[22,] <- c(13000000, 16000000)
        centromere[23,] <- c(58632012, 61632012)
        # centromere[24,] <- c(10104553, 13104553)
    }
    centromere <- apply(centromere, 1, mean);
    return(centromere);
}

.convertChromosomeToArm <- function(dataframe, build) { #changed 22/06/2009; 
    cat("Dividing chromosomes into arms using centromere positions from", build, "\n\n");
    centromere  <- .getCentromere(build);
    chr <- dataframe[,2]
    bp <- dataframe[,3]
    chrlev <- unique(chr)
    a<-1
    chrarms <- c()
    for(i in chrlev){
    print(i)
    chri <- which(chr==i)
    bpi <- bp[chri]
    wbpi <- length(which(bpi<=centromere[i]))
    wbpil <- length(which(bpi>centromere[i]))
    if(wbpi>0) {chrarms <- c(chrarms,rep(a,wbpi));a<-a+1}
    if(wbpil>0) {chrarms <- c(chrarms,rep(a,wbpil));a<-a+1}  
    }  
    dataframe[,2] <- chrarms
    return(dataframe)
}

