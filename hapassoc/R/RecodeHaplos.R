# HapAssoc- Inference of trait-haplotype associations in the presence of uncertain phase
# Copyright (C) 2003  K.Burkett, B.McNeney, J.Graham

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

########################################################################

RecodeHaplos<-function(dat,numSNPs,allelic,maxMissingGenos=1, verbose=TRUE) {

  # Split dat into nonSNP and SNP data - also report marker names if verbose=TRUE
  # so the user can see if pre.hapassoc is separating genetic and non-genetic 
  # as expected.
  ncols<-ncol(dat)
  if(!allelic) {    
      if(verbose) {
        cat(paste("Haplotypes will be based on the following SNPs (genotypic format): \n",
                   paste(names(dat)[(ncols-numSNPs+1):ncols],collapse=", "),"\n"))
        cat(paste("Remaining variables are: \n",paste(names(dat)[1:(ncols-numSNPs)],collapse=", "),"\n"))
      }
      # turn genetic d.f. to allelic d.f.
      nhdmnames <- names(dat)[1:(ncols-numSNPs)] # names are lost if ncols-numSNPs is a single column
      dat1<-dat[,1:(ncols-numSNPs)]
      for(i in (ncols-numSNPs+1):ncols) {
          dat1<-data.frame(dat1, substr(as.character(dat[,i]),1,1),
                                 substr(as.character(dat[,i]),2,2))

          names(dat1)[(length(dat1)-1):length(dat1)] <-
             paste(names(dat)[i],".",1:2, sep="")
      }
      dat <- dat1
      ncols <- ncol(dat)
      names(dat)[1:(ncols-2*numSNPs)] <- nhdmnames
  }  else {
      if(verbose) {
       mnms<-names(dat)[(ncols-2*numSNPs+1):ncols]
       cat("Haplotypes will be based on the following SNPs (allelic format): \n")
       for(ii in 1:(length(mnms)/2)) 
          cat(paste(" SNP ",ii,": ",mnms[2*ii-1],"/",mnms[2*ii],"\n",sep=""))
       cat(paste("Remaining variables are: \n",paste(names(dat)[1:(ncols-2*numSNPs)],collapse=", "),"\n"))
      }
  }
  # Substitute empty strings "" with NA:
  for(i in (ncols-2*numSNPs+1):ncols) {
      dat[!nzchar(as.vector(dat[,i])),i]<-NA
  }
  nonsnpcols<-ncols-2*numSNPs
  snpcols<-ncols-nonsnpcols
  nonSNPdat<-dat[,1:(ncols-2*numSNPs)]
  nhdmnames<-names(dat)[1:(ncols-2*numSNPs)] #save the col names for the output - they will disappear when we typecast to matrix below

  nonSNPdat <- as.matrix(nonSNPdat)#as.data.frame(nonSNPdat)
  SNPdat<-as.matrix(dat[,(ncols-2*numSNPs+1):ncols])

  for (i in seq(length=ncol(SNPdat)/2, from=1, by=2))
       if (length(na.omit(unique(c(SNPdat[,i],SNPdat[,i+1])))) > 2)
           stop(paste("Alleles", paste(na.omit(unique(c(SNPdat[,i],SNPdat[,i+1]))),collapse=", "),
                "at marker",i,"-- no more than 2 alleles allowed. \n"))

  # if SNP data are given as nucleotide letters, create SNPkey,
  # a vector of the more frequent nucleotide in each SNP
  # and encode SNPdat into 0's and 1's

  if (mode(SNPdat[,1])=="character") {
      SNPdatAlpha <- SNPdat
      SNPkey <- vector("character", 0)  # most frequent nucleotide in SNP
      SNPkey1 <- vector("character", 0) # least frequent nucleotide in SNP  

      for (i in seq(length=ncol(SNPdat)/2, from=1, by=2)) {
          nucleotide1 <- unique(na.omit(c(SNPdat[,i],SNPdat[,i+1])))[1]
          nucleotide2 <- unique(na.omit(c(SNPdat[,i],SNPdat[,i+1])))[2]   
          if (sum(na.omit(SNPdat[,i:(i+1)])==nucleotide1) >= sum(na.omit(SNPdat[,i:(i+1)])==nucleotide2)) {
              SNPkey <- c(SNPkey, nucleotide1); SNPkey1<- c(SNPkey1, nucleotide2)
          }
          else {
              SNPkey <- c(SNPkey, nucleotide2); SNPkey1<- c(SNPkey1, nucleotide1)
          }

          SNPdat[SNPdat[,i]==tail(SNPkey, n=1),i]<-0
          SNPdat[SNPdat[,i+1]==tail(SNPkey, n=1),i+1]<-0
          SNPdat[!SNPdat[,i]==0,i]<-1
          SNPdat[!SNPdat[,i+1]==0,i+1]<-1          
      }      
  }

  # Pre-process input to deal with missing data
  preProcDat<-handleMissings(SNPdat,nonSNPdat,numSNPs,maxMissingGenos)
  # Note: we could get rid of these 2 as.matrix(..) typecasts if we fix up handleMissings sometime in the future - MP.nov.2003
  SNPdat<-as.matrix(preProcDat$SNPdat)
  nonSNPdat<-as.matrix(preProcDat$nonSNPdat)
  ID<-preProcDat$ID

  # Initialization and setup:
  row.names(nonSNPdat)<-as.character(1:nrow(nonSNPdat))
  haploLabs<-makeHaploLabN(0:(2^numSNPs-1),numSNPs=numSNPs)
  numHaplos<-nrow(haploLabs)

  # Build data frames to hold nonSNP design matrix, haplotype data design
  # matrix and initial weights. The haplo design matrix will have a column
  # corresponding to each possible haplotype.
  wt<-NULL
  nonHaploDM<-matrix(nrow=(nrow(SNPdat)*2^(numSNPs-1)),ncol=nonsnpcols) #part of design matrix corresponding to non-haplo data
  nhdmidx<-1 #row index initializer for nonHaploDM matrix
  haploDM<-matrix(nrow=(nrow(SNPdat)*2^(numSNPs-1)),ncol=numHaplos) #part of design matrix corresponding to haplotypes
  hdmidx<-1 #row index initializer for haploDM matrix
  haploMat<-matrix(nrow=(nrow(SNPdat)*2^(numSNPs-1)),ncol=snpcols) #the actual haplotypes for each (pseudo-) individual
  hmatidx<-1 #row index initalizer for haploMat matrix
  ID.vec<-matrix(nrow=(nrow(SNPdat)*2^(numSNPs-1)),ncol=1) #to keep track of pseudo-individuals original ID
  ididx<-1 #row index for ID.vec

  heteroVec<-rep(NA,numSNPs) #Initializer

  # Main loop to construct design matrix-- do SNP and non-SNP data separately
  # Loop over subjects and if a subject's geno data does not determine
  # haplos, enumerate all consistent haplos and add pseudo-individuals
  # to the design matrices for each possible haplo specification  
  for(i in 1:nrow(SNPdat)){

    # Function isHetero returns a vector heterovec of T's and F's
    # describing whether the person is heterozygous at each locus
    # This has now been inlined below (function isHetero to be removed??)
    for(j in 1:numSNPs){heteroVec[j]<-(.subset(SNPdat,i,2*j-1)!=.subset(SNPdat,i,2*j))}
    numHetero<-sum(heteroVec)

    # The rows of matrix myhaplos are possible haplotype combinations
    # for the ith subject. The first column of the matrix has a character
    # string of 0's and 1's for a binary number describing the first
    # haplotype and the second column has a string for a binary number
    # describing the second haplotype.

    myhaplos<-getHaplos(SNPdat[i,],heteroVec)
    numHaploComb<-nrow(myhaplos)

    for(j in 1:numHaploComb) { #loop over haplo combos consistent w/ obs data
      haploDM[hdmidx,]<-codeHaploDM(myhaplos[j,],haploLabs)
      hdmidx<-hdmidx+1
    }
    for(j in 1:numHaploComb){
      nonHaploDM[nhdmidx,]<-nonSNPdat[i,]
      nhdmidx<-nhdmidx+1
      }
    for(j in 1:numHaploComb){
      ID.vec[ididx]<-ID[i]
      ididx<-ididx+1
    }
    for(j in 1:nrow(myhaplos)) {
    	haploMat[hmatidx,]<-myhaplos[j,]
	hmatidx<-hmatidx+1
    }

  }  #end for loop over subjects

  haploDM<-haploDM[1:(hdmidx-1),]
  haploMat<-haploMat[1:(hmatidx-1),]
  ID.vec<-ID.vec[1:(ididx-1),]
  nonHaploDM<-nonHaploDM[1:(nhdmidx-1),]

  #Now renormalize the weights to sum to 1 before returning. Only necessary if 
  #there was missing SNP data. With missing SNPs handleMissings creates 
  #copies of a person, one for each possible complete set of single-locus
  #genotypes. E.g. if there is one allele of one SNP missing for a person,
  #that person will be copied into two "people" each with weight 1.
  for(i in 1:(ididx-1)) {
    wt[i]<-1/sum(ID.vec==ID.vec[i])
  }

  #Only need to return the columns of haploDM that have non-zero column sums:
  myColSums<-colSums(haploDM)
  haploDM<-haploDM[,myColSums>0, drop=FALSE]
  haploDM<-data.frame(haploDM)

  #Here we re-format our  haploMat to be compatible with the old output format and store it
  #in HaploMat2.  It *may* be that in the future, haploMat2 could be dropped if the
  #computationally-friendly format of haploMat is preferred over the older format in whatever
  #code is calling RecodeHaplos and making use of the output. -Matt
  n<-ncol(haploMat)
  haploMat2<-matrix(nrow=nrow(haploMat),ncol=2)
  haploMat2[,1]<-paste("h",haploMat[,1],sep="")
  haploMat2[,2]<-paste("h",haploMat[,n/2+1],sep="")
  for(i in 2:(n/2)) {
	haploMat2[,1]<-paste(haploMat2[,1],haploMat[,i],sep="")
	haploMat2[,2]<-paste(haploMat2[,2],haploMat[,i+n/2],sep="")
  }
  #Need to protect columns of haploMat2 from being coerced into factors
  #with the I() function. BM Jan/04
  haploMat2<-data.frame(haplo1=I(haploMat2[,1]),haplo2=I(haploMat2[,2]))

  #Put column names just like the old format. -Matt
  hdmnames<-makeHaploLab(0:(2^numSNPs-1),numSNPs)

  if (exists("SNPkey")) {   # data file SNPs encoded with nucleotide letters
     # fix names(haploDM):
      for (i in 1:length(hdmnames)) {
          for (j in 1:nchar(hdmnames[i])) {
              if(substr(hdmnames[i],j,j)==0)
                  substr(hdmnames[i],j,j) <- SNPkey[j]
              else
                  substr(hdmnames[i],j,j) <- SNPkey1[j]    
          }
      }

      # fix haploMat:
      for (h in 1:length(haploMat2)) {
          for (i in 1:length(haploMat2[,h])) {
              for (j in 2:nchar(haploMat2[i,h])) {
                  if(substr(haploMat2[i,h],j,j)==0)
                      substr(haploMat2[i,h],j,j) <- SNPkey[j-1]
                  else
                      substr(haploMat2[i,h],j,j) <- SNPkey1[j-1]    
              }
          }
      }
  }                

  #Only the labels with >0 column sums though
  names(haploDM)<-paste("h",hdmnames[myColSums>0],sep="")

  nonHaploDM<-data.frame(nonHaploDM)
  #Put column names just like the old format. -Matt
  names(nonHaploDM)<-nhdmnames
  return(list(nonHaploDM=nonHaploDM,haploDM=haploDM,haploMat=haploMat2,
              wt=wt, ID=ID.vec))
}


## Other functions called in RecodeHaplos:

########################################################################
handleMissings<-function(SNPdat,nonSNPdat,numSNPs,maxMissingGenos)
{
  temnonSNPdat<-na.omit(nonSNPdat)
  omittedRows<-attr(temnonSNPdat,"na.action") # Which rows were removed
  nonSNPdat<-data.frame(temnonSNPdat)
  if(!is.null(omittedRows)) {
    warning(paste(length(omittedRows),
                  "subjects removed because of missing nongenetic data\n"))
    SNPdat<-SNPdat[-omittedRows,] # Remove these from SNPdat too
  }

  numMissingGenos<-rep(0,nrow(SNPdat)) # Count SNP genos with missing data
  for(i in 1:numSNPs) {
    numMissingGenos <- numMissingGenos+
                         (is.na(SNPdat[,(2*i-1)])|is.na(SNPdat[,2*i]))
  }

  # Remove people with too many missing genotypes
  ind<-numMissingGenos<=maxMissingGenos

  if(sum(!ind)>0) {
    warning(paste(sum(!ind),
            "subjects with missing data in more than",maxMissingGenos,
            "genotype(s) removed\n"))

    nonSNPdat<-nonSNPdat[ind,,drop=FALSE]
    SNPdat<-data.frame(SNPdat[ind,])
    numMissingGenos<-numMissingGenos[ind]
  }

  ID <- c(1:nrow(SNPdat)) # initial ID's
  missingGenos<-(numMissingGenos>0)
  if(any(missingGenos)) {
    # First save copies of those with missing data
    temSNPdat<-SNPdat[missingGenos,,drop=FALSE]

    temnonSNPdat<-nonSNPdat[missingGenos,,drop=FALSE]
    temID<-ID[missingGenos]

    # Reduce SNPdat and nonSNPdat to people who have no missing data
    SNPdat<-SNPdat[!missingGenos,,drop=FALSE]
    nonSNPdat<-nonSNPdat[!missingGenos,,drop=FALSE]
    ID<-ID[!missingGenos]

# data.frame casting turns numeric characters to factors.
# That created a problem when a rare SNP has only one factor level
# while completePhenos has 2 factor levels.  -SB
for(ii in 1:ncol(SNPdat)) {
SNPdat[[ii]]<-as.numeric(as.character(SNPdat[[ii]]))
}
for(ii in 1:ncol(temSNPdat)) {
temSNPdat[[ii]]<-as.numeric(as.character(temSNPdat[[ii]]))
}

    # Now augment the SNP and nonSNP data by enumerating sets of complete
    # SNP genos (call these "complete phenos") consistent with the
    # observed phenotypes.
    for(i in 1:sum(missingGenos)) {
      missingVec<-is.na(temSNPdat[i,])
      completePhenos<-getPhenos(temSNPdat[i,],numSNPs,missingVec)
      numPhenos<-nrow(completePhenos)
      for(j in 1:numPhenos) { #loop over complete phenos consistent w/ obs data
        SNPdat<-rbind(SNPdat,completePhenos[j,])
        nonSNPdat<-data.frame(rbind(as.matrix(nonSNPdat),
                   as.matrix(temnonSNPdat[i,])),
                   row.names=as.character(1:nrow(SNPdat)))
        ID <- c(ID,temID[i])
      }
    } # end loop over subjects with missing data
  } #end if(any(missingGenos))

  return(list(SNPdat=SNPdat,nonSNPdat=nonSNPdat,ID=ID))
}

getPhenos<-function(snps,numSNPs,missingVec) {
  # Inefficient but simple approach: consider both a 0 or 1 for each NA,
  # e.g for one locus with NA/NA --> 0/0, 0/1, 1/0, 1/1.
  # Then order the alleles       --> 0/0, 0/1, 0/1, 1/1
  # Then remove duplicates       --> 0/0, 0/1, 1/1

  # Enumerate all possible values for missings
  k<-sum(missingVec)
  # Can use makeHaploLab to enumerate all possible alleles for missing vals
  misAlleles<-makeHaploLab(0:(2^(k)-1),numSNPs=k)


  misAlleles<-matrix(as.numeric(unlist(strsplit(misAlleles,split=""))),
                         ncol=k,byrow=TRUE) #turn labels into a numeric matrix
  numPhenos<-nrow(misAlleles)

  myPhenos<-matrix(NA,nrow=numPhenos,ncol=length(snps))

  if(any(missingVec==FALSE)) {
      knownAlleles <- snps[!missingVec]
      myPhenos[, !missingVec] <- matrix(rep(knownAlleles, numPhenos),
        ncol = length(knownAlleles), byrow = TRUE)
  }

  knownAlleles<-snps[!missingVec]
  myPhenos[,!missingVec]<-matrix(rep(knownAlleles,numPhenos),
                                 ncol=length(knownAlleles),
                                 byrow=TRUE)
  myPhenos[,missingVec]<-misAlleles

  # Now order alleles
  for(i in 1:numSNPs) {
    ind<-myPhenos[,(2*i-1)]>myPhenos[,2*i] #these are the 1/0 genos
    myPhenos[ind,(2*i-1)]<-0; myPhenos[ind,2*i]<-1
  }

  # now reduce to unique rows with built-in unique.array function
  myPhenos<-unique.array(myPhenos)

  return(myPhenos)
}

makeHaploLab<-function(x,numSNPs=2) {

  #Function used to construct labels for SNP haplos; e.g. with 3 SNPs we want
  #labels 000, 001, ..., 111 which in each case is n plus one of
  #the numbers 0,1,...,(2^3-1) represented in base2.
  #Takes x as a base 10 number and returns haplo label
  #For example if numSNPs is 3 then x=0 would lead to the string "000",
  #x=1 would lead to "001" and so on. Function can take x as a vector so
  #an example of usage is: mylabs<-makeHaploLab(0:(2^3-1),numSNPs=3)
  #
  #Basic idea is to fill in digits of the base2 numbers left to right.
  #Example: numSNPs=3, x=7. Then x= 1*2^2 + 1*2^1 + 1*2^0 = 111 in base2.
  #Start by filling in the first digit, then second , then third.

  ans<-"" #start label as blank

  for(i in (numSNPs-1):0) {
    digit<-floor(x/2^i)
    ans<-paste(ans,as.character(digit),sep="") #update answer
    x<-x-digit*2^i
  }
  return(ans)
}

# A version of makeHaploLab that returns a numeric vector instead of a string of haplotype pairs
makeHaploLabN<-function(x,numSNPs=2) {
	len<-length(x)
	ans<-matrix(0,nrow=len,ncol=numSNPs)

	for(i in (numSNPs-1):0){
		digit<-floor(x/2^i)
		ans[,numSNPs-i]<-digit
		x<-x-digit*2^i
	}
	return(ans)
}


########################################################################

isHetero<-function(SNPvec,numSNPs) {

  #Function to take a vector of SNP data for a person and figure out
  #which loci person is heterozygous for. Returns logical vector.
  #Assumes genotypes are in pairs in the vector, e.g.
  #SNPvec = (M1.allele1,M1.allele2,M2.allele1,M2.allele2,...)

  if(length(SNPvec)/2 != numSNPs)
    stop("SNPvec not compatible with numSNPs\n")

  ans<-rep(NA,numSNPs)

  for(i in 1:numSNPs) {
    #marker i's data are in elements 2*i-1 and 2*i of the SNPvec
    ans[i]<-(SNPvec[2*i-1]!=SNPvec[2*i])
  }
  return(ans)
}
########################################################################


getHaplos<-function(SNPvec,heteroVec){
	nloci<-length(heteroVec)
	k<-sum(heteroVec)

	if(k<=1) {
		haplo<-matrix(nrow=1,ncol=length(SNPvec))
		mid<-length(SNPvec)/2
		haplo[1,(1:mid)]<-SNPvec[2*(1:nloci)-1]
		haplo[1,(mid+1):length(SNPvec)]<-SNPvec[2*(1:nloci)]
	}
	else {
		heteroStates1<-makeHaploLabN(0:(2^(k-1)-1),numSNPs=k)
		heteroStates2<-1-heteroStates1
		haplo<-matrix(NA,ncol=2*nloci,nrow=2^(k-1))
		hvec1<-rep(FALSE,2*nloci)
		hvec2<-rep(FALSE,2*nloci)
		hvec1[1:nloci]<-heteroVec
		hvec2[(nloci+1):(2*nloci)]<-heteroVec
		haplo[,hvec1]<-heteroStates1
		haplo[,hvec2]<-heteroStates2
		if((nloci-k)>0) {
			homoStates<-SNPvec[2*(1:nloci)][!heteroVec]
			homoStates<-matrix(rep(homoStates,2^(k-1)),ncol=(nloci-k),byrow=TRUE)
			hvec1<-rep(FALSE,2*nloci) ##was len
			hvec2<-rep(FALSE,2*nloci) ##was len
			hvec1[1:nloci]<-!heteroVec
			hvec2[(nloci+1):(2*nloci)]<-!heteroVec
			haplo[,hvec1]<-homoStates
			haplo[,hvec2]<-homoStates
		}
	}
	return(haplo)
}

##########################################################################


codeHaploDM<-function(haplos,haploLabs){

  n=length(haplos)
  nsnp=ncol(haploLabs)

  ans1<-t(haplos[1:(n/2)]==t(haploLabs))
  ans2<-t(haplos[(n/2+1):n]==t(haploLabs))
  ans11<-ans1[,1]
  ans22<-ans2[,1]
  for(i in 2:nsnp) {
	ans11<-ans11&ans1[,i]
	ans22<-ans22&ans2[,i]
  }
  ans=ans11+ans22

  return(ans)
}

