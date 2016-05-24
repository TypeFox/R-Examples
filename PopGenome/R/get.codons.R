# ------------------------------------------------------------
# get Codon Informations 
# ------------------------------------------------------------
setGeneric("get.codons", function(object, regionID) standardGeneric("get.codons"))
 setMethod("get.codons", "GENOME",
 function(object, regionID){

#if(!SNP.DATA){
CodonInfo <- codontable()
Triplets  <- CodonInfo$Triplets
Protein   <- CodonInfo$Protein

# Generate polarity Codon Vector
POLARITY <- character(64)
POLARITY[c(1:4,11,17:24,33:44,49:52)] <- "nonpolar"
POLARITY[c(5:10,13:14,53:58,61:62)]   <- "polar"
POLARITY[c(12,15:16)]                 <- "stop codon"
POLARITY[c(25:32,59:60,63:64)]        <- "basic"
POLARITY[45:48]                       <- "acidic"
# --------------------------------

# Get the codons from the GENOME slot
minor  <- sapply(object@region.data@codons[[regionID]],function(x){return(x[1])}) # FIXME really minor ? done 
mayor  <- sapply(object@region.data@codons[[regionID]],function(x){return(x[2])}) # FIXME really mayor ? done

#print(minor)
#print(mayor)

# ReCheck minor/mayor Codons
if(!object@snp.data){

minor.alleles <- object@region.data@minor.alleles[[regionID]]
subst         <- object@region.data@biallelic.substitutions[[regionID]]

 for (xx in 1:length(minor)){

	if(minor.alleles[xx]==subst[2,xx]){
        save.minor <- minor[xx]
        minor[xx]  <- mayor[xx]
	mayor[xx]  <- save.minor
	}

 }
}

# Create Polarity
ALTPOL <- POLARITY[minor]
REFPOL <- POLARITY[mayor]


# Get Triplets
ALT    <- Triplets[minor, ,drop=FALSE]
REF    <- Triplets[mayor, ,drop=FALSE]

# Character Codons
MATCHNUCS <- c("T","C","G","A","N")
NucALT    <- matrix(MATCHNUCS[ALT], ncol=3) 
NucREF    <- matrix(MATCHNUCS[REF], ncol=3) 
NucALT    <- apply(NucALT,1,function(x){paste(x,collapse="")})
NucREF    <- apply(NucREF,1,function(x){paste(x,collapse="")})

# Create Protein Vector 
ALT     <- Protein[1,minor]
REF     <- Protein[1,mayor]

# Get coding SNP positions
codeids  <- which(!is.na(object@region.data@synonymous[[regionID]]))

# create syn/nonsyn vector
erg      <- as.logical(object@region.data@synonymous[[regionID]][codeids])


# Biallelic Positions
bial.pos <- object@region.data@biallelic.sites[[regionID]][codeids]

# Collect information and store in data.frame 
XX    <- c("Position","Codons (minor)","Codons (major)","Protein (minor)","Protein (major)","synonymous", "Polarity (major)","Polarity (minor)")
DATA  <- data.frame(bial.pos,NucALT, NucREF, ALT, REF, as.logical(erg), ALTPOL, REFPOL)
colnames(DATA) <- XX

return(DATA)

#}











#############################################################################################











if(FALSE){#FIXME
CodonInfo <- codontable()
Triplets  <- CodonInfo$Triplets
Protein   <- CodonInfo$Protein

# Generate polarity Codon Vector
POLARITY <- character(64)
POLARITY[c(1:4,11,17:24,33:44,49:52)] <- "nonpolar"
POLARITY[c(5:10,13:14,53:58,61:62)]   <- "polar"
POLARITY[c(12,15:16)]                 <- "stop codon"
POLARITY[c(25:32,59:60,63:64)]        <- "basic"
POLARITY[45:48]                       <- "acidic"
# --------------------------------


# Generate Codon Tables 
bial.pos     <- object@region.data@biallelic.sites[[regionID]]
minor        <- object@region.data@biallelic.substitutions[[regionID]][1,]
major	     <- object@region.data@biallelic.substitutions[[regionID]][2,]

if(reverse.strand){
cod.pos      <- (reading.start.pos-bial.pos)   %%3
}else{
cod.pos      <- (bial.pos - reading.start.pos) %%3
}

# Generate Codon-Positions # Rows are SNP-positions
codons <- matrix(,length(cod.pos),3)
for (yy in 1:length(cod.pos)){
    if(reverse.strand){

     if(cod.pos[yy]==0){codons[yy,]=c(bial.pos[yy],bial.pos[yy]-1,bial.pos[yy]-2);next}
     if(cod.pos[yy]==1){codons[yy,]=c(bial.pos[yy]+1,bial.pos[yy],bial.pos[yy]-1);next}
     if(cod.pos[yy]==2){codons[yy,]=c(bial.pos[yy]+2,bial.pos[yy]+1,bial.pos[yy]-2);next}	

    }else{

     if(cod.pos[yy]==0){codons[yy,]=c(bial.pos[yy],bial.pos[yy]+1,bial.pos[yy]+2);next}
     if(cod.pos[yy]==1){codons[yy,]=c(bial.pos[yy]-1,bial.pos[yy],bial.pos[yy]+1);next}
     if(cod.pos[yy]==2){codons[yy,]=c(bial.pos[yy]-2,bial.pos[yy]-1,bial.pos[yy]);next}
    
    }
}


## Reading the reference chromosome
file.info <- .Call("get_dim_fasta",ref.chr)
CHR       <- .Call("get_ind_fasta",ref.chr,1,file.info[[1]][2])

# Create codons with nucleotides
Nuc.codons    <- CHR[codons]
Nuc.codons    <- matrix(Nuc.codons,ncol=3)

ALT           <- Nuc.codons
REF           <- Nuc.codons

komplement <- c(4,3,2,1,5)

for(yy in 1: dim(Nuc.codons)[1]){

 if(reverse.strand){

  #convert to komplement nucleotides
  REF[yy,]  <- komplement[REF[yy,]]
  ALT[yy,]  <- komplement[ALT[yy,]]
  minor[yy] <- komplement[minor[yy]]
  major[yy] <- komplement[major[yy]]
  ###########

 if(cod.pos[yy]==0){REF[yy,1] <- minor[yy];ALT[yy,1] <- major[yy];next}
 if(cod.pos[yy]==1){REF[yy,2] <- minor[yy];ALT[yy,2] <- major[yy];next}
 if(cod.pos[yy]==2){REF[yy,3] <- minor[yy];ALT[yy,3] <- major[yy];next}

  
 }else{ 

 if(cod.pos[yy]==0){REF[yy,1] <- minor[yy];ALT[yy,1] <- major[yy];next}
 if(cod.pos[yy]==1){REF[yy,2] <- minor[yy];ALT[yy,2] <- major[yy];next}
 if(cod.pos[yy]==2){REF[yy,3] <- minor[yy];ALT[yy,3] <- major[yy];next}

 }

} 

# Character Codons
MATCHNUCS <- c("T","C","G","A","N")
NucALT    <- matrix(MATCHNUCS[ALT], ncol=3) 
NucREF    <- matrix(MATCHNUCS[REF], ncol=3) 

NucALT <- apply(NucALT,1,function(x){paste(x,collapse="")})
NucREF <- apply(NucREF,1,function(x){paste(x,collapse="")})

# Coding Codons ...
ALT <- codonise64(ALT)
REF <- codonise64(REF)

# Polarity
ALTPOL <- POLARITY[ALT]
REFPOL <- POLARITY[REF]

# Proteins
ALT <- Protein[1,ALT]
REF <- Protein[1,REF]

# ALT and REF Proteins
CHECK <- cbind(ALT,REF)

# Check syn/nonsyn
erg         <- apply(CHECK,1,function(x){return(length(unique(x)))})
erg[erg==2] <- FALSE #nonsyn
erg[erg==1] <- TRUE  #syn


# Collect information and store in data.frame 
XX    <- c("Position","Codons (major)","Codons (minor)","Protein (major)","Protein (minor)","synonymous", "Polarity (major)","Polarity (minor)")
DATA  <- data.frame(bial.pos,NucALT, NucREF, ALT, REF, as.logical(erg), ALTPOL, REFPOL)
colnames(DATA) <- XX

return(DATA)
}# END of if SNP.DATA

})
