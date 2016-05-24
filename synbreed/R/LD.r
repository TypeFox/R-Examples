pairwiseLD <- function(gpData,chr=NULL,type=c("data.frame","matrix"),use.plink=FALSE,ld.threshold=0,ld.window=99999,rm.unmapped=TRUE){

    # catch errors
    if(is.null(gpData$geno)) stop("no genotypic data available")
    if(!gpData$info$codeGeno) stop("use function 'codeGeno' before")
    if(is.null(gpData$map))  stop("no map information available")
    
    type <- match.arg(type)
    if(ld.threshold>0 & type=="matrix") warning("'ld.threshold not used for type='matrix'")
    if(any(is.na(gpData$geno))) stop("no missing values allowed, try to impute using 'codeGeno'")

    # extract information from gpData  (only use mapped markers)
    if(rm.unmapped){
      mapped <- !(is.na(gpData$map$chr) | is.na(gpData$map$pos)) 
    }
    else mapped  <- rep(TRUE,ncol(gpData$geno))
    
    # restrict data to mapped markers
    gpData$geno <- gpData$geno[,mapped]
    gpData$map <- gpData$map[mapped,]
    
    linkageGroup <- as.character(gpData$map$chr)
    if(!rm.unmapped){
      gpData$map$chr <- as.character(gpData$map$chr)
      gpData$map$chr[is.na(linkageGroup)] <- "NA"
      gpData$map$chr <- as.factor(gpData$map$chr)
      linkageGroup[is.na(linkageGroup)] <- "NA"
    }
    pos <- gpData$map$pos
    names(pos) <- rownames(gpData$map)
   
    
    # select chromosomes if 'chr' is specified
    lg <- unique(linkageGroup)
    if(!is.null(chr)){ 
        lg <- chr
        if(any(chr=="all")) stop("option chr='all' not yet possible") #linkageGroup <- rep("all",length(linkageGroup))
        # NOTE: positions must be ordered consequtively within chromosomes
    }                                                
  
    # initialize return LD value data.frame list
    retList <- list()
    # initialize return LD value matrix list
    retMat <- list()

  
    # loop over all chromosomes (linkage groups)
    for (i in 1:length(lg)){
       
       if(use.plink){ # i.e. if there are 3 genotypes
        # call PLINK to compute the LD as r2
        sel <- rownames(gpData$map[gpData$map$chr!= lg[i],])
        gpTEMP <- discard.markers(gpData,which=sel)
        pre <- paste("chr",lg[i],sep="")
        write.plink(gpTEMP,type=type,ld.threshold=ld.threshold,ld.window=ld.window,prefix=pre)
        system(paste("plink --script ",pre,"plinkScript.txt",sep="")) 
       
        # distances between markers
        if(type=="matrix"){
           distance <- as.matrix(dist(pos[linkageGroup == lg[i]],diag=FALSE,upper=FALSE))
        }
        
        # read data from PLINK
        if(type=="matrix") {
          ld.r2 <- as.matrix(read.table(paste(pre,".ld",sep="")))
          colnames(distance) <- rownames(distance) <- colnames(ld.r2 ) <- rownames(ld.r2 ) <- names(pos)[linkageGroup == lg[i]]
        }
        if(type=="data.frame"){
          ld.r2.df.plink <- read.table(paste(pre,".ld",sep=""),header=TRUE,stringsAsFactors=FALSE)
          #distance <- abs(pos[ld.r2.df.plink$SNP_A]-pos[ld.r2.df.plink$SNP_B])
          ld.r2.df <- with(ld.r2.df.plink,data.frame(marker1=SNP_A,marker2=SNP_B,r2=R2,dist=abs(BP_A-BP_B),stringsAsFactors=FALSE))
        }  
       
        
       
       } # end if(use.plink)
       
       else{  # i.e. if there are 2 genotypes (e.g. DH lines)
           # read information from data
           markeri <- gpData$geno[,linkageGroup==lg[i]]
           p <- ncol(markeri)
           mn <- colnames(markeri)
           posi <- pos[linkageGroup==lg[i]]
       
           ld.r2 <- cor(markeri,method="spearman",use="pairwise.complete.obs")^2
 
           if(type=="data.frame"){
              ld.r2i <- ld.r2[lower.tri(ld.r2)]
              # index vectors for LD data.frame
              rowi <- rep(1:p,times=(p:1)-1)
              coli <- p+1 - sequence(1:(p-1))
              coli <- coli[length(coli):1]
              # distance between markers
              disti <- abs(posi[rowi] - posi[coli])
              ld.r2.df <- data.frame(marker1=mn[rowi],
                                     marker2=mn[coli],
                                     r2=ld.r2i,
                                     dist=disti,
                                     stringsAsFactors=FALSE)
              if(i==length(lg)) if(!rm.unmapped) {
                ld.r2ini <- cor(gpData$geno[,linkageGroup!=lg[i]],markeri,method="spearman",use="pairwise.complete.obs")^2
                ld.r2.dfini <- data.frame(marker1=rep(colnames(ld.r2ini), each=nrow(ld.r2ini)), 
                                          marker2=rep(rownames(ld.r2ini), ncol(ld.r2ini)), 
                                          r2=as.numeric(ld.r2ini), 
                                          dist=rep(NA,ncol(ld.r2ini)*nrow(ld.r2ini)), 
                                          stringsAsFactors=FALSE)
                ld.r2.df <- rbind(ld.r2.df, ld.r2.dfini)
              }
          }
          if(type=="matrix"){
              # matrix of distances 
              distance <- as.matrix(dist(pos[linkageGroup == lg[i]],diag=FALSE,upper=FALSE))
              colnames(distance) <- rownames(distance) <- colnames(gpData$geno)[linkageGroup == lg[i]]
          }
       }

       
       
       # create dataset with information from above in a data.frame
       if(type=="data.frame") retList[[i]] <- ld.r2.df 
       # and as a matrix
       if(type=="matrix") retMat$LD[[i]] <- ld.r2           # omit lower/upper triangle?
       if(type=="matrix") retMat$distance[[i]] <- distance  
    }  
    
    if(type=="data.frame") names(retList) <- paste("chr", lg, sep="_")
    if(type=="matrix") names(retMat$LD) <- names(retMat$distance) <- paste("chr", lg, sep="_")


      # return values 
      #if(type=="both") return(list(dataFrame=retList,matrix=retMat))
      if(type=="data.frame") {
        class(retList) <- "LDdf"
        return(retList)
      }  
      if(type=="matrix"){
         class(retMat) <- "LDmat"
         return(retMat)
      }   
}
