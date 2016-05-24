#' Alignment Filter
#'
#' Filter Protein Exon Alignment File for Specific Alleles.
#' @param Align Protein Alignment Object.
#' @param Alleles to be pulled.
#' @param Locus Locus to be filtered against.
#' @note This function is for internal BIGDAWG use only.
AlignmentFilter <- function(Align, Alleles, Locus) { 
  
  getCols <- c(match(c("Trimmed","P group","Unknowns","NullPositions"),colnames(Align)),
               grep("Position\\.",colnames(Align)))              
  Align.sub <- unique(Align[,getCols])
  Align.sub <- Align.sub[Align.sub[,'Trimmed'] %in% Alleles,]
  
  if(!is.null(nrow(Align.sub))) {
  
    Alleles.S <- names(which(table(Align.sub[,'Trimmed'])==1))
    Alleles.M <- names(which(table(Align.sub[,'Trimmed'])>1))
    
    if(length(Alleles.M>0)) {
      # Removing Duplicates at 2-Field Level
      Align.tmp <- list()
      
      for(m in 1:length(Alleles.M)) {
        Allele <- Alleles.M[m]
        Alignsub.Grp <- Align.sub[which(Align.sub[,"Trimmed"]==Allele),]
        
        #for multiple matches, check first for P Group match
        PgroupAlleles <- grep("P",Alignsub.Grp[,'P group'],fixed=T)
        
        if(length(PgroupAlleles)==0){
          #if multiple matches exist, use the match with the least unknowns in sequence
          Unknowns.Grp <- which(Alignsub.Grp[,'Unknowns']==min(Alignsub.Grp[,'Unknowns']))
          if(length(Unknowns.Grp)>1) {Unknowns.Grp <- Unknowns.Grp[1]}
          Align.tmp[[Allele]] <- Alignsub.Grp[Unknowns.Grp,]
          
        } else if(length(PgroupAlleles)==1){
          Align.tmp[[Allele]] <- Alignsub.Grp[PgroupAlleles,]
          
        } else if(length(PgroupAlleles)>1) {
          #if multiple matches exist, use the match with the least unknowns in sequence
          Alignsub.Grp <- Alignsub.Grp[PgroupAlleles,] 
          Unknowns.Grp <- which(Alignsub.Grp[,'Unknowns']==min(Alignsub.Grp[,'Unknowns']))
          if(length(Unknowns.Grp)>1) {Unknowns.Grp <- Unknowns.Grp[1]}
          Align.tmp[[Allele]] <- Alignsub.Grp[Unknowns.Grp,]
        }
        
      }; rm(m)
      
      Align.tmp <- do.call(rbind,Align.tmp)
      AlignMatrix <- rbind(Align.tmp,Align.sub[Align.sub[,'Trimmed'] %in% Alleles.S,])
      
    } else {
      
      AlignMatrix <- Align.sub[Align.sub[,'Trimmed'] %in% Alleles.S,]
      
    }
    
    AlignMatrix <- cbind(rep(Locus,nrow(AlignMatrix)),AlignMatrix)
    rownames(AlignMatrix) <- NULL
    colnames(AlignMatrix)[1] <- "Locus"
    colnames(AlignMatrix)[which(colnames(AlignMatrix)=="Trimmed")] <- "Allele.2D"  
    AlignMatrix <- AlignMatrix[order(AlignMatrix[,'Allele.2D']),]
  
  } else {
    
    AlignMatrix <- Align.sub
    
  }
    
  return(AlignMatrix)
}

#' Amino Acid Contingency Table Build
#'
#' Build Contingency Tables for Amino Acid Analysis.
#' @param x Filtered alignmnet list element.
#' @param y Phenotype groupings.
#' @note This function is for internal BIGDAWG use only.
AAtable.builder <- function(x,y) {
  #x = list element for filtered alignment
  #y = HLA_grp
  
  # Create count Grp 0 v Grp 1 (Control v Case)
  x[,2] <- gsub(" ","Null",x[,2])
  x[,2] <- gsub("\\*","Unknown",x[,2])
  x[,2] <- gsub("\\.","InDel",x[,2])
  
  Residues <- unique(x[,2])
  AminoAcid.df <- mat.or.vec(nr=length(Residues),2)
  colnames(AminoAcid.df) <- c("Group.0", "Group.1")
  rownames(AminoAcid.df) <- Residues
  
  y[,2:3] <- apply(y[,2:3],MARGIN=c(1,2),GetField,Res=2)
  
  # Grp 0 (Control)
  Grp0 <- y[which(y[,'grp']==0),]
  Grp0cnt <- table(c(x[match(Grp0[,2],x[,1]),2],
                     x[match(Grp0[,3],x[,1]),2]))
  PutRange <- match(rownames(AminoAcid.df),names(Grp0cnt))
  AminoAcid.df[,'Group.0'] <- Grp0cnt[PutRange]
  
  # Grp 1 (Case)
  Grp1 <- y[which(y[,'grp']==1),]
  Grp1cnt <- table(c(x[match(Grp1[,2],x[,1]),2],
                     x[match(Grp1[,3],x[,1]),2]))
  PutRange <- match(rownames(AminoAcid.df),names(Grp1cnt))
  AminoAcid.df[,'Group.1'] <- Grp1cnt[PutRange]
  
  AminoAcid.df[is.na(AminoAcid.df)] <- 0
  
  return(AminoAcid.df)
  
}

#' Contingency Table Check
#'
#' Checks amino acid contingency table data frame to ensure required variation exists.
#' @param x contingency table.
#' @note This function is for internal BIGDAWG use only.
AA.df.check <- function(x) {
  # Returns true if insufficient variations exists
  
  # Replace NAs
  x[is.na(x)] <- 0
  
  if(is.null(nrow(x))){
    return(T)
  } else if(nrow(x)<2){
    return(T)
  } else if (nrow(x)==2) {
    ExpCnts <- chisq.test(x)$expected
    if( min(ExpCnts)!=0 & length(which(ExpCnts<5))==0 ){
      return(F)
    } else {
      return(T)
    }
  } else {
    ExpCnts <- chisq.test(x)$expected
    #unbinned
    OK.rows <- as.numeric(which(apply(ExpCnts,min,MARGIN=1)>=5))
    if(length(OK.rows)==0) {
      return(T)
    } else if(length(OK.rows)>=2) {
      unbinned <- x[OK.rows,]
    } else {
      unbinned <- do.call(cbind,as.list(x[OK.rows,]))
      rownames(unbinned) <- rownames(x)[OK.rows]
    }
    #binned
    Rare.rows <- as.numeric(which(apply(ExpCnts,min,MARGIN=1)<5))
    if(length(Rare.rows)>=2) {
      binned <- colSums(x[Rare.rows,])
      New.df <- rbind(unbinned,binned) 
    } else {
      New.df <- x
    }
    ExpCnts.df <- chisq.test(New.df)$expected
    if (length(which(ExpCnts.df<5))/length(which(ExpCnts.df>=1))<0.2){
      return(F)
    } else {
      return(T)
    }
  } 
}
