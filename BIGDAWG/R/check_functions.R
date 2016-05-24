#' HLA formatting check
#'
#' Checks data to see if HLA data is properly formatted.
#' @param x All columns of HLA genotyping data.
#' @note This function is for internal BIGDAWG use only.
CheckHLA <- function(x) {
  #Return TRUE if properly formatted HLA
  
  x[is.na(x)] <- "00:00" # NA cells
  x[x=="^"] <- "00:00" # absent cells
  test <- apply(x,MARGIN=c(1,2),FUN=function(z) length(unlist(strsplit(as.character(z),split=":"))))
  test <- apply(test,MARGIN=2,FUN=min)
  Flag <- as.logical(min(test)==2)
  
  return(Flag)
}

#' Loci presence check
#'
#' Checks available loci against data specific to ensure complete overlap.
#' @param x Loci available in exon protein list alignment object.
#' @param y Unique column names
#' @note This function is for internal BIGDAWG use only.
CheckLoci <- function(x,y) {
  #Returns TRUE if absent locus(loci) encountered
  #x=Loci available in ExonPtnList
  #y=Loci.Set from data
  
  Output <- list()
  
  y <- unique(unlist(y))
  
  Flag <- ( !sum(y %in% x) == length(y) )
  Output[['Flag']] <- Flag
  if(Flag) {
    Output[['Loci']] <- paste(y[!y %in% x],collapse=",")
  } else {
    Output[['Loci']] <- NA
  }
  
  return(Output)
}

#' Allele presence check
#'
#' Checks available alleles against data to ensure complete overlap.
#' @param x Exon protein list alignment object.
#' @param y Genotypes from data file
#' @param z1 loci in data file
#' @param z2 Genotype column names
#' @note This function is for internal BIGDAWG use only.
CheckAlleles <- function(x,y,z1,z2) {
  #Returns TRUE if unknown allele(s) encountered
  #Checks at 3 levels of resolution: Full, 2-Field, and 1-Field
  #x=ExonPtnList
  #y=genos
  #z1=loci
  #z2=loci.ColNames
  
  Output <- list()
  for(i in 1:length(z1)) {
    
    Locus.tmp <- z1[i]
    
    # Available Alleles
    x.locus <- x[[Locus.tmp]]
    
    # Data Alleles
    y.locus <- y[which(z2==Locus.tmp)]
    y.locus <- matrix(c(unlist(y.locus[,1]),y.locus[,2]),ncol=1)
    y.locus <- na.omit(y.locus)
    
    # Remove Absent Allele Calls
    Allele.rm <- which(y.locus=="^")
    if(length(Allele.rm)>0) { y.locus <- y.locus[-Allele.rm] }
    
    #identify minimum resolution
    tmp <- sapply(y.locus,FUN=strsplit,split=":")
    tmp <- lapply(tmp,unlist)
    Res <- min(unlist(lapply(tmp,length)))
    
    y.alleles <- sort(unique(unlist(y.locus)))
    x.alleles <- unique(x.locus[,'Allele'])
    
    #format according to minimum resolution
    if (Res>=2) {
      y.alleles <- unique(sapply(y.alleles,GetField,Res=2))
      x.alleles <- unique(x.locus[,'Trimmed'])
    } else if (Res==1) {
      y.alleles <- unique(sapply(y.alleles,GetField,Res=1))
      x.alleles <- unique(sapply(unique(x.locus[,'Trimmed']),GetField,Res=1))
    }
    
    Output[[Locus.tmp]] <- list(Flag=ifelse(!sum(y.alleles %in% x.alleles)==length(y.alleles),T,F),
                                Allele=paste(Locus.tmp,paste(y.alleles[!y.alleles %in% x.alleles],collapse=","),sep="*"))
    
    i = i + 1
    Output
    
  }
  return(Output)
}

#' Data summary function
#'
#' Summary function for sample population within data file.
#' @param Tab Loci available in exon protein list alignment object.
#' @param All.ColNames Column names from genotype data.
#' @param rescall HLA resolution set for analysis.
#' @param HLA HLA bigdawg argument passed to function
#' @note This function is for internal BIGDAWG use only.
PreCheck <- function(Tab,All.ColNames,rescall,HLA) {
  
  Grp0 <- which(Tab[,2]==0)
  Grp1 <- which(Tab[,2]==1)
  nGrp0 <- length(Tab[Grp0,2])
  nGrp1 <- length(Tab[Grp1,2])
  
  if(min(nGrp0,nGrp1)==0) {
    Err.Log("Case.Con")
    stop("Analysis Stopped.",call. = F)
  }
  
  Loci <- as.list(unique(All.ColNames[3:length(All.ColNames)]))
  nLoci <- length(Loci)
  GTYPE <- Tab[,3:ncol(Tab)]
  colnames(GTYPE) <- All.ColNames[3:length(All.ColNames)]
  nGTYPE <- unlist(lapply(Loci,function(x) length(unique(unlist(GTYPE[,which(colnames(GTYPE)==x)])))))
  Grp0un <- unlist(lapply(Loci,function(x) length(unique(unlist(GTYPE[Grp0,which(colnames(GTYPE)==x)])))))
  Grp1un <- unlist(lapply(Loci,function(x) length(unique(unlist(GTYPE[Grp1,which(colnames(GTYPE)==x)])))))
  
  nMissing <- unlist(lapply(Loci,function(x) sum(is.na(GTYPE[,which(colnames(GTYPE)==x)]))))
  Grp0miss <- unlist(lapply(Loci,function(x) sum(is.na(GTYPE[Grp0,which(colnames(GTYPE)==x)]))))
  Grp1miss <- unlist(lapply(Loci,function(x) sum(is.na(GTYPE[Grp1,which(colnames(GTYPE)==x)]))))
  
  cat("  Sample Summary\n")
  cat("    Sample Size (n):",nrow(Tab),"\n")
  cat("    ...Number of Controls/Cases:",paste(paste(nGrp0,nGrp1,sep="/"),collapse=", "),"\n")
  cat("    Allele Count (2n):",nrow(Tab)*2,"\n")
  cat("    Total loci in file:",nLoci,"\n")
  cat("    Unique loci:",paste(Loci,collapse=", "),"\n") 
  cat("    Unique alleles per locus:",paste(nGTYPE,collapse=", "),"\n")
  cat("    ...Unique in Controls/Cases:",paste(paste(Grp0un,Grp1un,sep="/"),collapse=", "),"\n")
  cat("    Missing alleles per locus:",paste(nMissing,collapse=", "),"\n")
  cat("    ...Missing in Controls/Cases:",paste(paste(Grp0miss,Grp1miss,sep="/"),collapse=", "),"\n")
  cat("\n")
  
  if(HLA) {
    
    Grp0res <- max(unlist(lapply(Loci,function(x) max(unlist(lapply(strsplit(unlist(GTYPE[Grp0,which(colnames(GTYPE)==x)]),split=":"),length))))))
    Grp1res <- max(unlist(lapply(Loci,function(x) max(unlist(lapply(strsplit(unlist(GTYPE[Grp1,which(colnames(GTYPE)==x)]),split=":"),length))))))
    
    if(max(Grp0res,Grp1res)>4) {
      Err.Log("High.Res")
      stop("Analysis Stopped.",call. = F)
    }
    
    cat("  Observed Allele Resolution\n")
    cat("    Max Resolution Controls:",paste(Grp0res,"-Field",sep=""),"\n")
    cat("    Max Resolution Cases:",paste(Grp1res,"-Field",sep=""),"\n")
    cat("    Defined Resolution:",rescall,"\n")
    if(Grp0res!=Grp1res){ cat("  ***** Warning. \n") }
    if(Grp0res!=Grp1res){ cat("  ***** There may exist a Case-Control field resolution imbalance.\n") }
    if(Grp0res!=Grp1res){ cat("  ***** Considering trimming to",paste(min(Grp0res,Grp1res),"-Field resolution.",sep=""),"\n") }
    cat("\n")
  }
  
  if(HLA) {
    Output <- list(SampleSize=nrow(Tab),
                   No.Controls=nGrp0,
                   No.Cases=nGrp1,
                   AlleleCount=nrow(Tab)*2,
                   TotalLoci=nLoci,
                   Loci=paste(Loci,collapse=", "),
                   AllelePerLocus=paste(nGTYPE,collapse=", "),
                   MissingPerLocus=paste(nMissing,collapse=", "),
                   MaxResGrp0=paste(Grp0res,"-Field",sep=""),
                   MaxResGrp1=paste(Grp1res,"-Field",sep=""),
                   SuggestedRes=paste(min(Grp0res,Grp1res),"-Field",sep=""),
                   SetRes=rescall)
  } else {
    Output <- list(SampleSize=nrow(Tab),
                   No.Controls=nGrp0,
                   No.Cases=nGrp1,
                   AlleleCount=nrow(Tab)*2,
                   TotalLoci=nLoci,
                   Loci=paste(Loci,collapse=", "),
                   AllelePerLocus=paste(nGTYPE,collapse=", "),
                   MissingPerLocus=paste(nMissing,collapse=", "),
                   SetRes=rescall)
  }
  
  return(do.call(rbind,Output))
  
}
