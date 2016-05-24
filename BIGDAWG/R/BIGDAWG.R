#' BIGDAWG wrapper function
#'
#' This is the main wrapper function for each analysis.
#' @param Data Name of the genotype data file.
#' @param HLA Logical indicating whether data is HLA class I/II genotyping data only.
#' @param Run.Tests Specifics which tests to run.
#' @param Loci.Set Input list defining which loci to use for analyses (combinations permitted). 
#' @param All.Pairwise Logical indicating whether all pairwise.
#' @param Trim Logical indicating if HLA alleles should be trimmed to a set resolution.
#' @param Res Numeric setting what desired resolution to trim HLA alleles.
#' @param EVS.rm Logical indicating if expression variant suffixes should be removed.
#' @param Missing Numeric setting allowable missing data for running analysis (may use "ignore").
#' @param Cores.Lim Interger setting the number of cores accessible to BIGDAWG (Windows, Cores = 1L).
#' @param Results.Dir Optional, full directory name for BIGDAWG output.
#' @param Output Logical indicating if data should be written to output directory (mainly used for testing).
#' @examples
#' ### The following examples use the synthetic data set bundled with BIGDAWG
#' 
#' # Haplotype analysis with no missing genotypes for two loci
#' # Significant haplotypic association with phenotype 
#' BIGDAWG(Data="HLA_data", Run.Tests="H", Missing=0, Loci.Set=list(c("DRB1","DQB1")))
#' 
#' # Hardy-Weinberg and Locus analysis ignoring missing data
#' # Significant associations with phenotype at all but DQB1
#' BIGDAWG(Data="HLA_data", Run.Tests=c("HWE","L"), Missing="ignore")
#' 
#' # Hardy-Weinberg analysis trimming data to 2-Field resolution
#' # Significant deviation at DQB1
#' BIGDAWG(Data="HLA_data", Run.Tests="HWE", Trim=TRUE, Res=2)
BIGDAWG <- function(Data, HLA=TRUE, Run.Tests, Loci.Set, All.Pairwise=FALSE, Trim=FALSE, Res=2, EVS.rm=FALSE, Missing=2, Cores.Lim=1L, Results.Dir, Output=TRUE) {
  
  options(warn=-1)
  
  MainDir <- getwd()
  
  #Define Output Directory
  if (Output) {
    if(missing(Results.Dir)) {
      OutDir <- paste(MainDir,"/output ",format(Sys.time(), "%d%m%y %H%M%S"),sep="")
      dir.create(OutDir)
    } else {
      OutDir <- Results.Dir
    }
  }
  
  on.exit(setwd(MainDir), add = TRUE)
  
  if(Sys.info()['sysname']=="Windows") {Cores=1L} else {Cores=parallel::detectCores()}
  
  #Check for updated ExonPtnList ... UpdatePtnList. Use if found.
  UpdatePtnList <- NULL; rm(UpdatePtnList)
  UPL <- paste(path.package('BIGDAWG'),"/data/UpdatePtnAlign.RData",sep="")
  if( file.exists(UPL) ) { load(UPL) ; EPL <- UpdatePtnList ; rm(UpdatePtnList); UPL.flag=T } else { EPL <- ExonPtnList ; UPL.flag=F}
  
  cat("\n>>>>>>>>>>>>>>>>>>>>>>>>> BIGDAWG Analysis <<<<<<<<<<<<<<<<<<<<<<<<<\n\n")
  
############################################################################################################################################  
### Read in Data ###########################################################################################################################
  
  NAstrings=c("NA","","****","-","na","Na")
  
  if (Data=="HLA_data") {
    Tab <- BIGDAWG::HLA_data
    colnames(Tab) <- toupper(colnames(Tab))
    All.ColNames <- gsub(".1","",colnames(Tab),fixed=T)
    rownames(Tab) <- NULL
    setwd(OutDir)
    DRBFLAG <- NULL
    
  } else {
    
    # Read in data and Pre-process
    if(!file.exists(Data)) {
      Err.Log("Bad.Filename",Data)
      stop("Analysis stopped.",call.=F) }
    Tab <- read.table(Data, header = T, sep="\t", stringsAsFactors = F, na.strings=NAstrings, fill=T, comment.char = "#", strip.white=T)
    if(HLA==T) { if(sum(grepl("DRB3.4.5",colnames(Tab)))>0) { colnames(Tab) <- gsub("DRB3.4.5","DRB345",colnames(Tab))  } }
    colnames(Tab) <- toupper(colnames(Tab))
    All.ColNames <- gsub(".1","",colnames(Tab),fixed=T)
    rownames(Tab) <- NULL
    Tab <- rmABstrings(Tab)
    setwd(OutDir)
    
    # Separate DRB345 if exists as single column pair and check zygosity
    if(HLA==T) {
      if(sum(grepl("DRB345",colnames(Tab)))>0) {
        
        cat("Processing DRB3/4/5 column data.\n")
        
        DRBFLAG <- T
        getCol <- grep("DRB345",colnames(Tab))
        
        # Stop if not Locus*Allele formatting
        if( sum(grepl("\\*",Tab[,getCol]))==0 ) {
          Err.Log("Bad.DRB345.format")
          stop("Analysis Stopped.",call. = F)
        }
        
        # Expand DRB3/4/5 to separate column pairs
        Tab <- DRB345.parser(Tab)
        
        # Check for locus hemizygosity by DR haplotype
        Tab.list <- lapply(seq_len(nrow(Tab)),FUN=function(z) Tab[z,])
        tmp <- lapply(Tab.list,FUN=DRB345.zygosity)
        tmp.df <- lapply(tmp,"[[",1)
        Tab <- as.data.frame(do.call(rbind,tmp.df))
        All.ColNames <- gsub(".1","",colnames(Tab),fixed=T)
        
        #Identify DR345 flagged haplotypes
        DR.Flags <- list()
        
        tmp.DR3 <- do.call(rbind,lapply(tmp,"[[",2))
          if(length(which(tmp.DR3==T)>0)) {
            DR.Flags[['DRB3']] <- c("DRB3",paste(Tab[which(tmp.DR3==T),1],collapse=","))
          } else { DR.Flags[['DRB3']] <- NULL }
        
        tmp.DR4 <- do.call(rbind,lapply(tmp,"[[",3))
          if(length(which(tmp.DR4==T)>0)) { 
            DR.Flags[['DRB4']] <- c("DRB4",paste(Tab[which(tmp.DR4==T),1],collapse=",")) 
          } else { DR.Flags[['DRB4']] <- NULL }
        
        tmp.DR5 <- do.call(rbind,lapply(tmp,"[[",4))
          if(length(which(tmp.DR5==T)>0)) { 
            DR.Flags[['DRB5']] <- c("DRB5",paste(Tab[which(tmp.DR5==T),1],collapse=","))
          } else { DR.Flags[['DRB5']] <- NULL }
        
        DR.Flags <- do.call(rbind,DR.Flags)
        if(Output) {
          if(!is.null(DR.Flags)) {
            colnames(DR.Flags) <- c("Flagged.Locus","Sample.ID") ; rownames(DR.Flags) <- NULL
            Err.Log("Bad.DRB345.hap") ; cat("\n")
            write.table(DR.Flags,file="Flagged_DRB345_Haplotypes.txt",sep="\t",quote=F,row.names=F,col.names=T)
          }
        }
        cat("\n")
      } else { DRBFLAG <- F }
    } else { DRBFLAG <- NULL }
    
    # Separate locus and allele names if data is formatted as Loci*Allele
    if(HLA==T) {
      if(sum(regexpr("\\*",Tab[,3:ncol(Tab)]))>0) {
        Tab[,3:ncol(Tab)] <- apply(Tab[,3:ncol(Tab)],MARGIN=c(1,2),FUN=function(x) unlist(lapply(strsplit(x,split="\\*"),"[",2)))
      }
    }
    
  }

############################################################################################################################################    
### Case-Control Pre-Check #################################################################################################################
  
  cat(">>>> CASE - CONTROL SUMMARY STATISTICS\n")
  cat(paste(rep("_",50),collapse=""),"\n")
  if (Trim) { rescall <- paste(Res,"-Field",sep="") } else { rescall <- "Not Defined" }
  Check <- PreCheck(Tab,All.ColNames,rescall,HLA)
  if(Output) { write.table(Check,file="Data_Summary.txt",sep=": ",col.names=F,row.names=T,quote=F); rm(Check,rescall) }

############################################################################################################################################
### Data Processing and Sanity Checks ######################################################################################################
  
  cat(">>>> DATA PROCESSING AND CHECKS.\n")
  cat(paste(rep("_",50),collapse=""),"\n")
  
  ## __________________ General processing and checks for any data
  
  # MISSING DATA
  if(Missing == "ignore") {
    cat("Ignoring any missing data.\n")
    cat("Consider setting a missing threshold or running without the haplotype ('H') analysis.\n")
    rows.rm <- NULL
  } else {
    if (Missing > 2) {
      cat("The number of allowable missing may affect performance.\nConsider running with a smaller 'Missing' value or without the haplotype ('H') analysis.\ncontinuing......")
    }
    cat("Removing any missing data. This will affect Hardy-Weinberg Equilibrium test.\n")
    geno.desc <- haplo.stats::summaryGeno(Tab[,3:ncol(Tab)], miss.val=NAstrings)
    test <- geno.desc[,2] + 2*geno.desc[,3]
    rows.rm <- which(test > Missing)
    if( length(rows.rm) > 0 ) {
      rows.rm <- which(test > Missing)
      ID.rm <- Tab[rows.rm,1]
      Tab <- Tab[-rows.rm,]
      write.table(ID.rm, file="Removed_SampleIDs.txt", sep="\t", row.names=F, col.names=F, quote=F)
    }
    rm(geno.desc,test,ID.rm)
    if(nrow(Tab)==0) { Err.Log("Bad.Format.HLA") ; stop("Analysis Stopped.",call. = F) }
  }

  
  # LOCI SET DEFINITIONS
  if (missing(Loci.Set)) {
    Set <- list(c(3:ncol(Tab)))
  } else {
    if(CheckLoci(unique(All.ColNames[3:ncol(Tab)]),Loci.Set)$Flag) {
      Err.Log("Bad.Locus.NA")
      stop("Analysis Stopped.",call. = F)
    } else {
      if(sum(grepl("All",Loci.Set))>0) { Loci.Set[[which(Loci.Set=="All")]] <- unique(All.ColNames[3:ncol(Tab)])  } 
      Set <- lapply(Loci.Set,FUN=function(x)seq(1,ncol(Tab),1)[All.ColNames %in% x]) }
  }
  
  # MULTICORE LIMITATIONS
  if (Cores.Lim!=1L) { 
    if(Sys.info()['sysname']=="Windows" & as.numeric(Cores.Lim)>1) {
      cat("You seem to be using Windows and specified a number of Cores > 1. Please see vignette.\n") }
    Cores = as.integer(Cores.Lim)
  }
  
  # RUN TESTS DEFINITIONS
  if (missing(Run.Tests)) { Run <- c("HWE","H","L","A") } else { Run <- Run.Tests }
  
  ## __________________ HLA specific checks
  if(Trim & !HLA) { cat("Trimming only relevant to HLA data, no trimming performed.\n") }
  if(EVS.rm & !HLA) { cat("Expression variant suffix stripping only relevant to HLA data, no stripping performed.\n") }
  
  if(HLA) {
    
    if(Trim | EVS.rm | "A" %in% Run ) { cat("Running HLA specific functions...\n") }
    
    # Sanity Check for Resoltion if Trim="T" and Trim Data
    if(Trim & CheckHLA(Tab[,3:ncol(Tab)])) {
      cat("--Trimming Data.\n")
      Tab.untrim <- Tab
      Tab[,3:ncol(Tab)] <- apply(Tab[,3:ncol(Tab)],MARGIN=c(1,2),GetField,Res=Res)
      rownames(Tab) <- NULL
    } else if (Trim) {
      Err.Log("Bad.Format.Trim")
      stop("Analysis Stopped.",call. = F)
    }
    
    # Sanity Check for Expresion Variant Suffix Stripping
    if(EVS.rm & CheckHLA(Tab[,3:ncol(Tab)])) {
      cat("--Stripping Expression Variants Suffixes.\n")
      Tab[,3:ncol(Tab)] <- apply(Tab[,3:ncol(Tab)],MARGIN=c(1,2),gsub,pattern="[[:alpha:]]",replacement="")
      EVS.loci <- as.list(names(EPL))
      EPL <- lapply(EVS.loci,EVSremoval,EPList=EPL)
      names(EPL) <- EVS.loci ; rm(EVS.loci)
    } else if (EVS.rm) {
      Err.Log("Bad.Format.EVS")
      stop("Analysis Stopped.",call. = F)
    }
    
    if ("A" %in% Run) {
      
      Release <- as.character(unlist(EPL[['Release']]))
      
      # Sanity Check for Known HLA loci
      cat(paste("--Checking loci against ",Release,".\n",sep=""))
      test <- CheckLoci(names(EPL),unique(All.ColNames[3:ncol(Tab)]))
      if( test$Flag ) {
        Err.Log("Bad.Locus.HLA")
        cat("Problem loci:",test$Loci,"\n")
        stop("Analysis stopped.",call. = F)
      }
      
      # Sanity Check for Known HLA alleles
      cat(paste("--Checking alleles against ",Release,".\n",sep=""))
      test <- CheckAlleles(EPL, Tab[,3:ncol(Tab)], unique(All.ColNames[3:ncol(Tab)]), All.ColNames[3:ncol(Tab)])
      if(sum(unlist(lapply(test,"[[",1)))>0) {
        Err.Log("Bad.Allele.HLA")
        tmp <- as.character(unlist(lapply(test[which(lapply(test,"[[",1)==T)],"[",2)))
        cat("Problem alleles:",tmp,"\n")
        stop("Analysis stopped.",call. = F)
      }
    }
    
  } # End HLA if statement and HLA specific functionalities

############################################################################################################################################    
### Write to Parameter File ################################################################################################################
  
  if(Output) {
    
    Params.Run <- list(Time = format(Sys.time(), "%a %b %d %X %Y"),
                       BD.Version = as.character(packageVersion("BIGDAWG")),
                       Cores.Available = Cores,
                       File = Data,
                       HLA.Data = HLA,
                       DRB345.Parsed = DRBFLAG,
                       Tests = paste(Run,collapse=","),
                       All.Pairwise = All.Pairwise,
                       Trim = Trim,
                       Resolution = ifelse(Trim,Res,NA),
                       Suffix.Stripping = EVS.rm,
                       Missing = Missing,
                       Samples.Removed = length(rows.rm))
    
    Params.Run <- do.call(rbind,Params.Run)
    write.table(Params.Run,file="Run_Parameters.txt",sep=": ", row.names=T, col.names=F, quote=F)
  }

############################################################################################################################################    
### Hardy Weignberg Equilibrium ############################################################################################################
  
  if ("HWE" %in% Run) {
    
    cat("\n>>>> STARTING HARDY-WEINBERG ANALYSIS...\n")
    cat(paste(rep("_",50),collapse=""),"\n")
    if(Trim) { cat("HWE performed at user resolution on controls.\n") } else { cat("HWE performed on controls.\n") }
    HWE <- HWEChiSq(Tab,All.ColNames)
    if(Output) { write.table(HWE,file="HWE.txt",sep="\t",col.names=T,row.names=F,quote=F) }
    
    cat("\n> HARDY-WEINBERG ANALYSIS COMPLETED\n")
    HWE <- as.data.frame(HWE)
    print(HWE,row.names=F,quote=F)
    cat("\n")
    
    rm(HWE)
    
  } #END HARDY-WEINBERG

############################################################################################################################################
####Set Loop Begin (loop through each defined locus/loci set) ##############################################################################
  
  if ( sum( c("H","L","A") %in% Run ) > 0 ) {
  
    cat("\n>>>>>>>>>>>>>>>>>>>>>>>>> Begin Locus Sets <<<<<<<<<<<<<<<<<<<<<<<<<\n\n")
    cat(paste("Your analysis has ", length(Set), " set(s).", sep=""),"\n")
    
    for(k in 1:length(Set)) {
      cat("\n")
      cat(paste(rep(">",35),collapse=""),"Running Set",k,"\n")
      
      cols <- Set[[k]]
      Tabsub <- Tab[,c(1,2,cols)]

############################################################################################################################################  
### Per Set Global Variables ###############################################################################################################
      
      SID <- Tabsub[,1] # sample IDs
      genos <- Tabsub[,3:ncol(Tabsub)] # genotypes
      genos[genos==""] <- NA
      grp <- Tabsub[, 2] # phenotype
      nGrp0 <- length(which(grp==0))*2 #nalleles
      nGrp1 <- length(which(grp==1))*2 #nalleles
      loci <- unique(gsub(".1","",colnames(genos),fixed=T)) # name of loci
      loci.ColNames <- gsub(".1","",colnames(genos),fixed=T) # column names
      nloci <- as.numeric(length(loci)) # number of loci
      
      if(Output) {
        
        OutSetDir <- paste(OutDir,"/set",k,sep="")
        dir.create(OutSetDir)
        setwd(OutSetDir)
        
        Params.set <- list(Set = paste("Set",k),
                       Loci.Run = paste(loci,collapse=","))
        
        Params.set <- do.call(rbind,Params.set)
        write.table(Params.set,file="Set_Parameters.txt",sep=": ", row.names=T, col.names=F, quote=F)
      }
      
      SAFE <- c(ls(),"SAFE")

############################################################################################################################################    
### Haplotype Analysis #####################################################################################################################
      
      if ("H" %in% Run) {
        
        cat(">>>> STARTING HAPLOTYPE ANALYSIS...","\n")
        cat(paste(rep("_",50),collapse=""),"\n")
        
        # Sanity check for set length and All.Pairwise=T
        if (nloci<2) {
          Err.Log("Loci.No")
          stop("Analysis Stopped.", call. = F)
        } else if (All.Pairwise & nloci<=2)  {
          Err.Log("Loci.No.AP")
          stop("Analysis Stopped.", call. = F) }
        
        HAPsets <- list() ; HAPsets[['Set.AllLoci']] <- genos
        
        # Define Pairwise Combinations to Run When Selected
        if(All.Pairwise) {
          
          # Define Combinations
          Combos <- t(combn(loci,2))
          
          cat("\nYou have opted to include all pairwise combinations for the haplotype analysis.\n")
          cat("There are", nrow(Combos), "possible locus combinations to run.\n" )
          
          # Define Pairwise Sets
          for(s in 1:nrow(Combos)) {
            Set.H <- loci.ColNames %in% Combos[s,]
            HAPsets[[paste(Combos[s,1],Combos[s,2],sep="_")]] <- genos[,c(Set.H)]
          }; rm(s)
          
        }
        
        for(h in 1:length(HAPsets)) {
          
          # Get haplotype loci set for analysis
          genos.sub <- HAPsets[[h]]
          
          # Run Analysis
          H.list <- H(genos.sub,grp)
          
          if(Output) {
            
            # File names for output
            name1 <- paste("haplotype_freqs.",names(HAPsets)[h],".txt",sep="")
            name2 <- paste("haplotype_binned.",names(HAPsets)[h],".txt",sep="")
            name3 <- paste("haplotype_OR.",names(HAPsets)[h],".txt",sep="")
            name4 <- paste("haplotype_chisq.",names(HAPsets)[h],".txt",sep="")    
            name5 <- paste("haplotype_table.",names(HAPsets)[h],".txt",sep="") 
            name6 <- paste("haplotype_bySubject.",names(HAPsets)[h],".txt",sep="") 
            
            ## write to file
            write.table(H.list[['freq']], name1, sep="\t", quote = F, row.names=F, col.names=T)
            write.table(H.list[['binned']], name2, sep="\t", quote = F, row.names=F, col.names=T)
            write.table(H.list[['OR']], name3, sep="\t", quote = F, row.names=F, col.names=T)
            write.table(H.list[['chisq']], name4, sep="\t", row.names = F, quote = F)
            write.table(H.list[['table']], name5, sep="\t", row.names = F, quote = F)
            
            Haplotypes <- cbind(SID,H.list[['Haplotypes']])
            colnames(Haplotypes)[1] <- colnames(Tabsub)[1]
            write.table(Haplotypes, name6, sep="\t", row.names = F, quote = F)
            
          }
          
          cat("\n> HAPLOTYPE ANALYSIS COMPLETED:",names(HAPsets)[h],"\n")
          overall.chisq <- H.list[['chisq']]
          overall.chisq$X.square <- round(as.numeric(levels(overall.chisq$X.square)),digits=5)
          print(overall.chisq, row.names=F)
          cat("\n")
          
        }# END hapset Loop
        
        rm(list=ls()[!(ls() %in% SAFE)])
        
      } #END HAPLOTYPE

############################################################################################################################################  
### Locus Level ############################################################################################################################
      
      if ("L" %in% Run) {
        
        cat("\n>>>> STARTING LOCUS LEVEL ANALYSIS...\n")
        cat(paste(rep("_",50),collapse=""),"\n")
        
        Allele.binned <- list() # Alleles binned during chi-square test
        Allele.freq <- list() # Alleles Frequencies
        overall.chisq <- list() # Chi-Square Table
        ORtable <- list() # Odds Ratio Table
        Final_binned <- list() # Contingency Table
        
        for(j in 1:nloci) {
          
          # Get Locus
          Locus <- loci[j]
          
          # Run Locus Level Analysis
          L.list <- L(loci.ColNames,Locus,genos,grp,nGrp0,nGrp1)
          
          # Build Output Lists
          Allele.binned[[Locus]] <- L.list[['binned']]
          Allele.freq[[Locus]] <- L.list[['freq']]
          overall.chisq[[Locus]] <- L.list[['chisq']]
          ORtable[[Locus]] <- L.list[['OR']]
          Final_binned[[Locus]] <- L.list[['table']]
          
        }# END locus
        
        if(Output) {
          ## write to file
          write.table(do.call(rbind,Allele.freq), file = paste("Locus_freqs.txt",sep=""), sep="\t", row.names = F, col.names=T, quote = F)
          write.table(do.call(rbind,Final_binned), file = paste("Locus_table.txt",sep=""), sep="\t", row.names = F, col.names=T, quote = F)
          write.table(do.call(rbind,Allele.binned), file = paste("Locus_binned.txt",sep=""), sep="\t", row.names = F, col.names=T, quote = F)
          write.table(do.call(rbind,ORtable), file = paste("Locus_OR.txt",sep=""), sep="\t", row.names = F, col.names=T, quote = F)
          write.table(do.call(rbind,overall.chisq), file = paste("Locus_chisq.txt",sep=""), sep="\t", row.names = F, col.names=T, quote = F)  
        }
        
        cat("\n> LOCUS LEVEL ANALYSIS COMPLETED","\n")
        overall.chisq <- do.call(rbind,overall.chisq)
        print(overall.chisq,row.names=F)
        cat("\n")
        
        rm(list=ls()[!(ls() %in% SAFE)])
        
      } #END LOCUS

##########################################################################################################################################
### Amino Acid Level #######################################################################################################################
      
      if(HLA) {
        if ("A" %in% Run) {
          
          genos[genos=='^'] <- "00:00"
          
          cat("\n>>>> STARTING AMINO ACID LEVEL ANALYSIS...\n")
          cat(paste(rep("_",50),collapse=""),"\n")
          
          if( UPL.flag ) { cat("Using updated protein exon alignments.\n") }
          
          # Amino Acid Analysis Sanity Checks
          if(Res<2 | !CheckHLA(genos))  {
            cat("You have opted to run the amino acid analysis.\n")
            Err.Log("Low.Res")
            stop("Analysis stopped.",call. = F)
          }
          
          # Define Lists for Per Loci Running Tallies
          AAlog <-  list()
          AminoAcid.binned <- list()    
          AminoAcid.freq <- list()    
          overall.chisq <- list()    
          ORtable <- list()    
          Final_binned <- list()
          
          # Loop Through Loci
          for(x in 1:nloci){
            
            # Get Locus
            Locus <- loci[x]
            
            # Read in Locus Alignment file for Locus specific exons
            ExonAlign <- EPL[[Locus]]; rownames(ExonAlign) <- NULL
            
            # Run Amino Acid Analysis
            A.list <- A(loci.ColNames,Locus,genos,grp,nGrp0,nGrp1,ExonAlign,Cores)
            
            # Build Output Lists
            AAlog[[Locus]] <- A.list[['log']]
            AminoAcid.binned[[Locus]] <- A.list[['binned']]
            overall.chisq[[Locus]] <- A.list[['chisq']]
            ORtable[[Locus]] <- A.list[['OR']]
            Final_binned[[Locus]] <- A.list[['table']]
            AminoAcid.freq[[Locus]] <- A.list[['freq']]
            
          }; rm(x) #locus loop
          
          if(Output) {
            ## write to file
            write.table(do.call(rbind,AAlog), file = "AA_log.txt", sep="\t", row.names = F, col.names=T, quote = F)
            write.table(do.call(rbind,AminoAcid.freq), file = "AA_freqs.txt", sep="\t", row.names = F, col.names=T, quote = F)
            write.table(do.call(rbind,AminoAcid.binned), file = "AA_binned.txt", sep="\t", row.names = F, col.names=T, quote = F)
            write.table(do.call(rbind,ORtable), file = "AA_OR.txt", sep="\t", row.names = F, col.names=T, quote = F)
            write.table(do.call(rbind,overall.chisq), file = "AA_chisq.txt", sep="\t", row.names = F, col.names=T, quote = F)
            write.table(do.call(rbind,Final_binned), file = "AA_table.txt", sep="\t", row.names = F, col.names=T, quote = F)
            write.table(Release, file = "Set_Parameters.txt", sep="\t", row.names = F, col.names=F, quote = F, append=T)
          }
          
          cat("\n> AMINO ACID ANALYSIS COMPLETED\n")
          
          cat("Significant Amino Acid Position(s):","\n")
          tmp <- do.call(rbind,overall.chisq); rownames(tmp) <- NULL
          tmp.sig <- tmp[which(tmp[,'sig']=="*"),]; rownames(tmp.sig) <- NULL
          if(nrow(tmp.sig)>0) { print(tmp.sig,row.names=F,quote=F) }
          
          rm(list=ls()[!(ls() %in% SAFE)])
          
        } #END AMINO ACID
      }

############################################################################################################################################  
### End Analyses ###########################################################################################################################
      
    }; rm(k)
    
##  Set Loop End
############################################################################################################################################  
  
  }# END SET LOOP IF
  
  cat("\n>>>>>>>>>>>>>>>>>>>>>>>>>> End Analysis <<<<<<<<<<<<<<<<<<<<<<<<<<\n")
  
  options(warn=0)
  
}# END FUNCTION
