# Tumorboost normalization for normal-tumor study
# @title tumorboost normalization
# @param dataSetName name of the data set.
# @param normalTumorArray  A data.frame with 2 columns: "normal" and "tumor".
# The first column contains the names of normal files and the second the names of associated tumor files.
# @param plot if TRUE, the function will produced graphics
#
#
# This function is mainly based on the code from http://aroma-project.org/vignettes/tumorboost-highlevel from Pierre Neuvial
#
#
tumorboost<-function(dataSetName,normalTumorArray,plot=TRUE)
{
  require(aroma.cn)
  require(aroma.core)
  require(R.devices)
  require(R.filesets)
  
  log <- verbose <- Arguments$getVerbose(-1, timestamp=TRUE);
  rootPath <- "totalAndFracBData";
  rootPath <- Arguments$getReadablePath(rootPath);
  dataSet <- paste0(dataSetName,",ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY");
  
  
  ######################################################
  ####################   Data Import
  ######################################################
  
  ##########  fraction allele B
  ds <- aroma.core::AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
    
  # Identify all sample names by extracting all unique arrays name (ignoring tags)
  #sampleNAMES = getId(R.filesets::getNames(ds),normalTumorArray) 
  normalTumorMatrix=getNormalTumorMatrix(R.filesets::getNames(ds),normalTumorArray)
  #uniqueNAMES <- sort(unique(sampleNAMES))
  
  ########## import total copy number
  if(plot)
  {
    dsC <- aroma.core::AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath);
    #print(dsC);
    #sampleNAMESC <-getId(R.filesets::getNames(dsC),normalTumorArray);
    normalTumorMatrixC=getNormalTumorMatrix(R.filesets::getNames(dsC),normalTumorArray)
  }
  
  ##########  genotypage
  
  # Identify available genotype calls
  rootPath <- "callData";
  rootPath <- Arguments$getReadablePath(rootPath);
  
  genotypeTag <- "NGC";
  gsN <- aroma.core::AromaUnitGenotypeCallSet$byName(dataSet, tags=genotypeTag, chipType="*");
  
  # Keep only normal genotype files (not needed here, but could be needed in other situations)
  #types <- getStatus(R.filesets::getNames(gsN),normalTumorArray)
  
  #keep <- which(is.element(types, c("normal")));
  #gsN <- extract(gsN, normalSample);
  
  #get id of genotypae call file
  #genNAMES <- normalTumorMatrix[,1]
  
  
  #loop on all pairs (normal, tumor)
  for(indSample in 1:nrow(normalTumorMatrix))
  {
    normalSample=normalTumorMatrix[indSample,1]
    tumorSample=normalTumorMatrix[indSample,2]
    
    ########## frac allele B
    # Extract the two arrays with this name, which should be the tumor and the normal
        
    pair <- normalTumorMatrix[indSample,]
    stopifnot(length(pair) == 2);
      
    # Order as (tumor,normal)
    types <- getStatus(R.filesets::getNames(extract(ds,pair)),normalTumorArray)

    if(which(types=="normal")==1)
      o=2:1
    else
      o=1:2
        
    types <- types[o];
    pair <- pair[o];
      
    # Extract (tumor, normal) pair
    dsPair <- extract(ds, pair);
    dsT <- extract(dsPair, 1);
    #print(dsT);
    dsN <- extract(dsPair, 2);
    #print(dsN);

      
    ########## genotypage of normalSample
    gcN <- extract(gsN,which(gsN$Names==R.filesets::getNames(dsC)[normalSample]))
    #gcN2 <- extract(gsN,normalSample)
    
      
    ######################################################
    ####################   Apply TumorBoost Normalization
    ######################################################
      
    # Create a list of matched data sets
    dsList <- list(normal=dsN, tumor=dsT, callsN=gcN);
    # Make sure they are ordered the same way such that the
    # k:th array corresponds to the same sample in all sets.

    #print(dsList);
      
    #dummy <- lapply(dsList, FUN=function(ds) print(R.filesets::getFile(ds,1)));
      
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Normalize allele B fractions for tumors given matched normals
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tbn <- aroma.cn::TumorBoostNormalization(dsList$tumor, dsList$normal, gcN=dsList$callsN, tags=c("*", "NGC"));
    dsTN <- aroma.core::process(tbn, verbose=log);
    #print(dsTN);

    #if the user wants the different plot
    if(plot)
    {    
      #id <- getId(R.filesets::getNames(dsTN),normalTumorArray)
      dsTN=extract(dsTN,which(dsTN$Names==R.filesets::getNames(dsC[tumorSample])))
      
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Create a list of matched data sets
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dsList <- list(normal=dsN, tumor=dsT, tumorN=dsTN, callsN=gcN);
      rootPath <- "totalAndFracBData";
      rootPath <- Arguments$getReadablePath(rootPath);
      tumorboostPlot(ds,dsList,dataSetName,normalTumorArray,tumorSample,normalSample,dsC,normalTumorMatrixC)   
    }

    rm(dsN, dsT, dsTN, gcN);
    

  }##end loop files

}##end function

# plot part of the tumorboost function
#
# @param ds
# @param dsList list containing the fracB normal, tumor and the genotype call
# @param dataSetName name of the data set
# @param normalTumorArray  matrix containing information about about the name of files. The matrix has 3 columns.
# The first is the file name, the second the id to used for store the result and the third is the tag associated to the file (for example, "RC" or "TUMOR" ).
# @param tumorSample id of the current tumor file
# @param normalSample id of the current normal file
# @param dsC copynumber data
# @param sampleNAMESC names of the file of copy number
#
#
# This function is mainly based on the code from http://aroma-project.org/vignettes/tumorboost-highlevel from Pierre Neuvial
#
tumorboostPlot=function(ds,dsList,dataSetName,normalTumorArray,tumorSample,normalSample,dsC,normalTumorMatrixC)
{
  require(aroma.affymetrix)
  require(aroma.core)
  require(R.filesets)
  require(R.methodsS3)
  require(R.oo)
  
  ########## load the total copy number signal for a pair (normal,tumor) 
  pairC <- normalTumorMatrixC[normalTumorMatrixC[,2]==tumorSample,]
    
  stopifnot(length(pairC) == 2);
    
  # Order as (tumor,normal)
  types <- getStatus(R.filesets::getNames(extract(dsC,pairC)),normalTumorArray)
    
  if(which(types=="normal")==1)
    o=2:1
  else
    o=1:2
    
  types <- types[o];
  pairC <- pairC[o];
    
  # Extract (tumor, normal) pair for CN
  dsPairC <- extract(dsC, pairC);      

  
  ##########
  #path for save the plot
  figPath <- Arguments$getWritablePath(paste0("figures/",dataSetName,"/signal/"));
  siteTag <- R.filesets::getTags(ds);
  siteTag <- paste(siteTag[-1], collapse=","); 
  
  #
  ugp <- aroma.core::getAromaUgpFile(dsList$tumor);
  unf <- aroma.core::getUnitNamesFile(ugp);
  
  hg=R.filesets::getTags(ugp)[grep("hg",R.filesets::getTags(ugp))] #human genome reference 
  
  # prefix of SNP
  platform <- aroma.core::getPlatform(ugp);
  if (platform == "Affymetrix") 
  {
    require("aroma.affymetrix") || R.methodsS3::throw("Package not loaded: aroma.affymetrix");
    snpPattern <- "^SNP|^S-";
  } 
  else if (platform == "Illumina") 
  {
    snpPattern <- "^rs[0-9]";
  } 
  else 
  { 
    R.methodsS3::throw("Unknown platform: ", platform);
  }
  
  cat("Saving graphics for sample ",R.filesets::getNames(dsC[tumorSample]),"\n")
  
  #find ploidy for chromosome 23 and 24
  gender=findGender(R.oo::getName(dsC),normalSample,ugp)
  
  #loop on chr
  for(chromosome in 1:25)
  {
    chrTag <- sprintf("Chr%02d", chromosome);

    #check for the existence of the file
    figName <- sprintf("%s,%s", R.filesets::getNames(dsC[tumorSample]), chrTag);
    pathname <- filePath(figPath, sprintf("%s.png", figName));
    
    
    units <- aroma.core::getUnitsOnChromosome(ugp, chromosome=chromosome);
    unitNames <- aroma.core::getUnitNames(unf,units=units);##names of the SNP and CN probes
    posChr <- aroma.core::getPositions(ugp, units=units);#position of the probes on the chromosome
    
    ##########################################################
    ################  FRACTION ALLELE B  #####################
    ##########################################################
    # Identify SNP units
    keep <- (regexpr(snpPattern, unitNames) != -1);
    SNPunits <- units[keep];
    posSNP <- aroma.core::getPositions(ugp, units=SNPunits);#position des SNPs sur le chromosome      
    
    
    # Extract Allele B fractions (defined only for SNP probes)
    kk <- 1;
    dfList <- lapply(dsList, FUN=R.filesets::getFile, kk);
    beta <- lapply(dfList, FUN=function(df) df[SNPunits,1,drop=TRUE]);
    beta <- as.data.frame(beta);
    beta <- as.matrix(beta);
    names <- colnames(beta);
    names[names == "tumorN"] <- "normalized tumor";
    
    
    ##########################################################
    #####################  COPY NUMBER  ######################
    ##########################################################
    
    #tumor CN
    C <- R.filesets::extractMatrix(dsPairC, units=units);
    
    #the ploidy depends of the gender for the cromosome 23 and 24
    if(chromosome==23 || chromosome==24)
    {
      ploidy=NA
      if(gender=="XY")
        ploidy=1#male:chr23=chrX, chr24=chrY, the ploidy is 1
      else if(gender=="XX" && chromosome == 23)
      {
        ploidy=2#female:chr23=chrXX, the ploidy is 2
      }
      #female:chr24=nothing, the ploidy is NA
      
      #normalize by the normal sample        
      C <- ploidy*C[,1]/C[,2]
    }
    else
    {
      #normalization of CN
      C <- 2*C[,1]/C[,2];
    }
    
    print(head(C))
    print(head(posSNP))
    ####################### SORTIE GRAPHIQUE
    
    #dimensions
    x <- posChr/1e6;
    xlim <- range(x, na.rm=TRUE);
    xlab <- "Position (Mb)";    
    width <- 1280;
    aspect <- 0.6*1/3;
    ylim <- c(0,6);
    ylab <- "Copy number";
    
    # Plot total CNs
    fig <- R.devices::devNew("png", pathname, label=figName, width=width, height=4*aspect*width);
    par(mfrow=c(4,1))
    par(mar=c(2.7,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE);
    axis(side=1);
    axis(side=2, at=c(0,2,4,6));
    points(x, C, pch=".");
    label <- sprintf("%s", R.filesets::getNames(dsC[tumorSample]));
    stext(side=3, pos=0, label);
    stext(side=3, pos=1, chrTag);
    
    
    x <- posSNP/1e6;
    ylim <- c(-0.05,1.05);
    ylim <- c(-0.1,1.1);
    ylab <- "Allele B Fraction";
    cols <- as.integer(beta[,"callsN"] != 1) + 1L;
    #print(head(beta))
    
    # Plot Allele B fractions for normal, tumor and normalized tumor
    for (cc in 1:3)
    {
      tag <- colnames(beta)[cc];
      name <- names[cc];
      plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE);
      axis(side=1);
      axis(side=2, at=c(0,1/2,1));
      points(x, beta[,cc], pch=".", col=cols);
      label <- sprintf("%s (%s)", R.filesets::getNames(dsC[tumorSample]), name);
      stext(side=3, pos=0, label);
      stext(side=3, pos=1, chrTag); 
    }
    R.devices::devDone();
    cat("*")
    
    
  }##end boucle chromosome
}

