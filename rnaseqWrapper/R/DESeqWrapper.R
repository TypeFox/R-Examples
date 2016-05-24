DESeqWrapper <-
function(countData, # count file as numeric matrix/data.frame
                         conditions, # Character vector with the treatment groups
                         whichGeneNames = 0, # Which column has the gene names, 0 is rownames (preferred)
                         outNamePrefix = "DESeqOutputs", # name to prepend to written files
                         comps = "allPairwise", # character vector of which contrasts to run in the format 1vs2 with numbers matching order of conditions above
                         ## I should fix this to make it more flexible and include by name
                         conds = NULL, # character vector of the conditions for each column. If "NULL" will use grep on conditions
                         colorSet = NULL, # Colors to use for groups
                         makePDFs = TRUE, # should plots be saved as pdfs?
                         writeScaled = FALSE, # Should the corrected data be written?
                         writeDE = TRUE, # Should the DE data be written?
                         pCut = 0.05, ## What corrected Pvalue should be used for significance?
                         dispMethod = "pooled",
                         dispSharingMode="maximum"
         ){

  if(!require(DESeq)){
    stop("The Bioconductor package 'DESeq' is required to run this function. ",
         "Please install it now using:\n\n",
         'source("http://bioconductor.org/biocLite.R")\n',
         'biocLite("DESeq")\n')
  }
  
  ## Returns a list with:
  ##    list of DE results (each a data.frame)
  ##    the normalized reads
  ##    a data.frame of whether each gene is sig for each comparison
  
  if(is.null(colorSet)){
    colorSet <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF")  
  }
  
  ## Prep the output folder, if needed:
  dirsToMake <- unlist(strsplit(outNamePrefix,"/"))
  if(length(dirsToMake)>1){
    for(k in 1:(length(dirsToMake)-1)){
      suppressWarnings(try(dir.create(paste(dirsToMake[1:k],collapse="/"))))
    }
  }
  if(length(grep("/$",outNamePrefix))==1){
    suppressWarnings(try(dir.create(paste(dirsToMake[1:length(dirsToMake)],collapse="/"))))
  }
  
  conditions <- as.character(conditions)
  
  # Set the conds for each sample
  if(is.null(conds)){
    conds <- character(length=length(colnames(countData)))
    for(k in 1:length(conditions)){
      conds[grep(conditions[k],colnames(countData))] <- conditions[k]
    }
  }
  conds <- factor(conds,levels=conditions)
  
  ## Set the colors for plotting
  condCol <- conds
  levels(condCol) <- colorSet[1:length(levels(condCol))]
  
  ## Set rownames if not set
  if (whichGeneNames == 0 | whichGeneNames == "rownames" | whichGeneNames == "row.names" | whichGeneNames == "" | is.null(whichGeneNames)){
    
  } else if(is.numeric(whichGeneNames) & whichGeneNames > 0){
    rownames( countData ) <- countData[,whichGeneNames]  
    countData <- countData[,-whichGeneNames]
  } else if(is.character(whichGeneNames)){
    whichCol <- grep(whichGeneNames,colnames(countData))
    if(length(whichCol)==1){
      rownames( countData ) <- countData[,whichCol]
      countData <- countData[,-whichCol]
    } else if(length(whichCol)>1){
      stop("whichGeneNames matches more than one column")
    } else if(length(whichCol)<1){
      stop("whichGeneNames does not match any columns")
    }
  } else {
    stop("whichGeneNames must be numeric or character")
  }
  
  ## Prepare the countFile
  message("preparing count table\n")
  cds <- newCountDataSet( countData, conds )
  cds <- estimateSizeFactors( cds )
  normcds <- round(as.data.frame(counts( cds, normalized=TRUE )))
  
  ## Write the file, if desired:
  if(writeScaled){
    normCountFileName <- paste(outNamePrefix,"normalizedReadCounts.txt",sep="") # make a new file, saves the folder information
    write.table(normcds, file=normCountFileName, quote=FALSE, sep="\t", row.names=TRUE)    
  }

  ## Final step to prep
  #cds <- estimateDispersions( cds, method="pooled-CR" )
  cds <- estimateDispersions( cds , dispMethod, dispSharingMode)
  
  
  ## catch extra things that likely mean allPairwise:
  if(comps==0|comps==""|is.null(comps)){
    comps <- "allPairwise"
  }
  
  # make all comparisons
  if(comps=="allPairwise"){
    comps <- character()
    for(k in 1:length(conditions)){
      for(j in 2:length(conditions)){
        if(j>k){
          comps <- c(comps,paste(k,"vs",j,sep="")  )
        }
      }
    }
  }
  
  
  ## Run the Differential expression tests of interest
  message("\n\nRunning Differential Expression Tests\n")
  deOuts <- list()
  ## Initalize a data.frame to tell if sig
  isSig <- data.frame(row.names=rownames( countData ))
  for (k in 1:length(comps)){
    
    condsInt <- as.numeric(unlist(strsplit(comps[k], "vs")))
    cond1 <- conditions[condsInt[1]]
    cond2 <- conditions[condsInt[2]]
    
    outname <- paste(outNamePrefix,cond1,"vs",cond2,sep="")
    message(paste("\n\tTesting: ",cond1,"vs",cond2,"\n",sep=" "))
    
    res <- nbinomTest( cds, cond1, cond2 )
    names(res)[3] <- paste('baseMean',cond1,sep='_')
    names(res)[4] <- paste('baseMean',cond2,sep='_')
    isSig[[paste(cond1,"vs",cond2,"_isSig",sep="")]] <- (res$padj<pCut & !is.na(res$padj) )
    
    
    deglist <- res[ order(res$pval), ]

    # 5% FDR
    resSig <- res[ res$padj < pCut & !is.na(res$padj), ]
    resSig <- resSig[ order(resSig$pval), ]
    nSig <- length(resSig[,1])
    
    deOuts[[k]] <- res
    names(deOuts)[k] <- paste(cond1,cond2,sep="vs")
    if(writeDE){
      message("\t\tWriting DEgeneList")
      write.table(res, file=paste(outname,"DESeq_output.txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE)
      write.table(resSig, file=paste(outname,"DEG_output_sigGenes.txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE)
    }
    
    ## Do visualizations
    if(makePDFs){
      message("\t\tPreparing visualization Outputs\n")
     
      pdf(paste(outname,"_DEseqVis.pdf",sep=""),title=paste(cond1,cond2,sep=" vs "))
      
#       plotDispEsts( cds ,cond1,cond2)
      plotDE( res ,cond1,cond2)
      hist(res$pval, breaks=100, col="skyblue", border="slateblue", main=paste("Histogram of p-values: ",cond1," vs ",cond2,sep=""), xlab="p-value")
      
      cdsBlind <- estimateDispersions( cds, method="blind" )
      vsd <- getVarianceStabilizedData( cdsBlind )
      
      ## Plots of  'direct (lfc) versus moderated log-ratios (mod_lfc)'
#       mod_lfc <- (rowMeans( vsd[, conditions(cds)==cond1, drop=FALSE] ) - rowMeans( vsd[, conditions(cds)==cond2, drop=FALSE] ))
#       lfc <- res$log2FoldChange
#       finite <- is.finite(lfc)
#       table(as.character(lfc[!finite]), useNA="always")
#       largeNumber <- 10
#       lfc <- ifelse(finite, lfc, sign(lfc) * largeNumber)
#       
#       logdecade <- 1 + round( log10( 1+rowMeans(counts(cdsBlind, normalized=TRUE)) ) )
#       colors <- colorRampPalette( c( "gray", "blue" ) )(6)[logdecade]
#       plot( lfc, mod_lfc, pch=20, cex=.4, asp=1,
#             main=paste("Scatterplot of direct (lfc) versus moderated log-ratios (mod_lfc): ",cond1," vs ",cond2,sep=""), 
#             col = ifelse( finite, colors, "purple" ) )
#       abline( a=0, b=1, col="#40C04040" )
#       
      #plot( lfc, mod_lfc, pch=20, cex=.3, col = ifelse( finite, "#80808040", "red" ) )
      #abline( a=0, b=1, col="#40404040" )
      
      select <- order(res$pval)[1:100]
      colors <- colorRampPalette(c("white","darkblue"))(100)
      #colors <- colorpanel(75,"medium blue","black","yellow")
      heatmap.mark(vsd[select,], col=colors, cexCol = 0.75, trace="none", key=TRUE, symkey=FALSE, density.info="none", scale = "row", labRow = FALSE,ColSideColors=as.character(condCol),scaleLabel="",main='Top 100 DE' )
      legend(x="topleft",inset=c(-.01,.13),bty="n",legend=levels(as.factor(conds)),fill=levels(condCol),cex=.5,title="Conditions")
      
      if(nSig > 3){
        select2 <- order(res$pval)[1:nSig]
        heatmap.mark(vsd[select2,], col=colors, cexCol = 0.75, trace="none", key=TRUE, symkey=FALSE, density.info="none", scale = "row", labRow = FALSE,ColSideColors=as.character(condCol),scaleLabel="",main='All Significant differences' )
        legend(x="topleft",inset=c(-.01,.13),bty="n",legend=levels(as.factor(conds)),fill=levels(condCol),cex=.5,title="Conditions")
      }
      
      dists <- dist( t( vsd ) )
      heatmap( as.matrix( dists ), symm=TRUE, scale="none", margins=c(10,10), col = colors, labRow = paste( pData(cdsBlind)$condition, pData(cdsBlind)$type ) )
      
      dev.off()
      
    }
    
    
  }
  
  if(writeDE | writeScaled){
    message("Files written to: ",outNamePrefix, sep=" ")
  }
  
  out <- list(deOutputs = deOuts, normalizedReads = normcds, isSignificant = isSig)
  return(out)
}
