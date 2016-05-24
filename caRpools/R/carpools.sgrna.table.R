carpools.sgrna.table = function(wilcox=NULL, deseq=NULL, mageck=NULL, dataset=NULL, dataset.names = NULL, namecolumn=1, fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), type="enriched", cutoff.deseq = 0.05, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, plot.genes="overlapping", cutoff.hits=NULL, sgrna.file=NULL, labelgenes=NULL, write=FALSE, datapath=getwd(), analysis.name="Screen")
{
  #requireNamespace(xlsx)
  # Call function to gene hit list
  if(!is.null(labelgenes))
  {
    df.output = labelgenes
  }
  else
  {
    df.output = generate.hits(wilcox=wilcox, deseq=deseq, mageck=mageck, type=type, cutoff.deseq = cutoff.deseq, cutoff.wilcox = cutoff.wilcox, cutoff.mageck = cutoff.mageck, cutoff.override=cutoff.override, cutoff.hits=cutoff.hits, plot.genes=plot.genes)
    
  }
 
  #print(df.output)
  # Call plot.read.count.vs to create paired plots
  untreated.list = list(dataset[[1]], dataset[[2]])
  treated.list = list(dataset[[3]],dataset[[4]])
  
  
  if(length(df.output) >=1)
  {
    
    # Sorting 
    df.output = sort(df.output[!is.na(df.output)], na.last=TRUE, decreasing=FALSE) 
    #print(df.output)
    # Open file for writing
    if(identical(write,TRUE))
    {
      xlsx::write.xlsx(df.output, file=paste(datapath, paste(analysis.name, paste("HITS-sgRNA", paste(type, "xls", sep="."), sep="-"), sep="_"), sep="/"), sheetName="List of Genes")
    }
      
    
    for(i in 1:length(df.output))
    {
      
      # Output will be: Table with ALL sgRNAs for that gene including log2 Foldchange information AND Sequence of the sgRNA
      if(is.null(sgrna.file))
      {
        stop("No sgRNA file provided. sgrna.file is empty.")
      }

      # Normalize readcount using the normnalising function
      # Then MEAN the replicates
      untreated.list.mean=apply(
        do.call(
          "cbind",
          lapply(
            untreated.list,
            function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1  
            
          )
        ),
        1,
        mean)
      treated.list.mean=apply(
        do.call(
          "cbind",
          lapply(
            treated.list,
            function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
          )
        ),
        1,
        mean
      )
      
      # get gene names
      gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
      # get design name
      designs=treated.list[[1]][,namecolumn]
      
      
      dataset.combined <- data.frame( 
        untreated = untreated.list.mean,
        treated = treated.list.mean,
        genes=gene.names,  
        log2foldchange=log2(as.numeric(treated.list.mean)/as.numeric(untreated.list.mean)),  
        row.names = designs,
        sgRNA = as.character(designs),
        stringsAsFactors=FALSE 
      )
      
      # add sgRNA Target Sequence
      #dataset.combined$sequence = apply(dataset.combined, 1, function(z)
      #{
      # get sequence from sgrna.file data frame and paste it to dataset.combined
      #str(sgrna.file)
      # return(sgrna.file[as.character(sgrna.file[,"names"]) == as.character(z["sgRNA"]),"sequence"])
      #})
      #str(dataset.combined)
      # plot data
      
        dataset.table = data.frame(
          designs = dataset.combined[dataset.combined$genes==df.output[i],"sgRNA"],
          genes = dataset.combined[dataset.combined$genes==df.output[i],"genes"],
          stringsAsFactors=FALSE)
        
        dataset.table$sequence = apply(dataset.combined[dataset.combined$genes==df.output[i],],1, function(y)
        {
          toreturn = sgrna.file[sgrna.file$names == y["sgRNA"] ,"sequence"]
          if(toreturn == "")
          {
            stop(paste("Reference file and Read Count files do not match for", y["sgRNA"] , sep=" "))
          }
        else { return(toreturn)}
        })
        
        dataset.table$log2Foldchange = apply(dataset.combined[dataset.combined$genes==df.output[i],],1, function(y)
        {
          return(log2(dataset.combined$treated[dataset.combined$sgRNA==y["sgRNA"]]/dataset.combined$untreated[dataset.combined$sgRNA==y["sgRNA"]]))
        })
      
      # Sort according to fold change
      dataset.table = dataset.table[order(dataset.table$log2Foldchange, decreasing = TRUE),]
      
      if(identical(write,TRUE))
      {
        xlsx::write.xlsx(dataset.table, file=paste(datapath, paste(analysis.name, paste("HITS-sgRNA", paste(type, "xls", sep="."), sep="-"), sep="_"), sep="/"), sheetName=as.character(df.output[i]), append=TRUE )
        
      }
     
       # write.xlsx(dataset.table, file=paste(datapath, paste(analysis.name, paste("HITS-sgRNA", paste(type, "xls", sep="."), sep="-"), sep="_"), sep="/"), sheetName=df.output[i], append=TRUE)
        return(dataset.table)
        
      
      #dev.off()
      #par(mfrow=c(1,1),las=1) 
      
      # return data table for this gene

    }
    

  }

}