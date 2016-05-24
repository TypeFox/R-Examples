compare.analysis = function(wilcox=NULL, deseq=NULL, mageck=NULL, type="enriched", cutoff.deseq = NULL, cutoff.wilcox = NULL, cutoff.mageck = NULL, cutoff.override=TRUE, cutoff.hits=5, output="list", sort.by=c("mageck","pval","fdr"),plot.method=c("wilcox","mageck", "deseq"), plot.feature=c("pval","fdr","pval"), pch=16)
{
  #type = enriched | depleted
  #output = list | rank | venn
  #sort.by[1] = riger | wilcox | deseq 
  #sort.by[3] = pval | fc OR if sorting.output.method=DSS AUC | sgRNA


# input:
#   data: data from at least 2 data analysis -> check how many!
#   sortby: pval | fc -> to determine sorting for any of those
#   separate: TRUE | FALSE -> whether comparison is determined for enriched and depleted separately
#   cutoff: NUMERIC -> how many top candidates to use
#   sort.output: wilcox | deseq | mageck -> to which scoring the final output table will be sorted

# final returned output shall look like this

# status                gene        riger.rank  riger.pval  wilcox.rank   wilcox.pval   deseq.rank  deseq.pval    dss.rank  dss.AUC
# enriched/depleted     gene name   rank        pval        rank            pval            rank        pval          rank      AUC
#
# rank means at which position from most enriched/depleted did it occur in the screening set
# will be NA if not present within cutoff

# sort data separately for ENRICHED and DEPLETED

# Then sort by P.VAL (mageck, wilcox, Deseq)

###################################################

# first check which data sets are provided and put them into a new df,
# also make separate if desired before sorting
if(type=="enriched")
{
  
  
  # separate by enrichment 
  
  
  if(!is.null(wilcox))
  {
  # wilcox
  df.wilcox = wilcox[wilcox$foldchange >= 1 ,c("foldchange","p.value")]
  
  df.wilcox = df.wilcox[order(
    if(sort.by[2] == "pval")
    {df.wilcox$p.value}
    else { df.wilcox$foldchange } ,
    decreasing = 
      if(sort.by[2] == "pval")
      {FALSE}
    else { TRUE }
  ),]
  }
  
  if(!is.null(deseq))
  {
  
  #DESeq2
    # Deseq make sure it is -1 * foldchange, as deseq compare UNTREATED vs TREATED< all other treated vs untreated!
  df.deseq = as.data.frame(deseq$genes[deseq$genes$log2FoldChange >= 0 ,c("log2FoldChange","padj")])
  
  df.deseq = df.deseq[order(
    if(sort.by[2] == "pval")
    {df.deseq$padj}
    else { df.deseq$log2FoldChange } ,
    decreasing = 
      if(sort.by[2] == "pval")
      {FALSE}
    else { TRUE }
  ),]
  }
  

  
  if(!is.null(mageck))
  {
  # MAGeCK
    df.mageck = mageck$genes[, c("pos","rank.pos","sgrna.pos.good")]
    
    colnames(df.mageck) = c("fdr","rank","sgrna")
    df.mageck = df.mageck[order(
      df.mageck$rank,
      decreasing=TRUE,
      na.last=TRUE
    ),]
    
  }
  
  # Loading of data complete
  
  # Combine datasets using the cutoff
  # put in the rank of that gene in the data analysis table, if not present, put NA
  
  
} else if(type=="depleted")
{
  
  if(!is.null(wilcox))
  {
  
  # wilcox
  df.wilcox = wilcox[wilcox$foldchange < 1 ,c("foldchange","p.value")]
  df.wilcox = df.wilcox[order(
    if(sort.by[2] == "pval")
    {df.wilcox$p.value}
    else { df.wilcox$foldchange } ,
    decreasing = FALSE
  ),]
  }
  
  if(!is.null(deseq))
  {
  
  #DESeq2
  df.deseq = as.data.frame(deseq$genes[deseq$genes$log2FoldChange < 0 ,c("log2FoldChange","padj")])
  
  df.deseq = df.deseq[order(
    if(sort.by[2] == "pval")
    {df.deseq$padj}
    else { df.deseq$log2FoldChange } ,
    decreasing = FALSE
  ),]
  }
  
  

  
  if(!is.null(mageck))
  {
    # MAGeCK
    df.mageck = mageck$genes[, c("neg","rank.neg","sgrna.neg.good")]
    
    colnames(df.mageck) = c("fdr","rank","sgrna")
    
    df.mageck = df.mageck[order(
      df.mageck$rank,
      decreasing=TRUE,
      na.last=TRUE
    ),]
    
  }
   
} else {
  stop("No type selected. Please select enriched or depleted.")
}
#### Data loading complete

# add ranking just by length of dataset
#if(!is.null(riger)) {df.riger$rank = c(1:nrow(df.riger))}
if(!is.null(wilcox)) {df.wilcox$rank = c(1:nrow(df.wilcox))}
if(!is.null(deseq)) {df.deseq$rank = c(1:nrow(df.deseq))}
# create joining data.frame that has all gene names, use any of the above dataset
                    
if(!is.null(wilcox)) { v.join = rownames(wilcox)
                 
} else if(!is.null(deseq)) { v.join = rownames(deseq) 

                           
} else if(!is.null(mageck)) { v.join = rownames(mageck) 
                    
} else {stop("No datasets provided")}

# Create output data.frame
# rownames = gene names
# then
# EITHER p.value, foldchange of RIGER, wilcox, DESEQ + AUC of DSS
# OR rank of RIGER, wilcox, DESEQ, DSS that are within cutoff

# remove all genes that have NA on all of them
# remove columns that only have NAs (means not data was provided)

if(output=="list")
  {
  
  if(is.null(sort.by[1]))
  {stop("No sorting method applied!")}
  
    df.output = data.frame(
      genes = v.join,
      stringsAsFactors=FALSE)
    
    rownames(df.output) = df.output$genes


      # wilcox
      if(!is.null(wilcox)) {
        
        df.output$wilcox.log2fc=apply(df.output,1,function(i) return(log2(df.wilcox[i["genes"],"foldchange"])))
        df.output$wilcox.pval=apply(df.output,1,function(i) return(df.wilcox[i["genes"],"p.value"]))
        
      }

      # DESEQ2
      if(!is.null(deseq)) {
        
        df.output$deseq.log2fc=apply(df.output,1,function(i) return(df.deseq[i["genes"],"log2FoldChange"]))
        df.output$deseq.pval=apply(df.output,1,function(i) return(df.deseq[i["genes"],"padj"]))
        
      }
    
      # MAGeCK
      if(!is.null(mageck)) {
        
        df.output$mageck.fdr=apply(df.output,1,function(i) return(df.mageck[i["genes"],"fdr"]))
        df.output$mageck.rank=apply(df.output,1,function(i) return(df.mageck[i["genes"],"rank"]))
        
      }

    
    # Sorting 
  if (sort.by[3]!="genes" && sort.by[1] != "mageck")
    {
    df.output = df.output[order(df.output[,paste(sort.by[1],".",sort.by[3], sep="")], na.last=TRUE, decreasing=if(sort.by[3] == "fc") {TRUE} else {FALSE}),]
  } else if (sort.by[3]!="genes" && sort.by[1] == "mageck")
    {
    df.output = df.output[order(df.output[,paste(sort.by[1],".",sort.by[3], sep="")], na.last=TRUE, decreasing=FALSE),]
  } else
  {
  df.output = df.output[order(df.output[,sort.by[3]], na.last=TRUE, decreasing=FALSE),]
  }
      
  
 
    #df.output = df.output[-todelete,]
      
      # remove gene name column, we have rownames
    
  
    # Cutoff
    if(!is.null(cutoff.hits))
    {
      df.output = df.output[1:cutoff.hits,]
      
    }
    
    return(df.output)
    #return data

} else if(output == "rank") 
{

  if(is.null(sort.by[1]))
  {stop("No sorting method applied!")}
  
  df.output = data.frame(
    genes = v.join,
    stringsAsFactors=FALSE)
  
  rownames(df.output) = df.output$genes
  
  
  
  # wilcox
  if(!is.null(wilcox)) {
    
    df.output$wilcox=apply(df.output,1,function(i) return(df.wilcox[i["genes"],"rank"]))
    #df.output$wilcox.pval=apply(df.output,1,function(i) return(df.wilcox[i["genes"],"p.value"]))
    
  }
  
  # DESEQ2
  if(!is.null(deseq)) {
    
    df.output$deseq=apply(df.output,1,function(i) return(df.deseq[i["genes"],"rank"]))
    #df.output$deseq.pval=apply(df.output,1,function(i) return(df.deseq[i["genes"],"padj"]))
    
  }
  
  # MAGeCK
  if(!is.null(mageck)) {

    df.output$mageck=apply(df.output,1,function(i) return(df.mageck[i["genes"],"rank"]))
    
  }
  
  
  # Sorting
  
  df.output = df.output[order(df.output[,sort.by[1]], na.last=TRUE) ,]

      
      # remove gene name column, we have rownames
      df.output$genes=NULL

  # Cutoff
  if(!is.null(cutoff.hits))
  {
    df.output = df.output[1:cutoff.hits,]
    
  }
  #return data
  return(df.output)
  
#### FOR VENN DIAGRAMM  

}  else if(output == "venn")
  # CREATE VENN DIAGRAM
{  
  
  if(identical(cutoff.override,TRUE))
  {
    cutoff.wilcox = cutoff.hits
    cutoff.deseq = cutoff.hits
    cutoff.mageck = cutoff.hits
  }
  
  list.return=list()
 
  # wilcox
  if(!is.null(wilcox)) {
    
    venn.wilcox=df.wilcox[order(df.wilcox[,"p.value"], na.last=TRUE, decreasing=FALSE),]
    #df.output$wilcox.pval=apply(df.output,1,function(i) return(df.wilcox[i["genes"],"p.value"]))
  
    
    # Cutoff
    if(!is.null(cutoff.wilcox))
    {
      if(identical(cutoff.override,TRUE))
      {
        venn.wilcox = rownames(venn.wilcox[1:cutoff.hits,])
      }
      else
      {
        venn.wilcox = rownames(venn.wilcox[venn.wilcox$p.value <= cutoff.wilcox,])
      }
      
    }
    # return list
    list.return=c(list.return,list("Wilcox" = venn.wilcox))
  }
  
  # DESEQ2
  if(!is.null(deseq)) {
    
    venn.deseq=df.deseq[order(df.deseq[,"padj"], na.last=TRUE, decreasing=FALSE),]
    #df.output$deseq.pval=apply(df.output,1,function(i) return(df.deseq[i["genes"],"padj"]))
    # Cutoff
    if(!is.null(cutoff.deseq))
    {
      if(identical(cutoff.override,TRUE))
      {
        venn.deseq = rownames(venn.deseq[1:cutoff.hits,])
      }
      else
      {
        venn.deseq = rownames(venn.deseq[venn.deseq$padj <= cutoff.deseq,])
      }
    }
    # return list
    list.return=c(list.return,list("DEseq2" = venn.deseq))
  }
  
  
  # MAGeCK
  if(!is.null(mageck)) {
    
    venn.mageck=df.mageck[order(df.mageck[,"rank"], na.last=TRUE,decreasing=FALSE),]

    
    # Cutoff
    if(!is.null(cutoff.mageck))
    {
      if(identical(cutoff.override,TRUE))
      {
        venn.mageck = rownames(venn.mageck[1:cutoff.hits,])
      }
      else
      {
        venn.mageck = rownames(venn.mageck[venn.mageck$fdr <= cutoff.mageck,])
      }
    }
    # return list
    list.return=c(list.return,list("MAGeCK" = venn.mageck))
  }
  
  ## Now we have each gene name in the list, so we need to count
 
  #return data
  return(list.return)
  
} else if(output == "3dplot")
  
{
  # IDEA:
  
  # 3d scatter plot of 3 analysis methods
  # Scatter 1: p-value wilcox/riger/deseq vs. FDR/sgRNA Mageck vs. AUC/sgRNA dss
  # Scatter 2: log2FC wilcox/riger/deseqvs. FDR/sgRNA Mageck vs. AUC/sgRNA dss
  df.output.list = generate.hits(wilcox=wilcox, deseq=deseq, mageck=mageck, type=type, cutoff.deseq = cutoff.deseq, cutoff.wilcox = cutoff.wilcox, cutoff.mageck = cutoff.mageck, cutoff.override=cutoff.override, cutoff.hits=cutoff.hits)
  #print(cutoff.override)
  #print(cutoff.hits)
  #print(df.output.list)
  if(length(df.output.list) > 0)
  {
    
    df.output = data.frame(
      genes = df.output.list,
      stringsAsFactors=FALSE)
 
    rownames(df.output) = df.output$genes
    
    # what to plot
    # plot.method = x,y,z method
    # plot.feature= x,y,z output of method
    #  -> MUST MATCH
    
    
    
    # wilcox
    if(!is.null(wilcox)) {
      
      df.output$wilcox.log2fc=apply(df.output,1,function(i) return(log2(df.wilcox[i["genes"],"foldchange"])))
      df.output$wilcox.pval=apply(df.output,1,function(i) return(df.wilcox[i["genes"],"p.value"]))
      
    }
   
    # DESEQ2
    if(!is.null(deseq)) {
      
      df.output$deseq.log2fc=apply(df.output,1,function(i) return(df.deseq[i["genes"],"log2FoldChange"]))
      df.output$deseq.pval=apply(df.output,1,function(i) return(df.deseq[i["genes"],"padj"]))
      
    }
    
    # MAGeCK
    if(!is.null(mageck)) {
      
      df.output$mageck.fdr=apply(df.output,1,function(i) return(df.mageck[i["genes"],"fdr"]))
      df.output$mageck.sgrna=apply(df.output,1,function(i) return(df.mageck[i["genes"],"sgrna"]))
      
    }
    
    # add color table
    if (type=="enriched")
    { df.output$color=rgb(217,35,35, 255, maxColorValue=255) }
    else
    {
      df.output$color=rgb(46,98,166, 255, maxColorValue=255)
    }
    
    # Cutoff
    # the cutoff will be INDIVIDUALLY, which means this number indicates TOP HITS of each being plotted
    # sorting is done according to the feature selected for EACH set of data.

    
    # if there is any NA since no prediction was made, set it to zero for a respective thing
    
    
    # now prepare plotting
    # plot.method + plot.feature
    # what to plot
    # X-axis
    plot.xaxis = paste(plot.method[1],plot.feature[1], sep=".")
    
    if(plot.feature[1]=="sgrna" || plot.feature[1]=="log2fc")
    {
      df.output.x = df.output[order(df.output[,plot.xaxis], decreasing = TRUE),]
    }
    else
    {  
      df.output.x = df.output[order(df.output[,plot.xaxis], decreasing = FALSE),]
    }
    # Y-axis
    plot.yaxis = paste(plot.method[2],plot.feature[2], sep=".")
    
    
    if(plot.feature[2]=="sgrna" || plot.feature[2]=="AUC" || plot.feature[2]=="log2fc")
    {df.output.y = df.output[order(df.output[,plot.yaxis], decreasing=TRUE),]
    }
    else
    {
      df.output.y = df.output[order(df.output[,plot.yaxis], decreasing=FALSE),]
    }
    
    # Z-Axis
    plot.zaxis = paste(plot.method[3],plot.feature[3], sep=".")
    
    
    if(plot.feature[3]=="sgrna" || plot.feature[3]=="AUC" || plot.feature[3]=="log2fc")
    {
      df.output.z = df.output[order(df.output[,plot.zaxis], decreasing=TRUE),]
    }
    else
    {
      df.output.z = df.output[order(df.output[,plot.zaxis], decreasing=TRUE),]
    }
    
    # make -log10 for p-values
    # get to know using plot.feature 
    
    xlab=plot.xaxis
    ylab=plot.yaxis
    zlab=plot.zaxis
    
    
    if(plot.feature[1]=="pval" || plot.feature[1]=="fdr")
    {
      df.output.x[,plot.xaxis] = as.numeric(-log10(df.output.x[,plot.xaxis]))
      df.output.y[,plot.xaxis] = as.numeric(-log10(df.output.y[,plot.xaxis]))
      df.output.z[,plot.xaxis] = as.numeric(-log10(df.output.z[,plot.xaxis]))
      # remove NAs
      df.output.x[,plot.xaxis] = as.numeric(apply(df.output.x,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(i[plot.xaxis])} ) )
      df.output.y[,plot.xaxis] = as.numeric(apply(df.output.y,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(i[plot.xaxis])} ) )
      df.output.z[,plot.xaxis] = as.numeric(apply(df.output.z,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(i[plot.xaxis])} ) )
      
      xlab = paste("-log10 of", plot.xaxis, sep=" ")
      
    }
    if(plot.feature[2]=="pval" || plot.feature[2]=="fdr")
    {
      df.output.x[,plot.yaxis] = as.numeric(-log10(df.output.x[,plot.yaxis]))
      df.output.y[,plot.yaxis] = as.numeric(-log10(df.output.y[,plot.yaxis]))
      df.output.z[,plot.yaxis] = as.numeric(-log10(df.output.z[,plot.yaxis]))
      # remove NAs
      df.output.x[,plot.yaxis] = as.numeric(apply(df.output.x,1, function(i) if(is.na(i[plot.yaxis]) || is.infinite(i[plot.yaxis])) {return(0)} else {return(i[plot.yaxis])} ) )
      df.output.y[,plot.yaxis] = as.numeric(apply(df.output.y,1, function(i) if(is.na(i[plot.yaxis]) || is.infinite(i[plot.yaxis])) {return(0)} else {return(i[plot.yaxis])} ) )
      df.output.z[,plot.yaxis] = as.numeric(apply(df.output.z,1, function(i) if(is.na(i[plot.yaxis]) || is.infinite(i[plot.yaxis])) {return(0)} else {return(i[plot.yaxis])} ) )
      
      ylab = paste("-log10 of", plot.yaxis, sep=" ")
    }
    if(plot.feature[3]=="pval" || plot.feature[3]=="fdr")
    {
      df.output.x[,plot.zaxis] = as.numeric(-log10(df.output.x[,plot.zaxis]))
      df.output.y[,plot.zaxis] = as.numeric(-log10(df.output.y[,plot.zaxis]))
      df.output.z[,plot.zaxis] = as.numeric(-log10(df.output.z[,plot.zaxis]))
      # remove NAs
      df.output.x[,plot.xaxis] = as.numeric(apply(df.output.x,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(i[plot.xaxis])} ) )
      df.output.y[,plot.xaxis] = as.numeric(apply(df.output.y,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(i[plot.xaxis])} ) )
      df.output.z[,plot.xaxis] = as.numeric(apply(df.output.z,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(i[plot.xaxis])} )) 
      
      zlab = paste("-log10 of", plot.zaxis, sep=" ")
    }
    
    #print(nrow(df.output.x))
    #print(nrow(df.output.y))
    #print(nrow(df.output.z))
    
    # remove NAs -> set to ZERO 
   # Go through each column and replce NA, or INF with 0 if numeric
    
    #print(apply(sapply(df.output.x,is.numeric), 2, function(x) sapply(x, is.numeric)) #sapply(x, function(i)
      #print(class(i))
      
      #if(is.numeric(i))  
      #{
      #  if(is.na(x) || is.infinite(x))
      #  {print("NO!")}
      #  else
       # { print("YES!")}
      #}
      
      #)
    #  )
   
   #str(df.output.x)
   
    # remove NAs
    df.output.x[,plot.xaxis] = apply(df.output.x,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(as.numeric(i[plot.xaxis]))} ) 
    df.output.y[,plot.xaxis] = apply(df.output.y,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(as.numeric(i[plot.xaxis]))} ) 
    df.output.z[,plot.xaxis] = apply(df.output.z,1, function(i) if(is.na(i[plot.xaxis]) || is.infinite(i[plot.xaxis])) {return(0)} else {return(as.numeric(i[plot.xaxis]))} ) 
    df.output.x[,plot.yaxis] = apply(df.output.x,1, function(i) if(is.na(i[plot.yaxis]) || is.infinite(i[plot.yaxis])) {return(0)} else {return(as.numeric(i[plot.yaxis]))} ) 
    df.output.y[,plot.yaxis] = apply(df.output.y,1, function(i) if(is.na(i[plot.yaxis]) || is.infinite(i[plot.yaxis])) {return(0)} else {return(as.numeric(i[plot.yaxis]))} ) 
    df.output.z[,plot.yaxis] = apply(df.output.z,1, function(i) if(is.na(i[plot.yaxis]) || is.infinite(i[plot.yaxis])) {return(0)} else {return(as.numeric(i[plot.yaxis]))} ) 
    df.output.x[,plot.zaxis] = apply(df.output.x,1, function(i) if(is.na(i[plot.zaxis]) || is.infinite(i[plot.zaxis])) {return(0)} else {return(as.numeric(i[plot.zaxis]))} ) 
    df.output.y[,plot.zaxis] = apply(df.output.y,1, function(i) if(is.na(i[plot.zaxis]) || is.infinite(i[plot.zaxis])) {return(0)} else {return(as.numeric(i[plot.zaxis]))} ) 
    df.output.z[,plot.zaxis] = apply(df.output.z,1, function(i) if(is.na(i[plot.zaxis]) || is.infinite(i[plot.zaxis])) {return(0)} else {return(as.numeric(i[plot.zaxis]))} ) 
    
    # start plotting
    #requireNamespace("scatterplot3d")
    
    #check if enough data is available
    
    
    
    # get min/max of ALL
    min.x = min(min(df.output.x[,plot.xaxis]),min(df.output.y[,plot.xaxis]),min(df.output.z[,plot.xaxis]))
    max.x = max(max(df.output.x[,plot.xaxis]),max(df.output.y[,plot.xaxis]),max(df.output.z[,plot.xaxis]))
    
    xlim = c(min.x, max.x)
    
    min.y = min(min(df.output.x[,plot.yaxis]),min(df.output.y[,plot.yaxis]),min(df.output.z[,plot.yaxis]))
    max.y = max(max(df.output.x[,plot.yaxis]),max(df.output.y[,plot.yaxis]),max(df.output.z[,plot.yaxis]))
    
    #print(c(min(df.output.x[,plot.yaxis]),min(df.output.y[,plot.yaxis]),min(df.output.z[,plot.yaxis])))
    #print(c(max(df.output.x[,plot.yaxis]),max(df.output.y[,plot.yaxis]),max(df.output.z[,plot.yaxis])))
    ylim = c(min.y, max.y)
    
   # print(ylim)
    min.z = min(min(df.output.x[,plot.zaxis]),min(df.output.y[,plot.zaxis]),min(df.output.z[,plot.zaxis]))
    max.z = max(max(df.output.x[,plot.zaxis]),max(df.output.y[,plot.zaxis]),max(df.output.z[,plot.zaxis]))
    
    zlim = c(min.z, max.z)
    #print(zlim)
    
    s3d.x <- scatterplot3d::scatterplot3d(df.output.x[,plot.xaxis],   # x axis
                           df.output.x[,plot.yaxis],     # y axis
                           df.output.x[,plot.zaxis],    # z axis
                           main=paste("Hit Comparison of\n", plot.xaxis, plot.yaxis, plot.zaxis, sep=" "),
                           xlab=xlab,
                           ylab=ylab,
                           zlab=zlab,
                           type="h",
                           pch=pch,
                           xlim = xlim,
                           ylim = ylim,
                           zlim = zlim,
                           grid=TRUE,
                           color=df.output.x$color,
    )
    s3d.coords <- s3d.x$xyz.convert(df.output.x[,plot.xaxis], df.output.x[,plot.yaxis], df.output.x[,plot.zaxis]) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=df.output.x$genes,               # text to plot
         cex=.7, pos=4)
    
    par(new=TRUE)
    
    s3d.y <- scatterplot3d::scatterplot3d(df.output.y[,plot.xaxis],   # x axis
                           df.output.y[,plot.yaxis],     # y axis
                           df.output.y[,plot.zaxis],    # z axis
                           main=paste("Hit Comparison of\n", plot.xaxis, plot.yaxis, plot.zaxis, sep=" "),
                           xlab=xlab,
                           ylab=ylab,
                           zlab=zlab,
                           xlim = xlim,
                           ylim = ylim,
                           zlim = zlim,
                           type="h",
                           pch=pch,
                           grid=TRUE,
                           color=df.output.x$color,
    )
    s3d.coords <- s3d.y$xyz.convert(df.output.y[,plot.xaxis], df.output.y[,plot.yaxis], df.output.y[,plot.zaxis]) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=df.output.y$genes,               # text to plot
         cex=.7, pos=4)
    
    par(new=TRUE)
    s3d.z <- scatterplot3d::scatterplot3d(df.output.z[,plot.xaxis],   # x axis
                           df.output.z[,plot.yaxis],     # y axis
                           df.output.z[,plot.zaxis],    # z axis
                           main=paste("Hit Comparison of\n", plot.xaxis, plot.yaxis, plot.zaxis, sep=" "),
                           xlab=xlab,
                           ylab=ylab,
                           zlab=zlab,
                           xlim = xlim,
                           ylim = ylim,
                           zlim = zlim,
                           type="h",
                           pch=pch,
                           grid=TRUE,
                           
                           color=df.output.x$color,
    )
    s3d.coords <- s3d.z$xyz.convert(df.output.z[,plot.xaxis], df.output.z[,plot.yaxis], df.output.z[,plot.zaxis]) # convert 3D coords to 2D projection
    text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=df.output.z$genes,               # text to plot
         cex=.7, pos=4)
    
    
    #return(df.output)
    #return(df.output)
    df.output.x$color = NULL
    df.output.y$color = NULL
    df.output.z$color = NULL
    df.output.list = list(df.output.x,df.output.y, df.output.z)
    return(df.output.list)
    
  }
  
  
####################

} 


}