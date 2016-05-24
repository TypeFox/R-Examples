stats.data = function(dataset, namecolumn = 1, fullmatchcolumn = 2, extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = NA, controls.target = NULL, controls.nontarget = "random", type="stats")
{

  # type can be:
  #   mapping | for mapping information
  #   stats | for general statistics of the dataset
  #   dataset | for indepth dataset statistics for each gene
  #   controls | for indepth statistics per control (target and non-targeting)
  
## Load data, extract gene names
gene.names = sub(extractpattern,"\\1",dataset[,namecolumn],perl=TRUE)

if(type=="mapping")
{
    # first get total read count
    readcount.mapped.total = round(sum(dataset[,fullmatchcolumn]))
    
    # get unmapped readcount from input
    readcount.unmapped.total = round(as.numeric(readcount.unmapped.total))
    
    # percentage
    readcount.total.percentage = (readcount.mapped.total / readcount.unmapped.total)*100
    
    
    # Make output
    # Print individual stats
    
    data.return= data.frame(
      type = c("Mapped Read Counts", "Total Read Counts", "Percentage"),
      readcount = c(readcount.mapped.total,readcount.unmapped.total, paste(round(readcount.total.percentage, digits=3), "%", sep=" ")),
      stringsAsFactors=FALSE)
    return(data.return)
}
if(type=="stats")
{
    # Average Dataset
    readcount.mean = round(mean(dataset[,fullmatchcolumn]))
    readcount.median = round(median(dataset[,fullmatchcolumn]))
    readcount.min = round(min(dataset[,fullmatchcolumn]))
    readcount.max = round(max(dataset[,fullmatchcolumn]))
    readcount.sd = round(sd(dataset[,fullmatchcolumn], na.rm=TRUE))
    unmapped.sgRNA = length(dataset[dataset[,fullmatchcolumn] == 0,namecolumn])
    
    data.return = data.frame(
      readcount = c("Mean", "Median", "SD", "Min", "Max", "# sgRNA not present"),
      sgRNA = c(readcount.mean,readcount.median,readcount.sd,readcount.min,readcount.max, unmapped.sgRNA),
      stringsAsFactors = FALSE)
    
    # Averages per gene
    dataset.g = aggregatetogenes(dataset, namecolumn=namecolumn, countcolumn=fullmatchcolumn, extractpattern=extractpattern, agg.function=sum)
    
    readcount.gene.mean = round(mean(dataset.g[,fullmatchcolumn]))
    readcount.gene.median = round(median(dataset.g[,fullmatchcolumn]))
    readcount.gene.min = round(min(dataset.g[,fullmatchcolumn]))
    readcount.gene.max = round(max(dataset.g[,fullmatchcolumn]))
    readcount.gene.sd = round(sd(dataset.g[,fullmatchcolumn], na.rm=TRUE))
    unmapped.gene = length(dataset[dataset.g[,fullmatchcolumn] == 0,namecolumn])
    
    data.return$gene = c(readcount.gene.mean,readcount.gene.median,readcount.gene.sd,readcount.gene.min,readcount.gene.max, unmapped.gene)
    
    return(data.return)
}
    
if(type=="dataset" || type=="controls")
{
    dataset$genes = gene.names
    
    sgRNA.counts = tapply(dataset[,fullmatchcolumn],dataset$genes,function(x) length(x[x!=0]))
    
  if(type == "dataset")
  { 
    
    df.gene = tapply(dataset[,fullmatchcolumn], dataset$genes, function(z)
      {
      return(c(round(mean(z)),round(median(z)),round(sd(z)),round(min(z)),round(max(z))))
    })
    
    # make data.frame
    df.attributes = attributes(df.gene)
    
    df.gene = do.call(rbind.data.frame, df.gene)
    df.gene$Name = unlist(df.attributes[[2]])
    colnames(df.gene) = c("readcount.mean","readcount.median","readcount.sd","readcount.min","readcount.max","Name")
    

      # make output for return of the whole dataset
      return(df.gene)
    }
    else if (type == "controls")
    {
      df.gene = as.data.frame(df.gene[as.character(df.gene[,namecolumn]) %in% c(controls.target,controls.nontarget),"Name"])

     colnames(df.gene) = "name"
     
     df.gene$readcount.mean = apply(df.gene, 1, function(u) {
       value = dataset[dataset$genes == u[namecolumn] , fullmatchcolumn]
       return(round(mean(value)))
     })
     df.gene$readcount.median = apply(df.gene, 1, function(u) {
       value = dataset[dataset$genes == u[namecolumn] , fullmatchcolumn]
       return(round(median(value)))
     })
     df.gene$readcount.sd = apply(df.gene, 1, function(u) {
       value = dataset[dataset$genes == u[namecolumn] , fullmatchcolumn]
       return(round(sd(value)))
     })
     df.gene$readcount.min = apply(df.gene, 1, function(u) {
       value = dataset[dataset$genes == u[namecolumn] , fullmatchcolumn]
       return(round(min(value)))
     })
     df.gene$readcount.max = apply(df.gene, 1, function(u) {
       value = dataset[dataset$genes == u[namecolumn] , fullmatchcolumn]
       return(round(max(value)))
     })

      data.return = df.gene[rownames(df.gene) %in% controls.target,]
      data.return = rbind.data.frame(data.return, df.gene[rownames(df.gene) %in% controls.nontarget,])
      
      # make output for return
      return(df.gene)
    }
   
}

}