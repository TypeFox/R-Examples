get.gene.info = function(data, namecolumn=1, extractpattern=expression("^(.+?)(_.+)"), host="www.ensembl.org", database="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", filters="ensembl_gene_id", attributes = c("hgnc_symbol"), return.val = "dataset", controls=FALSE)
{
#   if("biomaRt" %in% rownames(installed.packages()) == FALSE) {
#     source("http://bioconductor.org/biocLite.R")
#     biocLite("biomaRt",suppressUpdates=update , ask=FALSE)
#   }
#requireNamespace(biomaRt)
  
if(!is.null(data) && nrow(data)>=1)
{
    
  
  # we apply this to an sgRNA dataset, not single genes
  if(return.val=="dataset")
  {
    # check for attributes and filters to be only one, attributes must be also the filter
    if (length(filters)!=1) stop("None or more than one single filter selected")
    
    # if dataset is control
    if(identical(controls, TRUE))
    {
      # get gene names from controls
      
      gene.names = data[,namecolumn]
    }
    else
    {
      # get gene names from design file
      gene.names = sub(extractpattern,"\\1",data[,namecolumn],perl=TRUE)
    }
    
    
    #start biomaRt interface and check if biomaRt is available
    handling = biomaRt::useMart(database, host=host)
    if(!exists("handling"))
    {stop("biomaRt connection is not working. This can be a connectivity issue (e.g. proxy settings, internet connection) or the biomaRt service is currently not avaible. \n
          You can skip any data annotation by setting it to FALSE in the MIACCS file.")}
    handling = biomaRt::useDataset(dataset,mart=handling)
    gene.info = biomaRt::getBM(
      filters=filters,
      attributes= c(filters,attributes),
      values= gene.names,
      mart = handling)
    
    #str(gene.info)
    if(nrow(gene.info) >= 1)
    {
      data$replace = as.character(gene.names)
      
      #print(data$replace)
      #gene.info.ext <-gene.info
      #print(gene.info)
      
      
      # Gene info has - original identifier - new identifier
      # data has - namecolumn, .... - original gene name
      # we want to replace the namecolumn set with the new identifier, if there is no new one, we keep the old one
      
      # First duplicate check on gene info to remove duplicates in NEW gene identifier
      # output to file for additional reference
      # Remove duplicates
      gene.info <- gene.info[!duplicated(gene.info[,2]),]
      # Then apply on data with namecolumn and put in there, what gene info 2 is giving for those where gene.info 1 is equal to extractd namecolumn, replace namecolumn 
      data[,namecolumn] = apply(data,1, function(z){
        # compare name in gene.info[,1] with data namecolumn and set replace = gene.info[,2]
          #return(gene.info[gene.info[,1] == z["replace"] ,2]) 
        sgRNAident = sub(extractpattern,"\\2",as.character(z[namecolumn]),perl=TRUE)
        
        if(z["replace"] %in% gene.info[,1])
        {replacement = TRUE
        } else {
          replacement = FALSE
        }
        if(identical(controls,TRUE))
        {
          if(identical(replacement, TRUE))
          { return(gene.info[gene.info[,1] == z["replace"] ,2])
          } else {
            return(z["replace"])
          }
          #return(gene.info[gene.info[,1] == z["replace"] ,2])
        } else {
          
          if(identical(replacement, TRUE))
          { 
            return(as.factor(paste(as.character(gene.info[gene.info[,1] == z["replace"] ,2][1]),sgRNAident,sep="")))
          } else {
            return(as.factor(paste(as.character(z["replace"]),sgRNAident,sep="")))
          }
          #return(as.factor(paste(as.character(gene.info[gene.info[,1] == z["replace"] ,2]),sub(extractpattern,"\\2",as.character(z[namecolumn]),perl=TRUE),sep="")) )
        }
      })

      # output to file for additional reference
      data$replace=NULL
      
      # check for duplicates and use original identifier instead which is stored in gene names
      
      
      
      return(data)
      
    }
    else
    {
      return(as.data.frame(data))
    }
    
     
  }
  if(return.val=="info")
  {
    # We can enrich gene information based on what is provided for biomaRt
    # get gene names from design file
    # if dataset is control
    if(identical(controls, TRUE))
    {
      # get gene names from controls
      gene.names = data[,namecolumn]
    }
    else
    {
      # get gene names from design file
      gene.names = sub(extractpattern,"\\1",data[,namecolumn],perl=TRUE)
    }
    #start biomaRt interface
    handling = biomaRt::useMart(database, host=host)
    if(!exists("handling"))
    {stop("biomaRt connection is not working. This can be a connectivity issue (e.g. proxy settings, internet connection) or the biomaRt service is currently not avaible. \n
          You can skip any data annotation by setting it to FALSE in the MIACCS file.")}
    handling = biomaRt::useDataset(dataset,mart=handling)
    gene.info = biomaRt::getBM(
      filters=filters,
      attributes= c(filters,attributes),
      values= gene.names,
      mart= handling)
    
    cols = ncol(gene.info)
    gene = aggregate.data.frame(gene.names,by=list(gene.names), function(x) return(x[1]))
    gene$Group.1=NULL
    data.return=data.frame(
      gene = as.character(gene$x),
      stringsAsFactors=FALSE)
    
    
    data.return[,colnames(gene.info)] = NA
    for(m in 1:cols)
    {
      for(i in 1:nrow(gene.info))
      {
        data.return[data.return[,namecolumn]== gene.info[i,1],m] = gene.info[i,m]
      }
      #data.return = cbind.data.frame(data.return, gene.info[,1:cols])
    }
  
    colnames(data.return) = colnames(gene.info)
    data.return[,is.na(colnames(data.return))] = NULL
    return(data.return)
    
  }
  
}

  
}