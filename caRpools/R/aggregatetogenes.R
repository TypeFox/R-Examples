aggregatetogenes=function(data.frame, namecolumn=1, countcolumn=2, agg.function=sum,extractpattern=expression("^(.+?)_.+"), type="aggregate"){
  
  # extract gene names
  if(is.null(data.frame$Gene))
  {data.frame$Gene  = as.character(sub(extractpattern,"\\1",data.frame[,namecolumn],perl=TRUE))}
 
  #str(data.frame)
 
  # make factor to character
  data.frame[, namecolumn] = sapply(data.frame[,namecolumn], FUN=function(x){as.character(x)})
  # to make aggregate working, we overwrite designIDs by its genes
  
  #str(data.frame)
  
  if(type=="aggregate")
  {
    # use aggregate function
    countdata = aggregate(data.frame, list(data.frame$Gene),
                          FUN=function(x){
                            if(is.numeric(x)){
                              agg.function(x)
                            }else{
                              x[1]
                            }
                          })
    #str(countdata)
    row.names(countdata) = countdata$Group.1
    countdata$Group.1 = NULL
    countdata[,namecolumn] = row.names(countdata)
    colnames(countdata) = colnames(data.frame)
    #str(countdata)
    
  }
  else if(type=="annotate")
  {
    countdata = data.frame
    colnames(countdata) = colnames(data.frame)
    
  }
  else
  {stop("No type set!")}
  
  return(countdata)
}