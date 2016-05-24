gene.remove=function(data,namecolumn=1,toremove=NULL, extractpattern=expression("^(.+?)_.+")){
  # dataset can be either gene or design
  # get gene names
  
  # dataset must be on design level
  # check for poisitions in which a certain gene is
  todelete = which(sub(extractpattern,"\\1",data[,namecolumn],perl=TRUE) %in% toremove)
  # delete those rows which did contain designs belonging to a gene of interest
  data = data[-todelete,]
  
  return(data)
}