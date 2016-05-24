addFieldsFromHeader  <- function(structure)
{
# addFieldsFromHeader 
#
#   structure = addFieldsFromHeader(structure)
#       take a header list on a structure and create a structure of header
#       fields. Each field is composed by a cell array per individual. 
#       The refered individuals are in id field  
#
  header            <- structure$header
  value             <- structure$value
  structure$header  <-  NULL
  structure$value   <-  NULL
  
  ## only if name is a field  
  couples <- vector(length=0)
  if (isfield(structure,'name')){
    names <- structure$name;
    
    for (i in 1: length(names)){
      for (j in 1: length(header)){
        if (names[[i]]==header[[j]])
          couples <- rbind(couples ,c(i, j))
        }
    }
  }
  
  ##create fields directly from header and value
  idx_id <- which(header=='id') 
  structure$id <- value[,idx_id];
  structure$value <- matrix(nrow=length(structure$id), ncol=dim(couples)[1])
  for (j in 1 : dim(couples)[1]){
    structure$value[,j] <- value[, couples[j,2]]
  }
  
  ## Only if "id" is a field
  if (isfield(structure,'id')){
      ids = unique(structure$id)
      structure <- splitField(structure,'time',ids)
      structure <- splitField(structure,'amount',ids)
      structure <- splitField(structure,'rate',ids)
      structure <- splitField(structure,'tinf',ids)
      structure <- splitField(structure,'type',ids)
  
      structure$id <- ids
      structure
  #
  }else{
      stop('id must be a field')
  }

}

splitField  <- function(structure,field,ids)
{
  # splitField 
  #
  #   structure = splitField(structure,field,ids)
  #   splitField split a field in a structure according to ids  
  #   field. It create a cell per individual. 
  #   
  
  if (isfield(structure,field)){
    temp  <-  list();
    for (i in 1:length(ids)){
      temp[[i]]  <-  structure[[field]][structure$id==ids[i]]
    }
    structure[[field]]  <-  temp
  }
  structure
}

