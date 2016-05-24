BilogData<-function(resp, data.name="mydata", location=NULL, ret.val=FALSE, idvar=NULL){
  if (is.data.frame(resp)==FALSE)
 	stop("Data must be a data frame")	
  
  if (is.null(location)==TRUE) 
  	{location=getwd()}
  
  if (is.null(idvar)==TRUE)  
  	{nit = ncol(resp)}
  else 
  	{nit = ncol(resp) -1}
 
  if (nit > 9999) 
        stop("cannot have more than 9999 items")
    d = paste(data.name, ".dat", sep = "")
    d.save = paste(location, "/", d, sep="")
    np = nrow(resp)
    if (np > 999999) 
        stop("cannot have more than 999999 observations")
  if (is.null(idvar)==TRUE) 
  		{ids = sprintf("%06d", 1:np)}
  else 
  		{ids = sprintf("%06d", as.integer(unlist(resp[idvar])))}
  		
   if (is.null(idvar)==TRUE) 
  		{resp = resp}
    else 
  		{vars = names(resp) %in% idvar; resp = resp[!vars]}

    resp = cbind(ids, resp)
    write.table(resp, file = d.save, append = FALSE, sep = "", row.names = FALSE, col.names = FALSE, na = ".", quote = FALSE)
        
  if (ret.val) 
    {write.table(resp, file = "", sep = "", row.names = FALSE, col.names = FALSE, na = ".", quote = FALSE)}    
}
