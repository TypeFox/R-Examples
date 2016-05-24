fextractassoc <-
function(dnndata)
{

  if(class(dnndata)  != "nb") 
    stop(paste(sQuote("dnndata"), "must be a", 
               sQuote("nb"), "object", sep=" "))
  
  animID <- attr(dnndata, "region.id")
  
  #Extract into list object
  new.assoc<- sapply(1:length(animID),function(x) paste(c(x,dnndata[[x]])))
  #Ensure list object contains interactions
  if(is.list(new.assoc)==FALSE) new.assoc <- lapply(1:ncol(new.assoc), function(x) c(new.assoc[,x]))
  #Convert factors into numeric objects
  new.assoc2 <- lapply(new.assoc,as.numeric)
  new.assoc3 <- lapply(lapply(new.assoc2, function(x){replace(x, x == 0, NA)}),
                       function (x) x[!is.na(x)])
  # if function contains more than one interaction, run merging function TWICE to include ALL higher order interactions  
  if(length(new.assoc3) >= length(animID))
    if(max(unlist(lapply(new.assoc3,length))) > 1)
      if(max(tapply(rep(1,length(unlist(new.assoc3))),as.factor(unlist(new.assoc3)),sum))>2)
        new.assoc3 <- fmergeassoc(fmergeassoc(new.assoc3))
  return(new.assoc3)
}
