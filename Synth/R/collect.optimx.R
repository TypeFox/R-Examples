collect.optimx <-
function(res,opt="min"){ 

# reorder dataframe for best solution
if(opt=="min"){
  res <- res[order(res$value,decreasing=FALSE),]
}
if(opt=="max"){
  res <- res[order(res$value,decreasing=TRUE),]
}
 
# index for value
val <- which(colnames(res)=="value") 
# index for pars
pars <- 1:(val-1)

# extract best solution
out <- list(out.list=res,
            par=res[1,pars],
            value=res[1,val]
            )
return(invisible(out))

}

