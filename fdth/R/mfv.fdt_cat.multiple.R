mfv.fdt_cat.multiple <- function(x, ...)
{
 
 if(class(x)[1]=='fdt_cat.multiple' | class(x)[1] == 'fdt_cat'){
 res <- list()

 for(i in 1:length(x)){

  y <- x[[i]][[1]]
  class(y) <- c('fdt_cat.multiple',
                'fdt_cat',
                'data.frame')
  res[[i]] <- mfv.fdt_cat(y)

 }


}
 else {

  res <- lapply(x,
                mfv.fdt_cat)

 }

return(res)
}
