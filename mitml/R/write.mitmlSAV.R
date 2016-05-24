write.mitmlSAV <- function(x, filename){
# write to native SPSS format

  if(!"mitml"%in%class(x) & !"mitml.list"%in%class(x)) stop("'x' must be of class 'mitml' or 'mitml.list'.")
  if(!grepl(".sav$",filename)) filename <- paste(filename,".sav",sep="")

  if("mitml"%in%class(x)){
    il <- mitmlComplete(x,0:x$iter$m)
  }else{ 
    il <- x
  }

  for(ii in 1:length(il)){
    il[[ii]] <- cbind(ii-1,il[[ii]])
    colnames(il[[ii]])[1] <- "Imputation_"
  }
  out <- do.call(rbind, il)
  haven::write_sav(out, filename)

  invisible() 

}

