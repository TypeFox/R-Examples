summary.cooccur <-
function(object, ...){
  ptab <- object$results
  cat("Call:\n")
  print(object$call)
  if("omitted" %in% names(object)){
    cat(paste("\nOf ",object$pot_pairs," species pair combinations, ",object$omitted," pairs (",round(object$omitted / object$pot_pairs * 100,2)," %) were removed from the analysis because expected co-occurrence was < 1 and ",sep=""))
    cat(paste(object$pairs," pairs were analyzed","\n",sep=""))
  }else{
  cat(paste("\n",object$pairs," pairs were analyzed","\n",sep=""))
  }
  cat("\nCooccurrence Summary:\n")
  cooccur_list <- c(round(object$species,0),
                    round(mean(object$sites,na.rm = T),0), # MODIFIED
                    round(object$positive,0),
                    round(object$negative,0),
                    round(object$random,0),
                    round(object$unclassifiable,0),
                    round(object$percent_sig,1))
  names(cooccur_list) <- c("Species",
                           "Sites",
                           "Positive",
                           "Negative",
                           "Random",
                           "Unclassifiable",
                           "Non-random (%)")
  
  #print(cooccur_list)
  class(cooccur_list) <- "summary.cooccur"
  return(cooccur_list)
  #NextMethod("summary")
}
