#'@rdname gmc
#'@method gmc CRSM
#'@export

gmc.CRSM <-
function(object, splitcrit="score", ...){

  if(splitcrit[1] == "score"){
    sc  <- rowSums(object$data)
    scm <- ifelse(sc > median(sc), 1,0)
  }
  else{
    if(!is.vector(splitcrit)){stop("Error: split criterium has to be a vector!", call. = FALSE)}
    if(length(table(splitcrit)) != 2) {"Error: split criterium must have exactly two values for plotting"}
    scm <- splitcrit
  }
  sp_dat <- split(as.data.frame(object$data_p, stringsAsFactors=F), as.factor(scm), drop=FALSE)
  sp_res <- lapply(sp_dat, function(dat) CRSM(dat, low=0, high=1))

  itempar_split <- sapply(sp_res, function(re) list(re$itempar))
  itemse_split <- sapply(sp_res, function(re) list(re$se.item.mean))
  lambda_split <- sapply(sp_res, function(re) re$disppar)

  plot(itempar_split[[1]], itempar_split[[2]], type="p", pch=20, xlab="item par. group 1", ylab="item par. group 2", xlim=c(-ceiling(max(abs(unlist(itempar_split)))),ceiling(max(abs(unlist(itempar_split))))), ylim=c(-ceiling(max(abs(unlist(itempar_split)))),ceiling(max(abs(unlist(itempar_split))))))
  abline(0,1)
  text(itempar_split[[1]], itempar_split[[2]], labels=colnames(sp_dat[[1]]), pos=1)

  cat(" \t Disppar Group 1: ", lambda_split[1], "\n",
      "\t Disppar Group 2: ", lambda_split[2])
}
