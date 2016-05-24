multistationsreport <- function(...,
                                indices = c("meanflow", "Q95", "MAM1", "MAM7", "MAM10", "MAM30", "MAM90", "baseflowindex", "recession"),
                                recessionmethod = "MRC",
                                recessionseglength = 7,
                                recessionthreshold = 70,
                                recessiontrimIRS = 0.1,
                                lflist = NULL
){

  ind <- match.arg(indices, several.ok = TRUE)
  dotargs <- list(...)
  if(!is.null(lflist)){
    args <- c(dotargs,lflist)} else {
      args <- dotargs}
  if(is.null(lflist)){
    names <- as.character(as.list(substitute(list(...)))[-1L])
  }else{
    names <- c(as.character(as.list(substitute(list(...)))[-1L]),names(lflist))
  }
  table <- data.frame(row.names = names)
  for(ii in args){
    lfcheck(ii)
  }
  if("meanflow" %in% ind){
    table$meanflow <- sapply(args,meanflow)
  }
  if("Q95" %in% ind){
    table$Q95 <- sapply(args,Q95)
  }
  if("MAM1" %in% ind){
    table$MAM1 <- sapply(args,MAM, n=1)
  }
  if("MAM7" %in% ind){
    table$MAM7 <- sapply(args,MAM, n=7)
  }
  if("MAM10" %in% ind){
    table$MAM10 <- sapply(args,MAM, n=10)
  }
  if("MAM30" %in% ind){
    table$MAM30 <- sapply(args,MAM, n=30)
  }
  if("MAM90" %in% ind){
    table$MAM90 <- sapply(args,MAM, n=90)
  }
  if("baseflowindex" %in% ind){
    table$BFI <- sapply(args,BFI)
  }
  if("recession" %in% ind){
    table$recession <- sapply(args,recession,method = recessionmethod, seglength = recessionseglength, threshold = recessionthreshold,plotMRC = FALSE, trimIRS=recessiontrimIRS)
  }
  table
}
