#' makeMRVS
#'
#' This function will calculate matching random variant sets (MRVS) idential to AVS
#' @param avs A GRanges object which is outputted by makeAVS function
#' @param bgSize An integer for the number of MRVS. Default: 100
#' @param mc.cores Number of cores to use. Default: 8
#' @keywords VSE
#' @examples
#' \dontrun{
#' ld <- loadLd("ld.csv", type="raggr")
#' avs <- makeAVS(ld)
#' makeMRVS(avs, bgSize=100, mc.cores=8)
#' }
#' #As an example, we have added MRVS (size=200) for Breast Cancer AVS.
#' load(file.path(system.file("extdata", "bca.mrvs.200.Rda", package="VSE")))
#' @import GenomicRanges
#' @importFrom parallel mclapply
#' @importFrom IRanges IRanges
#' @export
makeMRVS <- function(avs, bgSize=100, mc.cores=6){
  if (!exists("nullblocks.08")){
      nullblocks.08 <- list()
    tmpfile <- tempfile(fileext = ".rda")
    url <- "http://www.hansenhelab.org/null0.8.rda"
    download.file(url, destfile = tmpfile, method = "curl")
    load(tmpfile)
  }
  if (!exists("nullblocks.08")){
    stop("Downloading null failed. Are you connected to internet?")
  }
  no_of_tags <- length(avs)
  tally <- nullblocks.08$tally
  null_id_list <- matrix(NA, nrow = no_of_tags, ncol = bgSize+1)
  for (i in 1:no_of_tags){
    null_id_list[i,1] <- as.character(elementMetadata(avs[[i]])[1,2])
     ld_tally <- length(avs[[i]])
     tlist <- tally[tally$X0.8 == ld_tally, 1]
     null_id_list[i, c(2:ncol(null_id_list))] <- as.character(tlist[sample(length(tlist), bgSize, replace = ifelse(length(tlist)<bgSize, TRUE, FALSE))])
  }
  message(paste0("Using ", mc.cores, " cores"))
  mrvs <- list()
  for (i in 2:ncol(null_id_list)){
    message(paste0("Computing MRVS no. ", i-1))
    null_glist <- GRangesList(
                    parallel::mclapply(null_id_list[,i],
                             function(x){
                                chr <- strsplit(x, ":")[[1]][1]
                                pos <- strsplit(x, ":")[[1]][2]
                                var <- paste0("nullblocks.08", "$", chr)
                                chr <- as.factor(unlist(lapply(chr, function(z) gsub("chr","",z))))
                                lds <- eval(parse(text=var))[eval(parse(text=var))$BP_A %in% pos, 1]
                                gr <- GRanges(seqnames = chr, ranges = IRanges::IRanges(start=lds, width=1))
                            },
                            mc.cores = mc.cores
                    )
                  )
    mrvs <- c(mrvs, null_glist)
  }
  return(mrvs)
}
