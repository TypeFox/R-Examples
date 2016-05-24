#' loadLd
#'
#' This function will load a raggr output csv or user defined bed file to GRanges object. It is recommended to generate a LD file using http://raggr.usc.edu. If you prefer to use other softwares to calculate LD SNPs (e.g., plink), please format the output to bed files. Please make sure that you use a cutoff of 0.8 for r2 value.
#' @param data An input file. Must be a raggr generated csv file or a bed file. If a bed file, it must contain at least five columns: chr, start, end, LD_snp_id, tag_snp_id.
#' @param type The type of input file. Must be "bed" or "raggr".
#' @keywords LD,bed,raggr,GRanges
#' @examples
#' ld<-loadLd(file.path(system.file("extdata", "ld_BCa_raggr.csv", package="VSE")), type="raggr")
#' @import GenomicRanges
#' @export
loadLd <- function(data, type){
  if (type == "bed"){
    df <- read.table(data, sep="\t")
    if (is.numeric(df[1,2])){
      cnames <- c("chr","start","end", "ldId","tagId")
    } else {
      cnames <- colnames(df)
    }
    colnames(df) <- cnames
    if (!("chr" %in% cnames)){
      stop("No chr column found in df");
    }
    if (!("start" %in% cnames)){
      stop("No start column found in df");
    }
    if (!("end" %in% cnames | "stop" %in% cnames)){
      stop("No end column found in df");
    }
    df <- df[,-c(6,ncol(df))]
  } else if (type == "raggr") {
    ld.df <- read.csv(data, header=T)
    df <- data.frame(chr=ld.df$SNP1.Chr,
                    start=ld.df$SNP2.Pos,
                    end=ld.df$SNP2.Pos+1,
                    idLd=sapply(as.character(ld.df$SNP2.Name),
                                function(x) strsplit(x, ":")[[1]][1]),
                    idTag=sapply(as.character(ld.df$SNP1.Name),
                                 function(x) strsplit(x, ":")[[1]][1]))
  } else {
    df <- read.table(data, sep="\t")
    if (ncol(df) != 7){
      stop("Does not have 7 columns");
    } else {
      if (!is.numeric(df[1,2])){
        df <- df[-1,]
      }
      cnames <- c("chrTag","startTag","idTag","chr","start","idLd","r")
      colnames(df) <- cnames;
    }
    df$end <- df$start + 1;
    df <- df[df$r >= 0.8,c(4,5,8,6,3)]
  }
  makeGRangesFromDataFrame(df, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)
}
