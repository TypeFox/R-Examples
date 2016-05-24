###NAMESPACE ADDITIONS###
#' @importFrom BiocInstaller  biocVersion
#' @importFrom stats family pnorm pt qnorm rchisq rnorm runif
#' @importFrom reader cat.path reader shift.rownames
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline lines points rect text plot
#' @importFrom methods as callNextMethod is new prototype representation setAs setClass setGeneric setMethod setValidity
#' @importClassesFrom GenomicRanges GNCList GRanges GenomicRanges GenomicRangesORmissing
#' @importFrom GenomicRanges GRanges GRangesList
#' @importMethodsFrom GenomicRanges "names<-" length names start end 
#' @importMethodsFrom GenomicRanges width strand mcols  "mcols<-" show findOverlaps subsetByOverlaps
#' @importMethodsFrom GenomicRanges length names  "names<-" "dimnames<-" "["  "[<-"  "[["  "[[<-"  "$"  "$<-" cbind rbind
#' @importMethodsFrom GenomeInfoDb "seqlevels"  "seqlevels<-" "genome<-" "genome"  seqinfo  "seqinfo<-" seqnames  "seqnames<-" 
#' @importClassesFrom IRanges RangedData 
#' @importFrom IRanges "%over%" IRanges  RangedData  showAsCell 
#' @importMethodsFrom IRanges "colnames<-" "rownames<-" "universe<-" showAsCell
#' @importMethodsFrom IRanges as.data.frame as.list as.matrix cbind rbind colnames elementLengths
#' @importMethodsFrom IRanges end findOverlaps subsetByOverlaps gsub intersect lapply
#' @importMethodsFrom IRanges mean nrow ncol order as.list subjectHits queryHits 
#' @importMethodsFrom IRanges ranges rownames runLength space  flank  reduce resize
#' @importMethodsFrom IRanges start universe unlist width  "start<-"  "width<-"  "end<-" ranges "ranges<-"
#' @importFrom "GenomicFeatures" makeTxDbFromUCSC  exonsBy  transcriptsBy
#' @importMethodsFrom "GenomicFeatures"  exonsBy  transcriptsBy  as.list
#' @importClassesFrom "rtracklayer"  ChainFile
#' @importMethodsFrom "rtracklayer"  liftOver  import.chain
#' @importMethodsFrom "genoset"  genome "genome<-"
#' @importFrom "genoset"  chr  chrIndices  chrInfo  chrNames  chrOrder   "chrNames<-"
#' @importFrom "genoset"  "isGenomeOrder"  "locData"  "toGenomeOrder"  "locData<-"
#' @importFrom "biomaRt"  useMart  useDataset  getBM
#' @importClassesFrom "biomaRt"  Mart
#' @importFrom parallel  mclapply
#' @import Rcpp BiocGenerics NCmisc S4Vectors
#' @importFrom utils capture.output download.file read.table write.table read.delim
#' @importFrom stats median 
#' @importFrom graphics par
###END NAMESPACE###

#DataFrame
#seqlevels
#seqlevels<-
#genome<-

# doNotimportFrom utils capture.output download.file read.table write.table head tail data  
# dontimportClassesFrom "genoset" RangedDataOrGenomicRanges
# importNoClassesFrom("GenomicRanges", GRanges)
# importNoClassesFrom("IRanges", Rle, RangedData)
# doNotimportFrom AnnotationDbi head tail ncol as.list colnames get exists sample 
# doNotimportFrom(BiocGenerics,strand, "strand<-", colnames, cbind, rbind, unlist, order, rownames, ncol, as.vector, paste, as.data.frame)
# doNotimportFrom GenomicRanges "seqlevels"  "seqlevels<-" Seqinfo seqlengths
# importNoClassesFrom "GenomicFeatures" TranscriptDb

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("humarray version 1.0.0\n")
}

.onLoad <- function(libname, pkgname) {
  # library.dynam("humarray", pkgname, libname)
  #~/github/iChip/ImmunoChip_ChipInfo_New.RData
  options(chip.info="") # if you can access the file, you won't need to change this path
  options(ucsc="hg19") # depends on which analysis, need to set something though
  #data("iChipRegionsB36", "egSymb", "ImmunoChipB37", "hg18ToHg19","hg38ToHg19","hg19ToHg18","hg19ToHg38",
  #     package=pkgname, envir=parent.env(environment()))
  options(save.annot.in.current=1)  # 1 = TRUE, stores annotation in current folder to speed up subsequent lookups
}



#require(GenomicRanges); require(IRanges); require(reader); require(genoset)


########################
## internal functions ##
########################

finitize <- function(X) {
  if(is.data.frame(X)) { X <- as.matrix(X) }
  return(X[is.finite(X)])
}

minna <- function(...) {
  if(length(list(...))==1) { 
    min(finitize(...),na.rm=TRUE)
  } else {
    min(...,na.rm=TRUE)
  }
}

maxna <- function(...) {
  if(length(list(...))==1) { 
    max(finitize(...),na.rm=TRUE)
  } else {
    max(...,na.rm=TRUE)
  }
}

meanna <- function(...) {
  if(length(list(...))==1) { 
    mean(finitize(...),na.rm=TRUE)
  } else {
    mean(...,na.rm=TRUE)
  }
}

medianna <- function(...) {
  if(length(list(...))==1) { 
    median(finitize(...),na.rm=TRUE)
  } else {
    median(...,na.rm=TRUE)
  }
}

sdna <- function(...) {
  if(length(list(...))==1) { 
    sd(finitize(...),na.rm=TRUE)
  } else {
    sd(...,na.rm=TRUE)
  }
}

sumna <- function(...) {
  if(length(list(...))==1) { 
    sum(finitize(...),na.rm=TRUE)
  } else {
    sum(...,na.rm=TRUE)
  }
}

sortna <- function(...) {
  sort(..., na.last=TRUE)
}


# internal
l10 <- function(x) { O <- log10(x); O[!is.finite(O)] <- NA; return(O) }
# internal
Cor <- function(...) { cor(...,use="pairwise.complete") }
# internal
pt2 <- function(q, df, log.p=FALSE) {  2*pt(-abs(q), df, log.p=log.p) }




#' Manage flexible input for the build parameter
#' 
#' The genome annotation version for internals in this package should always be
#' of the form 'hgXX', where XX can be 15,16,17,18,19,38. However most functions
#' allow flexible entry of this parameter as a build number, e.g, 36,37,38, or as
#' 'build36', 'b36', etc. This function sanitizes various forms of input to the 
#' correct format for internal operations. 
#' @param build the input to be sanitized. 
#' @param allow.multiple logical, whether to force a single value, or allow a vector
#' of build strings as input
#' @param show.valid logical, if TRUE, show a list of supported values.
#' @return build string in the correct 'hgXX' format.
#' @export
#' @examples
#' ucsc.sanitizer(36)
#' ucsc.sanitizer("build38")
#' ucsc.sanitizer("b37")
#' ucsc.sanitizer(show.valid=TRUE)
ucsc.sanitizer <- function(build,allow.multiple=FALSE,show.valid=FALSE) {
  build.alt <- c("hg15","hg20","hg17","hg18","hg19","hg38",17,18,19,20,35,36,37,38,
                 "build35","build36","build37","build38","b35","b36","b37","b38")
  build.new <- c("hg15","hg38",rep(c("hg17","hg18","hg19","hg38"),times=5))
  if(show.valid) { return(cbind(valid=build.alt,mapsTo=build.new)) }
  build <- build.new[match(tolower(build),build.alt)]
  if(is.null(build)) { build <- "hg19"; warning("build was NULL (see getOption('ucsc')), set to hg19") }
  if(any(is.na(build))) { 
    warning("Illegal build parameter '",build[1],"', defaulting to hg19") 
    build[is.na(build)] <- "hg19" 
  }
  if(allow.multiple) {
    return(build)
  } else {
    return(build[1])
  }
}


# ok as long as at least one non-missing snp in the summary
#' See snpStats::col.summary. Same in every way, except for the undesirable
#' behaviour of snpStats when a SNP has 100% missing values it is ignored
#' in the call-rate summary (rather than given a zero). This can unintentionally
#' mean that call-rate filters do not filter SNPs with 100% missing values.
#' This function is simply a wrapper that cleans up this problem.
# col.summary2 <- function(object,...) {
#   if(!is(object)[1]=="SnpMatrix")   { stop("'object' must be a SnpMatrix (snpStats package)") } 
#   if(any(!names(list(...)) %in% c("rules","uncertain"))) { 
#     stop("... contained invalid arguments to snpStats::col.summary") }
#   
# }

#internal
pduplicated <- function(X) {
  if(length(Dim(X))>1) {  stop("can only enter a vector into this function") }
  return((duplicated(X,fromLast=T) | duplicated(X,fromLast=F)))
}


#internal
comma <- function(...) {
  paste(...,collapse=",")
}






#internal function to properly sort chromosome labels as text
order.chr <- function(chrs) {
  # sort chr nms
  if(is.numeric(chrs)) { chrs <- paste(chrs) }
  if(!is.character(chrs)) { stop("chrs should be a character or integer vector") }
  asn <- function(X) { suppressWarnings(as.numeric(X)) }
  textz <- is.na(asn(chrs))
  nums <- chrs[!textz]
  txts <- chrs[textz]
  ns <- which(!textz)[order(asn(nums))]
  #print(max(ns,na.rm=TRUE))
  #print(order(txts)); print(txts)
  #print(which(textz))
  ts <- which(textz)[order(txts)]
  out <- c(ns,ts)
  return(out)
}

#internal
sort.chr <- function(chr) { chr[order.chr(chr)] }


#internal
# standardize snp ids so they would always appear the same, all _,.;, etc replaced with _
# all names leading with a number preceeded with X. mainly 'make.names' standard R-conventions
clean.snp.ids <- function(snpid.list) {
  snpid.list <- make.names(snpid.list)
  snpid.list <- gsub(".","_",snpid.list,fixed=T)
  return(snpid.list)
}


## internal function with extra mapping hits for immunochip that aren't in the chain file for 36-37
hard.coded.conv <- function() {
  chrzM <- c("7","7","9","5","7","14","17","4","8","8","15","7","6","6","2","2","4","17","19")
  pos36M <- c("142154515","142160115","132183222","17767156","141943232","27591752","41560151",
              "103951975","17510484","17501697","81350958","141911612","74644736","74644390",
              "1203295","21043693","4020119","78644427","52569727")
  pos37M <- c("142474939","142480539","135153668","17731427","142224511","28521898","44204373",
              "103732866","17466212","17457420","83559954","142108941","74588007","74587661",
              "1213294","21190209","3969218","81051007","47877928")
  rsidM <- c("rs10952532","rs10952534","rs11243704","rs11953245","rs17274","rs1952843",
             "rs2016730","rs223413","rs2427715","rs2517168","rs2621228",
             "rs2855938","rs2917890","rs2917891","rs4971417","rs6547409","rs6842556",
             "rs7502442","rs755327")
  chrzI <- c("3","3","3","3","3","6","7","7","8","8","8","8","17","17","X")
  pos36I <- c("50875374","50882163","50885514","50908888","195567372","119257505",
              "50323690","67383261","10961083","10961130","10975096","10975127",
              "21628754","59781521","75211826")
  pos37I <- c("50900354","50907147","50910499","50908888","194086083","119150813",
              "50353144","67745402","10923673","10923720","10937686","10937717",
              "21704627","59781521","75295444")
  rsidI <- c("rs12639243","rs62717061","rs4346541","imm_3_50908888","rs4974514","rs284919",
             "rs7804185","rs3113138","rs2898255","rs2409687","rs7827367","rs6601557",
             "rs17052332","rs1131012","rs929032")
  chrz <- c(chrzI,chrzM)
  pos36 <- c(pos36I,pos36M)
  pos37 <- c(pos37I,pos37M)
  rsid <- c(rsidI,rsidM)
  return(list(chr=chrz,pos36=pos36,pos37=pos37,rs.id=rsid))
}

substitute.36s <- function(granges=NULL) {
 pos36 <- matrix(c(c(3,"imm_3_50875337",50875337),
 c(3,"imm_3_50882163",50882163),
 c(3,"imm_3_50908888",50908888),
 c(7,"rs1574660",141711704),
 c(17,"rs1131012",59781521),
 c("X","rs16994803",147800037),
 c("X","rs12013571",148519953)),ncol=3,byrow=T)
 if(is.null(granges)) { return(pos36) }
 typ <- is(granges)[1]
 X <- as(granges,"GRanges")
 if(!is(X)[1]=="GRanges") { stop("conversion to GRanges failed") }
 ii <- narm(match(pos36[,2],rownames(X)))
 start(X)[ii] <- rep(1,length(ii))
 end(X)[ii] <- as.numeric(pos36[,3])
 start(X)[ii] <- as.numeric(pos36[,3])
 X <- as(X,typ)
 return(X)
}




# internal
# Remove trailing letter from non-unique rs-ids
#
# #@examples
# snp.ids <- rsnpid(25)
# snp.ids[1:2] <- paste0(snp.ids[1:2],"b")
# snp.ids[19:20] <- paste0(snp.ids[19:20],"c")
# snp.ids[6:7] <- paste0(snp.ids[6:7],"d")
# snp.ids[11:12] <- paste0(snp.ids[11:12],"a")
# snp.ids
# rmv.trail(snp.ids)
rmv.trail <- function(rs.ids,suffix=c("b","c","d","a")) {
  if(!is.character(suffix)) { stop("suffix must a character vector") }
  ind <- NULL
  for (cc in 1:length(suffix)) {
    ind <- c(ind,grep(suffix[cc],rs.ids))
  }
  ind <- unique(ind)
  X <- rs.ids[ind]
  nX <- nchar(X)
  last.chars <- substr(X,nX,nX)
  sufz <- (last.chars %in% c("a","b","c","d"))
  X[sufz] <- substr(X[sufz],1,nX[sufz]-1)
  rs.ids[ind] <- X
  return(rs.ids)
}

# internal
# Add trailing letter(s) to non-unique rs-ids
# #@examples
# snp.ids <- rsnpid(15)
# snp.ids
# add.trail(snp.ids)
# snp.ids <- snp.ids[sample(15,30,replace=TRUE)]
# snp.ids
# add.trail(snp.ids)
add.trail <- function(rs.ids,suffix=c("b","c","d","a")) {
  rs.ids <- rmv.trail(rs.ids)
  for (txt in suffix) {
    dupz <- duplicated(rs.ids)
    if(length(which(dupz))>0) {
      rs.ids[dupz] <- paste0(rmv.trail(rs.ids[dupz]),txt)
    }
  }
  dupz <- duplicated(rs.ids)
  if(length(which(dupz))>0) { warning("more than ",length(suffix),
                                      " duplications of at least 1 individial rs-id, suffixes exhausted, duplicates remain")
  }
  return(rs.ids)
}

# internal
rsampid <- function(n,pref="ID0") { paste0(pref,pad.left(1:n,"0")) }


# internal
rsnpid <- function(n) { 
  id.len <- sample(c(3:8),n,replace=T,prob=c(0.01, 0.01, 0.01, 0.10, 0.50, 0.37))
  each.id <- function(l) { sapply(l,function(n) { paste(replicate(n,sample(1:9,1)),collapse="",sep="") }) }
  sufz <- each.id(id.len)
  ids <- paste0("rs",sufz)
  return(ids)
}


# internal
# Remove that pesky 'elementmetadata.' prefix from column names that have been converted from GRanges
emd.rmv <- function(X, rmv.genome=TRUE) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(has.method("mcols",X, where=environment(emd.rmv))) {
    colnames(mcols(X)) <- gsub("elementMetadata.","",colnames(mcols(X)))
    ii <- which(colnames(mcols(X))=="genome")
    if(length(ii)>0) {
      gn <- mcols(X)[,ii[1]]
      if(length(unique(gn))==1) {
        mcols(X) <- mcols(X)[,-ii[1]]
      }
    }
  } else {
    if(has.method("colnames",X, where=environment(emd.rmv))) {
      colnames(X) <- gsub("elementMetadata.","",colnames(X))
      ii <- which(colnames(X)=="genome")
      if(length(ii)>0) {
        gn <- X[,ii[1]]
        if(length(unique(gn))==1) {
          X <- X[,-ii[1]]
        }
      }
    } else {
      stop("X must have column names, expecting GRanges, RangedData or data.frame")
    }
  }
  return(X)
}



chrOrder2 <- function (chr.names) {
  if(!is.character(chr.names)) { warning("expecting character() type for chr.names argument") }
  simple.names = gsub("^chr", "", chr.names)
  name.is.numeric = grepl("^[0-9]+$", simple.names, perl = T)
  numeric.names = chr.names[name.is.numeric][order(as.numeric(simple.names[name.is.numeric]))]
  non.numeric.names = chr.names[!name.is.numeric][order(chr.names[!name.is.numeric])]
  all.names = c(numeric.names, non.numeric.names)
  return(all.names)
}


# internal# iFunctions
chrNames2 <- function(X) {
  #requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(nrow(X)==0) { return(character(0)) }
  out <- as.character(unique(seqnames(X)))
  return(out)
}


# internal from genoset
TGORD <- function (ds, strict = TRUE) {
  if (strict == TRUE) {
    if (!isTRUE(all.equal(chrOrder2(chrNames2(ds)), chrNames2(ds)))) {
      ds = ds[chrOrder2(chrNames2(ds))]
    }
  }
  row.order = order(as.integer(space(ds)), start(ds))
  if (is.unsorted(row.order)) {
    return(ds[row.order, , drop = FALSE])
  }
  else {
    return(ds)
  }
}


# internal from genoset
TGOGR <- function (ds, strict = TRUE) {
  if (strict == TRUE) {
    if (!isTRUE(all.equal(chrOrder(seqlevels(ds)), seqlevels(ds)))) {
      seqlevels(ds) = chrOrder(seqlevels(ds))
    }
  }
  row.order = order(as.integer(seqnames(ds)), start(ds))
  if (is.unsorted(row.order)) {
    ds = ds[row.order, , drop = FALSE]
  }
  return(ds)
}



# internal # iFunctions
# version of toGenomeOrder() that is guaranteed to work for IRanges or GRanges
toGenomeOrder2 <- function(X,...) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges"); requireNamespace("genoset")
  if(is(X)[1] %in% c("GRanges","RangedData","ChipInfo")) {
    if(is(X)[1]=="RangedData") {
      return(TGORD(X))
    } else {
      return(TGOGR(X))
    }
  } else {
    typ <- is(X)[1]
    if(!typ %in% c("GRanges","RangedData","ChipInfo")) { warning("unsupported type '",typ,"' for toGenomeOrder(), failure likely") }
    alreadyThere <-("strand" %in% colnames(X))
    out <- TGOGR(as(X,"GRanges"),strict=T) #genoset::
    X <- as(out,"RangedData")
    if(("strand" %in% colnames(X)) & !alreadyThere) {
      X <- X[,-which(colnames(X) %in% "strand")]
    }
    return(X)
  }
}

#internal # iFunctions
# version of chrInfo() that is guaranteed to work for IRanges or GRanges
chrInfo2 <- function(X) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(is(X)[1] %in% c("GRanges","ChipInfo")) {
    return(genoset::chrInfo(X))
  } else {
    typ <- is(X)[1]
    if(!typ %in% c("GRanges","RangedData","ChipInfo")) { warning("unacceptable type '",typ,"' for chrInfo2(), failure likely") }
    out <- chrInfo(as(X,"GRanges"))
    return(out)
  }
}

#internal # iFunctions
# version of chrIndices() that is guaranteed to work for IRanges or GRanges
chrIndices2 <- function(X,...) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(is(X)[1] %in% c("GRanges","ChipInfo")) {
    return(genoset::chrIndices(X,...))
  } else {
    typ <- is(X)[1]
    if(!typ %in% c("GRanges","RangedData","ChipInfo")) { warning("unacceptable type '",typ,"' for chrIndices2(), failure likely") }
    out <- chrIndices(as(X,"GRanges"))
    return(out)
  }
}

#internal # iFunctions
# version of chr() that is guaranteed to work for IRanges or GRanges
chr2 <- function(X) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(is(X)[1] %in% c("GRanges","ChipInfo")) {
    return(genoset::chr(X))
  } else {
    if(is(X)[1]=="RangedData") {
      return(space(X))
    } else {
      if(is.null(X)) { warning("X was NULL, expecting RangedData/GRanges"); return(NULL) }
      warning("chr2() function applies only to RangedData objects, attempting to pass ",is(X)[1]," to chr()")
      return(genoset::chr(X))
    }
  }
}



#internal
make.divisor <- function(unit=c("b","kb","mb","gb"), par.name="scale (scl)") {
  valid.units <- c("k","m","g","b")
  unit <- tolower(unit[1]);
  unit <- substr(unit,1,1)
  if(!unit %in% valid.units) { warning("invalid entry to ",par.name," defaulting to base-pairs") ; unit <- "b" }
  divisor <- switch(unit,k=1000,m=10^6, g=10^9, b=1)
  return(divisor)
}

#internal
plotdf <- function(expr,fn="myTempPlot.pdf") {
  pdf(fn)
{ expr }
dev.off()
cat("wrote plot to",cat.path(getwd(),fn),"\n")
}


# gets limits of a plot space on current device
plot.get.area <- function() {
  success <- tryCatch(mine <- par("usr"),error=function(e) { F } )
  if(all(!success)) { warning("could not get plot limits - no plot open perhaps?"); return(NULL) }
  xlim=mine[1:2]
  ylim=mine[3:4]
  return(list(xlim=xlim,ylim=ylim))
}


#internal function
chrnums.to.txt <- function(X,do.x.y=TRUE) {
  cond <- paste(X) %in% paste(1:22)
  if(any(cond)) { X[cond] <-  paste0("chr",X[cond]) }
  if(do.x.y) {
    X <- gsub("X","chrX",X)
    X <- gsub("Y","chrY",X)
    X <- gsub("23","chrX",X)
    X <- gsub("24","chrY",X)
    X <- gsub("25","chrXY",X)
    X <- gsub("26","chrM",X)
    X <- gsub("chrXchrY","XY",X)
    X <- gsub("chrYchrX","YX",X) 
    X <- gsub("M","chrM",X)
    X <- gsub("XY","chrXY",X)
    X <- gsub("chrchr","chr",X)
  } else {
    X[X %in% paste(23:100)] <- paste0("chr",X[X %in% paste(23:100)])
  }
  return(X)
}

#internal function
chrnames.to.num <- function(X,keep.let=TRUE) {
  X <- tolower(X)
  if(!keep.let) {
    X <- gsub("chrM","26",X)
    X <- gsub("chrXY","25",X) 
    X <- gsub("chrY","24",X)
    X <- gsub("chrX","23",X)
  } else {
    X <- gsub("chrM","M",X)
    X <- gsub("chrXY","XY",X) 
    X <- gsub("chrY","Y",X)
    X <- gsub("chrX","X",X)
  }
  X <- gsub("chrchr","",X) 
  X <- gsub("chr","",X)
  X <- toupper(X)
  return(X)
}


# iFunctions
# internal, tidy chromosome names using extra chromosomal annotation into rough chromosomes
tidy.extra.chr <- function(chr,select=FALSE) {
  # most relevant to hg18
  chr <- paste(chr)
  SEL_c6 <- grep("c6",chr,ignore.case=T)
  SEL_c5 <- grep("c5",chr,ignore.case=T)
  SEL_NT <- grep("NT",chr,ignore.case=T)
  # most relevant to hg19
  SEL_LRG <- grep("LRG",chr,ignore.case=T)
  SEL_HG <- grep("HG",chr,ignore.case=T)
  SEL_GL <- grep("GL",chr,ignore.case=T)
  SEL_HS <- grep("HSCHR",chr,ignore.case=T)
  if(select) {
    # create TRUE/FALSE as to whether list elements have weird chromosome codes
    all <- unique(c(SEL_c6,SEL_c5,SEL_NT,SEL_LRG,SEL_HG,SEL_GL,SEL_HS))
    return(!1:length(chr) %in% all)
  } else {
    # transform weird chromosomes into more palatable codes
    chr[SEL_c6] <- 6  # prevent issues with c6_COX, c6_QBL  
    chr[SEL_c5] <- 5  # prevent issues with c5_H2  
    chr[SEL_NT] <- "Z_NT"  # merge all NT regions to one label
    chr[SEL_LRG] <- "Z_LRG"  # merge all NT regions to one label
    chr[SEL_HG] <- "Z_HG"  # merge all NT regions to one label
    chr[SEL_GL] <- "Z_GL"  # merge all NT regions to one label
    X <- names(table(chr))
    X <- X[grep("HSCHR",X)]
    if(length(X)>0) {
      HSC <- gsub("_","",substr(gsub("HSCHR","",X),1,2))
      for(cc in 1:length(X)) {
        #cat("replacing ",X[cc]," with ",HSC[cc],"\n",sep="")
        chr[chr==X[cc]] <- HSC[cc]
      }
    }
    return(chr)
  }  
}


#internal
gene.duplicate.report <- function(ga,full.listing=F,colname="gene",silent=FALSE) {
  # for a RangedData object, report on any multiple listings for the same gene
  if(is(ga)[1]!="RangedData") { warning("not a RangedData object") ; return(NULL) }
  if(colname=="gene") {
    if("gene" %in% tolower(colnames(ga)))
    { 
      gene.col <- (which(tolower(colnames(ga)) %in% c("gene","genes","geneid")))
    } else {
      gene.col <- 0
    }
  } else {
    if(colname %in% colnames(ga)) { 
      gene.col <- which(colnames(ga)==colname) 
    } else { 
      stop("colname not found in ga") 
    } 
  }
  if(length(gene.col)>0) { gene.col <- gene.col[1] } else { warning("no 'gene' column"); return(NULL) }
  colnames(ga)[gene.col] <- "gene" #force this colname
  duplicate.report <- T  ### when would this be FALSE???
  if(duplicate.report) {
    culprits <- unique(ga$gene[which(duplicated(ga$gene))])
    n.gene.multi.row <- length(culprits)
    culprit.ranges <- ga[ga$gene %in% culprits,]
    total.culprit.rows <- nrow(culprit.ranges)
    start.same.ct <- end.same.ct <- 0; which.ss <- NULL
    for (cc in 1:length(culprits)) { 
      mini <- (ga[ga$gene %in% culprits[cc],]) 
      if(full.listing) {
        cat(colname,":",culprits[cc],"# same start:",anyDuplicated(start(mini)),
            "# same end:",anyDuplicated(end(mini)),"\n") }
      start.same.ct <- start.same.ct+anyDuplicated(start(mini))
      end.same.ct <- end.same.ct+anyDuplicated(end(mini))
      if(anyDuplicated(start(mini)) | anyDuplicated(end(mini))) { which.ss <- c(which.ss,cc) }
    }
    if(!silent) {
      cat(" ",colname,"s with split ranges:\n"); print(culprits,quote=F); cat("\n")
      cat(" ",colname,"s with same start or end:\n"); print(culprits[which.ss],quote=F); cat("\n")
      cat(" total ",colname,"-segments with same start",start.same.ct,"; total with same end:",end.same.ct,"\n")
    }
  }
  return(culprits)
}



# internal function
validate.dir.for <- function(dir,elements,warn=F) {
  # in case the 'dir' input list object is not the standardised list form, convert
  # allows flexible use of list or regular directory specifications in plumbCNV functions
  if(is.null(dir)) { cat("directory empty\n"); return(NULL) }
  if(!is.list(dir)) {
    if(warn) { cat(elements[cc],"'dir' object wasn't a list\n")}
    dir <- as.list(dir); names(dir)[1:length(dir)] <- elements[1:length(dir)] 
  }
  for (cc in 1:length(elements)) {
    if(!elements[cc] %in% names(dir)) { 
      dir[[paste(elements[cc])]] <- "" ;
      if(warn) { stop(paste("dir$",elements[cc]," was empty.. set to current\n",sep="")) } 
    }
  }
  return(dir)
}



################# end internals #################
#if(getwd()!= "/home/ncooper"){
#  require(genoset)
#}

# examples for package
# 
# all.support <- chip.support()
# snp.info <- ChipInfo(chr=all.support[,"Chr"],pos=all.support[,"Pos"],ids=rownames(all.support),chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                      A1=all.support[,"A1"], A2=all.support[,"A2"])
# 
# gr.snp.info <- with(all.support,makeGRanges(chr=Chr,pos=Pos,row.names=rownames(all.support)))
# 
# snp.info <- ChipInfo(GRanges=gr.snp.info,chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                      A1=all.support[,"A1"], A2=all.support[,"A2"])
# 
# snp.info[chr(snp.info)=="MT",] # look at the mitochondrial SNP
# QCcode(snp.info)[chr(snp.info)=="MT"] <- 1 # exclude it by changing the fail code
# snp.info[["MT"]] # revisit and see it now registers as 'fail'
# QCfail(snp.info)
# xx <- conv.36.37(chr=all.support[,"Chr"],pos=all.support[,"Pos"],ids=rownames(all.support))
# ucsc(si36[["XY"]])
# si37 <- convTo37(snp.info)
# si36 <- convTo36(si37)
# snp.info <- ChipInfo(chr=all.support[,"Chr"],pos=all.support[,"Pos37"],ids=rownames(all.support),chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                     A1=all.support[,"A1"], A2=all.support[,"A2"],build=37)


#' Class to represent SNP annotation for a microarray
#' 
#' This class annotates a microarray SNP chip with data for each SNP including chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18, hg19 or hg38 (build 36/37/38) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, 
#' and you can also load your own annotation if using a different chip. The class
#' is basically a GRanges object, modified to always have columns for A1, A2 (alleles), 
#' rs.id, and a quality control flag. The default display is tidier than GRanges, it has
#' nice coersion to and frame data.frame and subsetting by chromosome using [[n]] has been
#' added, in addition to normal [i,j] indexing native to GRanges.
#' Note that with this package the first time annotation is used it might be slow, but
#' subsequent calls should be fast.
#' METHODS
#'  "[[", show, print, length, dim, rownames, initialize
#'  build, chip, rs.id, A1, A2, QCcode, QCcode<-, QCpass, QCfail
#'  convTo36, convTo37, convTo38
#' COERCION
#'  can use 'as' to convert to and from: GRanges, RangedData, data.frame
#'@section Fields: 
#'  \describe{
#'    \item{\code{seqnames}:}{Object of class \code{"Rle"}, containing chromosomes for each range, see GRanges.}
#'    \item{\code{ranges}:}{Object of class \code{"IRanges"}, containing genomic start and end, see GRanges.}
#'    \item{\code{strand}:}{Object of class \code{"Rle"}, containing plus or minus coding for forward or reverse strand, see GRanges.}
#'    \item{\code{seqinfo}:}{Object of class \code{"Seqinfo"}, containing chromosome listing, see GRanges.}
#'    \item{\code{chip}:}{Name, class \code{"character"}, containing user description of the chip, e.g, 'immunoChip'.}
#'    \item{\code{build}:}{Object of class \code{"character"}, annotation version, e.g, hg18, hg19, hg38, etc.}
#'    \item{\code{elementMetaData}:}{Object of class \code{"DataFrame"}, see GRanges, but with specific column names:
#'  A1, A2, QCcode and rs.id.}
#'  }
# @name ChipInfo-class
#' @aliases ChipInfo-method
#' @rdname ChipInfo-class
#' @exportClass ChipInfo
#' @author Nick Cooper
setClass("ChipInfo",
         contains="GRanges",
         slots=list(chip="character", 
                    build="character"),
         prototype=prototype(
              seqnames=Rle(factor()),
              ranges=IRanges::IRanges(),
              strand=Rle(GenomicRanges::strand()),
              elementMetadata=DataFrame(A1=NULL, A2=NULL,  QCcode=integer(), rs.id=NULL, chip.id=NULL),
              seqinfo=Seqinfo(),
              metadata=list(),
              chip=character(), 
              build=character()
            )
)



#' rownames method for GRanges objects
#' 
#' rownames: Returns the row names.
# @name rownames
#' @param x a GRanges object
#' @return rownames: Character vector of row names (e.g, SNP IDs).
#' @rdname GRanges-methods
#' @exportMethod rownames
setMethod("rownames", "GRanges", function(x) names(x@ranges))
setMethod("rownames<-", "GRanges", function(x,value) { names(x@ranges) <- value; return(x) })



# seqnames="Rle",
# ranges="IRanges",
# strand="Rle",
# elementMetadata="DataFrame",
# seqinfo="Seqinfo",
# metadata="list",
# ,
# prototype=prototype(
#   seqnames=IRanges::Rle(factor()),
#   ranges=IRanges::IRanges(),
#   strand=IRanges::Rle(GenomicRanges::strand()),
#   elementMetadata=IRanges::DataFrame(A1=NULL, A2=NULL,  QCcode=integer(), rs.id=NULL),
#   seqinfo=GenomicRanges::Seqinfo(),
#   chip=character(), 
#   build=character(),
# )

#' rownames method for ChipInfo objects
#' 
#' rownames: Returns the row names.
# @name rownames
#' @param x a ChipInfo object
#' @return rownames: Character vector of row names (SNP IDs).
#' @rdname ChipInfo-methods
#' @exportMethod rownames
setMethod("rownames", "ChipInfo", function(x) names(x@ranges))


#' dim method for ChipInfo objects
#' 
#' dim: Returns the dimension
# @name dim
#' @return dim: same as length
#' @rdname ChipInfo-methods
#' @exportMethod dim
setMethod("dim", "ChipInfo", function(x) dim(mcols(x)))


#' Length method for ChipInfo objects
#' 
#' length: Returns the number of rows
# @name length
#' @return length: integer, number of rows, same as inherited nrow()
#' @rdname ChipInfo-methods
#' @exportMethod length
setMethod("length", "ChipInfo", function(x) length(x@ranges))


#' Retrieve the Chip name for ChipInfo
#' 
#' Simply returns the name of the chip, e.g, 'ImmunoChip'
# @name chip
#' @param x a ChipInfo object
#' @return character string
#' @rdname chip-methods
#' @export
setGeneric("chip", function(x) standardGeneric("chip"))


#' @exportMethod chip
#' @rdname chip-methods
setMethod("chip", "ChipInfo", function(x) x@chip)


#' Retrieve the UCSC build for a ChipInfo object
#' 
#' Returns the UCSC build of the chip object, e.g, 'hg18', 'hg19', or 'hg38'
# @name ucsc
#' @param x a ChipInfo object
#' @return character, 'hg18', 'hg19', or 'hg38'
#' @rdname ucsc-methods
#' @export
setGeneric("ucsc", function(x) standardGeneric("ucsc") )

#' @rdname ucsc-methods
#' @exportMethod ucsc
setMethod("ucsc", "ChipInfo", function(x) x@build)


#' Access rs-ids for ChipInfo
#' 
#' Returns the rs-ids for the chip object, e.g, "rs689", etc
#' Only if these are annotated internally, or else a vector of NAs
# @name rs.id
#' @param x a ChipInfo object
#' @param b logical, whether to show 'b' suffixes on rs.ids which
#' are created in the background to allow duplicate ids to be uniquely
#' represented for lookup and reference purposes.
#' @return rs-ids: character vector of IDs (or NAs)
#' @rdname rs.id-methods
#' @export
setGeneric("rs.id", function(x,b=TRUE) standardGeneric("rs.id") )

#' @rdname rs.id-methods
#' @exportMethod rs.id
setMethod("rs.id", "ChipInfo", function(x,b=TRUE) { 
  u <- mcols(x) ;  
  if("rs.id" %in% colnames(u)) { 
    U <- u[,"rs.id"] 
    if(!b) { U <- gsub("b","",U) }
    return(U)
  } else { return(NULL) } 
})


#' Access chip-ids for ChipInfo
#' 
#' Returns the chip-ids for the chip object, e.g, "imm_1_898835", etc
#' Only if these are annotated internally, or else a vector of NAs
#' Note that the main purpose of this is because sometimes chip-ids
#' do not satisfy conditions to be an R column/row name, e.g, start
#' with a number, illegal characters, etc. So this allows certain
#' functions to return the actual chip names that would match the 
#' official manifest. These will largely be the same as the rownames,
#' but the rownames will always be valid R column names, converted from
#' the original using clean.snp.ids() [internal function]
# @name chip-id
#' @param x a ChipInfo object
#' @return chip ids: character vector of IDs (or NAs)
#' @rdname chipId-methods
#' @export
setGeneric("chipId", function(x) standardGeneric("chipId") )

#' @rdname chipId-methods
#' @exportMethod chipId
setMethod("chipId", "ChipInfo", function(x) { 
  u <- mcols(x) ;  
  if("chip.id" %in% colnames(u)) { 
    U <- u[,"chip.id"] 
    return(U)
  } else { return(NULL) } 
})

#' Access alleles for ChipInfo
#' 
#' A1/A2: Returns the letter for the A1/A2 alleles for the chip object, 
#' e.g, 'A','C','G','T', etc
#' Only if these are annotated internally, or else a vector of NAs
#' @param x a ChipInfo object
#' @rdname allele-methods
#' @family Alleles
#' @export
setGeneric("A1", function(x) standardGeneric("A1") )

#' @rdname allele-methods
#' @family Alleles
#' @exportMethod A1
setMethod("A1", "ChipInfo", function(x) { u <- mcols(x) ;  if("A1" %in% colnames(u)) { u[,"A1"] } else { NULL } })


#' @return character vector of allele codes (or NAs)
#' @rdname allele-methods
#' @family Alleles
#' @export
setGeneric("A2", function(x) standardGeneric("A2") )

#' @rdname allele-methods
#' @family Alleles
#' @exportMethod A2
setMethod("A2", "ChipInfo", function(x) { u <- mcols(x) ;  if("A2" %in% colnames(u)) { u[,"A2"] } else { NULL } })


#' Set quality control pass or fail codes for ChipInfo
#' 
#' A1<-/A2<-: Allows user to set the allele codes for each SNP of the chip object, 
#' e.g, A,C,G,T,K, etc. If you are using allele codes this is likely to 
#' necessary as each genotyping produces a different set of allele codes.
#' If using in conjunction with snpStats, remember that allele codes are
#' always flipped to be alphabetical, so the reference allele is the later
#' letter in the alphabet. Note, assignment to A2 needs to be done separately.
#' @param value new allele codes, e.g, A,C,G,T
#' @return A1<-: updates the ChipInfo object specified with new allele codes for the 'A1' slot
#' @rdname allele-methods
#' @family Alleles
#' @export
setGeneric("A1<-", function(x,value) standardGeneric("A1<-") )


#' @rdname allele-methods
#' @family Alleles
#' @exportMethod "A1<-"
setMethod("A1<-", "ChipInfo", function(x,value) {
  return(.updateAllele(x,value,"A1"))
} )




#' @return A2<-: updates the ChipInfo object specified with new allele codes for the 'A2' slot
#' @rdname allele-methods
#' @family Alleles
#' @export
setGeneric("A2<-", function(x,value) standardGeneric("A2<-") )


#' @rdname allele-methods
#' @family Alleles
#' @exportMethod "A2<-"
setMethod("A2<-", "ChipInfo", function(x,value) {
  return(.updateAllele(x,value,"A2"))
} )


#internal
.updateAllele <- function(x,value, allele="A1") {
  if(length(Dim(value))!=1) { stop("value must be a vector") }
  if(length(x)==length(value)) {
    if(is.character(value)) {
      mcols(x)[,allele] <- paste(value)
    } else {
      stop("only character values can be inserted into the ",allele," column, A,C,G,T, etc")
    }
  } else {
    stop("mismatching lengths, tried to insert ",length(value),"new allele codes into ChipInfo with ",length(x)," rows")
  }
  return(x)
}


#' Access quality control pass or fail codes for ChipInfo
#' 
#' Returns the pass or fail codes for each SNP of the chip object, 
#' e.g, 0,1,..,n etc
#' Only if these are added manually, or else all will be 'pass' (=0)
# @name QCcode
# why not? param x a ChipInfo object
#' @return integer vector of pass/fail codes
#' @rdname QC-methods
#' @family QC
#' @export
setGeneric("QCcode", function(x) standardGeneric("QCcode") )

#' @rdname QC-methods
#' @family QC
#' @exportMethod QCcode
setMethod("QCcode", "ChipInfo", function(x) { 
  u <- mcols(x) ;  if("QCcode" %in% colnames(u)) { u[,"QCcode"] } else { NULL } 
})


#' Set quality control pass or fail codes for ChipInfo
#' 
#' QCcode<-: Allows user to set the pass or fail codes for each SNP of the chip object, 
#' e.g, 0,1,..,n etc. 0 is always pass, >0 is always fail, but each integer
#' can be used to represent a different failure type, or for simplicity, stick
#' to 0 and 1, ie, just pass and fail.
#' @param x a ChipInfo object
#' @param value new pass/fail codes, e.g, 0,1,...,n
#' @return QCcode<-: updates the object specified with new pass/fail codes for the 'QCcode' slot
#' @rdname QC-methods
#' @family QC
#' @export
setGeneric("QCcode<-", function(x,value) standardGeneric("QCcode<-") )


#' @rdname QC-methods
#' @family QC
#' @exportMethod QCcode
setMethod("QCcode<-", "ChipInfo", function(x,value) {
  if(length(x)==length(value)) {
    if(is.numeric(value)) {
      mcols(x)[,"QCcode"] <- as.integer(value)
    } else {
      stop("only numeric values can be inserted into the QCcode column, 0=pass, higher integers are failure codes")
    }
  } else {
    stop("mismatching lengths, tried to insert ",length(value),"new values into ChipInfo with ",length(x)," rows")
  }
  return(x)
} )


#' Filter ChipInfo to for only SNPs passing QC
#' 
#' QCpass: Returns the subset of the ChipInfo object for which SNPs pass quality
#' control, according to the QCcodes() slot == 0.
# @name QCpass
# @param x a ChipInfo object
#' @return QCpass: ChipInfo object for which SNPs pass quality control
#' @rdname QC-methods
#' @family QC
#' @export
setGeneric("QCpass", function(x) standardGeneric("QCpass") )


#' Filter ChipInfo to for only SNPs failing QC
#' 
#' QCfail: Returns the subset of the ChipInfo object for which SNPs fail quality
#' control, according to the QCcodes() slot > 0.
# @name QCfail
# @param x a ChipInfo object
#' @param type integer between 1 and 100, failure type (user can assign own coding scheme)
#' @return QCfail: ChipInfo object for which SNPs fail quality control
#' @rdname QC-methods
#' @family QC
#' @export
setGeneric("QCfail", function(x,type=NA) standardGeneric("QCfail") )


#' @rdname QC-methods
#' @family QC
#' @exportMethod QCpass
setMethod("QCpass", "ChipInfo", function(x) { 
  ii <- which(QCcode(x)==0)
  if(length(ii)>0) { return(x[ii,]) } else { warning("No SNPs passed QC"); return(NULL) } })


#' @rdname QC-methods
#' @family QC
#' @exportMethod QCfail
setMethod("QCfail", "ChipInfo", function(x,type=NA) { 
  ii <- which(QCcode(x)!=0)
  if(is.numeric(type)) { 
    if(type %in% 1:100) {
      ii <- which(QCcode(x)==type) 
    } else { 
      warning("type must be an integer between 1 and 100, returning all failures") 
    }
  }
  if(length(ii)>0) { return(x[ii,]) } else { warning("All SNPs passed QC"); return(NULL) } })


#' Subset ChipInfo by chromosome
#' 
#' Returns the subset of the ChipInfo object for which SNPs are on
#' the chromosome specified, by either number or character.
#' @param x a ChipInfo object
#' @param i a chromosome number or letter, i.e, one of seqlevels(x)
#' @param j always leave missing, not applicable for this method.
#' @param ... further arguments - again there should not be any
#' @return ChipInfo object for the subset of SNPs on chromosome i
#' @rdname ChipInfo-subset
#' @exportMethod "[["
setMethod("[[", "ChipInfo", function(x,i,j,...) { 
  dotArgs <- list(...)
  if (length(dotArgs) > 0)
    dotArgs <- dotArgs[names(dotArgs) != "exact"]
  if (!missing(j) || length(dotArgs) > 0)
    stop("invalid subsetting")
  if (missing(i))
    stop("subscript is missing")
  if (!is.character(i) && !is.numeric(i)) 
    stop("invalid subscript type")
  if (length(i) < 1L)
    stop("attempt to select less than one element")
  if (length(i) > 1L)
    stop("attempt to select more than one element")
  cn <- chrNames(x)
  if (is.numeric(i) && !is.na(i) && (i < 1L || i > length(cn)))
    stop("subscript out of bounds")
  # do the selection #
  if(i %in% paste(chr(x))) {
    out <- x[chr(x)==i,]
  } else {
    if(is.numeric(i)) {
      out <- x[match(chr(x),chrNames(x))==i,]
    } else {
      stop("unknown index")
    }
  }
  out@build <- x@build
  out@chip <- x@chip
  return(out)
} )


#' Convert ChipInfo between build 36/37/38 coordinates
#' 
#' Returns the a ChipInfo object with positions updated to build 36/37/38
#' coordinates, assuming that the build() slot was entered correctly. 
#' Ensure that the value of ucsc(x) is correct before
#' running this function for conversion; for instance, if the coordinates 
#' are already build 37/hg19, but ucsc(x)!="hg19" (incorrect value), then
#' these coordinates will be transformed in a relative manner rendering the
#' result meaningless.
# @name convTo37
#' @param x a ChipInfo object
#' @return convTo37: Returns a ChipInfo object with the build updated to hg19 coordinates
#' @family conversion
#' @rdname conv-methods
#' @export
setGeneric("convTo37", function(x) standardGeneric("convTo37"))
          
#' @family conversion
#' @aliases convTo37
#' @rdname conv-methods
#' @exportMethod convTo37
setMethod("convTo37", "ChipInfo", function(x) {
  if(ucsc.sanitizer(ucsc(x)) %in% c("hg18","hg38")) {
    if(ucsc.sanitizer(ucsc(x)) == c("hg18")) {
      u <- conv.36.37(ranges=as(x,"GRanges"))
    } else {
      u <- conv.38.37(ranges=as(x,"GRanges"))
    }
    if(length(u)==length(x)) { 
      x@ranges <- u@ranges
      all.eq <- TRUE
      if(all.eq) { all.eq <- length(seqlevels(x))==length(seqlevels(u)) }
      if(all.eq) { all.eq <- all(sort(seqlevels(x))==sort(seqlevels(u))) }
      if(!all.eq) {
        #print(seqlevels(x)); print(seqlevels(u))
        #warning("conversion altered the chromosomes"); 
        seqlevels(x) <- c(seqlevels(x),seqlevels(u)[!seqlevels(u) %in% seqlevels(x)])  #x@seqinfo <- u@seqinfo 
      }
      xx <- as(x@seqnames,"character")
      uu <- as(u@seqnames,"character")
      if(any(xx!=uu)) { x@seqnames <- u@seqnames }
      x@build <- "hg19"
    } else { 
      stop("conversion to build37/hg19 failed, input had length ",length(x),"; output had length ",length(u)) 
    } 
  } else {
    if(ucsc.sanitizer(ucsc(x))!="hg19") { 
      warning("input object was not tagged as hg18/build36/hg38/build38 [@build], left unchanged") 
    } else {
      warning("object is already using hg19/build37, no change")
    }
  }
  return(x)
})



# @name convTo36
# @param x a ChipInfo object
#' @return convTo36: Returns a ChipInfo object with the build updated to hg18 coordinates
#' @family conversion
#' @rdname conv-methods
#' @export
setGeneric("convTo36", function(x) standardGeneric("convTo36"))


#' @family conversion
#' @aliases convTo36
#' @rdname conv-methods
#' @exportMethod convTo36
setMethod("convTo36", "ChipInfo", function(x) {
  if(ucsc.sanitizer(ucsc(x))=="hg19") {
    u <- conv.37.36(ranges=as(x,"GRanges"))
    if(length(u)==length(x)) { 
      x@ranges <- u@ranges
      all.eq <- TRUE
      if(all.eq) { all.eq <- length(seqlevels(x))==length(seqlevels(u)) }
      if(all.eq) { all.eq <- any(sort(seqlevels(x))==sort(seqlevels(u))) }
      if(!all.eq) {
        #warning("conversion altered the chromosomes"); 
        seqlevels(x) <- c(seqlevels(x),seqlevels(u)[!seqlevels(u) %in% seqlevels(x)])  #x@seqinfo <- u@seqinfo 
      }
      xx <- as(x@seqnames,"character")
      uu <- as(u@seqnames,"character")
      if(any(xx!=uu)) { x@seqnames <- u@seqnames }
      x@build <- "hg18"
    } else { 
      stop("conversion to build36/hg18 failed") 
    } 
  } else {
    if(ucsc.sanitizer(ucsc(x))!="hg18") { 
      warning("input object was not tagged as hg19/build37 [@build], left unchanged") 
    } else {
      warning("object is already using hg18/build36, no change")
    }
  }
  return(x)
})



# @name convTo38
# @param x a ChipInfo object
#' @return convTo38: Returns a ChipInfo object with the build updated to hg38 coordinates
#' @family conversion
#' @rdname conv-methods
#' @export
setGeneric("convTo38", function(x) standardGeneric("convTo38"))


#' @family conversion
#' @aliases convTo38
#' @rdname conv-methods
#' @exportMethod convTo38
setMethod("convTo38", "ChipInfo", function(x) {
  if(ucsc.sanitizer(ucsc(x))=="hg18") { stop("can't convert 36 to 38; must convert from 36 to 37, then convert from 37 to 38") }
  if(ucsc.sanitizer(ucsc(x))=="hg19") {
    u <- conv.37.38(ranges=as(x,"GRanges"))
    if(length(u)==length(x)) { 
      x@ranges <- u@ranges
      all.eq <- TRUE
      if(all.eq) { all.eq <- length(seqlevels(x))==length(seqlevels(u)) }
      if(all.eq) { all.eq <- any(sort(seqlevels(x))==sort(seqlevels(u))) }
      if(!all.eq) {
        #warning("conversion altered the chromosomes"); 
        seqlevels(x) <- c(seqlevels(x),seqlevels(u)[!seqlevels(u) %in% seqlevels(x)])  #x@seqinfo <- u@seqinfo 
      }
      xx <- as(x@seqnames,"character")
      uu <- as(u@seqnames,"character")
      if(any(xx!=uu)) { x@seqnames <- u@seqnames }
      x@build <- "hg38"
    } else { 
      stop("conversion to build38/hg38 failed") 
    } 
  } else {
    if(ucsc.sanitizer(ucsc(x))!="hg38") { 
      warning("input object was not tagged as hg19/build37 [@build], left unchanged") 
    } else {
      warning("object is already using hg18/build36, no change")
    }
  }
  return(x)
})


#' Display method for ChipInfo objects
#' 
#' show: Displays a preview of the object
# @name show
#' @param object a ChipInfo object
#' @return show: Displays a preview of the object
#' @rdname ChipInfo-methods
#' @exportMethod show
setMethod("show", "ChipInfo", 
          function(object) { showChipInfo(object,up.to=10,head.tail=5,show.strand=FALSE) } )


#' Print a ChipInfo object to the console
#' 
#' print: See 'show' as the behaviour is very similar and ... are just arguments of 'show'.
#' The key difference with 'print' instead of 'show' is that by default the parameter
#' 'up.to' is set to 50, so that any ChipInfo object (or subset) of less than or equal
#' to 50 rows will be displayed in its entirety, rather than just the top/bottom 5 rows. 
# @name print
# @param x a ChipInfo object
#' @param ... further arguments to showChipInfo()
#' @return print: Prints the object to terminal using 'showChipInfo()'.
#' @rdname ChipInfo-methods
#' @exportMethod print
setMethod("print", "ChipInfo", 
          function(x,...) { showChipInfo(x,...) } )


#' Constructor (wrapper) for ChipInfo annotation object
#' 
#' This class annotates a microarray SNP chip with data for each SNP including chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18, hg19 or hg38 (build 36/37/38) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' is basically a GRanges object, modified to always have columns for A1, A2 (alleles), 
#' rs.id, and a quality control flag. The default display is tidier than GRanges, it has
#' nice coersion to and frame data.frame and subsetting by chromosome using [[n]] has been
#' added, in addition to normal [i,j] indexing native to GRanges.
# @name ChipInfo
#' @param GRanges a GRanges object containing chromosome, start/end = position, and strand
#' information for the chip object to be created, also rownames should be used to code
#' the chip-ids for each SNP.
#' @param chr optional, alternative to using 'GRanges' to input SNP locations, enter here 
#' a vector of chromosome numbers/letters for each SNP. The recommended coding is: 
#' 1:22, X, Y, XY, MT
#' @param pos optional, vector of positions (integers), use in conjunction with 'chr' and
#'  'ids' as an alternative way to input SNP position information instead of GRanges.
#' @param ids optional, vector of SNP chip-ids, use in conjunction with 'chr' and
#'  'pos' as an alternative way to input SNP position information instead of GRanges.
#' @param chip character, name of the chip you are making this annotation for (only used
#' for labelling purposes)
#' @param build character, either "hg18" or "hg19". Will also accept build number, 36 or 37.
#' This indicates what coordinates the object is using, and will be taken into account by
#' conversion functions, and annotation lookup functions throughout this package.
#' @param rs.id 'rs' ids are standardized ids for SNPs, these usually differ from each chips'
#' own IDs for each snp. If you don't know these, or can't find them, they can be left blank,
#' but will render the functions 'rs.to.id()' and 'id.to.rs()' useless for this ChipInfo object.
#' @param chip.id chip ids are the chip-specific ids for SNPs, these usually differ between chips'
#' even for the same snp. If you don't know these, or can't find them, they can be left blank,
#' but will render the function 'chip.id()' useless for this ChipInfo object. The main purpose
#' of this parameter is for when the real chip ids are not valid R row/column name strings, and
#' by using this column, some functions can return the real chip ids instead of the sanitized 
#' version
#' @param A1 the first allele letter code for each SNP, e.g, usually "A","C","G", or "T", but
#' you can use any scheme you like. Can be left blank.
#' @param A2, as for A1, but for allele 2.
#' @param QCcode optional column to keep track of SNPs passing and failing QC. You can completely
#' ignore this column. It works based on integer codes, 0,1,2, you may wish to use simple 0 and 1,
#' for pass and fail respectively, or else 0 can be pass, and 1,2,... can indicate failure for 
#' different criteria. 0 will always be treated as a pass and anything else as a fail, so you
#' can code fails however you wish.
#' @export
ChipInfo <- function(GRanges=NULL, chr=NULL, pos=NULL, ids=NULL, chip="unknown chip", build="",
                     rs.id=NULL, chip.id=NULL, A1=NULL, A2=NULL, QCcode=NULL) {
  if(build!="") { build <- ucsc.sanitizer(build) }
  LL <- max(c(length(chr),length(GRanges)),na.rm=T)
  if(length(A1)!=LL | length(A2)!=LL) { A1 <- A2 <- rep(NA,times=LL) }
  if(length(rs.id)!=LL) { rs.id <- rep(NA,times=LL) } else { 
    if(any(duplicated(rs.id))) { rs.id <- add.trail(rs.id) } # appends letters to stop duplicates
  }
  if(length(chip.id)!=LL) { chip.id <- rep(NA,times=LL) } else { 
    if(any(duplicated(narm(chip.id)))) { stop("chip.id shouldn't contain duplicates") } # appends letters to stop duplicates
  }
  if(length(QCcode)!=LL) { QCcode <- rep(0,LL) }
  if(is.null(GRanges)) {
    GRanges <- makeGRanges(chr=chr,pos=pos,row.names=ids)
  } else {
    if(is(GRanges)[1]!="GRanges") { GRanges <- as(GRanges,"GRanges") }
  }
  df <- DataFrame(A1=A1,A2=A2,QCcode=QCcode,rs.id=rs.id,chip.id=chip.id)
  #print(build)
  return(new("ChipInfo", seqnames=GRanges@seqnames, ranges=GRanges@ranges,  strand=GRanges@strand,
            elementMetadata=df, seqinfo=GRanges@seqinfo,
             chip=chip, build=build))
}


#' Initialize (constructor) method for ChipInfo
#' 
#' Use the 'ChipInfo()' wrapper to construct ChipInfo objects from scratch
#  @name ChipInfo
#' @param .Object An object generated from the ChipInfo class prototype,
#'  see methods:initialize
#' @param ... Additional arguments to initialize. None recommended.
#' @rdname ChipInfo-class
#' @exportMethod initialize
setMethod("initialize", "ChipInfo",
              function(.Object, ...){
          		  callNextMethod(.Object, ...)
          	  })


#' As("ChipInfo", "GRanges")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("ChipInfo", "GRanges",
      function(from) { 
        #print(is(from)); print(from@seqnames)
        out <- GRanges(from@seqnames,ranges=from@ranges,strand=from@strand,
                       seqinfo=from@seqinfo,elementMetadata=from@elementMetadata,genome=ucsc(from))
        return(out)
      }
)

#' As("ChipInfo", "GRanges")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("ChipInfo", "RangedData",
      function(from) { 
        out <- as(as(from,"GRanges"),"RangedData")
        if("strand" %in% colnames(out)) { out <- out[,-which(colnames(out) %in% "strand")] }
        colnames(out) <- gsub("elementMetadata.","",colnames(out),fixed=TRUE)
        colnames(out) <- gsub("elementMetadata","",colnames(out),fixed=TRUE)
        return(out)
      }
)


#' As("ChipInfo", "GRanges")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("ChipInfo", "data.frame", function(from) { ranged.to.data.frame(as(from,"GRanges"),include.cols=TRUE) })


#' As("GRanges", "ChipInfo")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("GRanges", "ChipInfo", 
      function(from) { 
        bb <- genome(from)
        if(all(is.na(bb))) { build <- "" } else {
          if(length(unique(bb))==1) { build <- ucsc.sanitizer(bb[1]) } else { build <- "" }  
        }
        cN <- colnames(mcols(from))
        ii <- match(toupper("A1"), toupper(cN))
        if(!is.na(ii)) { a1 <- mcols(from)[,ii] } else { a1 <- NULL }
        ii <- match(toupper("A2"), toupper(cN))
        if(!is.na(ii)) { a2 <- mcols(from)[,ii] } else { a2 <- NULL }
        ii <- match(toupper("rs.id"), toupper(cN))
        if(!is.na(ii)) { rss <- mcols(from)[,ii] } else { rss <- NULL }
        ii <- match(toupper("chip.id"), toupper(cN))
        if(!is.na(ii)) { cid <- mcols(from)[,ii] } else { cid <- NULL }
        ii <- match(toupper("QCcode"), toupper(cN))
        if(!is.na(ii)) { qcc <- mcols(from)[,ii] } else { qcc <- NULL }
        ChipInfo("ChipInfo",GRanges=from,chip="unknown chip",build=build,rs.id=rss,chip.id=cid,A1=a1,A2=a2,QCcode=qcc)
      }
)

#' As("RangedData", "ChipInfo")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("RangedData", "ChipInfo", function(from) { as(as(from,"GRanges"),"ChipInfo") } )

#' As("data.frame", "ChipInfo")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("data.frame", "ChipInfo", 
      function(from) { 
        rr <- df.to.GRanges(from,chr="seqnames") 
        return(as(as(rr,"GRanges"),"ChipInfo"))
      } 
)



# improved coersion functions for data.frames to RangedData/GRanges

#' As("data.frame", "RangedData")
#'
#' @name as
#' @export
setAs("data.frame", "RangedData", function(from) { return(df.to.ranged(from,GRanges=FALSE))  } )

#' As("data.frame", "GRanges")
#'
#' @name as
#' @export
setAs("data.frame", "GRanges", function(from) { return(df.to.ranged(from,GRanges=TRUE))  } )

#' As("RangedData", "data.frame")
#'
#' Note that for automatic conversion of a data.frame to RangedData/GRanges, a column named 'chr' 
#' or 'seqnames' in the data.frame is expected/required to make the conversion effectively. 
#' Otherwise use 'ranged.to.data.frame()'
#' @name as
#' @export
setAs("RangedData", "data.frame", function(from) { return(ranged.to.data.frame(from))  } )

#' As("GRanges", "data.frame")
#'
#' @name as
#' @export
setAs("GRanges", "data.frame", function(from) { return(ranged.to.data.frame(from))  } )



# No roxygen required i think?
setValidity("ChipInfo",
            function(object) {
              if (!is.character(chip(object)) || length(chip(object)) != 1 || is.na(chip(object))) {
                return("'chip' slot must be a single string") 
              }
              if (!is.character(ucsc(object)) || length(ucsc(object)) != 1 || is.na(ucsc(object))) {
                return("'build' slot must be a single string") 
              } else {
               # requireNamespace(humarray)
                if(!paste(ucsc(object)) %in% c("",paste(as.vector(ucsc.sanitizer(show.valid=T)[,1])))) {
                  return("'build' must be a string, 36/37/38 or hg18/hg19/hg38") 
                }
              }
            }
)

#' Display a ChipInfo object
#' 
#' Returns a preview of a ChipInfo object to the console. This
#' is similar to a GRanges preview, but the seqlevels are hidden, the UCSC
#' build and chip name are displayed, start and end are merged to the virtual
#' label 'position' (as it's assume we are dealing with SNPs, not ranges), the strand
#' by default is hidden, and the integer codes for pass/fail in QCcodes() are 
#' displayed as 'pass' or 'fail', even though this is not how they are represented internally.
#' This is called by the default 'show' method for ChipInfo objects. 
#' @param x a ChipInfo object
#' @param margin margin for display, usually ""
#' @param with.classinfo logical, whether to display class information
#' @param print.seqlengths logical, whether to display sequence lengths below
#' the main output listing (e.g, chromsomes). Usually tidier when this is FALSE.
#' @param ... hidden arguments including: 'head.tail'; number of SNPs to display 
#' at start/end (only the head and tail are shown as these objects are generally
#' very large with >100K SNPs); 'up.to'; only SNPs at the start and end are generally
#' displayed, however this parameter specifies that when there are <= 'up.to' SNPs,
#' then all SNPs will be displayed; 'show.strand'; logical, by default the strand is 
#' hidden, particularly given that the strand can vary between different datasets 
#' of the same chip. Setting to TRUE will display the strand.
#' @return print compact preview of the object to the standard output (terminal)
#' @seealso \code{\link{ChipInfo}}
#' @export
showChipInfo <- function (x, margin = "", with.classinfo = FALSE, print.seqlengths = FALSE,...) 
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  qc <- QCcode(x)
  bb <- ucsc(x)
  if(bb=="") { bb <- "unknown" }
  if(length(qc)==lx) { QC <- rep("pass",lx); QC[qc>0] <- paste0("fail",QC[qc>0]) ; x$QCcode <- QC }
  cat("ChipInfo for ",chip(x)," with ", lx, " ", ifelse(lx == 1L, "SNP", 
                   "SNPs")," using ",bb," coordinates",":\n", sep = "")
  out <- makePrettyMatrixForCompactPrinting2(x, .makeNakedMatFromChipInfo,...)
  if (nrow(out) != 0L) 
    rownames(out) <- paste0(margin, rownames(out))
  print(out, quote = FALSE, right = TRUE)
}


#internal
extraColumnSlots2 <- function(x) {
  sapply(extraColumnSlotNames2(x), slot, object = x, simplify = FALSE)
}

#internal
setGeneric("extraColumnSlotNames2",
           function(x) standardGeneric("extraColumnSlotNames2"))

setMethod("extraColumnSlotNames2", "ANY", function(x) character())


#internal
.makeNakedMatFromChipInfo <- function (x,show.strand=TRUE) 
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  if(!show.strand) {
    ans <- cbind(seqnames = as.character(seqnames(x)), ranges = showAsCell(ranges(x)))
  } else { 
    ans <- cbind(seqnames = as.character(seqnames(x)), ranges = showAsCell(ranges(x)),strand = as.character(strand(x)))
  }
  extraColumnNames <- extraColumnSlotNames2(x)
  if (length(extraColumnNames) > 0L) {
    ans <- do.call(cbind, c(list(ans), lapply(extraColumnSlots2(x), 
                                              showAsCell)))
  }
  if (nc > 0L) {
    df <- mcols(x)
    if(tail(colnames(df),1)=="chip.id") { df <- df[,-ncol(df)] }
    tmp <- do.call(data.frame, c(lapply(df, showAsCell), 
                                 list(check.names = FALSE)))  ### hide chip.id here!
    ans <- cbind(ans, `|` = rep.int("|", lx), as.matrix(tmp))
    if(all(colnames(ans)[1:2]==c("seqnames","ranges"))) { colnames(ans)[1:2] <- c("chr","pos") }
  }
  ans
}


#internal
makePrettyMatrixForCompactPrinting2 <- function (x, makeNakedMat.FUN,head.tail=6,up.to=50, show.strand=TRUE) 
{
  lx <- NROW(x)
  if(lx <= up.to) { head.tail <- up.to }
  nhead <- head.tail
  ntail <- head.tail
  if (lx < (nhead + ntail + 1L)) {
    ans <- makeNakedMat.FUN(x,show.strand=show.strand)
    ans_rownames <- .rownames3(names(x), lx)
  }
  else {
    top_idx <- 1:nhead
    if (nhead == 0) 
      top_idx <- 0
    bottom_idx = (lx - ntail + 1L):lx
    if (ntail == 0) 
      bottom_idx <- 0
    ans_top <- makeNakedMat.FUN(x[top_idx, , drop = FALSE],show.strand=show.strand)
    ans_bottom <- makeNakedMat.FUN(x[bottom_idx, , drop = FALSE],show.strand=show.strand)
    ans <- rbind(ans_top, matrix(rep.int("...", ncol(ans_top)), 
                                 nrow = 1L), ans_bottom)
    ans_rownames <- .rownames3(names(x), lx, top_idx, bottom_idx)
  }
  rownames(ans) <- format(ans_rownames, justify = "right")
  ans
}

#internal
.rownames3 <- function (names = NULL, len = NULL, tindex = NULL, bindex = NULL) 
{
  if (is.null(tindex) && is.null(bindex)) {
    if (len == 0L) 
      character(0)
    else if (is.null(names)) 
      paste0("[", seq_len(len), "]")
    else names
  }
  else {
    if (!is.null(names)) {
      c(names[tindex], "...", names[bindex])
    }
    else {
      s1 <- paste0("[", tindex, "]")
      s2 <- paste0("[", bindex, "]")
      if (all(tindex == 0)) 
        s1 <- character(0)
      if (all(bindex == 0)) 
        s2 <- character(0)
      c(s1, "...", s2)
    }
  }
}



#' Select chromosome subset for ranged objects
#' 
#' Returns the object filtered for specific chromosomes for a ranged object
#' @param object a ChipInfo, GRanges or RangedData object
#' @param chr vector, string or numeric of which chromosome(s) to select
#' @return vector of chromosome values for each range/SNP
#' @rdname chrSel-methods
#' @export
setGeneric("chrSel", function(object,chr) standardGeneric("chrSel"))


#' Select chromosome subset for RangedData objects
#' 
#' Returns the object filtered for specific chromosomes for a RangedData object
#' @rdname chrSel-methods
#' @exportMethod chrSel
setMethod("chrSel", "RangedData", function(object,chr) {
  return(humarray::chrSelect(object,chr))
})


#' Select chromosome subset for GRanges objects
#' 
#' Returns the object filtered for specific chromosomes for a GRanges object
#' @rdname chrSel-methods
#' @exportMethod chrSel
setMethod("chrSel", "GRanges", function(object,chr) {
  return(humarray::chrSelect(object,chr))
})

#' Select chromosome subset for ChipInfo objects
#' 
#' Returns the object filtered for specific chromosomes for a GRanges object
#' @rdname chrSel-methods
#' @exportMethod chrSel
setMethod("chrSel", "ChipInfo", function(object,chr) {
  return(humarray::chrSelect(object,chr))
})



#' Chromosome method for RangedData objects
#' 
#' Return the list of chromosome values from a RangedData object
#' @param object RangedData object
#' @return vector of chromosome values for each range/SNP
#' @rdname chrm-methods
#' @export
setGeneric("chrm",function(object) standardGeneric("chrm"))

#' @rdname chrm-methods
#' @importMethodsFrom genoset  chr  
#' @exportMethod chrm
setMethod("chrm", "RangedData", function(object) {
  return(chr2(object))
})

#' @rdname chrm-methods
#' @importMethodsFrom genoset  chr  
#' @importFrom genoset  chr  
#' @exportMethod chrm
setMethod("chrm", "GRanges", function(object) {
  return(genoset::chr(object))
})

#' @rdname chrm-methods
#' @importMethodsFrom genoset  chr  
#' @importFrom genoset  chr  
#' @exportMethod chrm
setMethod("chrm", "ChipInfo", function(object) {
  return(genoset::chr(object))
})



#' @importMethodsFrom "genoset"  toGenomeOrder  
#' @exportMethod toGenomeOrder
setMethod("toGenomeOrder", "RangedData", function(ds) {
  return(TGORD(ds))
})



#' @importMethodsFrom "genoset"   chrIndices  
#' @exportMethod chrIndices
setMethod("chrIndices", "RangedData", function(object) {
  return(chrIndices2(object))
})



#' @importMethodsFrom "genoset"   chrInfo 
#' @exportMethod chrInfo
setMethod("chrInfo", "RangedData", function(object) {
  return(chrInfo2(object))
})




#' @importMethodsFrom "genoset"  chrNames  
#' @exportMethod chrNames
setMethod("chrNames", "RangedData", function(object) {
  return(chrNames2(object))
})


#' Plot method for GRanges objects
#' 
#' See plotRanges()
# @name plot
#' @param x a GRanges or RangedData object
#' @param y not used for plotRanges
#' @param ... further arguments, see plotRanges()
#' @rdname plot-methods
#' @aliases GRanges GRanges-method
#' @seealso \code{\link{plotRanges}}
#' @exportMethod plot
setMethod("plot", "GRanges", function(x,y,...) {
  plotRanges(ranged=x,...)
})


#' Plot method for RangedData objects
#' 
# @name plot
#' @rdname plot-methods
#' @aliases RangedData RangedData-method
#' @exportMethod plot
setMethod("plot", "RangedData", function(x,y,...) {
  plotRanges(ranged=x,...)
})

##################################
## General Annotation Functions ##
##################################


#' Download GWAS hits from t1dbase.org
#' 
#' Retrieve human disease top GWAS hits from t1dbase in build hg19 coords (37).
#' 28 Diseases currently available
#' @param disease integer (1-28), or character (abbreviation), or full name of one of the listed
#' diseases. A full list of options can be obtained by setting show.codes=TRUE.
#' @param snps.only logical, default is just to return a list of rs-ids. Setting FALSE gives a table
#' @param show.codes logical, if set to TRUE, instead of looking up t1dbase, will simply return
#' a table of available diseases with their index numbers and abbreviations.
#' @return A character vector of SNP rs-ids
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references PMID: 20937630
#' @examples
#' get.immunobase.snps(disease="CEL") # get SNP ids for celiac disease
#' get.immunobase.snps(disease="AS") # get SNP ids for Ankylosing Spondylitis in build-37/hg19
#' get.immunobase.snps(show.codes=TRUE) # show codes/diseases available to download
#' get.immunobase.snps(disease=27) # get SNP ids for Alopecia Areata
#' get.immunobase.snps("Vitiligo")
get.immunobase.snps <- function(disease="T1D",snps.only=TRUE,show.codes=FALSE) {
  disease.codes <- c("Type 1 Diabetes", "Crohns Disease","Rheumatoid Arthritis",
                     "Systemic Scleroderma",  "Ulcerative Colitis","Inflammatory Bowel Disease",  "Multiple Sclerosis",
                     "Bipolar Disorder",  "Diabetes Mellitus",  "Coronary Artery Disease",  "Hypertension",  "Celiac Disease",
                     "Systemic Lupus Erythematosus",  "Ankylosing Spondylitis",  "Type 2 Diabetes",  "Sjogren Syndrome",
                     "Graves' Disease",  "Juvenile Rheumatoid Arthritis",  "Vitiligo",  "Primary Biliary Cirrhosis",
                     "Psoriasis",  "Idiopathic Membranous Nephropathy",  "Immunoglobulin A Deficiency",
                     "Autoimmune Thyroid Disease",  "Juvenile Idiopathic Arthritis",  "Narcolepsy",  "Alopecia Areata",
                     "Alzheimer's Disease")
  abbr <- c("T1D","CD","RA","SCL","UC","IBD","MS","BD","DM","CAD","HYP",
            "CEL","SLE","AS","T2D","SS","GD","JRA","VIT",
            "PBC","PSO","IMN","IGA","ATD","JIA","NAR","AA","AD")
  code.table <- cbind(Abbreviation=abbr,FullNames=disease.codes)
  if(show.codes) { cat("values for the 'disease' parameter can be specified by the following index numbers or abbreviations:\n")
                   print(code.table,quote=F) ; return() }
  disease <- disease[1]
  if(toupper(disease) %in% abbr) { 
    disN <- match(toupper(disease),toupper(abbr)) 
  } else {
    if(toupper(disease) %in% toupper(disease.codes)) {
      disN <- match(toupper(disease),toupper(disease.codes))
    } else {
      if(as.numeric(disease) %in% 1:length(disease.codes)) { 
        disN <- as.numeric(disease)
      } else {
        stop("Invalid input for 'disease', use show.codes=TRUE to see list of codes/abbreviations")
      }
    }
  }
  #if(is.null(build)) { build <- getOption("ucsc") }
  build <- "hg19" # ucsc.sanitizer(build)
  if(!build %in% c("hg18","hg19")) { stop("only hg18 and hg19 are supported for this function") }
  filenm <- cat.path(dir=getwd(),pref=tolower(abbr[disN]),"hits",suf=build,ext="tab")
  cat("attempting to download",abbr[disN],"hits from t1dbase\n")
  url36 <- paste("http://www.immunobase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=",disN,"&type=assoc&build=GRCh36",sep="")
  url37 <- paste("http://www.immunobase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=",disN,"&type=assoc&build=GRCh37",sep="")
#  url36 <- paste("http://www.t1dbase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=",disN,"&type=assoc&build=GRCh36",sep="")
#  url37 <- paste("http://www.t1dbase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=",disN,"&type=assoc&build=GRCh37",sep="")
  urL <- switch(build, hg18=url36,  hg19=url37)
  success <- T
  success <- tryCatch(download.file(urL ,filenm ,quiet=T),error=function(e) { F } )
  #prv(filenm,urL,success)
  if(!is.logical(success)) { success <- T }
  if(success) {
    t1dh <- readLines(filenm)
    firsts <- substr(t1dh,1,2)
    t1dh <- t1dh[firsts!="##"]
    len.lst <- strsplit(t1dh,"\t")
    rsids <- sapply(len.lst,"[",3)
    if(substr(rsids[1],1,2)!="rs") { rsids <- rsids[-1] }
    if(length(rsids)<1) { 
      cat("download successful but the list of hits for",disease.codes[disN],"was empty\n") 
    } else {
      cat("download successful for",disease.codes[disN],"\n")
    }
    if(snps.only) {
      return(unique(rsids))
    } else {
      return(read.delim(filenm,comment.char="#"))
    } 
  } else {
    stop("couldn't reach t1dbase website at: ",urL)
  }
}




#' Retreive GO terms from biomart for a given gene list
#' 
#' Gene-ontology terms (GO-terms) are commonly used for testing for simple functional
#' enrichment for pathways, etc. This function can retrieve biological function, 
#' cellular component, or molecular description, depending on the parameters chosen.
#' @param gene.list a list of gene, use HGNC names, like COMT, HLA-C, CTLA4, etc.
#' @param bio logical, whether to return biological process GO terms
#' @param cel logical, whether to return cellular component GO terms
#' @param mol logical, whether to return molecular function GO terms
#' @param host.txt character, the argument to pass to biomaRt::useMart(). Default is 
#' 'may2009.archive.ensembl.org', but more recently the recommended link is 'www.ensembl.org'
#' @return data.frame containing the gene name in the first column, chromosome in the
#' second column, and the GO terms in the third column, where one gene has multiple
#' GO terms, this will produce multiple rows, so there will usually be more rows
#' than genes entered. The data.frame can have 3,4 or 5 columns depending on
#' how many GO terms are selected.
#' @export
#' @examples
#' get.GO.for.genes(c("CTLA4","PTPN2","PTPN22")) # biological terms (default)
#' get.GO.for.genes(c("CTLA4","PTPN2","PTPN22"),cel=TRUE) # add cellular GO terms
get.GO.for.genes <- function(gene.list,bio=T,cel=F,mol=F,host.txt="may2009.archive.ensembl.org") {
 # must.use.package(c("biomaRt","genoset","gage"),T)
  mart.txt <- "ENSEMBL_MART_ENSEMBL"
  ens <- biomaRt::useMart(mart.txt,
                 dataset="hsapiens_gene_ensembl",
                 host= host.txt,
                 path="/biomart/martservice",
                 archive=FALSE)
  ens <- biomaRt::useDataset("hsapiens_gene_ensembl",mart=ens)
  #egSymb <- reader("~/github/iChip/egSymb.rda")
  egSymb <- humarray::egSymb
  base.attr <- c("hgnc_symbol", "chromosome_name")
  if(bio) { base.attr <- c(base.attr,"go_biological_process_description") }
  if(cel) { base.attr <- c(base.attr,"go_cellular_component_description") }
  if(mol) { base.attr <- c(base.attr,"go_molecular_function_description") }
  #  dat <- getBM(attributes = c("hgnc_symbol", "chromosome_name",
  #                              "start_position", "end_position", "band"), filters = "hgnc_symbol",
  #               values = egSymb[,2], mart = ens)
  results <- biomaRt::getBM(attributes = base.attr, filters = "hgnc_symbol",
                   values = c(gene.list), mart = ens)
  return(results)
}




#' Convert ensembl ids to HGNC gene ids 
#' 
#' Retrieve the gene IDs (HGNC) corresponding to a list of ensembl gene ids.
#' Note that this will not find all IDs found on ensembl.org, as it uses bioMart which
#' seems to be incomplete, but this only pertains to a small minority of genes, so this
#' function should have general utility for most applications. This is of course the case
#' at the time of writing - bioMart is likely to be updated at some point.
#' @param ens character, a list of ensembl gene ids, of the form ENSG00xxxxxxxxx
#' @param ... further arguments to get.gene.annot()
#' @param dir character, 'dir' is the location to download gene and cytoband information; if
#' left as NULL, depending on the value of getOption("save.annot.in.current"), the annotation
#' will either be saved in the working directory to speed-up subsequent lookups, or deleted 
#' after use.
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param name.dups logical, if TRUE then duplicates will have a suffix appended to force the
#' list to be unique (e.g, so it would be usable as rownames, or in a lookup table). Otherwise
#' duplicate entries will just appear in the list multiple times
#' @param name.missing logical, if TRUE then missing values will be named as MISSING_n (n=1
#'  to # of missing), ensuring a valid unique name if the results are to be used as rownames,
#' etc. If FALSE then these will be left as NA. 
#' @return Returns a vector of HGNC gene ids corresponding to the 'ens' ensembl ids entered,
#' any ids not found will be returned as MISSING_n (n=1 to # of missing), if name.missing=TRUE.
#' If name.missing is FALSE then missing will be set to NA. Similarly with 'name.dups', if
#' duplicates are found and name.dups is true, each will be appended with suffix _n; else
#' their names will be left as is.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{GENE.to.ENS}}, \code{\link{rs.to.id}}, \code{\link{id.to.rs}}; eg2sym, sym2eg from package 'gage'
#' @examples
#' \donttest{
#' setwd(tempdir())
#' ENS.ids <- c("ENSG00000183214", "ENSG00000163599", "ENSG00000175354", "ENSG00000134460")
#' ENS.to.GENE(ENS.ids)
#' gene.ids <- c("HLA-B","IFIH1","fake_gene!","FUT2")
#' ENS.to.GENE(GENE.to.ENS(gene.ids)) # lookup fails for the fake id, gives warning
#' }
ENS.to.GENE <- function(ens,dir=NULL,build=NULL,name.dups=FALSE,name.missing=TRUE,...) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
#  must.use.package(c("biomaRt","genoset","gage"),T)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.gene.annot(...,dir=dir,bioC=FALSE,ens.id=TRUE,GRanges=FALSE)
  ### now have the gene data with the ensembl ids ##
  #print(head(ga))
  indx <- match(ens,ga$ens.id)
  missin <- length(which(is.na(indx))); valid <- length(indx)-missin
  if(valid<1) { warning("did not find any ENSEMBL ids from 'ens' in the bioMart human gene reference"); return(NULL) }
  if(missin>0) { warning(out.of(missin,(valid+missin))," of 'ens' did not match any ENSEMBL ids in the bioMart human gene reference") }
  outData <- ga$gene[indx]
  if(name.missing & any(is.na(outData))) {
    outData[is.na(outData)] <- paste("MISSING",pad.left(1:length(which(is.na(outData))),"0"),sep="_")
  }
  if(any(duplicated(outData))) { 
    if(name.dups) { 
      cnt <- 2
      while(any(duplicated(outData))) { 
        if(cnt==2) {
          outData[duplicated(outData)] <- paste(outData[duplicated(outData)],cnt,sep="_")
        } else {
          outData[duplicated(outData)] <- gsub(paste("_",cnt-1,sep=""),paste("_",cnt,sep=""),outData[duplicated(outData)])
        }
        cnt <- cnt + 1
      }
    } else { 
      warning("duplicated gene names produced, select 'name.dups=TRUE' to append numbers to make these unique")
    }
  }
  return(outData)
}


#' Convert gene ids to ensembl ids
#' 
#' Retrieve the ensembl IDs corresponding to a list of common gene names (HGNC format).
#' @param genes character, gene labels, e.g, "APOE"
#' @param ... further arguments to get.gene.annot()
#' @param dir character, 'dir' is the location to download gene and cytoband information; if
#' left as NULL, depending on the value of getOption("save.annot.in.current"), the annotation
#' will either be saved in the working directory to speed-up subsequent lookups, or deleted 
#' after use.
#' @return Returns a vector of HGNC gene ids corresponding to the 'ens' ensembl ids entered,
#' any ids not found will be returned as NA.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{GENE.to.ENS}}, \code{\link{rs.to.id}}, \code{\link{id.to.rs}}; eg2sym, sym2eg from package 'gage'
#' @examples
#' \donttest{
#' setwd(tempdir())
#' gene.ids <- c("MYC","PTPN2","IL2RA","APOE")
#' GENE.to.ENS(gene.ids)
#' }
GENE.to.ENS <- function(genes,dir=NULL,...) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.gene.annot(...,dir=dir,bioC=FALSE,ens.id=TRUE,GRanges=FALSE)
  ### now have the gene data with the ensembl ids ##
  indx <- match(genes,ga$gene)
  missin <- length(which(is.na(indx))); valid <- length(indx)-missin
  if(valid<1) { warning("did not find any gene ids from 'genes' in the bioMart human gene reference"); return(NULL) }
  if(missin>0) { warning("at least one of 'genes' did not match any gene ids in the bioMart human gene reference") }
  outData <- ga$ens.id[indx]
  return(outData)
}



#' Retrieve the 'n' closest GENE labels or positions near specified locus
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' @param pos integer, genomic position, should be between 1 and the length of the chromosome 'chr'
#' @param n integer, the number of nearest GENEs to seek, if there aren't enough in the annotation
#' then NAs will fill the gaps to force the return value length to equal 'n'
#' @param side character, can be 'either', 'left' or 'right' and specifies which side of the 'pos'
#' to look for nearest genes (where left is decreasing genomic position and right is increasing)
#' @param ids logical, if TRUE will return GENE labels, 
#' or if FALSE will return the chromosome positions of the genes
#' @param limit integer, a limit on the maximum distance from the position 'pos' can be specified
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param ga RangedData object, e.g, result of get.gene.annot(); gene annotation to save download
#' time if repeatedly calling this function
#' @export
#' @seealso \code{\link{expand.nsnp}}, \code{\link{nearest.snp}}, \code{\link{get.gene.annot}}
#' @return Set of GENE ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr'.
#' If the number of gemes on the chromosome or the bounds of the 'side' and 'limit' parameters
#' restrict the number returned to less than 'n' then the return value will be padded with NAs.
#' @examples
#' \donttest{
#' nearest.gene(1,159000000,n=10) # return ids
#' nearest.gene(1,159000000,n=10,build=37)
#' nearest.gene(1,159000000,n=10,build=36,ids=FALSE) # return positions
#' nearest.gene(1,159000000,n=10,build=37,ids=FALSE)
#' nearest.gene(6,25000000,n=10,build=37,ids=FALSE,side="left")  # only genes to the left of the locus
#' nearest.gene(6,25000000,n=10,build=37,ids=FALSE,side="right") # only genes to the right of the locus
#' }
nearest.gene <- function(chr, pos, n=1, side=c("either","left","right"),ids=TRUE,limit=NULL,build=NULL, ga=NULL) { 
  # ids - whether to return ichip SNP ids or positions
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(pos)>1) { warning("pos should be length 1, using only first entry"); pos <- pos[1] }
  if(is.null(build)) { build <- getOption("ucsc") }
  chrom <- paste(chr)
  build <- ucsc.sanitizer(build)
  if(is(get.gene.annot)[1]=="RangedData") { 
    if(!"gene" %in% colnames(ga)) { ga <- NULL }
  } else { ga <- NULL }
  if(is.null(ga)) {  
    ga <- get.gene.annot(build=build,GRanges=F) 
    if(!exists("ga")) { stop("couldn't find gene annotation") }  ## load object: ga [gene database]
    ga <- ga[ga$gene!="",]
  }
  side <- tolower(side[1]); 
  if(!side %in% c("either","left","right")) {
    side <- "either"; warning("invalid side argument, defaulting to 'either'") }
  if(!is.null(limit)) { if(!is.numeric(limit)) { limit <- NULL; warning("invalid limit argument, defaulting to NULL") } }
  all.chr <- paste(chr2(ga))
  all.st <- start(ga)[all.chr %in% chrom]
  all.en <- end(ga)[all.chr %in% chrom]
  #prv(all.st,all.en)
  if(length(all.st)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  difzS <- pos-all.st
  difzE <- pos-all.en
  all.true <- difzS==difzS
  if(is.null(limit)) { lfilt <- all.true } else { lfilt <- abs(difzS)<=limit | abs(difzE)<=limit }
  within <- difzS>0 & difzE<0
  if(side=="left") { filt <- ( difzE>0 & lfilt ) | within }
  if(side=="right") { filt <- ( difzS<0 & lfilt ) | within }
  if(side=="either") { filt <- ( all.true & lfilt ) | within }
  #print(length(which(filt)))
  tab <- rbind(abs(difzS),abs(difzE))
  minz <- apply(tab,2,min,na.rm=T)
  Difz <- abs(minz[filt])
  if(length(Difz)<n)  { warning("fewer than ",n," genes found for 'chr' specified (within 'limit'), NAs returned") }
  indx <- order(Difz)[1:n]
  # prv(minz,Difz,filt,indx)
  subi <- ga[["gene"]][all.chr %in% chrom][filt]
  #prv(ga,subi,all.chr,chrom)
  if(ids) {
    out <- ga[["gene"]][all.chr %in% chrom][filt][indx]
  } else {
    out <- start(ga)[(all.chr %in% chrom)][filt][indx]
  }
  return(out)
}




#' Plot genes to annotate figures with genomic axes
#' 
#' Quite often it is helpful to visualize genomic locations in the context of Genes
#' in the same region. This function makes it simple to overlay genes on plots
#' where the x-axis is chromosomal location.
#' @param chr chromosome number/name that the plot-range lies on
#' @param scl character, the scale that the x axis uses, ie, "b","kb","mb", or "gb", meaning
#' base-pairs, kilobases, megabases or gigabase-pairs.
#' @param y.ofs numeric, y-axis-offset, depending on what units are on your y-axis,
#' you may prefer to specify an offset so that the gene annotation is drawn at an appropriate
#' level on the vertical axis, this value should be the centre of annotation
#' @param width depending on the range of your y-axis, you might want to expand or reduce
#' the vertical width of the gene annotation (in normal graph units), default
#' when width=NA is 10 percent of the y-axis size.
#' @param txt logical, TRUE to include the names of genes on top of their representation
#' on the plot, or if FALSE, genes are drawn without labels.
#' @param chr.pos.offset if for some reason zero on the x-axis is not equal to 'zero' on
#' the chromsome, then this offset can correct the offset. For instance if you were using
#' a graph of the whole genome and you were plotting genes on chromosome 10, you would
#' set this offset to the combined lengths of chromosomes 1-9 to get the start point
#' in the correct place.
#' @param gs GRanges or RangedData object, this is annotation for the location of genes.
#' This will be retrieved using get.gene.annot() if 'gs' is NULL. THere may be several reasons
#' for passing an object directly to 'gs'; firstly speed, if making many calls then you won't
#' need to load the annotation every time; secondly, if you want to use an alternative annotation
#' you can create your own so long as it is a GRanges/RangedData object and contains a column
#' called 'gene' (which doesn't strictly have to contain gene labels, it could be any feature
#' you require, eg., transcript names, etc).
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param box.col genes are drawn as boxes, this sets the colour of the boxes
#' @param txt.col this sets the colour of the label text (Gene names)
#' @param join.col for exons, or multipart genes, joins are made between the sections with
#' a central line, this sets the colour of that line.
#' @param ... further arguments to 'rect', the graphics function used to plot the 'genes'.
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome and start and end positions, other information depends on 
#' specific parameters documented above
#' @export
#' @examples
#' # EXAMPLE PLOT OF SOME SIMULATED SNPS on chr21-p11.1 #
#' # do we need to require(GenomicRanges)? #
#' setwd(tempdir())
#' loc <- c(9.9,10.2)
#' Band(chr=21,pos=loc*10^6)
#' rr <- in.window(rranges(50000),chr=21,pos=loc,unit="mb") # make some random MHC ranges
#' # create some SNPs and plot
#' rr3 <- rr; end(rr3) <- start(rr3) 
#' rownames(rr3) <- paste0("rs",sample(10^6,nrow(rr3)))
#' plotRanges(rr3,col="blue",scl="mb",xlim=loc,xlab="Chr21 position (Mb)",ylab="")
#' # NOW add UCSC hg18 GENE annotation to the plot #
#' \donttest{ plotGeneAnnot(chr=21,pos=c(9.95,10.1),scl="mb",y.ofs=1,build=36) }
plotGeneAnnot <- function(chr=1, scl=c("b","kb","mb","gb"), y.ofs=0, width=NA, txt=T, chr.pos.offset=0,
                            gs=NULL, build=NULL, dir=NULL, box.col="green", txt.col="black", join.col="red", ...)
{
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  dir <- validate.dir.for(dir,"ano")
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(!is(gs)[1] %in% c("RangedData","GRanges")) { gs <- get.gene.annot(dir=dir,build=build,GRanges=FALSE) }
  if(is(gs)[1]=="GRanges") { gs <- as(gs,"RangedData") }
  if(!"gene" %in% colnames(gs)) { warning("didn't find 'gene' column in annotation") ; return(NULL) }
  Col <- c("green", "darkgreen")

  # get set of genes in range of the graph section + remove duplicate genes/exons
  plot.area <- plot.get.area();
  if(all(is.null(plot.area))) { pos <- c(1,Inf) } 
  x.lim <- plot.area$xlim;
  y.lim <- plot.area$ylim;
  pos <- x.lim*make.divisor(scl) # overrides pos
  if(any(is.na(pos))) { pos <- c(1,Inf) } 
  rng.genez <- in.window(gs,chr,pos,full.overlap=F, rmv.dup=T,unit=scl)
  if(nrow(rng.genez)<1) { warning("no genes found in range") ; return(NULL) }
  old.cc <- 1 ; tp <- 2
  # set vertical alignments for annotation
  old.auto <- F
  if(y.ofs==0) { y.ofs <- min(y.lim) + .1*(max(y.lim)-min(y.lim)) }
  if(is.na(width)) { width <- .1*(max(y.lim)-min(y.lim)) } 
  y.cent <- y.ofs
  y.bot <- y.ofs-(width/2)
  y.top <- y.ofs+(width/2)
  # text alignment
  tps <- y.bot + c(.18,.29,.46,.64,.82)[c(1,3,5)]*width
  # x position with scaling (e.g, Mb units = 10^6)
  #unit <- tolower(unit[1]) ; mult <- switch(unit,b=0,kb=3,mb=6,gb=9); pos <- pos*10^mult
  x.scl <- make.divisor(scl)
  cnrlo <- (start(rng.genez)/x.scl)+chr.pos.offset
  cnrhi <- (end(rng.genez)/x.scl)+chr.pos.offset
  gnnm <- (rng.genez$gene)
  n.genes <- length(cnrlo)
  txt.cex <- .75; if(n.genes>10) { txt.cex <- .5 } ; if(n.genes>100) { txt.cex <- .35 } # more genes = smaller labels
  #print(n.genes)
  for (cc in 1:n.genes) {
    # draw rectangle and label for each gene
    #prv(cnrlo[cc],y.top,cnrhi[cc], y.bot)
    rect(cnrlo[cc],y.top,cnrhi[cc], y.bot,border=box.col,...)
  }
  for (cc in 1:n.genes) {
    if (gnnm[old.cc]==gnnm[cc] & cc!=1)
    {
      link <- c(cnrhi[old.cc],cnrlo[cc])
      if(link[1]<link[2]) {
        lines(link,y=rep(y.cent,2),lwd=1,col=join.col,lty="dotted")
      }
    } else {
      if(txt) {
        if(cnrlo[cc] < min(pos/x.scl)) { 
          txt.x <- mean(c(min(pos/x.scl),min(cnrhi[cc],max(pos/x.scl))),na.rm=T)  
        } else { txt.x <- cnrlo[cc] }
        text(txt.x,tps[tp],gnnm[cc],col=txt.col,cex=txt.cex,pos=4,las=2,offset=0)
      }
    }
    old.cc <- cc  
    tp <- tp+1; if(tp==4) {tp <- 1}; #if(n.genes <=5) { tp <- 2 }
  }
  return(rng.genez)
}



#' Retrieve locations of Immunoglobin regions across the genome
#' 
#' Returns the locations of immunoglobin regions in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' For instance, for CNV research, these regions are known to be highly structurally complex
#' and can lead to false positive CNV-calls, so are often excluded.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be immunoglobin chromosome, start and end positions.
#' @examples
#' get.immunog.locs()
#' get.immunog.locs(bioC=FALSE)
#' get.immunog.locs(text=TRUE,build=37)
get.immunog.locs <- function(build=NULL,bioC=TRUE,text=FALSE,GRanges=TRUE) {
  # http://atlasgeneticsoncology.org/Genes/GC_IGH.html
  # ^ has list of CELL line breakpoints... for future filtering?
  nchr <- 22
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(build[1]=="hg18") {
    # hg18
    chr <- c(22,14,2,14)
    stz <- c(20715572,105065301,88937989,21159897)
    enz <- c(21595082,106352275,89411302,22090937)
  } else {
    if(build[1]=="hg38") {
      #hg38
      chr <- c(22,14,2,14)
      stz <- c(20668232,105594256,88857361,21621904)
      enz <- c(22922910,107281230,89917421,22552944)
    } else {
      # hg19
      chr <- c(22,14,2,14)
      stz <- c(22385572,105994256,89156874,22090057)
      enz <- c(23265082,107281230,89630187,23021097)
    }
  }
  nmz <- c("ig_c22","ig_c14_a","ig_c2","ig_c14_b")
  reg.dat <- rep("immunoglobin",length(chr))
  if(bioC | text) {
    #must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chr,
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(text) { outData <- ranged.to.txt(outData) } else {
      if(GRanges) { outData <- as(outData,"GRanges") }
    }
  } else {
    outData <- vector("list",nchr); names(outData) <- paste("chr",1:nchr,sep="")
    for (cc in 1:nchr) {
      if(cc %in% chr) {
        outData[[cc]] <- list(start=stz[chr==cc],end=enz[chr==cc])
      }
    }
  }
  return(outData)
}


#' Return Centromere locations across the genome
#' 
#' Returns the locations of centromeres in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame.
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param autosomes logical, if TRUE, only return results for autosomes, if FALSE, also include
#' X and Y.
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be centromere chromosome and start and end positions.
#' @examples
#' setwd(tempdir())
#' get.centromere.locs()
#' get.centromere.locs(bioC=FALSE,autosomes=TRUE)
#' get.centromere.locs(text=TRUE)
get.centromere.locs <- function(dir=NULL,build=NULL,
                                bioC=TRUE,GRanges=TRUE,text=FALSE,autosomes=FALSE)
{
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  dir <- validate.dir.for(dir,c("ano"),warn=FALSE); success <- TRUE
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  local.file <- cat.path(dir$ano,"cyto")
  tt <- get.cyto(build=build,bioC=FALSE,dir=dir,GRanges=FALSE)
  chrn <- paste(1:22)
  if(!autosomes) { 
    chrn <- c(chrn,c("X","Y"))
  }
  nchr <- length(chrn)
  my.chr.range <- vector("list",nchr)
  names(my.chr.range) <- paste("chr",chrn,sep="")
  for (cc in 1:nchr) {
    just.centros <- tt[paste(tt[,5])=="acen",]
    just.chr <- just.centros[which(paste(just.centros[,1])==names(my.chr.range)[cc]),]
    my.chr.range[[cc]] <- list(start=min(just.chr[,2]), end=max(just.chr[,3]))
  }
  reg.dat <- rep("centromere",nchr)
  nmz <- paste(reg.dat,chrn,sep="_")
  stz <- sapply(my.chr.range,"[[",1)
  enz <- sapply(my.chr.range,"[[",2)
  if(bioC | text) {
    #must.use.package(c("genoset","IRanges"),bioC=TRUE)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=gsub("chr","",chrn),
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=TRUE)
    if(text) { 
      outData <- ranged.to.txt(outData) 
    } else {
      if(GRanges){
        outData <- as(outData,"GRanges")
      }
    }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}


#' Return Cytoband/Karyotype locations across the genome
#' 
#' Returns the locations of cytobands/karyotype-bands in the human genome, for a given build, as
#' a data.frame, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param refresh logical, whether to re-download the file if the existing file has become corrupted
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be centromere chromosome and start and end positions.
#' @examples
#' require(BiocInstaller)
#' setwd(tempdir())
#' get.cyto()
#' cyto.frame <- get.cyto(bioC=FALSE)
#' prv(cyto.frame)
#' get.cyto(build=36)
get.cyto <- function(build=NULL,dir=NULL,bioC=TRUE,GRanges=TRUE,refresh=FALSE) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  local.file="cyto"
  if(is.null(dir)) {
    local.file <- cat.path("",local.file,suf=build,ext="tar.gz")
  } else { 
    local.file <- cat.path(dir,local.file,suf=build,ext="tar.gz")
  }
  if(!file.exists(local.file) | refresh) {
    golden.path <- paste("http://hgdownload.cse.ucsc.edu/goldenPath/",build,"/database/cytoBand.txt.gz",sep="")
    success <- tryCatch( download.file(url=golden.path,local.file,quiet=T),error=function(e) { F } )
    if(is.logical(success)) {
      if(!success) { warning("couldn't reach ucsc website! try sourcing cytoband data elsewhere"); return(NULL) } }
    tt <- reader(local.file,header=FALSE)
    if(is.null(dir)) { unlink(local.file) }
  } else {
    tt <- reader(local.file)
  }
  colnames(tt) <- c("chr","start","end","band","negpos")
  write.table(tt,file=local.file,col.names=T,row.names=F,sep="\t",quote=F)
  mychr <- gsub("chr","",tt$chr,fixed=T)
  fullbands <- paste(mychr,tt$band,sep="")
  if(bioC ) {
    st <- as.numeric(tt$start)
    en <- as.numeric(tt$end)
   # must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=st,end=en,names=fullbands),space=mychr,
                          negpos=tt$negpos,universe=build[1])
    #prv(outData)
    outData <- toGenomeOrder2(outData) ##,strict=T)
    if(GRanges) { outData <- as(outData,"GRanges") }
  } else {
    outData <- tt 
    if("band" %in% colnames(outData)) {
      ## make 'chr-band' rownames to be consistent with the RangedData object if bioC=T
      rownames(outData) <- fullbands
      #outData <- outData[,-which(colnames(outData) %in% "band")]
    }
  }
  return(outData)
}



#' Get HapMap recombination rates for hg18 (build 36)
#' 
#' Recombination rate files can be used to calculate recombination distances
#' for genome locations, in centimorgans. This function downloads these reference
#' files from the hapmap NCBI website. At the time of writing they were only 
#' availble for build 36. If using a more recent build I suggest using the
#' conversion function conv.37.36(), then recomWindow(), then conv.36.37() to 
#' get recombination distances for other builds. If getOption("save.annot.in.current")
#' is <=0 then no files will be kept. Otherwise an object containing this mapping data
#' will be saved in the local directory if dir=NULL, or else in the directory specified.
#' Allowing this reference to be saved will greatly increase the speed of this function
#' for subsequent lookups
#' @param dir character, location to store binary file with the recombination maps for
#' chromosomes 1-22. If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param verbose logical, if the binary file is not already downloaded, when verbose
#' is TRUE, there will be some output to the console indicating the progress of the
#' download. If FALSE, all output is suppressed.
#' @param refresh logical, if you already have the binary file in the current directory,
#' this argument will let you re-download and re-generate this file, e.g, if the file
#' is modified or corrupted this will make a new one without having to manually delete it
#' @param compress logical, this argument is passed to 'save' and will result in a larger
#' binary file size, but quicker loading times, so 'FALSE' is recommended for faster retrieval.
#' @export
#' @return Returns a list object of length 22, containing the recombination map files
#' as 22 separate data.frame's.
#' @examples
#' \donttest{
#' ## not run as it takes roughly 2 minutes to download and read-in ##
#' setwd(tempdir())
#' rec.map <- get.recombination.map(getwd())
#' file.on.disk <- "rrates_genetic_map_chr_1_22_b36.RData"
#' if(file.exists(file.on.disk)) { unlink(file.on.disk) } # remove the downloaded file
#' }
get.recombination.map <- function(dir=NULL,verbose=TRUE,refresh=FALSE, compress=FALSE) {
  n.chr <- 22
  hap.dir <- "http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/latest/rates/"
  temp.dir <- "recombinationratesGF13fDR1er119"
  local.file <- "rrates_genetic_map_chr_1_22_b36.RData"
  if(!file.exists(temp.dir)) { dir.create(temp.dir) } 
  local.files=paste0(temp.dir,"/genetic_map_chr",1:n.chr,"_b36.txt")
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(!is.null(dir)) {
    local.files <- cat.path(dir,local.files,ext="txt")
    local.file <- cat.path(dir,local.file,ext="RData")
  }
  if(!file.exists(local.file) | refresh) {
    if(verbose) { cat("Downloading recombination data from: ",hap.dir,"\n") }
    hapmap.urls <- cat.path(dir=hap.dir,fn=basename(local.files))
    success <- TRUE
    for (cc in 1:n.chr) {
      #print(hapmap.urls[cc])
      success <- tryCatch( download.file(url=hapmap.urls[cc],local.files[cc],quiet=T),error=function(e) { F } )
      if(verbose) { loop.tracker(cc,n.chr*2) }
    }
    if(is.logical(success)) {
      if(!success) { warning("couldn't download at least one of the files from: ",hap.dir); return(NULL) } }
    map.files.list <- vector("list",n.chr)
    for (cc in 1:n.chr) {
      map.files.list[[cc]] <- read.table(local.files[cc],header=TRUE)
      if(is.data.frame(map.files.list[[cc]])) { 
        unlink(local.files[cc]) 
      } else { warning("downloaded map file was corrupt for chr",cc) }
      if(verbose) { loop.tracker(n.chr+cc,n.chr*2) }
    }
    if(file.exists(temp.dir)) { file.remove(temp.dir) }   # delete the temporary directory
  } else {
    map.files.list <- reader(local.file)
  }
  if(!is.null(dir)) { save(map.files.list,file=local.file,compress=compress) }
  if(length(map.files.list)!=n.chr) { stop("Unfortunately the object derived seems corrupted") }
  names(map.files.list) <- paste0("chr",1:n.chr)
  return(map.files.list)
}


#' Get exon names and locations from UCSC
#' 
#' Various R packages assist in downloading exonic information but often the input required is 
#' complex, or several lines of code are required to initiate, returning an object that
#' might require some manipulation to be useful. This function simplifies the job 
#' considerably, not necessarily requiring any arguments. The object returned can be
#' a standard data.frame or a bioconductor GRanges/RangedData object. The raw annotation
#' file downloaded will be kept in the working directory so that subsequent calls to
#' this function run very quickly, and also allow use offline.
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param transcripts logical, if TRUE, return transcripts rather than exons
#' @param GRanges logical, if TRUE and bioC is also TRUE, then returned object will be GRanges, otherwise
#' it will be RangedData
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome, start and end positions, transcript id number and name
#' @examples
#' \donttest{
#' setwd(tempdir())
#' get.exon.annot()
#' }
get.exon.annot <- function(dir=NULL,build=NULL,bioC=T, transcripts=FALSE, GRanges=TRUE) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  ## load exon annotation (store locally if not there already)
  from.scr <- T
  txt <- if(transcripts) { "trans" } else { "exon" }
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,"ano")
    ex.fn <- cat.path(dir$ano,pref=txt,"Annot",suf=build,ext="RData")
    if(file.exists(ex.fn)) {
      tS <- reader(ex.fn)
      if(transcripts) {
        if(is(tS)[1]=="data.frame") { from.scr <- F }
      } else {
        if(is(tS)[1]=="GRanges") { from.scr <- F }
      }
    }
  }
 # must.use.package("GenomicFeatures",T)
  if(from.scr) {
    #must.use.package("gage",T) this is where egSymb came from, but is now not needed
    # get transcripts from build table 'knownGene'
    success <- tryCatch(txdb <- suppressWarnings(makeTxDbFromUCSC(genome=build,
                                                         tablename="knownGene"))  ,error=function(e) { F } )
    if(is.logical(success)) { 
      if(!success) {
        warning("Couldn't reach build website! try again later or, \n",
                "if in europe/uk, there may still be a bug in rtracklayer; \n",
                "Installing the latest version of R and bioconductor\n",
                "and running biocLite('rtracklayer'), should fix this")
        return(NULL) }
    }
    if(!transcripts) {
      ex = exonsBy(txdb, by="gene")
      tS <- toGenomeOrder2(unlist(ex),strict=T)
    } else {  
      tS = transcriptsBy(txdb, by="gene")
      #egSymb <- egSymb <- reader("~/github/iChip/egSymb.rda")
      egSymb <- humarray::egSymb
      select <- match(names(tS),egSymb[,1])
      names(tS)[!is.na(select)] <- egSymb[,2][select[!is.na(select)]]
      tS <- as.data.frame(tS)
    }
  }
  if(exists("ex.fn")) { save(tS,file=ex.fn) }
  if(transcripts) {
    chrs <- paste(tS$seqnames); chrs <- gsub("chr","",chrs)
    tS$seqnames <- chrs
    if(all(colnames(tS)==c("element","seqnames","start","end","width","strand","tx_id","tx_name"))) {
      colnames(tS) <- c("gene","chr","start","end","width","strand","txid","txname")
    } else {
      cat(" unexpected colnames found using makeTxDbFrombuild()\n")
      if(bioC) { cat(" therefore returning data.frame instead of RangedData object\n") ;
                 bioC <- F }
    }
    if(bioC) {
      tS <- RangedData(ranges=IRanges(start=tS$start,end=tS$end),
                       space=tS$chr,gene=tS$gene, strand=tS$strand,
                       txid=tS$txid, txname=tS$txname,universe=build)
      tS <- toGenomeOrder2(tS,strict=T)
      if(GRanges) { tS <- as(tS,"GRanges") }
    }
    return(tS)
  } else {
    rownames(tS) <- paste(1:nrow(tS))
    ei <- mcols(tS)[["exon_id"]]
    ei2 <- add.trail(ei,suffix=strsplit("abcdefghijklmnopqrstuvwxyz","")[[1]])
    mcols(tS)[["exon_name"]] <- ei2
    if(bioC) {
      if(GRanges) {
        return(tS)
      } else {
        return(as(tS,"RangedData"))
      }
    } else {
      return(ranged.to.data.frame(tS))
    }
  }
}



#' Get human gene names and locations from biomart
#' 
#' Various R packages assist in downloading genomic information but often the input required is 
#' complex, or several lines of code are required to initiate, returning an object that
#' might require some manipulation to be useful. This function simplifies the job 
#' considerably, not necessarily requiring any arguments. The object returned can be
#' a standard data.frame or a bioconductor GRanges/RangedData object. The raw annotation
#' file downloaded will be kept in the working directory so that subsequent calls to
#' this function run very quickly, and also allow use offline.
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param duplicate.report logical, whether to provide a report on the genes labels that are listed
#' in more than 1 row - this is because some genes span ranges with substantial gaps within them
#' @param one.to.one logical, as per above, some genes have duplicate entries, sometimes for simplicity
#' you want just one range per gene, if this parameter is set TRUE, one range per gene is enforced,
#' and only the widest range will be kept by default for each unique gene label
#' @param remap.extra logical, whether to remap chromosome annotation for alternative builds and
#' unconnected segments to the closest regular chromosome, e.g, mapping MHC mappings to chromosome 6
#' @param discard.extra logical, similar to above, but if TRUE, then any non-standard chromosome
#' genes will just be discarded
#' @param only.named logical, biomart annotation contains some gene segments without names, if TRUE, then
#' such will not be included in the returned object (note that this will happen also if one.to.one is TRUE)
#' @param ens.id logical, whether to include the ensembl id in the dataframe
#' @param refresh logical, if you already have the file in the current directory,
#' this argument will let you re-download and re-generate this file, e.g, if the file
#' is modified or corrupted this will make a new one without having to manually delete it
#' @param GRanges logical, if TRUE and bioC is also TRUE, then returned object will be GRanges, otherwise
#' it will be RangedData
#' @param host.txt character, the argument to pass to biomaRt::useMart(). Default for build 36 is 
#' 'may2009.archive.ensembl.org', and for build 37, "feb2014.archive.ensembl.org" but for recent builds
#'  the recommended link is 'www.ensembl.org'
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome and start and end positions, other information depends on 
#' specific parameters documented above
#' @examples
#' \donttest{
#' setwd(tempdir())
#' get.gene.annot()
#' }
get.gene.annot <- function(dir=NULL,build=NULL,bioC=TRUE,duplicate.report=FALSE,
                           one.to.one=FALSE,remap.extra=FALSE,discard.extra=TRUE,only.named=FALSE,
                           ens.id=FALSE,refresh=FALSE,GRanges=TRUE, host.txt="") {
  # faster than exon, but only contains whole gene ranges, not transcripts
  # allows report on duplicates as some might be confused as to why some genes
  # have more than one row in the listing (split across ranges usually)
  # run with dir as NULL to refresh changes in COX
  #must.use.package(c("biomaRt","genoset","gage"),T)
  mart.txt <- "ENSEMBL_MART_ENSEMBL"
  verbose <- TRUE # hard coded at this stage!!
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(build %in% c("hg15","hg16","hg17")) { stop("older builds, prior to may 2009 are not supported by this function")}
  from.scr <- T
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,"ano")
    utxt <- ""; if(one.to.one) { utxt <- "_unq" }
    if(ens.id) { utxt <- paste(utxt,"ens",sep="_") }
    if(only.named) { utxt <- paste(utxt,"onm",sep="_") }
    gn.fn <- cat.path(dir$ano,"geneAnnot",pref=build,suf=utxt,ext="RData")
    if(file.exists(gn.fn) & !refresh) {
      dat <- get(paste(load(gn.fn)))
      from.scr <- F
    }
  }
  # colnames for output
  nm.list <- c("gene","chr","start","end","band")
  if(ens.id & bioC) { warning("ens.id=TRUE only has an effect when bioC=FALSE") }
  if(from.scr) {
    if(build=="hg18") {
      if(all(host.txt=="")) { host.txt <- "may2009.archive.ensembl.org" }
      ens <- biomaRt::useMart(mart.txt,
                              dataset="hsapiens_gene_ensembl",
                              host= host.txt,
                              path="/biomart/martservice",
                              archive=FALSE)
    } else {
      if(build=="hg19") {
        if(all(host.txt=="")) { host.txt <- "feb2014.archive.ensembl.org" }
        ens <- biomaRt::useMart(mart.txt,
                                dataset="hsapiens_gene_ensembl",
                                host= host.txt,
                                path="/biomart/martservice",
                                archive=FALSE)
      } else {
        # whatever the current mart is
        if(all(host.txt=="")) { host.txt <- "www.ensembl.org" }
        ens <- biomaRt::useMart(mart.txt,host=host.txt)
      }
    }
    ens <- biomaRt::useDataset("hsapiens_gene_ensembl",mart=ens)
    attr.list <- c("hgnc_symbol", "chromosome_name",
                   "start_position", "end_position", "band")
    if(ens.id) { attr.list <- c(attr.list,"ensembl_gene_id") }
    #    if(only.named & !ens.id & build=="hg18") {
    #      egSymb <- reader("~/github/iChip/egSymb.rda") # this is a hg18 list!
    #      dat <- biomaRt::getBM(attributes = attr.list, filters = "hgnc_symbol",
    #                   values = egSymb[,2], mart = ens)
    #    } else {
    dat <- biomaRt::getBM(attributes = attr.list, mart = ens)
    #    }
    if(exists("gn.fn")) { save(dat,file=gn.fn) }
  } 
  if(ens.id) { nm.list <- c(nm.list,"ens.id") }
  #return(dat)
  no.gene.names <- which(paste(dat[[1]])=="")
  #prv(no.gene.names)
  if((one.to.one | (only.named & !ens.id)) & length(no.gene.names)>0) { dat <- dat[-no.gene.names,] }
  if(remap.extra) {
    dat$chromosome_name <- tidy.extra.chr(dat$chromosome_name)
    #dat$chromosome_name[grep("c6",dat$chromosome_name,ignore.case=T)] <- 6  # prevent issues with c6_COX, c6_QBL  
    #dat$chromosome_name[grep("c5",dat$chromosome_name,ignore.case=T)] <- 5  # prevent issues with c5_H2  
    #dat$chromosome_name[grep("NT",dat$chromosome_name,ignore.case=T)] <- "Z_NT"  # merge all NT regions to one label
  }
  if(discard.extra) {
    ## http://www.lrg-sequence.org/ ##
    # note that if remapping is already done, then these won't be discarded unless remapping failed
    tt <- tidy.extra.chr(dat$chromosome_name,select=TRUE)
    badz <- which(!tt)
    #prv(badz)
    if(length(badz)>0) { dat <- dat[-badz,] } # remove LRG, GS, HG, NT, COX, etc annotation from set
  }
  if(bioC) {
    missin <- function(x) { is.na(x) | (x=="") | x=="NA" }
  #  keep <- rep(T,nrow(dat))
    stz <- dat$start_position; enz <- dat$end_position
    whichmis <- missin(stz) | missin(enz)
    if(any(whichmis)) { dat <- dat[!whichmis,]; warning("some start/end positions were missing") }
    stz <- dat$start_position; enz <- dat$end_position
    if(any(stz>enz)) { mm <- stz>enz; prv(cbind(stz,enz)[mm,]); tmp <- stz; stz[mm] <- enz[mm]; enz[mm] <- tmp[mm]; warning("some start positions were after end positions (swapped these)") }
    dat$start_position <- stz; dat$end_position <- enz
#    bad1s <- NULL
#    for (jj in 1:22) { ii <- ga[ga$chr==jj,"end"]>get.chr.lens()[jj]; if(any(ii)) { bad1s <- c(bad1s,which(ga$chr==jj)[ii]) } }
#    if(length(bad1s)>0) { dat <- dat[-bad1s,] ; warning(length(bad1s)," gene positions exceeded chromosome length") }
    outData <- RangedData(ranges=IRanges(start=dat$start_position,end=dat$end_position),
                          space=dat$chromosome_name,gene=dat$hgnc_symbol, band=dat$band, universe=build)
    outData <- toGenomeOrder2(outData,strict=T)
    if(duplicate.report | one.to.one) {
      genez <- outData$gene
      dG <- which(duplicated(genez))
      if(!duplicate.report) {
        dup.genes <- genez[dG]
      } else {
        dup.genes <- gene.duplicate.report(outData)
      }
      stz <- start(outData); enz <- end(outData); wdz <- width(outData)
      to.del <- to.ch <- ch.st <- ch.en <- NULL
      n.dup <- length(dup.genes)
      if(one.to.one & n.dup>0) {
        #return(dup.genes)
        indz <- sapply(as.list(unique(dup.genes)),function(X) { which(genez %in% X) })
        # keep the range that is widest
        st.en <- lapply(indz,function(X) { c(stz[X][wdz[X]==max(wdz[X])][1],enz[X][wdz[X]==max(wdz[X])][1]) } )
        to.del <- dG
        to.ch <- sapply(indz,min,na.rm=T)
        ch.st <- sapply(st.en,"[",1) 
        ch.en <- sapply(st.en,"[",2) 
        start(outData)[to.ch] <- 1
        end(outData)[to.ch] <- ch.en 
        start(outData)[to.ch] <- ch.st 
        outData <- outData[-to.del,]
        if(verbose) { warning(cat("kept widest ranges, merging",length(to.del)+length(to.ch),"duplicate gene labels to",length(to.ch),"labels\n")) }
      } else {
        #cat("no dups found")
      }
      # but haven't done anything about them or removed them! 
    }
    if(GRanges) { outData <- as(outData, "GRanges") }
  } else {
    outData <- dat; colnames(outData) <- nm.list
  }
  return(outData)
}



#' Derive Telomere locations across the genome
#' 
#' Returns the locations of telomeres in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param kb The number of base pairs at the start and end of a chromosome that are defined as
#' belonging to the telomere can be a little arbitrary. This argument allows specification
#' of whatever threshold is required.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param autosomes logical, if TRUE, only return results for autosomes, if FALSE, also include
#' X and Y.
#' @param mito.zeros logical, Mitochondria have no telomeres (are circular) but for some purposes you
#' might want zero values in order to match with other annotation that includes all chromosomes and MT.
#' TRUE adds zeros for chrMT, and FALSE excludes chrMT.
#' @export
#' @return Returns a text vector, GRanges or RangedData object, depending on input parameters. Contained
#' will be telomere chromosome and start and end positions.
#' @examples
#' setwd(tempdir())
#' get.telomere.locs()
#' get.telomere.locs(bioC=FALSE)
#' get.telomere.locs(text=TRUE)
get.telomere.locs <- function(dir=NULL,kb=10,build=NULL,bioC=TRUE,GRanges=TRUE,
                              text=FALSE,autosomes=FALSE,mito.zeros=FALSE)
{
  # the actual telomeres are typically about 10kb, but
  # for cnv-QC purposes want to exclude a larger region like 500kb
  # Mt have no telomeres, are circular, but for some purposes might want zero values in there
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  chr.lens <- get.chr.lens(dir=dir,build=build[1],autosomes=FALSE,mito=mito.zeros)
  n <- 1:22; if(!autosomes) { n <- c(n,"X","Y") } # Mt have no telomeres, are circular
  nchr <- length(n) #default
  if(mito.zeros) { n <- c(n,"M") }
  my.chr.range <- vector("list",length(n))
  names(my.chr.range) <- paste("chr",n,sep="")
  for (cc in 1:nchr) {
    one <- force.chr.pos(Pos=c(1,kb*1000),Chr=cc,build=build) # f..c..pos() makes sure is a valid range
    two <- force.chr.pos(Pos=chr.lens[cc]+c(-kb*1000,0),Chr=cc,build=build)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  if(mito.zeros) {
    # add null values for the Mitochondrial chromosome
    cc <- cc+1; one <- c(1,1);  two <- chr.lens[cc]+c(0,0)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  reg.dat <- rep("telomere",length(n)*2)
  chrz <- rep(n,each=2)
  nmz <- paste(reg.dat,chrz,rep(c("a","b"),times=length(n)),sep="_")
  stz <- as.vector(sapply(my.chr.range,"[[",1))
  enz <- as.vector(sapply(my.chr.range,"[[",2))
  if(bioC | text) {
    #must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chrz,
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(text) { outData <- ranged.to.txt(outData) } else { if(GRanges) { outData <- as(outData,"GRanges") } }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}



#' Get chromosome lengths from build database
#' 
#' Quick and easy way to retrieve human chromosome lengths. Can select from hg18/hg19 (ie, 
#'  build 36/37), or any future builds (hg20, etc) stored in the same location on the build website.
#'  Default is to return lengths for 22 autosomes, but can also retrieve X,Y 
#'  and Mitochondrial DNA lengths by 'autosomes=FALSE' or n=1:25. Even if not connected to 
#'  the internet can retrieve hard coded lengths for hg18 or hg19.
#'
#' @param dir directory to retrieve/download the annotation from/to (defaults to current getwd())
#'  if dir is NULL then will automatically delete the annotation text file from the local directory
#'   after downloading
#' @param build string, currently 'hg17','hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is getOption("ucsc"). Will also accept integers 17,18,19,35,36,37 as alternative arguments.
#' @param autosomes logical, if TRUE, only load the lengths for the 22 autosomes, else load X,Y,[MT] as well
#' @param len.fn optional file name to keep the lengths in
#' @param mito logical, whether to include the length of the mitochondrial DNA (will not include unless
#'  autosomes is also FALSE)
#' @param names logical, whether to name the chromosomes in the resulting vector
#' @param delete.after logical, if TRUE then delete the text file that these lengths were downloaded to. 
#' @param verbose logical, if TRUE display extra information on progress of chromsome retrieval
#'  If FALSE, then the file will be kept, meaning future lookups will be faster, and available offline.
#' @export
#' @examples
#'  setwd(tempdir())
#'  get.chr.lens(delete.after=TRUE) # delete.after simply deletes the downloaded txt file after reading
#'  get.chr.lens(build=35,autosomes=TRUE,delete.after=TRUE) # only for autosomes
#'  get.chr.lens(build="hg19",mito=TRUE,delete.after=TRUE) # include mitochondrial DNA length
get.chr.lens <- function(dir=NULL,build=NULL,autosomes=FALSE,len.fn="humanChrLens.txt",
                         mito=FALSE,names=FALSE, delete.after=FALSE, verbose=FALSE)
{
  # retrieve chromosome lengths from local annotation file, else download from build
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(dir)) { dir <- getwd() ; delete.after <- TRUE }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  chrlens.f <- cat.path(dir$ano,len.fn) # existing or future lengths file
  n <- 1:22; if(!autosomes) { n <- c(n,"X","Y","M") }
  hg18.backup <- c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
                   146274826,140273252,135374737,134452384,132349534,114142980,106368585,
                   100338915,88827254,78774742,76117153,63811651,62435964,46944323,
                   49691432,154913754,57772954,16571)
  hg19.backup <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
                   146364022,141213431,135534747,135006516,133851895,115169878,107349540,
                   102531392,90354753,81195210,78077248,59128983,63025520,48129895,
                   51304566,155270560,59373566,16571)
  hg38.backup <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,
                   145138636,138394717,133797422,135086622,133275309,114364328,107043718,
                   101991189,90338345,83257441,80373285,58617616,64444167,46709983,
                   50818468,156040895,57227415,16571)
  # backups for offline use
  if(build=="hg18") { 
    offline.backup <- hg18.backup 
  } else {
    if(build=="hg19") {
      offline.backup <- hg19.backup 
    } else {
      offline.backup <- hg38.backup 
    }
  }
  if(file.exists(chrlens.f))
  {
    # file seems to be in annotation directory already
    chrLens <- readLines(chrlens.f)
    if (length(chrLens)!=length(n))
    {
      #warning("Length of existing chromosome file didn't match expected:",length(n))
      notGot <- T
    } else {
      notGot <- F
      # we have the right length, but do we have the right version?
      if(build=="hg18" & ( length(which(chrLens %in% hg38.backup))>2 | length(which(chrLens %in% hg19.backup))>2) ) { notGot <- T }
      if(build=="hg19" & ( length(which(chrLens %in% hg18.backup))>2 | length(which(chrLens %in% hg38.backup))>2) ) { notGot <- T }
      if(build=="hg38" & ( length(which(chrLens %in% hg18.backup))>2 | length(which(chrLens %in% hg19.backup))>2) ) { notGot <- T }
      names(chrLens) <- paste0("chr",n)
    }
  } else { notGot <- T }
  if (notGot | (!build %in% c("hg18","hg19","hg38"))) {
    #download from build
    if(verbose) { cat("attempting to download chromosome lengths from genome build ... ") }
    urL <- switch(build,
                  hg17="http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/chromInfo.txt.gz",
                  hg18="http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz",
                  hg19="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz",
                  hg38="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz")
    success <- T
    success <- tryCatch(download.file(urL, chrlens.f,quiet=T),error=function(e) { F } )
    if(!is.logical(success)) { success <- T }
    if(success) {
      if(verbose) {  cat("download successful\n") }
      chrL.f <- readLines(chrlens.f)
      len.lst <- strsplit(chrL.f,"\t")
      nmz <- sapply(len.lst,"[",1)
      lnz <- sapply(len.lst,"[",2)
      nnn <- paste0("chr",n)
      want.chr.names <- match(nnn,nmz)
      want.chr.names <- want.chr.names[!is.na(want.chr.names)]
      #print(want.chr.names)
      chrLens <- lnz[want.chr.names]
      names(chrLens) <- nnn
    } else {
      warning("couldn't reach build website, so have used offline versions of chr lengths")
      if(!build %in% c("hg18","hg19","hg38")) { warning("no offline version for build version:",build) }
      chrLens <- paste(offline.backup)[1:length(n)]; names(chrLens) <- paste0("chr",n)
      #print(n)
      delete.after <- F
    }
    if(length(dir)==1 & dir[1]=="" & delete.after) {
      unlink(chrlens.f)
    } else {
      writeLines(chrLens,con=chrlens.f) # save file for future use
    }
  }
  if(!mito & length(chrLens)>22) { chrLens <- chrLens[-grep("M",n)] }
  if(names) {
    return(chrLens)
  } else {
    return(as.numeric(chrLens))
  }
}



#' Obtain a listing of known T1D associated genomic regions
#'
#' This function uses a full list of ichip dense regions combined with a list of t1d
#' SNPs to get the t1d regions. For type 1 diabetes researchers.
#' @param dense.reg GRanges or RangedData object, only use if you need to provide for a
#' build other than 36 or 37 (hg18/hg19).
#' @param build e.g, 36/hg18 or 37/hg19, if left as NULL current getOption('ucsc') will
#' be used.
#' @param invert logical, set to TRUE if you wish to get the set of NON-T1D regions.
#' @return a GRanges object with the specified type 1 diabetes (or inverse) ranges
#' @export
#' @examples
#' \donttest{
#' # t1d.reg <- get.t1d.regions()
#' # non.t1d <- get.t1d.regions(build=36,invert=TRUE)
#' }
get.t1d.regions <- function(dense.reg=NULL,build=NULL,invert=FALSE) {
  #source("~/github/iChip/iFunctions.R")
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  def.build <- getOption("ucsc")
  if(is.null(dense.reg)) { 
    #dense.reg <- reader("~/github/iChip/iChipFineMappingRegionsB36.RData") 
    dense.reg <- humarray::iChipRegionsB36
    if(build!="hg18") {
      if(build=="hg19") { 
        dense.reg <- conv.36.37(dense.reg) 
      } else { 
        if(build=="hg38") {
          dense.reg <- conv.37.38(conv.36.37(dense.reg))
        } else {
          stop("automatic loading for dense regions is only supported for builds 36, 37 or 38") 
        }
      }
    }
  }  
  if(is(dense.reg)[1]=="GRanges") { dense.reg <- as(dense.reg,"RangedData") }
  if(!is(dense.reg)[1]=="RangedData") { stop("dense.reg must be RangedData or GRanges") }
  ichip.regions <- dense.reg
  rs.ids <- get.immunobase.snps()
  locs <- Pos(rs.ids); chrs <- Chr(rs.ids)
  good <- !is.na(locs) & !is.na(chrs)
  locs <- locs[good]; chrs <- chrs[good]
  if(build!=def.build) {
    if(!all(c(build,def.build) %in% c("hg18","hg19"))) {
      stop("if build is not equal to getOption('ucsc') then build and this getOption('ucsc')",
           " parameter must both be equivalent to either hg18 or hg19. To use build 38,
           first set 'options(ucsc=38)'") 
    }
    if(build=="hg18") {
      #prv(chrs,locs)
      locs <- conv.37.36(chr=chrs,pos=locs)[,"start"]
    } else {
      locs <- conv.36.37(chr=chrs,pos=locs)[,"start"]
    }
  }
  t1dgr <- makeGRanges(chr=chrs,pos=locs,build=build)
  #prv(t1dgr,ichip.regions)
  t1d.regions <- suppressWarnings(subsetByOverlaps(set.chr.to.numeric(as(ichip.regions,"GRanges")),t1dgr))
  #  t1dgr <- as(makeGRanges(chr=chrs,pos=locs,build=build),"RangedData")
  #  t1d.regions <- find.overlaps(ichip.regions,ref=t1dgr,thresh=0.000000000001,ranges.out=TRUE)
  #prv(t1d.regions)
  if(invert) { t1d.regions <- invGRanges(t1d.regions,build=build) }
  return(t1d.regions)
}



#' Obtain subset of ranged object overlapping known T1D associated genomic regions
#'
#' Return subset of a ranged object that overlaps ichip dense mapped regions. For type 1 diabetes
#' and autoimmune disesase researchers.
#' @param X GRanges or RangedData object, ranged object for which you want the T1D subset
#' of ranges/SNPs.
#' @param ichip.regions GRanges or RangedData object, only use if you need to provide for a
#' build other than 36 or 37 (hg18/hg19), or for multiple lookups to avoid reloading each time
#' @param T1D.regions GRanges or RangedData object, only use if you need to provide for a
#' build other than 36 or 37 (hg18/hg19), or for multiple lookups to avoid reloading each time.
#' @param build e.g, 36/hg18 or 37/hg19, if left as NULL current getOption('ucsc') will
#' be used.
#' @param T1D.only logical, standard is to return type 1 diabetes (T1D) regions subset, but
#' if this parameter is set to FALSE, will return the subset for all 12 autoimmune diseases
#' mapped by the ImmunoChip consortium. (Cortes and Brown, 2010).
#' @param invert logical, set to TRUE if you wish to get the set of NON-T1D regions, or
#' non-immune dense regions when T1D.only=FALSE.
#' @return a GRanges object with the specified type 1 diabetes/autoimmune (or inverse) ranges
#' @export
#' @examples
#' \donttest{
#' # all.reg <- rranges(10000)
#' # t1d <- get.t1d.subset(all.reg) # T1D regions
#' # non.autoimmune <- get.t1d.subset(T1D.only=FALSE,build=36,invert=TRUE) # non-autoimmune regions
#' }
get.t1d.subset <- function(X,T1D.only=TRUE,build=NULL,ichip.regions=NULL,T1D.regions=NULL,invert=FALSE) {
  #source("~/github/iChip/iFunctions.R")
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(ichip.regions)) {
    #ichip.regions <- reader("~/github/iChip/iChipFineMappingRegionsB36.RData")
    ichip.regions <- humarray::iChipRegionsB36
  }
  if(T1D.only) {
    if(is.null(T1D.regions)) {
      T1D.regions <- get.t1d.regions(ichip.regions,build=build,invert=invert)
    }
    #T1D.regions <- as(T1D.regions,"RangedData")
    filt.sd <- suppressWarnings(subsetByOverlaps(set.chr.to.numeric(as(X,"GRanges")),set.chr.to.numeric(as(T1D.regions,"GRanges"))))
    #filt.sd <- find.overlaps(X,ref=T1D.regions,thresh=0.000000000000001,ranges.out=TRUE)
  } else {
#    filt.sd <- find.overlaps(X,ref=ichip.regions,thresh=0.000000000000001,ranges.out=TRUE)
    if(invert) { ichip.regions <- invGRanges(ichip.regions,build=build) }
    filt.sd <- suppressWarnings(subsetByOverlaps(set.chr.to.numeric(as(X,"GRanges")),set.chr.to.numeric(as(ichip.regions,"GRanges"))))
  }
  return(filt.sd)
}



#' Obtain subset of ranged object overlapping human genes
#'
#' Return subset of a ranged object that overlaps human genes (or custom ranges/exons).
#' A wrapper for subsetByOverlaps that tries to ensure equivalent builds and chromosome
#' labelling are taken care of automatically, and provides the default case of genic subsetting
#' where no explicit ref parameter is required.
#' @param X GRanges or RangedData object, ranged object for which you want the genic subset
#' of ranges/SNPs.
#' @param ref GRanges or RangedData object, only use if you need to provide for a
#' build other than 36 or 37 (hg18/hg19), or to use this function for another reference set,
#' for instance you could provide an object with exons
#' @param build e.g, 36/hg18 or 37/hg19, if left as NULL current getOption('ucsc') will
#' be used.
#' @return a GRanges object with the specified genic (or custom) ranges
#' @seealso \code{\link{get.gene.annot}}, \code{\link{get.exon.annot}}, \code{\link{subsetByOverlaps}}
#' @export
#' @examples
#' \donttest{
#' # all.reg <- rranges(1000) # random set of 1000 regions
#' # genic <- get.genic.subset(all.reg) # gene regions from the random set
#' # exonic <- get.genic.subset(all.reg,ref=get.exon.annot()) # exonic regions from the random set
#' }
get.genic.subset <- function(X,ref=NULL,build=NULL) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(!is(ref)[1] %in% c("RangedData","GRanges")) {
    if(!is.null(ref)) { warning("ref was not GRanges or RangedData, so reverted to get.gene.annot()") }
    ref <- get.gene.annot(build=build)
  }
  #if(!is.character(DB)) { stop() }
  #if(!DB %in% c("gene","exon"))
  #filt.sd <- find.overlaps(X,db=DB,thresh=0.00000000000001,ranges.out=TRUE,...)
  filt.sd <- suppressWarnings(subsetByOverlaps(set.chr.to.numeric(as(X,"GRanges")),set.chr.to.numeric(as(ref,"GRanges"))))
  return(filt.sd)
}


################## end annotation ##########################



######################
## Ranged Functions ##
######################



#' Wrapper to construct GRanges object from chr,pos or chr,start,end
#' 
#' Slightly simplifies the creation of a GRanges object, allowing flexible input of
#' chr, pos, or chr,start,end, and specification of rownames and the 'genome' parameter
#' for specifying the build/coordinate type, e.g, hg18, build 37, etc. Designed for
#' a simplified GRanges object without metadata, and where the 'strand' data is of
#' no interest, so if strand/metadata is to be used, use the original GRanges() constructor.
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions for the GRanges object
#' @param pos integer/numeric, for SNPs, can enter positions just once in 'pos' instead of 
#' entering the same value for start and end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param start integer/numeric, specify the start position of ranges to encode in the new
#' GRanges object, alongside 'end' (do not use 'pos' if using start+end)
#' @param end integer/numeric, specify the end position of ranges to encode in the new
#' GRanges object, alongside 'start' (do not use 'pos' if using start+end)
#' @param row.names character, rownames for the output object, e.g, unique IDs describing the 
#' ranges
#' @param ... further arguments to df.to.GRanges, such as 'fill.missing'
#' @return Returns a GRanges object with the ranges, build and rownames specified. Rownames
#' will be 1:nrow if the 'row.names' parameter is empty. The strand information will default
#' to '+' for all entries, and the metadata will be empty (this function is only for creation
#' of a very basic GRanges object).
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}, \code{\link{Pos}}, \code{\link{Pos.gene}}, \code{\link{Band}},
#'  \code{\link{Band.gene}}, \code{\link{Band.pos}}, \code{\link{Gene.pos}}, \code{zlink{df.to.GRanges}}
#' @examples
#' g1 <- makeGRanges(chr=c(1,4,"X"),pos=c(132432,434342,232222))
#' g2 <- makeGRanges(chr=c(22,21,21),start=c(1,1,1),end=c(1000,10000,100000),
#'                                               row.names=c("1K","10K","100K"))
#' g1 ; g2
makeGRanges <- function(chr,pos=NULL,start=NULL,end=NULL,row.names=NULL,build=NULL,...) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(start) & is.null(end) & !is.null(pos)) {
    dF <- cbind(paste(chr),round(as.numeric(pos)))
  } else {
    if(!is.null(start) & !is.null(end) & is.null(pos)) {
      pos <- cbind(round(as.numeric(start)),round(as.numeric(end)))
      dF <- cbind(paste(chr),pos) 
    } else {
      stop("must use either 'pos' or 'start' and 'end'")
    }
  }
  #prv(paste(row.names)); prv(dF)
  if(is.character(row.names)) { if(length(row.names)==nrow(dF)) { rownames(dF) <- row.names } else { warning("row.names had an incorrect length")} }
  if(!any(Dim(pos)==length(chr))) { stop("chr and pos must be of the same length") }
  if(length(Dim(pos))>1) { 
    if(!ncol(pos) %in% c(1,2)) { 
      stop("pos must be a vector of SNP locations, or a 2-column object with start and end coordinates") 
    } else {
      if(ncol(pos)==2) {
        if(any(pos[,2]<pos[,1])) { warning("end coordinates should be equal or greater than start coordinates") }
      }
    }
  }
  if(ncol(dF)==3) { 
    colnames(dF) <- c("chr","start","end") } else { colnames(dF) <- c("chr","pos") }
  #return(dF)
  ranged <- df.to.GRanges(dF,start=colnames(dF)[2],end=tail(colnames(dF),1),build=build,...) 
  if(!any(rownames(ranged) %in% row.names)) {
    if(is.character(row.names)){ if(length(row.names)==nrow(ranged)) { rownames(ranged) <- row.names }}
  }
  return(ranged)
}




#' Convert from build 37 to build 36 SNP coordinates
#' 
#' Convert range or SNP coordinates between builds using a chain file. Depending on the chain file
#' this can do any conversion, but the default will use the hg19 to hg18 (37-->36) chain file
#' built into this package. The positions to convert can be entered using using chr, pos vectors,
#'  or a RangedData or GRanges object. This function is a wrapper for liftOver() from rtracklayer,
#' providing more control of input and output and 'defensive' preservation of order and length
#' of the output versus the input ranges/SNPs. 
#' @param chr character, an optional vector of chromosomes to combine with 'pos' to describe
#'  positions to convert from build hg19 to hg18
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' a ranges object if this is provided along with 'chr'
#' @param ranges optional GRanges or RangedData object describing positions for which conversion
#' should be performed. No need to enter chr, pos if using ranges
#' @param ids if the ranges have ids (e.g, SNP ids, CNV ids), then by including this parameter
#' when using chr, pos input, the output object will have these ids as rownames. For ranges input
#' these ids would already be in the rownames of the GRanges or RangedData object, so use of
#' this parameter should be unnecessary
#' @param ... additional arguments to makeGRanges(), so in other words, can use 'start' and
#' 'end' to specify ranges instead of 'pos'.
#' @return Returns positions converted from build 37 to 36. If using the 'ranges' parameter 
#' for position input, the object returned will be of the same format. If using chr and pos 
#' to input, then the object returned will be a data.frame with columns, chr and pos with 
#' rownames 'ids'. Output will be the same length as the input, which is not necessarily the
#'  case for liftOver() which does the core part of this conversion. Using vector or GRanges 
#'  input will give a resulting data.frame or GRanges object respectively that has the same
#'  order of rownames as the original input. Using RangedData will result in an output that
#'   is sorted by genome order, regardless of the original order.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{conv.36.37}}, \code{\link{conv.37.38}}, 
#' \code{\link{conv.38.37}}, \code{\link{convTo37}}, \code{\link{convTo36}}
#' @examples
#' \donttest{
#' gene.labs <- c("CTLA4","IL2RA","HLA-C")
#' pp <- Pos.gene(gene.labs,build=37)
#' gg <- GRanges(ranges=IRanges(start=pp$start,end=pp$end),seqnames=pp$chr)
#' conv.37.36(gg) # order of output is preserved   ### HERE!!! ###
#' rr <- as(gg,"RangedData")
#' conv.37.36(rr) # note the result is same as GRanges, but in genome order
#' }
conv.37.36 <- function(ranges=NULL,chr=NULL,pos=NULL,...,ids=NULL) {
  #chain.file <- "~/github/iChip/hg19ToHg18.over.chain"
  #chain.file <- system.file("extdata", "hg19ToHg18.over.chain", package="humarray")
  chain.file <- humarray::hg19ToHg18
  return(conv.36.37(ranges=ranges,chr=chr,pos=pos,...,ids=ids,chain.file=chain.file))
}

#' Convert from build 38 to build 37 SNP coordinates
#' 
#' Convert range or SNP coordinates between builds using a chain file. Depending on the chain file
#' this can do any conversion, but the default will use the hg38 to hg19 (38-->37) chain file
#' built into this package. The positions to convert can be entered using using chr, pos vectors,
#'  or a RangedData or GRanges object. This function is a wrapper for liftOver() from rtracklayer,
#' providing more control of input and output and 'defensive' preservation of order and length
#' of the output versus the input ranges/SNPs. 
#' @param chr character, an optional vector of chromosomes to combine with 'pos' to describe
#'  positions to convert from build hg38 to hg19
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' a ranges object if this is provided along with 'chr'
#' @param ranges optional GRanges or RangedData object describing positions for which conversion
#' should be performed. No need to enter chr, pos if using ranges
#' @param ids if the ranges have ids (e.g, SNP ids, CNV ids), then by including this parameter
#' when using chr, pos input, the output object will have these ids as rownames. For ranges input
#' these ids would already be in the rownames of the GRanges or RangedData object, so use of
#' this parameter should be unnecessary
#' @param ... additional arguments to makeGRanges(), so in other words, can use 'start' and
#' 'end' to specify ranges instead of 'pos'.
#' @return Returns positions converted from build 38 to 37. If using the 'ranges' parameter 
#' for position input, the object returned will be of the same format. If using chr and pos 
#' to input, then the object returned will be a data.frame with columns, chr and pos with 
#' rownames 'ids'. Output will be the same length as the input, which is not necessarily the
#'  case for liftOver() which does the core part of this conversion. Using vector or GRanges 
#'  input will give a resulting data.frame or GRanges object respectively that has the same
#'  order of rownames as the original input. Using RangedData will result in an output that
#'   is sorted by genome order, regardless of the original order.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{conv.36.37}}, \code{\link{conv.37.36}}, 
#' \code{\link{conv.37.38}}, \code{\link{convTo37}}, \code{\link{convTo36}}
#' @examples
#' \donttest{
#' gene.labs <- c("CTLA4","IL2RA","HLA-C")
#' pp <- Pos.gene(gene.labs,build=38)
#' gg <- GRanges(ranges=IRanges(start=pp$start,end=pp$end),seqnames=pp$chr)
#' conv.38.37(gg) # order of output is preserved   ### HERE!!! ###
#' rr <- as(gg,"RangedData")
#' conv.38.37(rr) # note the result is same as GRanges, but in genome order
#' }
conv.38.37 <- function(ranges=NULL,chr=NULL,pos=NULL,...,ids=NULL) {
  #chain.file <- "~/github/iChip/hg38ToHg19.over.chain"
  #chain.file <- humarray::hg38ToHg19.over.chain
  #chain.file <- system.file("extdata", "hg38ToHg19.over.chain", package="humarray")
  chain.file <- humarray::hg38ToHg19
  return(conv.36.37(ranges=ranges,chr=chr,pos=pos,...,ids=ids,chain.file=chain.file))
}


#' Convert from build 37 to build 38 SNP coordinates
#' 
#' Convert range or SNP coordinates between builds using a chain file. Depending on the chain file
#' this can do any conversion, but the default will use the hg19 to hg38 (37-->38) chain file
#' built into this package. The positions to convert can be entered using using chr, pos vectors,
#'  or a RangedData or GRanges object. This function is a wrapper for liftOver() from rtracklayer,
#' providing more control of input and output and 'defensive' preservation of order and length
#' of the output versus the input ranges/SNPs. 
#' @param chr character, an optional vector of chromosomes to combine with 'pos' to describe
#'  positions to convert from build hg19 to hg38
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' a ranges object if this is provided along with 'chr'
#' @param ranges optional GRanges or RangedData object describing positions for which conversion
#' should be performed. No need to enter chr, pos if using ranges
#' @param ids if the ranges have ids (e.g, SNP ids, CNV ids), then by including this parameter
#' when using chr, pos input, the output object will have these ids as rownames. For ranges input
#' these ids would already be in the rownames of the GRanges or RangedData object, so use of
#' this parameter should be unnecessary
#' @param ... additional arguments to makeGRanges(), so in other words, can use 'start' and
#' 'end' to specify ranges instead of 'pos'.
#' @return Returns positions converted from build 37 to 38. If using the 'ranges' parameter 
#' for position input, the object returned will be of the same format. If using chr and pos 
#' to input, then the object returned will be a data.frame with columns, chr and pos with 
#' rownames 'ids'. Output will be the same length as the input, which is not necessarily the
#'  case for liftOver() which does the core part of this conversion. Using vector or GRanges 
#'  input will give a resulting data.frame or GRanges object respectively that has the same
#'  order of rownames as the original input. Using RangedData will result in an output that
#'   is sorted by genome order, regardless of the original order.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{conv.36.37}}, \code{\link{conv.37.36}}, 
#' \code{\link{conv.37.38}}, \code{\link{convTo37}}, \code{\link{convTo36}}
#' @examples
#' \donttest{
#' gene.labs <- c("CTLA4","IL2RA","HLA-C")
#' pp <- Pos.gene(gene.labs,build=37)
#' gg <- GRanges(ranges=IRanges(start=pp$start,end=pp$end),seqnames=pp$chr)
#' conv.37.38(gg) # order of output is preserved   ### HERE!!! ###
#' rr <- as(gg,"RangedData")
#' conv.37.38(rr) # note the result is same as GRanges, but in genome order
#' }
conv.37.38 <- function(ranges=NULL,chr=NULL,pos=NULL,...,ids=NULL) {
  #chain.file <- "~/github/iChip/hg19ToHg38.over.chain"
  #chain.file <- humarray::hg19ToHg38.over.chain
  #chain.file <- system.file("extdata", "hg19ToHg38.over.chain", package="humarray")
  chain.file <- humarray::hg19ToHg38
  return(conv.36.37(ranges=ranges,chr=chr,pos=pos,...,ids=ids,chain.file=chain.file))
}



#' Convert from build 36 to build 37 SNP coordinates
#' 
#' Convert range or SNP coordinates between builds using a chain file. Depending on the chain file
#' this can do any conversion, but the default will use the hg18 to hg19 (36-->37) chain file
#' built into this package. The positions to convert can be entered using using chr, pos vectors,
#'  or a RangedData or GRanges object. This function is a wrapper for liftOver() from rtracklayer,
#' providing more control of input and output and 'defensive' preservation of order and length
#' of the output versus the input ranges/SNPs. 
#' @param chr character, an optional vector of chromosomes to combine with 'pos' to describe
#'  positions to convert to an alternative build
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' a ranges object if this is provided along with 'chr'
#' @param ranges optional GRanges or RangedData object describing positions for which conversion
#' should be performed. No need to enter chr, pos if using ranges
#' @param ids if the ranges have ids (e.g, SNP ids, CNV ids), then by including this parameter
#' when using chr, pos input, the output object will have these ids as rownames. For ranges input
#' these ids would already be in the rownames of the GRanges or RangedData object, so use of
#' this parameter should be unnecessary
#' @param chain.file character, a file location for the liftOver chain file to use for the
#' conversion. If this argument is left NULL the default UCSC file that converts from hg18
#' to hg19 will be used. Can also use a 'Chain' object from rtracklayer created using
#' import.chain(). Alternate chain files for other conversions are available from
#' http://crossmap.sourceforge.net/, and you could also customize these or create your own.
#' So this function can be used for conversion between any in-out build combination, using
#' this argument, not just 36--37.
#' @param include.cols logical, whether to include any extra columns (e.g, in addition to positional
#' information) in the output object.
#' @param ... additional arguments to makeGRanges(), so in other words, can use 'start' and
#' 'end' to specify ranges instead of 'pos'.
#' @return Returns positions converted from build 36 to 37 (or equivalent for alternative chain 
#' files). If using the 'ranges' parameter for position input, the object returned will be of
#' the same format. If using chr and pos to input, then the object returned will be a data.frame
#' with columns, chr and pos with rownames 'ids'. Output will be the same length as the input,
#' which is not necessarily the case for liftOver() which does the core part of this conversion.
#' Using vector or GRanges input will give a resulting data.frame or GRanges object respectively
#' that has the same order of rownames as the original input. Using RangedData will result in an
#' output that is sorted by genome order, regardless of the original order. If ranges has no
#' rownames, or if 'ids' is blank when using chr, pos, ids of the form rngXXXX will be generated
#' in order to preserve the original ordering of locations.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{conv.37.36}}, \code{\link{conv.37.38}},
#'  \code{\link{conv.38.37}}, \code{\link{convTo37}}, \code{\link{convTo36}}
#' @references http://crossmap.sourceforge.net/
#' @examples
#' \donttest{
#' # various chain files downloadable from http://crossmap.sourceforge.net/ #
#' options(ucsc="hg18")
#' gene.labs <- c("CTLA4","IL2RA","HLA-C")
#' snp.ids <- c("rs3842724","rs9729550","rs1815606","rs114582555","rs1240708","rs6603785")
#' pp <- Pos(snp.ids); cc <- Chr(snp.ids)
#' conv.36.37(chr=cc,pos=pp,ids=snp.ids)
#' pp <- Pos(gene.labs)
#' gg <- GRanges(ranges=IRanges(start=pp$start,end=pp$end),seqnames=pp$chr)
#' conv.36.37(gg) # order of output is preserved
#' rr <- as(gg,"RangedData")
#' conv.36.37(rr) # note the result is same as GRanges, but in genome order
#' }
conv.36.37 <- function(ranges=NULL,chr=NULL,pos=NULL,...,ids=NULL,chain.file=NULL,include.cols=TRUE) {
 # require(rtracklayer); #require(genoset) #require(GenomicRanges); 
  if(!is.character(chain.file)) {
    #chain.file <- humarray::hg18ToHg19.over.chain 
    if(is(chain.file)[1]=="Chain") { 
      chn <- chain.file 
    } else {
      #chain.file <- system.file("extdata", "hg18ToHg19.over.chain", package="humarray")
      chn <- humarray::hg18ToHg19
    }
  } else {
    if(!file.exists(chain.file)) { stop("couldn't find chain file: ",chain.file) }
    chn <- import.chain(chain.file)
  }
  #toranged <- F
  outType <- is(ranges)[1]
  inr <- is.null(ranges)
  used.st.en <- all(c("start","end") %in% names(list(...)))
  if(!is.null(chr) & (!is.null(pos) | used.st.en)) {
    if(is.null(pos) & length(chr)==1) {
      if(length(list(...)$start)>1) {
        warning("when using start/end, 'chr' must have the same length as 'start'") 
      }
    }
    if(is.null(ids)) { ids <- paste0("rng",1:(max(length(chr),length(pos)))) } 
    ranges <- makeGRanges(chr=chr,pos=pos,row.names=ids,...)
    orn <- ids
  } else {
  	prv(ranges)
    if(is.null(rownames(ranges))) { rownames(ranges) <- paste0("rng",1:nrow(ranges)) }    
    orn <- rownames(ranges)
  }
  # return(ranges)
  if(is(ranges)[1]=="RangedData") { ranges <- as(ranges, "GRanges") }
  if(is(ranges)[1] %in% c("RangedData","GRanges")) {
    wd <- width(ranges)
    if(length(which(wd==1))/length(wd)>.95) { SNPs <- TRUE } else { SNPs <- FALSE }
    mcols(ranges)[["XMYINDEXX"]] <- rownames(ranges)
    mcols(ranges)[["XMYCHRXX"]] <- ocr <- chrm(ranges)
    if(include.cols) { 
      meta.data <- mcols(ranges); if(is.null(dim(meta.data))) { meta.data <- as.data.frame(meta.data) }
      rownames(meta.data) <- rownames(ranges)
      oo <- (colnames(meta.data) %in% c("XMYINDEXX","XMYCHRXX"))
      if(length(oo)>0) { meta.data <- meta.data[,-oo,drop=F] }
    }
    #prv(orn,ocr)
    opos <- start(ranges)
    ranges <- set.chr.to.char(ranges)
    #print(head(ranged))
    ranged.gr <- ranges # as(ranges,"GRanges"); #toranged <- T
  } else {
    stop("input specified resulted in an invalid GRanges/RangedData 'ranged' object, type ",is(ranges)[1]) 
  } 
  # change CHR-XY to CHR-X prior to liftOver, then change back #
#  xy.ind <- grep("XY",seqnames(ranged.gr))
  xy.ind <- grep("XY",as.character(seqnames(ranged.gr)))
  if(length(xy.ind)>0) {
    found.xy <- TRUE
    xy.id <- rownames(ranged.gr)[xy.ind]
    if(!"chrX" %in% seqlevels(ranged.gr)) { seqlevels(ranged.gr) <- c(seqlevels(ranged.gr),"chrX") }
    seqnames(ranged.gr)[xy.ind] <- "chrX"
  } else { found.xy <- FALSE }
  ranged.gr.37 <- liftOver(ranged.gr,chn)
  myfun <- function(x) { 
    data.frame(start=minna(start(x)),end=maxna(end(x))) 
  }
  if(!SNPs ) {
    new.coords.df <- do.call("rbind",lapply(ranged.gr.37,myfun))
    ranged.gr.37 <- ranged.gr
    if(!used.st.en & inr) { stop("need start and end arguments unless the dataset is all SNPs (width=1)" ) }
    ranges(ranged.gr.37) <- with(new.coords.df,IRanges(start=start,end=end))
    #seqlevels(ranged.gr.37) <- gsub("chr","",seqlevels(ranged.gr.37))
    seqlevels(ranged.gr.37) <- gsub("chr","",seqlevels(ranged.gr.37))
    out <- ranged.gr.37
  } else {
    seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    #seqnames(ranged.gr.37)<-gsub("chr","",seqnames(ranged.gr.37))
    out <- as(ranged.gr.37,"IRangesList")
    #seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    #new.coords.df <- as.data.frame(ranged.gr.37)
  }
  # seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
  # out <- as(ranged.gr.37,"IRangesList")
  out <- as(out,"RangedData")
  #return(ranged.gr.37)
  #ranged.gr.37 <- set.chr.to.numeric(ranged.gr.37)
  #if(!toranged | T) { return(ranged.gr.37) }
  ranged.gr.37 <- out #toGenomeOrder2(out)
  #return(ranged.gr.37)
  if(all(c("XMYINDEXX","XMYCHRXX") %in% colnames(ranged.gr.37))) {
    RN <- ranged.gr.37[["XMYINDEXX"]]
    nr <- nrow(ranged.gr.37)
    MAXDISPLAY <- 50
    if(length(orn)>length(RN)) { 
      cat("conversion failed for",length(orn[!orn %in% RN]),"rows, original positions kept:\n") ;  
      failz <- orn[!orn %in% RN]
      cat(comma(head(failz,MAXDISPLAY))) 
      if(length(failz)>MAXDISPLAY) { cat(", ... and",length(failz)-MAXDISPLAY,"more\n")  } else { cat("\n") }
      ln <- orn[!orn %in% RN]
      #return(ranges)
      newchr <- gsub("chr","",chrm(ranges[match(ln,ranges$XMYINDEXX),]))
      noopos <- start(ranges[match(ln,ranges$XMYINDEXX),])
      hcc <- hard.coded.conv()
      ifNAthen0 <- function(X) { X[is.na(X)] <- 0; return(X) }
      h36 <- which(noopos %in% hcc$pos36 & newchr==ifNAthen0(hcc$chr[match(noopos,hcc$pos36)])); l36 <- length(h36)
      h37 <- which(noopos %in% hcc$pos37 & newchr==ifNAthen0(hcc$chr[match(noopos,hcc$pos37)])); l37 <- length(h37)
      if(l36>=l37 & l36>0) {
        noopos[h36] <- hcc$pos37[match(noopos[h36],hcc$pos36)]
        cat("found",l36,"of the missing SNP hg18-hg19 lookups in an internal table\n")
      } else {
        if(l36<l37 & l37>0) {
          noopos[h37] <- hcc$pos36[match(noopos[h37],hcc$pos37)]
          cat("found",l37,"of the missing SNP hg19-hg18 lookups in an internal table\n")
        } else {
          ## no matches to extras table
        }
      }
      extra <- data.frame(Chr=newchr,Start=noopos,End=noopos)
      rownames(extra) <- ln
    } #else { cat("length is already the same\n") }
    Ind <- match(ranged.gr.37[["XMYINDEXX"]],orn)
    out <- data.frame(Chr=ranged.gr.37[["XMYCHRXX"]],Start=start(ranged.gr.37),End=end(ranged.gr.37),ind=Ind)
    rownames(out) <- RN
    if(length(orn)>length(RN)) {
      #prv(out,extra)
      out <- out[,-4] # 4 is the 'ind' column
      out <- rbind(out,extra)
      out <- out[orn,]
    } #else { cat("length is now the same\n") }
    #return(out) 
  } else { warning("missing key columns for chr, snp-name")  }
  #print(outType)
  #return(out)
  #prv(out)
  ranged.rd <- toGenomeOrder2(df.to.ranged(out))
  #print(colnames(ranged.rd))
  ranged.gr.37 <- as(ranged.rd,"GRanges")
  #print(colnames(ranged.gr.37))
  # if immunochip, these positions should be constant (only present in B36 really)
  if(all(c("imm_3_50875337","imm_3_50882163","imm_3_50908888") %in% rownames(ranged.gr.37))) { ranged.gr.37 <- substitute.36s(ranged.gr.37) }
  if(found.xy) {
    xy.ind <- match(xy.id,rownames(ranged.gr.37))
    lmis <- length(which(is.na(xy.ind)))
    if(lmis>0) { warning('liftOver function removed ",lmis," chrX/chrY ranges'); xy.ind <- narm(xy.ind) }
    if(!"XY" %in% seqlevels(ranged.gr.37)) { seqlevels(ranged.gr.37) <- c(seqlevels(ranged.gr.37),"XY") }
    seqnames(ranged.gr.37)[xy.ind] <- "XY"
  }
  if(outType=="GRanges") { 
    #return(ranged.gr.37)
    cn37 <- colnames(mcols(ranged.gr.37))
    if("ind" %in% cn37) { 
      mind <- as.numeric(mcols(ranged.gr.37)[["ind"]])
      #prv(ranged.gr.37,mind)
      ranged.gr.37 <- ranged.gr.37[order(mind),]
      mcols(ranged.gr.37) <- mcols(ranged.gr.37)[,-which(cn37 %in% "ind")] 
    } # else { warning("couldn't find index column, GRanges object not sorted in original order") }
    if(all(rownames(ranged.gr.37) %in% orn)) { ranged.gr.37 <- ranged.gr.37[orn,] }
    if(include.cols) {
      meta.data <- meta.data[rownames(ranged.gr.37),,drop=F]
      mcols(ranged.gr.37) <- cbind(mcols(ranged.gr.37),meta.data)
    }
    return(ranged.gr.37)
  } else {
    if(outType=="RangedData") {
      #if("ind" %in% colnames(ranged.gr.37)) { ranged.gr.37 <- ranged.gr.37[,-which(colnames(ranged.gr.37) %in% "ind")] }
      if(include.cols) {
        meta.data <- meta.data[rownames(ranged.gr.37),,drop=F]
        mcols(ranged.gr.37) <- cbind(mcols(ranged.gr.37),meta.data)
      }
      return(toGenomeOrder2(as(ranged.gr.37,"RangedData")))
    } else {
      #prv(ranged.gr.37)
      out <- ranged.to.data.frame(ranged.gr.37,include.cols=include.cols,use.names=TRUE)
      cn37 <- colnames(mcols(ranged.gr.37))
      if("ind" %in% cn37) { 
        mind <- as.numeric(mcols(ranged.gr.37)[["ind"]])
        out <- out[order(mind),,drop=FALSE]
        if("ind" %in% colnames(out)) { out <- out[,-which(colnames(out) %in% "ind")] }
      }
      if(is.null(dim(out))) { dim(out) <- c(length(out)/3,3) }
      if(all(rownames(out) %in% orn)) { out <- out[orn,] }
      return(out)
    }
  }
}




#' Extend an interval or SNP by distance in centimorgans (recombination distance)
#' 
#' It is straightforward to extend a genomic interval or position by a number of basepairs, or
#' a percentage, but extending by recombination units of centimorgans is more involved, requiring
#' annotation lookup. This function streamlines this process.
#' This function makes use of recombination rate hapmap reference files to calculate 
#' recombination distances for genome locations, in centimorgans. For a given position
#' (or vector), a window can be returned of a given extension on either side of the position,
#' for instance, 1 centimorgan to the left, and to the right of a SNP, giving a 2 centimorgan
#' range as a result. Warning - this function only uses build hg18/36, so please convert to
#' build 36 coordinates before using this function.
#' @param ranges optional GRanges or RangedData object describing positions for which we want to
#' generate windows, removing the need to enter chr, start and end
#' @param chr character, an optional vector of chromosomes to combine with 'start' and 'end'
#'  to describe positions for which to generate recombination windows
#' @param start integer, an vector of start points for chromosome ranges
#' @param end integer, an vector of end points for chromosome ranges
#' @param window numeric, number of centimorgans to extend the window either side of
#' the range or location (can be a fraction)
#' @param bp.ext numeric, optional number of base-pairs to extend the window by in addition
#' to the centimorgan extension
#' @param rec.map recombination map object (list of 22 data.frames) generated using 
#' 'get.recombination.map()'; if you are performing many of these operations, loading this 
#' object into your workspace and passing it on to this function will save loading it each 
#' time, and provide a speed advantage. Only use an object generated by get.recombination.map(),
#'  as otherwise the results will almost certainly be meaningless.
#' @param info logical, whether to display the derived window size and number of hapmap SNPs within
#' the window for each window derived
#' @seealso \code{\link{get.recombination.map}}, \code{\link{get.nearby.snp.lists}}, \code{\link{expand.nsnp}}
#' @author Chris Wallace and Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @export
#' @examples
#' \donttest{
#' # not run, as initial download of the recombination map takes nearly a minute #
#' recomWindow(chr=11,start=10000000,end=10000000,window=1,bp.ext=10000)
#' rd <- RangedData(ranges=IRanges(start=c(1.5,10.1)*10^7, end=c(1.55,10.1)*10^7),space=c(2,10))
#' rd # show original data
#' recomWindow(rd) # now extended by the interval
#' recomWindow(as(rd,"GRanges"),info=FALSE) # also works for GRanges
#' }
recomWindow <- function(ranges=NULL,chr=NA,start=NA,end=start,window=0.1,bp.ext=0, rec.map=NULL, info=TRUE) {
  #www <- window ; ccc <- chr; sss <- start; prv(ccc,sss,www,bp.ext)
  if(!is.numeric(bp.ext)) { warning("bp.ext must be numeric, setting to zero"); bp.ext <- 0 }
  if(!is.numeric(window)) { warning("window must be numeric, setting to 0.1 centimorgans"); window <- 0.1 }
  if(!all(is.na(chr))) { if(any(!paste(chr) %in% paste(1:22))) { 
    stop("this function only works for autosomes 1-22 [e.g, no X,Y or formatting like 'chr2', etc]") } }
  chr <- as.numeric(chr)
  typ <- is(ranges)[1]
  if(typ %in% c("RangedData","GRanges")) { 
    if(typ=="GRanges") { ranges <- as(ranges,"RangedData") }
    ranges <- toGenomeOrder2(ranges,strict=T)
    ss <- start(ranges); ee <- end(ranges); cc <- chr2(ranges)
    out <- recomWindow(chr=cc,start=ss,end=ee,window=window,bp.ext=bp.ext,info=info,rec.map=rec.map)
    if(length(out)==2) { out <- as.matrix(out); dim(out) <- c(1,2) } 
    outData <- RangedData(ranges=IRanges(start=out[,1],end=out[,2],names=rownames(ranges)),space=cc)
    outData <- toGenomeOrder2(outData,strict=TRUE)
    if(ncol(ranges)>0) {
      for (zz in 1:ncol(ranges)) { 
        if(typ=="GRanges") {
          outData[[colnames(ranges)[zz]]] <- mcols(ranges)[,colnames(ranges)[zz]]
        } else {
          outData[[colnames(ranges)[zz]]] <- ranges[[colnames(ranges)[zz]]]
        }
      }
    }
    if(is(ranges)[1]=="GRanges") { outData <- as(toGenomeOrder2(outData,strict=T),"GRanges") }
    return(outData)
  } else {
    if(all(!is.na(chr)) & all(!is.na(start)) & all(!is.na(end))) {
      if(length(chr)==length(start) & length(start)==length(end)) {
        if(length(chr)>1) {
          # run for a vector
          out <- matrix(ncol=2,nrow=length(chr)); colnames(out) <- c("start","end")
          for (dd in 1:length(chr)) {
            out[dd,] <- recomWindow(chr=chr[dd],start=start[dd],end=end[dd],
                                  window=window,bp.ext=bp.ext,info=info,rec.map=rec.map)
          }
          return(out)
        } else {
          ## continue as normal, just a single coordinate/range to process
        }
      } else {
        stop("invalid input, start, end and chr need to be the same length")
      }
    } else {
      stop("invalid input, either use a RangedData object, or else chr, start and end")
    }
  }
  
  #rate.fn <- sprintf("/dunwich/scratch/chrisw/HapMap/rates_rel22/genetic_map_chr%s_b36.txt.gz",chr)
  #print(rate.fn)
  #rates <- read.table(gzfile(rate.fn),header=TRUE)
  if(is.list(rec.map)) { if(length(rec.map)==22) { rates <- rec.map[[chr]] } } else {
    if(is.null(rec.map)) { rates <- get.recombination.map()[[chr]] } else {
      if(is.character(rec.map)) { rates <- get.recombination.map(dir=rec.map)[[chr]] } else {
        stop("invalid value for rec.map entered")
      }
    }
  }
  cm.st <- rates[which.min(abs(rates$position-start)),3]
  cm.en <- rates[which.min(abs(rates$position-end)),3]
  
  mx <- max(window,1)
  kk <- rates[which.min(abs(rates[,3]-(cm.st-window))) : which.min(abs(rates[,3]-(cm.en+window))),]
  if(info) { cat("n hapmap snps in window =",nrow(kk),"\n") }
  from <- min(kk[,1])
  to <- max(kk[,1])
  lft <- (start - from + bp.ext)
  rgt <- (to - end + bp.ext)
  if(lft<0) { lft <- abs(lft); from <- from-(2*(lft)) }
  if(rgt <0) { rgt <- abs(rgt); to <- to+(2*(rgt)) }
  ##
  if(info) {
    cat("new window size is\nleft: ",(start-from+bp.ext)/1000,"kb\tright: ",
      (to-end+bp.ext)/1000,"kb\ttotal: ",(to-from+(2*bp.ext))/1000,"kb\n",sep="")
  }
  if(info & bp.ext>0) { cat("in addition to cM distance, window was extended by",
                     bp.ext,"base pairs on either side\n")} 
  from <- max(c(0,(from-bp.ext)))
  to <- min(c((to+bp.ext),get.chr.lens()[chr][1]),na.rm=T)
  return(c(from,to))
}




#' Convert GRanges/RangedData to chr:pos1-pos2 vector
#' 
#' Takes a RangedData or GRanged object from some annotation lookup functions and converts to standard text
#' positions, such as what you might see on the UCSC genome browser, such as 
#' chr1:10,000,234-11,000,567 for a range, or chrX:234,432 for a SNP. Useful for printing
#' messages, concatenating positions to a single vector, or creating queries for databastes.
#' @param ranges A RangedData or GRanges object
#' @export
#' @return a text vector of the same length as 'ranges' with notation as described above
#' representing each position in the 'ranges' object
#' @seealso \code{\link{convert.textpos.to.data}}
#' @examples
#' ranged.to.txt(rranges())
ranged.to.txt <- function(ranges) {
  if(!is(ranges)[1] %in% c("RangedData","GRanges")) { stop("Not a GRanges or RangedData object") }
  text.out.a <- paste0("chr",chr2(ranges),":",format(start(ranges),scientific=F,trim=T))
  text.out.b <- paste0("-",format(end(ranges),scientific=F,trim=T))
  text.out.b[start(ranges)==end(ranges)] <- ""
  text.out <- paste0(text.out.a,text.out.b)
  return(text.out)
}






#' Select ranges only within the 22 autosomes in a ranged data object
#' 
#' Select only data from autosomes from a GRanges/RangedData object.
#' Will exclude X,Y, mitochondrial chromosome rows, and can automatically
#' detect whether chromosomes are coded as 'chr1' or just '1', etc.
#' @param ranges A RangedData or GRanges object
#' @param deselect logical, if TRUE, then will select non-autosomes
#' @export
#' @return an object of the same format as the input (ranges), except
#' with non-autosomal ranges removed.
#' @examples
#' rand.ranges <- rranges(chr.range=20:26)
#' rand.ranges # should include some non-autosomes
#' select.autosomes(rand.ranges) # only autosomes remain
select.autosomes <- function(ranges,deselect=FALSE) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) { 
    warning("not a RangedData or GRanges object"); return(ranges) 
  }
  if(length(unique(chr2(ranges))) < length(levels(chr2(ranges)))) {
    # this fixes the problem when a subset of a ranges object with less 
    #  chromosomes still has empty chr slots from previous object
    if(typ=="RangedData") { ranges <- ranges[as.numeric(unique(chr2(ranges)))] }
  } #else { cat("ok\n") }
  Chrz <- (rownames(chrInfo2(ranges)))
  chrz <- tolower(paste(Chrz))
  if(length(grep("chr",chrz))>0) {
    select1 <- which(chrz %in% paste("chr",1:22,sep=""))
    select2 <- which(chrz %in% paste(1:22))
    if(length(select2)>length(select1)) { select <- select2 } else { select <- select1 }
  } else {
    select <- which(chrz %in% paste(1:22))
  }
  if(deselect) {
    ok.chrs <- Chrz[!select]
  } else {
    ok.chrs <- Chrz[select]
  }
  if(typ=="RangedData") {
    return(ranges[ok.chrs])
  } else {
    return(ranges[chr2(ranges) %in% ok.chrs,])
  }
}






#' Convert RangedData/GRanges to a data.frame
#' 
#' Convert a RangedData/GRanges object to a data.frame with columns
#' chr, start and end. Default is to only translate the chromosome and
#' position information, which is faster. Using 'include.cols'=TRUE
#' allows all the columns from 'ranged' to be taken across to the resulting
#' data.frame.
#' @param ranged A RangedData or GRanges object
#' @param include.cols logical, whether to also bring across non-positional
#' columns to the resulting data.frame
#' @param use.names logical, whether to keep the rownames from the
#' original object for the output. Only has an effect when include.cols=FALSE,
#' otherwise original rownames are always kept.
#' @export
#' @seealso \code{\link{df.to.ranged}}, \code{\link{df.to.GRanges}}
#' @return A data.frame with columns chr, start and end, and depending on
#' chosen parameters, the same rownames as the input, and optionally the
#' same additional columns.
#' @examples
#' rd <- rranges(9,GRanges=FALSE, fakeids=TRUE)
#' rd[["fakecol"]] <- sample(nrow(rd))
#' rd[["rs.id"]] <- paste0("rs",sample(10000,9))
#' ranged.to.data.frame(rd)
#' ranged.to.data.frame(rd,,FALSE)
#' ranged.to.data.frame(rd,TRUE) # keep all the columns
#' df.to.GRanges(ranged.to.data.frame(rd,TRUE)) # inverse returns original
ranged.to.data.frame <- function(ranged,include.cols=FALSE,use.names=TRUE) {
  if(!include.cols) {
    u <- ranged.to.txt(ranged)
    v <- convert.textpos.to.data(u)
    if(!is.null(rownames(ranged)) & nrow(ranged)==nrow(v) & use.names) { rownames(v) <- rownames(ranged) }
    return(v)
  } else {
    u <- as.data.frame(ranged)
    cn <- tolower(colnames(u))
    if(is(ranged)[1]=="RangedData") {
      if("names" %in% cn) { 
        rownames(u) <- u[["names"]]
        u <- u[,-which(cn=="names")]
      } 
      if("space" %in% cn) { colnames(u)[which(cn=="space")] <- "chr" }
    } else {
      if(is(ranged)[1]=="GRanges") {
        if("seqnames" %in% cn) { colnames(u)[which(cn=="seqnames")] <- "chr" }
      } else {
        warning("'ranged' should be RangedData or GRanges, coercion could fail")
      }
    }
    return(u)
  }
}

#' Convert a data.frame with positional information to GRanges
#' 
#' Convert a data.frame containing chromosome and position information
#' to a GRanges object. Assumes the position information is contained in
#' columns named 'chr', 'start' and 'end' respectively (not case sensitive) 
#' although you can enter alternative column names for each as parameters. 
#' 'seqnames' will be automatically detected as an alternative to 'chr' if 
#' present. Column names that are default GRanges slot names such as 'seqnames',
#' 'ranges', 'strand', 'seqlevels', etc, will be removed during conversion, so
#' rename these if you want them to be translated into the resulting GRanges
#' objects' column metadata. If there is a column 'pos' but no columns 'start'
#'  and 'end' this will be detected automatically without needing to change
#'  the default parameters and start will equal end equals pos (ie., SNPs).
#' @param dat a data.frame with chromosome and position information 
#' @param ... additional arguments to df.to.ranged(), namely:
#' ids, start, end, width, chr, exclude and build
#' @export
#' @seealso \code{\link{ranged.to.data.frame}}, \code{\link{df.to.ranged}}
#' @return A RangedData or GRanges object. If 'dat' doesn't
#' use the default column names, specify these using parameters
#' ids, start, and end or width. Exclude will remove prevent any 
#' column names of 'dat' specified not to be translated to the 
#' returned GRanges object. 'build' specifies the 'genome'
#' slot of the resulting object. 'ids' allows specification of
#' a column to be converted to the rownames of the new object.
#' @examples
#' chr <- sample(1:22,10)
#' start <- end <- sample(1000000,10)
#' df1 <- cbind(chr,start,end)
#' df.to.GRanges(df1) # basic conversion
#' width <- rep(0,10)
#' df2 <- cbind(chr,start,width)
#' df.to.GRanges(df2,end=NULL,width="width") # define ranges with start and width
#' id.col <- paste0("ID",1:10)
#' rs.id <- paste0("rs",sample(10000,10))
#' df3 <- cbind(chr,start,end,id.col,rs.id)
#' df.to.GRanges(df3) # additional columns kept
#' df4 <- cbind(chr,start,end,id.col,rs.id, ranges=1:10)
#' df.to.GRanges(df4) # 'ranges' column excluded as illegal name
#' df.to.GRanges(df4, exclude="rs.id") # manually exclude column
#' df5 <- cbind(chr,start,end,rs.id)
#' rownames(df5) <- paste0("ID",1:10)
#' df.to.GRanges(df5) # rownames are kept
#' df.to.GRanges(df4,ids="id.col") # use column of 'dat' for rownames
df.to.GRanges <- function(dat,...) {
  return(df.to.ranged(dat=dat,...,GRanges=TRUE))
}


#' Convert a data.frame with positional information to RangedData/GRanges
#' 
#' Convert a data.frame containing chromosome and position information
#' to a RangedData or GRanges object. Assumes the position information is contained in
#' columns named 'chr', 'start' and 'end' respectively (not case sensitive) 
#' although you can enter alternative column names for each as parameters. 
#' 'seqnames' will be automatically detected as an alternative to 'chr' if 
#' present. If there is a column 'pos' but no columns 'start' and 'end' this
#' will be detected automatically without needing to change the default parameters
#' and start will equal end equals pos (ie., SNPs). Column names that are default 
#' GRanges slot names such as 'seqnames', 'ranges', 'strand', 'seqlevels', etc, will
#' be removed during conversion, so rename these if you want them to be translated 
#' into the resulting object.
#' @param dat a data.frame with chromosome and position information 
#' @param ids character string, an optional column name containing ids which
#' will be used for rownames in the new object, as long as the ids are unique.
#' If not, this option is overridden and the ids will simply be a normal column
#' in the new object.
#' @param start character, the name of a column in the data.frame contain
#' the start point of each range. Not case sensitive. In the case of SNP
#' data, a column called 'pos' will also be automatically detected without
#' modifying 'start' or 'end', and will be used for both start and end.
#' @param end character, the name of a column in the data.frame containing the
#' end point of each range, can also use 'width' as an alternative specifier,
#' in which case 'end' should be set to NULL. Not case sensitive. In the case of SNP
#' data, a column called 'pos' will also be automatically detected without
#' modifying 'start' or 'end', and will be used for both start and end.
#' @param width the name of a column in the data.frame containing 'width' of
#' ranges, e.g, SNPs would be width=0. This is optional, with 'start' and 'end'
#' being the default way to specify an interval. If using 'width' you must
#' also set 'end' to NULL.  Not case sensitive.
#' @param chr character, the name of the column in the data.frame containing
#' chromosome values. The default is 'chr' but 'seqnames' will also be
#' detected automatically even when chr='chr'. Not case sensitive.
#' @param exclude character string, and column names from the data.frame to 
#' NOT include in the resulting S4 object.
#' @param build the ucsc build for the result object which will apply to the
#' 'universe' (RangedData) or 'genome' slot (GRanges) of the new object.
#' @param GRanges logical, whether the resulting object should be GRanges (TRUE),
#' or RangedData (FALSE)
#' @param fill.missing logical, GRanges/RangedData objects cannot handle missing
#' chrs/positions, so if fill missing is selected, will insert values of chr99, and
#' start=end=1, and if FALSE, will exclude any row with a missing value from the
#' resulting object.
#' @export
#' @seealso \code{\link{ranged.to.data.frame}}, \code{\link{df.to.ranged}}
#' @return A RangedData or GRanges object. If 'dat' doesn't use the default 
#' column names 'chr', 'start'/'end' or 'pos', specify these using parameters 
#' 'ids', 'start', and 'end' or 'width'. Exclude will remove prevent any 
#' column names of 'dat' specified not to be translated to the returned GRanges
#' object. 'build' specifies the 'genome' slot of the resulting object. 'ids' 
#' allows specification of a column to be converted to the rownames of the new object.
#' @examples
#' chr <- sample(1:22,10)
#' start <- end <- sample(1000000,10)
#' df1 <- cbind(CHR=chr,Start=start,enD=end)
#' print(df1)
#' df.to.GRanges(df1) # not case sensitive!
#' width <- rep(0,10)
#' df2 <- cbind(chr,start,width)
#' df.to.GRanges(df2,end=NULL,width="width") # define ranges with start and width
#' id.col <- paste0("ID",1:10)
#' rs.id <- paste0("rs",sample(10000,10))
#' df3 <- cbind(chr,start,end,id.col,rs.id)
#' df.to.GRanges(df3) # additional columns kept
#' df4 <- cbind(chr,start,end,id.col,rs.id, ranges=1:10)
#' df.to.GRanges(df4) # 'ranges' column excluded as illegal name
#' df.to.GRanges(df4, exclude="rs.id") # manually exclude column
#' df5 <- cbind(chr,start,end,rs.id)
#' rownames(df5) <- paste0("ID",1:10)
#' df.to.GRanges(df5) # rownames are kept
#' df.to.GRanges(df4,ids="id.col") # use column of 'dat' for rownames
df.to.ranged <- function(dat, ids=NULL,start="start",end="end",width=NULL,
                                 chr="chr",exclude=NULL,build=NULL,GRanges=FALSE,
                                 fill.missing=TRUE) 
{
  ## abandon longer names as they clash with function names
  st <- paste(start); en <- paste(end); ch <- paste(chr); wd <- paste(width)
  if((!chr %in% colnames(dat)) & ("seqnames" %in% colnames(dat)) & GRanges) { ch <- "seqnames" }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  #must.use.package(c("genoset","IRanges"),T)
  g.illegal <- tolower(c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
                         "isCircular", "start", "end", "width", "element"))
  if(is.matrix(dat)) { dat <- as.data.frame(dat,stringsAsFactors=FALSE) }
  if(!is.data.frame(dat)) { stop("Error: not a dataframe")}
  key.nms <- c(ids,st,en,ch,wd)
  tries <- 0
  #print(key.nms); print(colnames(dat))
  if(st=="position" & is.null(en)) { en <- st } # if only 1 entered
  while(!all(key.nms %in% colnames(dat))) { 
    colnames(dat) <- tolower(colnames(dat)); key.nms <- tolower(key.nms)
    st <- tolower(st); en <- tolower(en); ch <- tolower(ch); wd <- tolower(wd)
    if(tries>2) {
      if((tolower(st)=="pos" | tolower(en)=="pos") & !(tolower(st)=="pos" & tolower(en)=="pos")) {
        st <- en <- "pos"
      } else {
        if(tolower(st)=="start" & tolower(en)=="end") { st <- en <- "pos" }
      }
    }
    key.nms <- c(ids,st,en,ch,wd)
    tries <- tries+1
    if(tries > 3) { if(!all(c(st,en,ch) %in% colnames(dat))) {
      warning("chromosome and position columns not found") } ; break }
  }
  if(!is.null(ids)) { 
    if(anyDuplicated(dat[[ids]])==0) { 
      id <- dat[[ids]] 
    } else { 
      key.nms <- key.nms[-match(ids,key.nms)] # allow non-unique ids as regular
      ids <- NULL
      warning("id must be unique to form rownames, will insert as a separate column") 
    }
  }
  if(is.null(ids)) { 
    if(!is.null(rownames(dat)) & all(rownames(dat)!=paste(1:nrow(dat)))) { 
      id <- rownames(dat)
    } else { 
      id <- paste(1:nrow(dat)) 
    }
  }
  ## not sure why here are adding 'chr' to X and Y?
  #this was here before? :  if(length(ch)>0) { ch1 <- gsub("Y","chrY",gsub("X","chrX",gsub("chr","",dat[[ch]],ignore.case=T))) } else { ch1 <- NULL }
  if(length(ch)>0) { ch1 <- gsub("chr","",dat[[ch]],ignore.case=T) } else { ch1 <- NULL }
  if(length(st)>0) { st1 <- as.numeric(dat[[st]]) } else { st1 <- NULL }
  if(length(en)>0) { en1 <- as.numeric(dat[[en]]) } else { en1 <- NULL }
  if(length(wd)>0) { en1 <- st1+as.numeric(dat[[wd]]) } # { en1 <- st1+dat[[wd]] }
  misser <- (is.na(ch1) | is.na(en1) | is.na(st1))
  if(any(misser)) {
    if(fill.missing) {
      ch1[is.na(ch1)] <- 99; st1[is.na(st1)] <- 1; en1[is.na(en1)] <- 1
    } else {
      ch1 <- ch1[!misser]; st1 <- st1[!misser]; en1 <- en1[!misser]; id <- id[!misser]
    }
  }
  #print(length(st1)); print(length(en1)); print(length(id)); print(length(ch1))
  outData <- GRanges(ranges=IRanges(start=st1,end=en1,names=id),seqnames=ch1); genome(outData) <- build[1]
  #outData <- RangedData(ranges=IRanges(start=st1,end=en1,names=id),space=ch1,universe=build[1])
  ###  ###  ###  outData <- toGenomeOrder2(outData,strict=T)
  # note when adding data subsequently that 'RangedData' sorts by genome order, so need
  # to re-sort any new data before adding.
  if(is.null(rownames(outData))) { rownames(outData) <- paste(1:nrow(outData)) }
  reorder <- match(rownames(outData),id)
  more.cols <- colnames(dat)[!colnames(dat) %in% key.nms]
  more.cols <- more.cols[!more.cols %in% exclude]
  if(is(outData)[1]=="GRanges") { more.cols <- more.cols[!more.cols %in% g.illegal] }
  if(length(more.cols)>0) {
    for (cc in 1:length(more.cols)) {
      u <- dat[[more.cols[cc]]][reorder]; #prv(u)
      if(is(outData)[1]=="GRanges") {
        mcols(outData)[[more.cols[cc]]] <- u
      } else {
        outData[[more.cols[cc]]] <- u
      }
    }
  }
  if(GRanges) {
    return(as(outData,"GRanges"))
  } else {
    cncn <- colnames(mcols(outData))
    outData <- as(outData,"RangedData")
    if(any(cncn %in% "strand")) {
      outData <- outData[,-which(cncn=="strand")]
    }
    #outData <- toGenomeOrder2(outData,strict=T)
    return(outData)
  }
}



#' Select chromosome subset of GRanges or RangedData object
#' 
#' One of the main differences between RangedData and GRanges is the way
#' of selecting the subset for a chromosome. RangedData just uses [n] where
#' 'n' is the chromosome name or number. Whereas GRanges, does not have a
#' method like this, so need to select using [chr(X)==chr.num,]
#' This wrapper allows selection of a chromosome or chromosomes regardless of
#' whether the object is RangedData or GRanges type.
#' @param X A GRanges or RangedData object
#' @param chr Vector, the chromosome(s) (number(s) or name(s)) to select
#' @param index logical, if FALSE, will assume 'chr' is a string, indicating the
#' chromosome name, if TRUE, if 'chr' is numeric, will assume it refers to the
#' chromosome index, which if there are some chromosomes not represented, may
#' be different to the name. E.g, an object with data for chromosomes 1,2,4,5
#' would select chromosome 5 with chr=4, if index=TRUE.
#' @export
#' @return returns an object of the same type as X, with only the chromosome
#' subset specified.
#' @examples
#' some.ranges <- rranges(100,chr.range=1:10)
#' chrSelect(some.ranges,6)
#' more.ranges <- rranges(10, chr.range=21:25)
#' chrSelect(more.ranges,1:22) # gives warning
#' select.autosomes(more.ranges)
chrSelect <- function(X,chr,index=FALSE) {
  typ <- is(X)[1]
  if(!typ %in% c("RangedData","GRanges","ChipInfo")) { stop("not a ChipInfo, GRanges or RangedData object") }
  if(nrow(X)==0) { warning("X has no ranges") ; return(X) }
  if(!(is.character(chr) | is.numeric(chr))) { stop("chr must be character or numeric type") }
  if(is.numeric(chr)) { if(!all(chr %in% 1:99)) { 
    stop("illegal chromosome index, valid range 1-99 [although 1-28 typical for human]") } }
  if(!any(paste(chr) %in% paste(unique(chrm(X))))) { stop("X did not have any chromosome from the list: ",paste(chr,collapse=",")) }
  if(typ=="RangedData") { if(index) { return(X[chr]) } else { return(X[paste(chr)]) } }
  all.chr <- chr2(X)
  if(!all(chr %in% unique(all.chr))) { 
    if(!any(chr %in% unique(all.chr))) { 
      warning("none of the specified chromosome indices were present in the GRanges object, returning NULL")
      return(NULL)
    } else { 
      warning("some of the specified chromosome indices were not present in the GRanges object") 
    }
  }
  return(X[all.chr %in% chr,])
}


#' Simulate a GRanges or RangedData object
#' 
#' For testing purposes, this function will generate a S4 ranged object
#' based on the human genome. The default is to produce ranges selected
#' from chromosomes, with probability of a position in each chromosome
#' equal to the length of that chromosome versus the whole genome. The
#' maximum position allocated within each chromosome will be within
#' the length bounds of that chromosome. You can specify SNPs (ie., start
#' =end), but the default is for random ranges. You can alter the UCSC
#' build to base the chromosome lengths on, and you can specify whether
#' chromosomes should appear as chr1,chr2,... versus 1,2,..
#' @param n integer, number of rows to simulate
#' @param SNP logical, whether to simulate SNPs (width 1, when SNPs=TRUE)
#'  or just ranges (when SNP=FALSE)
#' @param chr.range integer vector of values from 1 to 26, to specify which
#' chromosomes to include in the simulated object. 23-26 are X,Y,XY,MT 
#' respectively.
#' @param chr.pref logical, if TRUE chromosomes will be coded as chr1,chr2,...,
#' versus 1,2,.. when chr.pref=FALSE
#' @param order logical, if TRUE the object returned will be in genomic order,
#' otherwise the order will be randomized
#' @param equal.prob logical, when FALSE (default), random positions will be
#' selected on chromosomes chosen randomly according to the their length (i.e,
#' assuming every point on the genome has equal probability of being chosen.
#' If equal.prob=TRUE, then chromosomes will be selected with equal probability,
#' so you could expect just as many MT (mitochondrial) entries as Chr1 entries.
#' @param GRanges logical, if TRUE the returned object will be GRanges format,
#' or if FALSE, then RangedData format
#' @param build character, to specify the UCSC version to use, which has a small
#' effect on the chromosome lengths. Use either "hg18" or "hg19". Will also 
#' accept build number, e.g, 36 or 37.
#' @param fakeids logical, whether to add rownames with random IDs (TRUE) or 
#' leave rownames blank (FALSE). If SNP=TRUE, then ids will be fake rs-ids.
#' @export
#' @return returns a ranged object (GRanges or RangedData) containing data
#' for 'n' simulated genomic ranges, such as SNPs or CNVs across chromosomes in
#' 'chr.range', using UCSC 'build'.
#' @examples
#' rranges()
#' rr <- rranges(SNP=TRUE,chr.pref=TRUE,fakeids=TRUE)
#' width(rr) # note all have width 1
#' rr
#' tt <- table(chrm(rranges(1000)))
#' print(tt/sum(tt)) # shows frequencies at which the chr's were sampled
#' tt <- table(chrm(rranges(1000,equal.prob=TRUE)))
#' print(tt/sum(tt)) # shows frequencies at which the chr's were sampled
rranges <- function(n=10,SNP=FALSE,chr.range=1:26,chr.pref=FALSE,order=TRUE,equal.prob=FALSE,
                    GRanges=TRUE,build=NULL, fakeids=FALSE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(!is.numeric(chr.range)) { stop("chr.range must be a numeric integer vector rangeing from 1 to 26") }
  chr.range <- unique(chr.range[chr.range<=26 & chr.range>=1]) 
  cL <- get.chr.lens(mito=TRUE,build=build)[c(1:24,24,25)]
  if(equal.prob) {
    cP <- rep(1/length(chr.range),length(chr.range))
  } else {
    cP <- cL[chr.range]/sum(cL[chr.range]) # probabilities of a location being in each chromosome
  }
  if(!is.numeric(n)) { stop("'n' must be a numeric integer vector, representing the number of rows to simulate") }
  nn <- round(force.scalar(n,min=1,max=10^9,default=10))
  chrs <- sort(sample(chr.range,size=nn,replace=TRUE,prob=cP))
  cnts <- table(chrs)
  chr.lengths <- get.chr.lens(mito=TRUE)[c(1:24,24,25)]
  dubb <- if(SNP) { 1 } else { 2 }
  randoms <- starts <- ends <- vector("list",length(chr.range))
  for (cc in 1:length(cnts)){
    ii <- match(names(cnts[cc]),paste(chr.range))
    if(!is.na(ii)) {
      randoms[[cc]] <- oo <- sort(sample(chr.lengths[chr.range[ii]],size=cnts[cc]*dubb,replace=TRUE))
    } else { randoms[[cc]] <- NA }
    #prv(oo)
    if(!SNP) { 
      kk <- length(randoms[[cc]])/2
      #print(kk)
      starts[[cc]] <- randoms[[cc]][1+(2*(0:(kk-1)))]
      ends[[cc]] <- randoms[[cc]][2*(1:kk)]
    } 
  }
  if(SNP) { starts <- ends <- randoms }
  starts <- unlist(starts); ends <- unlist(ends)
  chrs <- paste0(if(chr.pref) { "chr" } else { "" },chrs)
  chrs <- gsub("23","X",chrs);  chrs <- gsub("24","Y",chrs)
  chrs <- gsub("25","XY",chrs);  chrs <- gsub("26","MT",chrs)
  gg <- GRanges(ranges=IRanges(start=starts,end=ends),seqnames=chrs)
  if(!order) {
    gg <- gg[order(rnorm(nrow(gg))),]
  } else {
    gg <- toGenomeOrder2(gg,strict=TRUE)
  }
  if(!GRanges) { gg <- as(gg,"RangedData"); gg <- gg[,-1] }
  if(n==0) { return(gg[-1,])} # returns empty ranges if that's what you really want
  if(fakeids) {
    if(SNP) { rownames(gg) <- rsnpid(nrow(gg)) } else { rownames(gg)  <- rsampid(nrow(gg)) }
  }
  return(gg)
}



#' Extract chromosome numbers from GRanges/RangedData 
#' 
#' Sometimes chromosomes are codeds as 1:22, sometimes there is also X,Y, etc, sometimes it's 
#' chr1, ch2, etc. This function extracts the set of chromosome labels used by a ranged object 
#' (ie, GRanges or RangedData) and converts the labels to numbers in a consistent way, so
#' 1:22, X, Y, XT, MT ==> 1:26, and optionally you can output the conversion table of codes to
#' numbers, then input this table for future conversions to ensure consistency.
#' @param ranged GRanges or RangedData object
#' @param warn logical, whether to display a warning when non autosomes are converted to numbers
#' @param table.out logical, whether to return a lookup table of how names matched to integers
#' @param table.in data.frame/matrix, col 1 is the raw text names, col 2 is the integer that should be assigned,
#'  col 3 is the cleaned text (of col 1) with 'chr' removed. the required form is outputted by this function if
#'  you set 'table.out=TRUE', so the idea is that to standardize coding amongst several RangedData objects you
#'  can save the table each time and ensure future coding is consistent with this. Note that chromosomes 1-22, X,
#'  Y, XY, and MT are always allocated the same integer, so table is only useful where there are extra NT, COX, HLA
#'  regions, etc.
#' @return a set of integers of length equal to the number of unique chromosomes in the ranged data.
#' @export
#' @examples
#' require(genoset)
#' gg <- rranges(1000)
#' chrNames(gg); chrNums(gg)
#' gg <- rranges(1000,chr.pref=TRUE) # example where chromosomes are chr1, chr2, ...
#' chrNames(gg); chrNums(gg)
#' lookup <- chrNums(gg,table.out=TRUE)
#' lookup
#' gg2 <- rranges(10)
#' chrNums(gg2,table.in=lookup) # make chromosome numbers using same table as above
chrNums <- function(ranged,warn=FALSE,table.out=FALSE,table.in=NULL) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranged)[1]
  if(!typ %in% c("RangedData","GRanges","ChipInfo")) { warning("not a GRanges, ChipInfo or RangedData object"); return(NULL) }
  lookup <- c("X","Y","XY","MT")
  txt1 <- chrNames2(ranged)
  txt <- gsub("chr","",txt1,fixed=T)
  nums <- suppressWarnings(as.numeric(txt))
  num.na <- length(nums[is.na(nums)])
  if(num.na>0) { 
    if(warn) { warning(paste("chromosome numbers requested for non-autosomes, will assign numbers >=23 to letters",
                             paste(txt[is.na(nums)],collapse=","))) }
    aux.ind <- match(txt,lookup)
    nums[!is.na(aux.ind)] <- 22+aux.ind[!is.na(aux.ind)]
    unmatched <- txt[is.na(nums)]
    if(!is.null(table.in)) {
      if((all(table.in[,1] %in% unmatched)) | (all(unmatched %in% table.in[,1]))) {
        if(all(unmatched %in% table.in[,1])) {
          out <- table.in[,2][match(unmatched,table.in[,1])]
          nums[is.na(nums)] <- as.numeric(out)
        } else {
          out <- table.in[,2][match(unmatched,table.in[,1])]
          nums[is.na(nums)][!is.na(out)] <- as.numeric(out)[!is.na(out)]
          st.num <- max(c(22+length(lookup),as.numeric(table.in[,2])),na.rm=T)+1
          nums[is.na(nums)][is.na(out)] <- st.num:(st.num+length(nums[is.na(nums)][is.na(out)])-1)
        }
      } else {
        out <- table.in[,2][match(unmatched,table.in[,1])]
        nums[is.na(nums)][!is.na(out)] <- as.numeric(out)[!is.na(out)]
        st.num <- max(c(22+length(lookup),as.numeric(table.in[,2])),na.rm=T)+1
        nums[is.na(nums)][is.na(out)] <- st.num:(st.num+length(nums[is.na(nums)][is.na(out)])-1)
      }
    } else {
      nums[is.na(nums)] <- 27:(27+length(nums[is.na(nums)])-1)
    }
  }
  if(table.out) {
    out <- cbind(txt1,nums,txt)
    return(out)
  } else {
    return(sortna(as.numeric(nums)))
  }
}



#' Expand genomic locations to the ranges covering the 'n' closest SNPs
#' 
#' Sometimes for chip data we want to create windows around some locus, and
#' fixed distance [see flank()], recombination distance [see recomWindow()] or a number of SNPs 
#' might be used. This function allows expansion of regions according to a set number of SNPs.
#' The result gives two regions for each row of a GRanges or RangedData object describing
#' the start and end of the left flanking 'nsnp' region, and right flanking 'nsnp' region
#' respectively.
#' @param ranged a GRanges or RangedData object describing the locations for
#' which we want to find regions encompassing 'nsnps' closest SNPs.
#' @param snp.info An object of type: ChipInfo, RangedData or GRanges, describing the set of SNPs
#' you are using (e.g, chip annotation). If left as null the ChipInfo object from chip.support() 
#' with default options() will be used
#' @param nsnp Number of nearest SNPs to return for each location
#' @param add.chr logical, whether to add a chromosome column for the output object
#' @seealso \code{\link{nearest.snp}}, \code{\link{chip.support}}, \code{\link{recomWindow}}
#' @export
#' @return Two regions for each row of a the 'ranged' object describing
#' the start and end of the left flanking 'nsnp' region, and right flanking 'nsnp' region
#' respectively. If 'ranged' has rownames these should stay in the same order in the resulting
#' object. Chromosome will be the final column if you set add.chr=TRUE.
#' @examples
#' rngs <- rranges()
#' # not run - slow ~5 seconds # expand.nsnp(rngs)
#' # not run - slow ~5 seconds # expand.nsnp(rngs,add.chr=TRUE)
expand.nsnp <- function(ranged,snp.info=NULL,nsnp=10, add.chr=FALSE) {
  if(is.null(snp.info)) { snp.info <- chip.support() }
  if(!is(snp.info)[1] %in% c("ChipInfo","RangedData","GRanges")) { 
    stop("snp.info must be of type: ChipInfo, RangedData or GRanges") }
  snp.info <- toGenomeOrder2(snp.info,strict=TRUE); rw.cnt <- 1
  all.fl <- matrix(ncol=4+(as.numeric(add.chr)), nrow=0)
  for(cc in chrNums(ranged)) {
    si <- chrSel(snp.info,paste(cc))
    nxt.nm <- rownames(si); pos <- start(si)
    rng <- chrSel(ranged,paste(cc)) # ranged[paste(cc)]
    st.en.snp <- rangeSnp(ranged=rng,snp.info=si)
    fl <- matrix(ncol=4, nrow=nrow(st.en.snp))
    fl[,2] <- start(rng); fl[,3] <- end(rng);
    for(dd in 1:nrow(st.en.snp)) {
      x1 <- pos[max(1,match(st.en.snp[dd,1],nxt.nm)-nsnp)]
      x2 <- pos[min(length(nxt.nm),match(st.en.snp[dd,2],nxt.nm)+nsnp)]
      #print(x1); print(x2); print(length(x1)); print(length(x2))
      fl[dd,1] <- x1[1]
      fl[dd,4] <- x2[1]
    }
    if(add.chr) { 
      chrz <- cc
      fl <- cbind(fl,chrz)
    }
    all.fl <- rbind(all.fl,fl)
  }
  fl <- (all.fl)
  fl[fl[,1]>fl[,2],1] <- fl[fl[,1]>fl[,2],2]
  fl[fl[,3]>fl[,4],1] <- fl[fl[,3]>fl[,4],4]
  cN <- c("left.start","left.end","right.start","right.end")
  if(add.chr) { cN <- c(cN,"chr") }
  colnames(fl) <- cN
  if(nrow(fl)==nrow(ranged) & !is.null(rownames(ranged))) { rownames(fl) <- rownames(ranged) }
  return(fl)
}




#' Find closest SNPs to the ends of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the end of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' ends of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by chip.support()
#' which will depend on your options() settings, see ?chip.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the ends of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the ends of.
#' @param nearest will preferably find an exact match but if nearest=TRUE, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output vector (character)
#' should be the same length as the number of ranges entered.
#' @examples
#' endSnp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' endSnp(rranges())
endSnp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(startSnp(ranged=ranged,snp.info=snp.info,chr=chr,pos=pos,start=F,end=T,nearest=nearest))
}


#' Find closest SNPs to the starts and ends of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the starts and ends
#' of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' starts/ends of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by chip.support()
#' which will depend on your options() settings, see ?chip.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the starts/ends of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the starts/ends of.
#' @param nearest will preferably find an exact match but if nearest=TRUE, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output will be a matrix
#' which should have the same number of rows as the number of ranges entered.
#' @examples
#' rangeSnp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' rangeSnp(rranges())
rangeSnp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(startSnp(ranged=ranged,snp.info=snp.info,chr=chr,pos=pos,start=T,end=T,nearest=nearest))
}


#' Find closest SNPs to the starts of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the starts 
#' of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' starts of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by chip.support()
#' which will depend on your options() settings, see ?chip.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the starts of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the starts of.
#' @param start logical whether to return the SNP nearest the range starts
#' @param end logical whether to return the SNP nearest the range ends
#' @param nearest will preferably find an exact match but if nearest=TRUE, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output will be a vector
#' which will have the same length as the input. Unless start=TRUE and end=TRUE, then will return a matrix
#' which should have the same number of rows as the number of ranges entered. Note that endSnp() is 
#' equivalent to using this function when end=TRUE and start=FALSE, and rangeSnp() is the same as setting
#' start=TRUE and end=TRUE.
#' @examples
#' startSnp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' startSnp(rranges())
startSnp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,start=T,end=F,nearest=T) {
  # will preferably find an exact match but if nearest=T, will fall-back on nearest match
  #must.use.package("genoset",T)
  nmz <- NULL
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) {
    if(!is.null(chr) & !is.null(pos)) {
      if(is.null(dim(pos))) { st <- pos[1]; en <- pos[2] } else {
        st <- pos[,1]; en <- pos[,2]
      }
      if(length(st)>length(chr)) { chr <- rep(chr[1],length(st)) } else { chr <- chr[1:length(st)] }
    } else {
      stop("if not using 'ranged' input, then chr and pos must be valid")
    }
  } else {
    st <- start(ranged); en <- end(ranged); chr <- chr2(ranged)
    nmz <- rownames(ranged)
  }
  if(is.null(snp.info)) { snp.info <- chip.support() }  # load default chip
  if(!is(snp.info)[1] %in% c("ChipInfo","RangedData","GRanges")) {
    stop("snp.info must be of type ChipInfo, RangedData or GRanges")
  } else {
    if(is.null(rownames(snp.info))) {  rownames(snp.info) <- paste(1:nrow(snp.info)) }
  }
  st.snps <- en.snps <- character(length(chr)) ; prch <- 0
  for (cc in 1:length(chr)) {
    if(chr[cc]!=prch) { 
      ref <- chrSel(snp.info,paste(chr[cc])) # snp.info[paste(chr[cc])]
      st.ref <- start(ref); rnref <- rownames(ref)
      if(is.null(ref)) { stop(paste("snp.info did not contain chr",chr[cc])) }
    }
    #exact
    if(start) { 
      ind <- match(st[cc],st.ref)
      if(any(is.na(ind)) & nearest) {
        difs <- abs(st[cc]-st.ref)
        ind <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind)==0) { st.snps[cc] <- NA } else {
        st.snps[cc] <- rnref[ind]
      }
    }
    if(end){
      ind2 <- match(en[cc],st.ref)
      if(any(is.na(ind2)) & nearest) {
        difs <- abs(en[cc]-st.ref)
        ind2 <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind2)==0) { en.snps[cc] <- NA } else {
        en.snps[cc] <- rnref[ind2]
      }
    }
    prch <- chr[cc]
  }
  if(start & !end) {
    return(st.snps)
  }
  if(!start & end) {
    return(en.snps)
  }
  #otherwise looks like want both
  out <- cbind(st.snps,en.snps)
  if(!is.null(nmz)) { if(length(nmz)==nrow(out)) { rownames(out) <- nmz } }
  return(out)
}


#' Force a valid genomic range, given the inputted coordinates
#'
#' Enter a pair of genomic locations representing a range for a given chromosome and this
#' function will ensure that no position is less than 1 or greater than the relevant chromosome
#'  lengths. Anything below will be coerced to 1, and anything above to the chromosome length.
#' @param Pos must be numeric, length 2, e.g, c(20321,30123)
#' @param Chr chromosome label
#' @param snp.info optional object to take boundaries from, the maxima and minima for each
#' chromosome within this object will take the place of the chromsome lengths / 1.
#' @param build ucsc build, only need to enter if this differs from getOption("ucsc")
#' @param dir directory to use for download of chromosome lengths (only if you wish to
#' keep the chromosome length file)
#' @export
#' @examples
#' pss <- ps <- c(345035,345035); ch <- 1
#' force.chr.pos(ps,ch)
#' pss[1] <- 0
#' force.chr.pos(pss,ch) # won't allow zero
#' pss[1] <- -1
#' force.chr.pos(pss,ch) # won't allow negative
#' pss[1] <- 645035012
#' force.chr.pos(pss,ch) # won't allow pos > chromosome length
force.chr.pos <- function(Pos,Chr,snp.info=NULL,build=NULL,dir=NULL) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  # convert any non autosomes to numbers:
  Chr <- paste(Chr)
  Chr[grep("c6",Chr,ignore.case=T)] <- 6  # prevent issues with c6_COX, c6_QBL  
  Chr[grep("X",Chr,ignore.case=T)] <- 23
  Chr[grep("Y",Chr,ignore.case=T)] <- 24
  Chr[grep("M",Chr,ignore.case=T)] <- 25
  Chr[grep("NT",Chr,ignore.case=T)] <- 26  # prevent issues with NT_11387, etc
  Chr <- as.numeric(Chr)
  if(any(!paste(Chr) %in% paste(c(1:26)))) { stop("invalid chromosome(s) entered") }
  if(any(paste(Chr) == paste(26))) { warning("'NT' chromosome(s) entered, not supported, NAs produced") }
  if(length(Pos)==2 & is.numeric(Pos)) {
    if(is(snp.info)[1]!="RangedData" & is(snp.info)[1]!="GRanges") { 
      maxln <- get.chr.lens(dir=dir,mito=T,autosomes=FALSE,build=build)[Chr] 
    } else { 
      maxln <- end(tail(snp.info[paste(Chr)],1)) # force start and end to be within 1:chr.len
    }
    mbs <- min(max(1,Pos[1]),(maxln-1)); mbe <- min(max(2,Pos[2]),maxln)
    return(c(mbs,mbe))
  } else {
    Pos <- NA; warning("Pos needs to be numeric length 2, min, max")
  }
  return(round(Pos))
}







#' Select all ranges lying within a chromosome window
#' 
#' Input a ranged object (ie., GRanges or RangedData) and this function will
#' return the subset from chromosome 'chr' and within the base-pair range specified
#' by 'pos', in units of 'unit'. By default ranges with ANY overlap are returned, but
#' it can be specified that it must be full overlap. Duplicates can be removed.
#' @param ranged GRanges or RangedData object
#' @param chr a chromosome, e.g, 1,2,3,...,22,X,Y,XY,MT or however chromosomes are 
#' annotated in 'ranged'
#' @param pos a numeric range (length 2), with a start (minima) and end (maxima), specifying
#' the window on the chromosome to select ranges from, base-pair units are specified by 'unit'.
#' @param full.overlap logical, the default is to return objects with ANY overlap with the window,
#' whereas setting this as TRUE, will only return those that fully overlap
#' @param unit the unit of base-pairs that 'pos' is using, eg, "b", "kb", "mb", "gb"
#' @param rmv.dup logical, whether to remove duplicate ranges from the return result. The default
#' is not to remove duplicates.
#' @return an object of the same type as 'ranged', but only containing the rows that
#' were within the specified bounds
#' @export
#' @examples
#' require(GenomicRanges)
#' iG <- get.immunog.locs()[2,] # select the 2nd iG region
#' ciG <- chrm(iG)  #  get the chromosome
#' posiG <- c(start(iG),end(iG)) # get the region start and end
#' rr <- rranges(10000) # create a large random GRanges object
#' in.window(rr,chr=ciG,pos=posiG) # set with ANY overlap of iG
#' in.window(rr,chr=ciG,pos=posiG,TRUE) # set with FULL overlap of iG
#' in.window(rr,chr=6,pos=c(25,35),unit="mb") # look between 25 - 35 MB on chr6 [ie, MHC]
in.window <- function(ranged,chr,pos,full.overlap=F, unit=c("b","kb","mb","gb"), rmv.dup=FALSE) {
  if(length(pos)>2 | !is.numeric(pos)) { warning("pos should be a start and end numeric range"); return(NULL) }
  if(length(pos)==1) { pos <- rep(pos,2) }
  all.chrz <- unique(chrNames2(ranged))
  if(length(chr)>1 | (!chr %in% all.chrz)) { warning("chr must be a value in",comma(all.chrz)); return(NULL) }
  typ <- is(ranged)[1]
  if(!any(typ %in% c("RangedData","IRanges","GRanges","RangesList"))) { 
    warning("'ranged' should be a RangedData type or similar"); return(NULL) }
  pos <- pos*make.divisor(unit,"unit")
  #unit <- tolower(unit[1]) ; mult <- switch(unit,b=0,kb=3,mb=6,gb=9); pos <- pos*10^mult
  # get set of genes in a position range for a chromosome
  chr.genez <- chrSel(ranged,paste(chr)) 
  if(full.overlap) {
    ranged <- chr.genez[which(start(chr.genez)>min(pos) & end(chr.genez)<max(pos)),]
  } else {
    # any overlap
    ranged <- chr.genez[which(end(chr.genez)>min(pos) & start(chr.genez)<max(pos)),]
  }
  if(rmv.dup) {
    # remove duplicate genes/exons
    ranged <- ranged[!(duplicated(start(ranged)) & duplicated(end(ranged)) & duplicated(width(ranged))),]
  }
  return(ranged)
}



#' Plot the locations specified in a GRanges or RangedData object
#' 
#' GRanges and RangedData objects are used in bioconductor to store genomic locations and
#' ranges, such as transcripts, genes, CNVs and SNPs. This function allows simple
#' plotting of this data directly from the ranged object. SNPs will be plotted as dots 
#' and ranges as lines. Either can be plotted using vertical bars at the start/end of each
#' range. There are options for labelling and other graphical parameters.
#' This package also creates a generic 'plot' method for GRanges and RangedData that
#' calls this function.
#' @param ranged GRanges or RangedData object with genomic ranges. Should only contain
#' one chromosome, but if not, the first will be used
#' @param labels by default labels for each range are taken from the rownames of 'ranged',
#' but if you want to use another column in the ranged object, specify the column name
#' or number to use to label these ranges on the plot. Or else input a character
#' vector the same length as ranged for custom labels.
#' @param do.labs logical, whether or not to display these labels
#' @param skip.plot.new logical, whether to append to an existing plot (TRUE), or start
#' a new plot (FALSE --> default)
#' @param lty line type to use, see '?lines()' - not used for SNP data when v.lines=FALSE
#' @param alt.y alternative y-axis values (other than the default ordering from the input)
#' This can be a vector of length 1 or length(ranged), or else a column name in ranged to 
#' take the values from
#' @param v.lines TRUE will plot the ranges as pairs of vertical lines, occupying the full
#' vertical extent of the plot, whereas FALSE will plot the ranges as individual horizontal lines
#' @param ylim numeric, length 2, the y-axis limits for the plot, same a 'ylim' for ?plot()
#' @param xlim numeric, length 2, the x-axis limits for the plot, same a 'xlim' for ?plot(),
#' This shouldn't usually be needed as the automatic x-limits should work well,
#'  however is here in case fine tuning is required.
#' @param scl character, the scale that the x axis uses, ie, 'b','kb','mb', or 'gb', meaning
#' base-pairs, kilobases, megabases or gigabase-pairs.
#' @param col character, colour, same as 'col' argument for plot(), etc.
#' @param srt integer, text rotation in degrees (see par) for labels
#' @param pos integer, values of '1', '2', '3' and '4', respectively indicate positions below, 
#' to the left of, above and to the right of the specified coordinates. See 'pos' in graphics:text()
#' @param lwd line width, see '?lines()' - not used for SNP data when v.lines=FALSE
#' @param pch point type, see '?points()' - not used for ranged data
#' @param cex font/symbol size, see '?plot()' - passed to plot, points if using SNP data 
#' @param ... further arguments to 'plot', so long as skip.plot.new==FALSE.
#' @export
#' @return Plots the ranges specified in 'ranged' to the current plot, or to a new plot
#' @examples
#' require(GenomicRanges)
#' rr <- in.window(rranges(5000),chr=6,pos=c(28,32),unit="mb") # make some random MHC ranges
#' rownames(rr) <- paste0("range",1:length(rr))
#' # plotRanges vertically 
#' #print(rr)
#' plotRanges(rr,v.lines=TRUE)
#' # make some labels and plot as horizontal lines #
#' rr2 <- rr[1:5,]; mcols(rr2)[["GENE"]] <- c("CTLA9","HLA-Z","BS-1","FAKr","teST")
#' plotRanges(rr2,label="GENE",scl="Mb",col="black",
#'             xlab="Chr6 position (megabases)",
#'             yaxt="n",ylab="",bty="n")
#' # create some SNPs and plot
#' rr3 <- rr; end(rr3) <- start(rr3) 
#' rownames(rr3) <- paste0("rs",sample(10^6,nrow(rr3)))
#' plotRanges(rr3,col="blue",yaxt="n",ylab="",bty="n")
plotRanges <- function(ranged,labels=NULL,do.labs=T,skip.plot.new=F,lty="solid", alt.y=NULL,
                        v.lines=FALSE,ylim=NULL,xlim=NULL,scl=c("b","Kb","Mb","Gb"),
                        col=NULL,srt=0,pos=4,pch=1,lwd=1,cex=1,...) {
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) { 
    warning("ranged needs to be a RangedData or GRanges object, plot likely to fail") ; return(NULL) }
  chk <- chrNums(ranged)
  typ <- is(ranged)[1]
  if(!is.null(alt.y)) {
    if(is.numeric(alt.y)) {
      if(length(alt.y)==1 | length(alt.y)==length(ranged)) {
        yy <- alt.y
      } else {
        warning("alt.y ignored, must be same length as ranged, or else length 1"); alt.y <- NULL
      }
    } else {
      if(is.character(alt.y)) {
        if(typ=="GRanges") { 
          cn <- colnames(mcols(ranged)); df <- mcols(ranged)
        } else {
          cn <- colnames(ranged); df <- ranged
        }
        if(!alt.y %in% cn) { stop("alternative y.axis column name ",alt.y," not found in 'ranged'") }
        yy <- df[,alt.y]; rm(df)
      } else { 
        warning("invalid value for alt.y, ignoring"); alt.y <- NULL
      }
    }
  }
  if(!is.null(labels)) {
    labels <- paste(labels)
    if(is.character(labels)) {
      if(length(labels)==1 | length(labels)==length(ranged)) {
        if(length(labels)==1) {
          if(typ=="GRanges") { 
            cn <- colnames(mcols(ranged)); df <- mcols(ranged)
          } else {
            cn <- colnames(ranged); df <- ranged
          }
          if(!labels %in% cn) { stop("labels column name ",labels," not found in 'ranged'") }
          lab <- df[,labels]; rm(df)
        } else {
          lab <- labels
        }
      } else {
        warning("labels ignored, must be same length as ranged, or else length 1"); labels <- NULL
      }
    } else {
      warning("invalid value for labels, ignoring"); labels <- NULL
    }
  } else {
    lab <- rownames(ranged) 
  } 
  if(length(chk)>1) { 
    warning(length(chk)," chromosomes in 'ranged', only using the first, chr",chk[1]) 
    ranged <- chrSel(ranged,1) 
  }
  if(all(width(ranged)<=1)) { theyAreSnps <- TRUE } else { theyAreSnps <- FALSE }
  scl <- make.divisor(scl)
  xl <- range(c(start(ranged),end(ranged)),na.rm=T)
  xl <- xl + ((diff(xl)*0.1)*c(-1,1))
  xl <- xl/scl
  nr <- nrow(ranged); if(is.null(nr)) { nr <- length(ranged) }
  if(is.null(alt.y)) {
    yl <- c(0,(nr+2))
  } else {
    yl <- range(yy,na.rm=T)
  }
  if(is.numeric(ylim) & length(ylim)==2) {
    ylim <- range(ylim,na.rm=T)
    ydif <- diff(ylim)
    yl <- ylim
  }
  if(is.numeric(xlim) & length(xlim)==2) {
    xlim <- range(xlim,na.rm=T)
    xdif <- diff(xlim)
    xl <- xlim
  }
  if(is.null(alt.y)) {
    YY <- seq(from=yl[1],to=yl[2],length.out=nr+2)[-1]
  } else {
    if(length(yy)==1) { YY <- rep(yy,length(nr)) } else { YY <- yy }
  }
  #print(YY)
  if(!is.null(col)) {
    if(length(col)==1) {
      col <- rep(col,times=nr) 
    } else {
      if(length(col)!=nr) { warning("col was not the same length as ranged, using first only"); col <- rep(col[1],nr) }
    }
  }
  if(is.null(col)) {
    if(nr>22) { colz <- rep("black",nr) } else { colz <- get.distinct.cols(nr) }
  } else { colz <- col[1:nr] }
  if(is.null(lab) & do.labs) { lab <- paste(1:nr) } # last resort
  if(!skip.plot.new) {
    position <- c(start(ranged[1,]),end(ranged[1,]))/scl
    Y <- YY[c(1,1)]
    #prv(position,Y)
    TY <- if(theyAreSnps) { "p" } else { "l" }
    if(v.lines) {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col="white", lty=lty, ...)
      abline(v=position,col=colz[1])
    } else {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col=colz[1], lty=lty, lwd=lwd, cex=cex, ...)
    }
    st <- 2
  } else {
    st <- 1
  }
  if(nr>1 | st==1) {
    for (cc in st:nr) {
      if(v.lines) {
        abline(v=c(start(ranged[cc,]),end(ranged[cc,]))/scl,col=colz[cc],lty=lty)
      } else {
        if(theyAreSnps) { 
          points(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], pch=pch, cex=cex)
        } else { 
          lines(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], lty=lty, lwd=lwd)
        }
      }
    }
  }
  if(do.labs) {
    for (cc in 1:nr) {
      if(v.lines) { YY <- rep(tail(YY,1),length(YY)) }
      V.scale <- (diff(head(YY,2))*0.5)
      if(length(V.scale)<1 | srt!=90) { V.scale <- 0 }
      text(x=start(ranged[cc,])/scl,y=YY[cc]+V.scale,labels=lab[cc],cex=0.6,pos=pos,offset=0,srt=srt)
    }
  }
}




#' Change the chromosome labels in a RangedData or GRanges object to string codes
#' 
#' @param ranged A GRanges or RangedData object
#' @param do.x.y logical, if TRUE then the usual numbers allocated to chromosomes, X,Y,XY, MT will
#' be allocated as 23,24,25,26 respectively. If false, these will just have 'chr' appended as a
#' prefix
#' @param keep logical, whether to keep additional metadata columns in the new object 
#' @seealso \code{\link{set.chr.to.numeric}}
#' @export
#' @return returns the 'ranged' object, but wherever a chromosome number was previously, a character
#' label, e.g, 'chr1', or 'X', will returned to replace the number, e.g, 1 or 23 respectively. 
#' If table.out is TRUE will return a list where the first element is the resulting object, and the second 
#' element is a table showing which numbers were converted to what label This table
#' can then be used for future conversions via the parameter 'table.in' to ensure consistency of
#' coding.
#' @examples
#' x <- rranges()
#' x
#' x <- set.chr.to.numeric(x) # make entirely numeric
#' x <- rranges(chr.range=20:26)
#' # next two will give warning about X, Y, etc
#' set.chr.to.char(x) # 23 = chrX, etc
#' set.chr.to.char(x,do.x.y=FALSE) # 23=chr23, etc
set.chr.to.char <- function(ranged,do.x.y=T,keep=T) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranged)[1]
  if(!typ %in% c("GRanges","RangedData")) { warning("not a GRanges or RangedData object"); return(NULL) }
  if(length(grep("chr",chrNames2(ranged)))<length(chrNames2(ranged))) {
    ranged <- toGenomeOrder2(ranged,strict=TRUE)
    #prv(ranged)
    mychr2 <- mychr <- paste(chr2(ranged))
    RN <- rownames(ranged)
    #all.nams <- chrNames2(ranged)
    #all.nums <- chrNums(ranged,table.in=table.in)
    if(length(grep("23",paste(mychr2)))>0) { 
      warning("use of arbitrary chromosome numbers for non-autosomes (i.e, >=23)",
              "can lead to annotation issues, try to use labels, X, Y, MT, and XY where possible") }
    sar <- select.autosomes(ranged)
    if(nrow(sar)>0) {
      all.nums.t <- chrNums(sar,table.in=NULL,table.out=T) 
      all.nams <- all.nums.t[,1]
      all.nums <- all.nums.t[,2]
      #mychr2 <- all.nums.t[,2][match(mychr,all.nums.t[,1])]
      for (cc in 1:length(all.nams)) { mychr2[which(mychr==all.nams[cc])] <- paste("chr",all.nums[cc],sep="") }
    } else {
      # no autosomes
    }
    if(do.x.y) {
      mychr2 <- gsub("X","chrX",mychr2)
      mychr2 <- gsub("Y","chrY",mychr2)
      mychr2 <- gsub("23","chrX",mychr2)
      mychr2 <- gsub("24","chrY",mychr2)
      mychr2 <- gsub("25","chrXY",mychr2)
      mychr2 <- gsub("26","chrM",mychr2)
      mychr2 <- gsub("chrXchrY","XY",mychr2)
      mychr2 <- gsub("chrYchrX","YX",mychr2) 
      mychr2 <- gsub("MT","chrM",mychr2)
      mychr2 <- gsub("XY","chrXY",mychr2)
      mychr2 <- gsub("chrchr","chr",mychr2)
    } else {
      mychr2[mychr2 %in% paste(23:100)] <- paste0("chr",mychr2[mychr2 %in% paste(23:100)])
    }
    #print(tail(mychr2)); print((all.nums))
    #prv(mychr2)
    if(any(is.na(mychr2))) { prv(mychr2[which(is.na(mychr2))]) }
    if(is.null(RN) | length(RN)!=nrow(ranged)) { RN <- 1:nrow(ranged) } #make sure RN's are valid
    if(is(ranged)[1]=="GRanges") {
      all.chr <- chr2(ranged)
      out <- GRanges(ranges=IRanges(start=start(ranged),end=end(ranged),names=RN),seqnames=mychr2)
    } else {
      out <- RangedData(ranges=IRanges(start=start(ranged),end=end(ranged),names=RN),space=mychr2)
    }
    out <- toGenomeOrder2(out,strict=TRUE)
    # prv(out)
    ###return(out)
    # need to allow for different indexing of dataframe part for GRanges
    ncr <- switch(typ,RangedData=ncol(ranged),GRanges=ncol(mcols(ranged)))
    if(is.null(ncr)) { ncr <- 0 }
    if(ncr>0 & keep) {
      cn <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
      for(cc in 1:ncr) {
        if(typ=="GRanges") {
          mcols(out)[[paste(cn[cc])]] <- mcols(ranged)[[paste(cn[cc])]]
        } else {
          out[[paste(cn[cc])]] <- ranged[[paste(cn[cc])]]
        }
      }
    }
    return(out)
  } else {
    #cat("no change\n")
    return(ranged)      # change not needed
  }
}




#' Change the chromosome labels in a RangedData or GRanges object to numbers
#' 
#' @param ranged A GRanges or RangedData object
#' @param keep logical, whether to keep additional metadata columns in the new object 
#' @param table.in matrix/data.frame object, usually a result of a prior run of 
#' set.chr.to.numeric(table.out=TRUE), which shows for each label (column 1), what
#' chromosome number should correspond. A way of ensuring consistent coding in different
#' sets.
#' @param table.out logical, if FALSE, the output will just be the object with updated 
#' chromosome labels. If TRUE, then the output will be a list, where the first element
#' is the updated object and the second object is a table describing the coding
#' scheme used to convert from labels to numeric indices.
#' @seealso \code{\link{set.chr.to.char}}
#' @export
#' @return returns the 'ranged' object, but wherever a chromosome label was previously a character
#' label, e.g, 'chr1', or 'X', will return as a number, e.g, 1 or 23 respectively. If table.out
#' is TRUE will return a list where the first element is the resulting object, and the second 
#' element is a table showing which labels were converted to what number. This table
#' can then be used for future conversions via the parameter 'table.in' to ensure consistency of
#' coding.
#' @examples
#' char <- rranges(chr.pref=TRUE)
#' char
#' set.chr.to.numeric(char)
#' # behaviour with X, Y, etc
#' char <- rranges(chr.range=c(20:26))
#' #' char
#' set.chr.to.numeric(char)
#' tab <- set.chr.to.numeric(char,table.out=TRUE)[[2]]
#' tab # codes used in conversion #
#' char <- rranges(chr.range=c(20:26))
#' set.chr.to.numeric(char, table.in=tab) # code using codes from 'tab'
set.chr.to.numeric <- function(ranged,keep=T,table.in=NULL,table.out=FALSE) {
  typ <- is(ranged)[1]
  if(!typ %in% c("GRanges","RangedData")) { warning("not a GRanges or RangedData object"); return(NULL) }
  if(table.out | suppressWarnings(any(is.na(as.numeric(paste(chr2(ranged))))))) {
    silly.name <- "adf89734t5b"
    ranged <- toGenomeOrder2(ranged,strict=T)
    if(typ=="GRanges") {
      mcols(ranged)[[silly.name]] <- paste(1:nrow(ranged))
    } else {
      ranged[[silly.name]] <- paste(1:nrow(ranged))
    }
    #prv(ranged)
    mychr2 <- mychr <- paste(chr2(ranged))
    #all.nams <- chrNames2(ranged)
    #all.nums <- chrNums(ranged,table.in=table.in)
    all.nums.t <- chrNums(ranged,table.in=table.in,table.out=T) 
    all.nams <- all.nums.t[,1]
    all.nums <- all.nums.t[,2]
    #mychr2 <- all.nums.t[,2][match(mychr,all.nums.t[,1])]
    for (cc in 1:length(all.nams)) { mychr2[which(mychr==all.nams[cc])] <- all.nums[cc] }
    #print(tail(mychr2)); print((all.nums))
    if(typ=="GRanges") {
      all.chr <- chr2(ranged)
      out <- GRanges(ranges=IRanges(start=start(ranged),end=end(ranged)),seqnames=mychr2,silly.name=mcols(ranged)[[silly.name]])
    } else {
      out <- RangedData(ranges=IRanges(start=start(ranged),end=end(ranged)),space=mychr2,silly.name=ranged[[silly.name]])
    }
    out <- toGenomeOrder2(out,strict=T)
    if(typ=="GRanges") {
      oo <- mcols(out)[["silly.name"]]
      rr <- mcols(ranged)[[silly.name]]
    } else {
      oo <- out[["silly.name"]]
      rr <- ranged[[silly.name]]
    }
    if(all(!is.na(oo))) {
      if(is.null(rownames(ranged))) { rownames(ranged) <- paste(1:nrow(ranged)) ; rmv.rn <- TRUE } else { rmv.rn <- FALSE }
      #iioo <- rownames(ranged)[match(oo,rr)]; print(iioo); print(is(iioo))
      rn <- narm(rownames(ranged)[match(oo,rr)])
      if(nrow(out)==length(rn) ) { rownames(out) <- rn } else { warning("rownames did not match number of rows") }
      if(rmv.rn) { rownames(out) <- NULL }
    } else {
      warning("index column was corrupted")
    }
    # prv(out)
    ncr <- switch(typ,RangedData=ncol(ranged),GRanges=ncol(mcols(ranged)))
    if(is.null(ncr)) { ncr <- 0 }
    if(ncr>0 & keep) {
      cn <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
      for(cc in 1:ncr) {
        if(typ=="GRanges") {
          mcols(out)[[paste(cn[cc])]] <- mcols(ranged)[[paste(cn[cc])]]
        } else {
          out[[paste(cn[cc])]] <- ranged[[paste(cn[cc])]]
        }
      }
    }
    cno <- switch(typ,RangedData=colnames(out),GRanges=colnames(mcols(out)))
    if(any(cno %in% "silly.name")) { out <- out[,-which(cno %in% "silly.name")] }
    cno <- switch(typ,RangedData=colnames(out),GRanges=colnames(mcols(out)))
    if(any(cno %in% silly.name)) { out <- out[,-which(cno %in% silly.name)] }
    if(table.out) {
      return(list(ranged=out,table.out=all.nums.t))
    } else {
      return(out)
    }
  } else {
    #cat("no change\n")
    cnr <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
    if(any(cnr %in% "silly.name")) { ranged <- ranged[,-which(cnr %in% "silly.name")] }
    return(ranged)      # change not needed
  }
}





#' Invert a ranged object
#' Select the empty space between ranges for the whole genome, for instance you may want
#' to overlap with everything NOT in a set of ranges.
#' @param X a ranged object, GRanges, RangedData or ChipInfo
#' @param inclusive logical, TRUE if the ends of ranges should be in the inverted object
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param pad.missing.autosomes logical, whether to add entire chromosomes to the inverted
#' range object when they are not contained within X
#' @return a ranged object of the same type as X, but with the inverse set of human genomic ranges selected
#' @export
#' @examples
#' \donttest{
#' X <- rranges()
#' invGRanges(X,inclusive=TRUE)
#' invGRanges(X)
#' invGRanges(X,pad.missing.autosomes=FALSE)
#' }
invGRanges <- function(X,inclusive=FALSE,build=NULL,pad.missing.autosomes=TRUE) {
  typ <- is(X)[1]
  if(!typ %in% c("GRanges","RangedData","ChipInfo")) { stop("invalid type for X; ",typ) }
  
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  X <- toGenomeOrder2(X)
  X <- set.chr.to.char(X)
  ch <- chrNames2(X)
  chrLs <- get.chr.lens(mito=T,names=T,build=build)
  chm <- ch
  chm[ch %in% c("chrXY","XY")] <- gsub("Y","",chm[ch %in% c("chrXY","XY")])
  ii <- match(chm,names(chrLs))
  if(any(is.na(ii))) { stop("contained chromosome name not in reference: ",comma(ch[is.na(ii)])) }
  chrL <- as.integer(chrLs[ii])
  #all.dat <- GRanges()
  if(length(ch)>0) { 
  	all.dat <- vector("list",length=length(ch)); names(all.dat) <- ch
  	offs <- if(inclusive) { 0 } else { 1 }
    for (cc in 1:length(ch)) {
      nxt.chr <- chrSel(X,ch[cc])
      st <- as.integer(start(nxt.chr)); en <- as.integer(end(nxt.chr))
      new.st <- as.integer(c(1,en+offs))
      new.en <- as.integer(c(st-offs,chrL[cc]))
      #prv(new.st,new.en)
      if(any(new.en<new.st)) {
        ind <- (rep(which(new.en<new.st),each=5)+rep(c(-2,-1,0,1,2),length(which(new.en<new.st))))
        ind <- ind[ind %in% 1:length(new.st)] #; prv(ind)
        cat("Found illegal start/end in ",ch[cc],"\n")
        print(head(cbind(chr=(rep(ch[cc],length(new.st))),start=new.st,end=new.en)[ind,]))
        if(length(grep("19",ch[cc]))>0) { cat("you may be using incorrect value of 'build' (current is '",build,"')\n",sep="") }
      }
      all.dat[[cc]] <- makeGRanges(chr=rep(ch[cc],length(new.st)),start=new.st,end=new.en)
    }
  } else {
  	all.dat <- NULL
  	warning("X was empty")
  	if(!pad.missing.autosomes) { return(NULL) }
  }
  if(pad.missing.autosomes) {
    autoz <- paste0("chr",1:22)
    misn <- (!autoz %in% ch)
    if(any(misn)) { 
      mis.list <- vector("list",length(which(misn))); names(mis.list) <- autoz[misn]
      for(dd in 1:length(which(misn))) {
        mis.list[[dd]] <- makeGRanges(chr=autoz[which(misn)[dd]],start=1,end=chrLs[which(misn)[dd]])
      }
      all.dat <- c(all.dat,mis.list)
    }
  }
  myDat <- do.call("rbind",args=lapply(all.dat,as,"RangedData"))
  myDat <- toGenomeOrder2(myDat)
  myDat <- as(myDat,typ)
  return(myDat)
}


################## end ranged ##########################

######################
## Simple Functions ##
######################


#' Normalize Lambda inflation factors to specific case-control count
#' 
#' Lambda inflation statistics are influenced by the size of the generating datasets. To facilitate
#' comparison to other studies, this function calculates then converts a given lambda from 
#' n cases and m controls, to be equivalent to 1000 cases and 1000 controls.
#' @param p.values numeric, a vector of analysis p.values, generated from n cases and m controls (although order switching n/m makes no difference to this function)
#' @param n integer, original number of cases that p.values were derived from
#' @param m integer, original number of controls that p.values were derived from
#' @return A normalized Lambda coefficient
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references Freedman M.L., et al. Assessing the impact of population stratification
#'  on genetic association studies. Nat. Genet. 2004;36:388-393.
#' @examples
#' # create some p-values with clear 'inflation' (divergence from uniform[0,1])
#' p.vec <- c(runif(3000)/200,runif(7000)) 
#' # let's imagine these p values come from 3000 cases and 5000 controls
#' L1000_a <- lambda_1000(p.vec,3000,5000)
#' # alternatively, imagine the sample sizes are 10 times larger
#' L1000_b <- lambda_1000(p.vec,30000,50000)
#' plot(sort(p.vec),type="l") 
#' L1000_a; L1000_b
lambda_1000 <- function(p.values,n=1000,m=1000) {
  if(!is.numeric(p.values)) { stop("p.values must be numeric") }
  if(!is.numeric(n)) { stop("n must be numeric") } else { n <- abs(round(n)) }
  if(!is.numeric(m)) { stop("m must be numeric") } else { m <- abs(round(m)) }
  Lnm <- median(p.to.Z(narm(p.values))^2,na.rm=T)/.454
  return(1 + ((Lnm-1)*(((1/n)+(1/m))/((1/1000)+(1/1000)))) )
}



#' Convert a chr:pos1-pos2 vector to a matrix
#' 
#' Takes standard text positions, such as what you might see on the UCSC genome browser, such as 
#' chr1:10,000,234-11,000,567 for a range, or chrX:234,432 for a SNP, and converts to 
#' with cols: chr, start, end.
#' @param text character vector, format like chr:pos1-pos2
#' @export
#' @return a matrix of the same length as 'ranges' with columns chr, start and end, and
#' rownames will be the same as the original text vector.
#' @seealso \code{\link{ranged.to.txt}}
#' @examples
#' txt <- ranged.to.txt(rranges())
#' convert.textpos.to.data(txt)
convert.textpos.to.data <- function(text) {
  do.one <- function(X) {
    chr.pos <- strsplit(X,":",fixed=T)[[1]]  
    chr.txt <- gsub("chr","",chr.pos[[1]],ignore.case=T)
    chr <- chr.txt #allow for X , as.integer(chr.txt); if(is.na(chr)) { print(chr.txt) }
    pos.txt <- strsplit(chr.pos[[2]],"-",fixed=T)[[1]]
    pos12 <- as.integer(pos.txt)
    out <- c(chr,pos12[1],pos12[2])
    if(is.na(out[3])) { out[3] <- out[2] }
    names(out) <- c("chr","start","end")
    return(out)
  }
  return(t(sapply(text,do.one)))
}




#' Make a compact version of gene annotation
#'
#' When adding gene annotation to genomic ranges, sometimes
#' there are many genes associated with a single feature, so
#' that compiling a table becomes awkward, if some rows contain
#' hundreds of genes. This function takes a character vector
#' of gene lists delimited by some separator and provides
#' a compact representation of the gene labels
#' 
#' @param x is a character vector of gene label listings, where multiple hits
#' are delimited by 'sep'
#' @param n number of genes to list before abbreviating
#' @param sep character, separator used to delimit genes in elements of x
#' @param others logical, TRUE to abbreviate with '+ # others' or FALSE to
#' append just the number of genes not listed.
#' @return a character vector with the form:
#' gene-1, gene-2, ..., gene-n, + length(gene-n) - n [others]
#' @export
#' @examples
#' my.genes <- c("ERAP1","HLA-C;CTLA4;IFIH","INS;MYC","AGAP1;APOE;DRDB1;FUT2;HCP5;BDNF;COMT")
#' compact.gene.list(my.genes)
#' compact.gene.list(my.genes,n=2,others=TRUE)
compact.gene.list <- function(x,n=3,sep=";",others=FALSE) {
  XX <- strsplit(x,sep,fixed=T)
  #prv(XX)
  XX <- lapply(XX,function(x) { x[x %in% c(""," ")] <- "unnamed gene"; if(length(x)==0) { x <- "unnamed gene" };return(x) })
 # prv(XX)
  lens <- sapply(XX,length)
  sel <- which(lens>n)
  all <- sapply(XX,function(x) { sel <- FALSE; if(length(x)>0) { sel <- 1:(min(length(x),n)) }; paste(x[sel],collapse=sep) })
  extrz <- lens[sel]-n
  if(others) { oth <- rep("others",length(extrz)); oth[extrz==1] <- "other" }
  all[sel] <- paste(all[sel],"+",extrz,if(others) { oth } else { "" } )
  return(all)  
}


#' Meta-analysis using odds ratio and standard error from 2 datasets
#' 
#' This function calculates meta analysis odds ratios, standard errors and p-values
#' using results from a table containing odds ratio and standard error data for analyses
#' of 2 different datasets (typically logistic regression, but other analyses can be
#' incorporated if an odds-ratio and SE can be derived, for instance one analysis might
#' be a case control logistic regression GWAS and the other a family TDT analysis).
#' @param X A data.frame with column names which should be entered in the parameters:
#' OR1, OR2, SE1, SE2, and optionally N1, N2. 
#' @param OR1 The column name of X containing odds ratios from the first analysis
#' @param OR2 Same as OR1 above but pertaining to the second analysis
#' @param SE1 The column name of X containing standard errors from the first analysis
#' @param SE2 Same as SE1 above but pertaining to the second analysis
#' @param N1 Only required if method="sample.size". Either the column name in X with the 
#' number of samples in the first analysis, of a vector of the same, or if N's is the
#'  same for all rows, a scalar' value can be entered
#' for each.
#' @param N2 Only required if method="sample.size". Same as N1 above but pertaining to analysis 2
#' @param Z1 Only use if method="sample.size" or "z.score". The column name in X with the 
#' z.scores in the first analysis.
#' @param Z2 Same as Z1 above but pertaining to analysis 2
#' @param method character, can be either 'beta', 'z.score' or 'sample.size', and upper/lower
#' case does not matter. 'Beta' is the default and will calculate meta-analysis weights using
#' the inverse variance method (based on standard errors), and will calculate the p-values
#' based on the weighted beta coefficients of the two analyses. 'Z.score' also uses inverse variance
#' but calculates p-values based on the weighted Z scores of the two analyses. 'Sample.size' uses
#' the sqrt of the sample sizes to weight the meta analysis and uses Z scores to calculate p values
#' like 'Z.score' does.#' 
#' @return The object returned should have the same number of rows and rownames as the data.frame
#'  X but columns are the meta analysis stastistics, namely:
#'   OR.meta, beta.meta, se.meta, z.meta, p.meta, which will contain the meta
#' analysis odds-ratio, beta-coefficient, standard error, z-score, and p-values respectively
#' for each row of X.
#' @export
#' @examples
#' X <- data.frame(OR_CC=c(1.8,1.15),OR_Fam=c(1.33,0.95),SE_CC=c(0.02,0.12),SE_Fam=c(0.07,0.5))
#' rownames(X) <- c("rs689","rs23444")
#' X
#' meta.me(X)
#' X <- data.frame(OR_CC=c(1.8,1.15),OR_CC2=c(1.33,0.95),
#'  SE_CC=c(0.02,0.12),SE_CC2=c(0.02,0.05),
#'  n1=c(5988,5844),n2=c(1907,1774))
#' # even with roughly the same number of samples the standard error will determine the influence of
#' # each analysis on the overall odds ratio, note here that the second SE for dataset goes
#' # from 0.5 to 0.05 and as a result the estimate of the odds ratio goes from 1.137 to 0.977,
#' # i.e, from very close to OR1, changing to very close to OR2.
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2") 
#' # sample size and z-score methods give similar (but distinct) results
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2",N1="n1",N2="n2",method="sample.size") 
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2",N1="n1",N2="n2",method="z.score")  # N's will be ignored
meta.me <- function(X,OR1="OR_CC",OR2="OR_Fam",SE1="SE_CC",SE2="SE_Fam",Z1=NA,Z2=NA,
                    N1=NA,N2=NA,method=c("beta","z.score","sample.size")) {
  #N1=18856,N2=7638
  validz <- c("beta","z.score","sample.size")
  method <- tolower(method[1])
  if(!method %in% validz) { method <- validz[1]; warning("invalid method entered, using 'beta' method") }
  if(!is(X)[1]=="data.frame") { stop("X must be a data.frame") }
  cnx <- colnames(X)
  if(is.null(rownames(X))) { rownames(X) <- paste(1:nrow(X)) }
  if(!all(c(OR1,OR2,SE1,SE2) %in% cnx)) { stop("X must contain column names specified by OR1,OR2,SE1,SE2") }
  ok <- FALSE
  if(method=="sample.size") {
    if(is.numeric(N1) & is.numeric(N2)) { if(length(N1)!=1 | length(N2)!=1) {
      stop("N1, N2 should either be scalar integers, or column names containing N's") } else {
        N.coln <- FALSE
      } }
    if(is.character(N1) & is.character(N2)) {
      if(!all(c(OR1,OR2,SE1,SE2) %in% cnx)) {
        stop("N1,N2 must contain either be scalar integers or column names with N's") } else {
          N.coln <- TRUE
        }
    }
  }
  OR.CC <- X[,OR1]
  beta.CC  <- log(X[,OR1])
  se.CC <- X[,SE1]
  OR.family <- X[,OR2]
  beta.family  <- log(X[,OR2])
  se.family <- X[,SE2]
  if(!is.na(Z1) & method!="beta") {
    if(Z1 %in% colnames(X)) {
      z.CC <- X[,Z1]
    } else { warnings("Z1 column not found, ignoring") }
  } else {
    z.CC <- beta.CC/se.CC
  }
  if(!is.na(Z2) & method!="beta") {
    if(Z2 %in% colnames(X)) {
      z.family <- X[,Z2]
    } else { warnings("Z2 column not found, ignoring") }
  } else {
    z.family <- beta.family/se.family
  }
  inv.CC <- 1 / (se.CC^2)
  inv.family <- 1 / (se.family^2)
  var.meta <- 1 / (inv.CC+inv.family)
  weight.CC <- inv.CC * var.meta
  weight.family <- inv.family * var.meta
  se.meta <- round(sqrt(var.meta), digits=3)
  if(method=="sample.size") {
    if(N.coln) {
      famN <- X[,N2]
      ccN <- X[,N1]
    } else {
      famN <- N2 # 3819*2  #3509*2   #  3819*2   #  10796
      ccN <- N1 # 6683+12173 # including CBR, or for UVA analyses use instead: 9416+6670
    }
    WeightFam = sqrt(famN)/(sqrt(famN)+sqrt(ccN))
    WeightCC <- 1-WeightFam
    # beta calculated the same way for sample.size method using sample size weights
    beta.meta <- round((WeightCC * beta.CC) + (WeightFam * beta.family),digits=3) # beta based
    z.meta <- round((WeightCC * z.CC) + (WeightFam * z.family),digits=6) # n-based, z-based
  } else {
 # beta calculated the same way for z.score method and beta method using inverse variance
    beta.meta <- round((weight.CC * beta.CC) + (weight.family * beta.family),digits=3) # beta based
    if(method=="z.score") {
      z.meta <- round((weight.CC * z.CC) + (weight.family * z.family),digits=6) # z-based
    } else {
      # default (beta) method
      z.meta <- beta.meta/se.meta
    }
  }
  OR.meta <- exp(beta.meta)
  p.meta <- 2*pnorm(-abs(z.meta))
  out <- (cbind(OR.meta,beta.meta,se.meta,z.meta,p.meta))
  colnames(out) <- c("OR.meta","beta.meta","se.meta","z.meta","p.meta")
  rownames(out) <- rownames(X)
  return(out)
}



################## end simple ##########################



################################
## ChipInfo Support Functions ##
################################


#' Retrieve current ChipInfo annotation object
#' 
#' This function returns the current 'ChipInfo' annotation object, containing chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18 or hg19 (build 36/37) coordinates. Can also be used to update the 
#' current object to a new object.
#' This package makes extensive use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' of the object used is 'ChipInfo' which is a GRanges object, modified to always
#' have columns for A1, A2 (alleles), rs.id, and a quality control flag. The
#' default display is tidier than GRanges, it has nice coersion to and frame data.frame
#' and indexing by chromosome using [[n]] has been added, in addition to normal [i,j]
#' indexing native to GRanges. A1 and A2 values are usually specific to each dataset
#' so for immunochip you may need to manually update these values to reflect
#' the allele coding in your own dataset.
#' @param build character, either "hg18", "hg19" or "hg38". Will also accept build numbers,
#' 36, 37 or 38.
#' @param refresh logical, FALSE to just load whatever object is already in memory (except
#' when first using a function in this package, there should be a ChipInfo object loaded),
#' or TRUE to reload from the original source. For instance you may wish to do this when 
#' you want to use a different chip, different build, or if the annotation has been modifed
#' via a manual correction).
#' @param alternate.file character, name of an alternative RData file containing a ChipInfo
#' object to use instead of the object found in getOption("chip.info"). This will replace
#' the current ChipInfo object.
#' @param warn.build logical, whether to warn if the 'build' argument does not match
#' the current value of getOption("ucsc"). The default is to display this warning,
#' but if you set this argument to FALSE this can be suppressed.
#' @return returns the current ChipInfo object [S4]. This may be slow first time, but
#' subsequent lookups should be much faster. Builds 36/38 are not stored explicitly so
#' will take a little while to convert the first time, but subsequent lookups should be
#' fast. To increase the speed save the object locally and use option(chip.info=<PATH>)
#' to set a custom path for future chip.support() calls [which are also made internally
#' by many of the function in this package]. This is also the option to set if you
#' want to add a ChipInfo object for a different chip, e.g, metabochip, exomechip, etc.
#' @seealso \code{\link{ChipInfo}}, \code{\link{build}}, \code{\link{rs.id}}, 
#'  \code{\link{QCfail}}, \code{\link{convTo36}}, 
#'  \code{\link{convTo37}}, \code{\link{A1}}, \code{\link{A2}}
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @export
#' @concept ImmunoChip MetaboChip microarray iSelect Illumina
#' @examples
#' chip.support() # shows the current ChipInfo object (default is 'ImmunoChip' build 37)
#' #/donttest{
#' chip.support(build=36) # gives warning as hg19 version is currently loaded
#' chip.support(build=36,refresh=TRUE)
#' getOption("chip.info") # shows the object is now saved in the tmp directory for subsequent calls
#' chip.support(build=38,refresh=TRUE)
#' #}
chip.support <- function(build=NULL,refresh=FALSE,alternate.file=NULL,warn.build=TRUE) {
  update.build <- T #refresh
  if(is.null(build)) { build <- getOption("ucsc") } else {  update.build <- TRUE  }
  build <- ucsc.sanitizer(build)
  ## NEED TO ADD SUPPORT HERE TO USE CORRECT OUT OF 36/37
  if(getOption("chip.info")=="") { refresh <- TRUE }
#  refresh <- T ###??????? why a param?
  #old# if(!exists("all.support",envir=globalenv())) { refresh <- T }  # change global environment to namespace of iChip package
  if(refresh) {
#    options(chip.info="") # refresh this?
    if(getOption("chip.info")=="ImmunoChip_37_builtIn" | update.build) {  options(chip.info="") } # refresh so erase current
    use.options <- TRUE
    if(is.character(alternate.file)) {
      if(file.exists(alternate.file)) {
        use.options <- FALSE
      }
    }
    if(use.options) {
      file <- getOption("chip.info")
      if(length(file)==0) { file <- "" }
    } else {
      file <- alternate.file
    }
    conv <- FALSE # whether to convert from specified build, e.g, compulsory for ichip if not hg19
    if(!file.exists(file)) {
      if(is(file)[1]=="character") { if(nchar(file)>1) { warning("file",file,"did not exist, using default (immunoChip)") }}
      all.support <- humarray::ImmunoChipB37
      if(!build %in% c("hg19","hg38","hg18")) {
            stop("unsupported build",build,
                 "for chip.support() function, please use hg18/19/38",
                  "or use option(chip.info=<PATH>) to specify a custom ChipInfo file")
      } else {
        if(!build=="hg19") {
          conv <- TRUE
        } else {
          options(chip.info="ImmunoChip_37_builtIn")
        }
      }
      #/home/ncooper/github/iChip/data/ImmunoChipB37.rda
      #file <- system.file("extdata", fname, package="humarray")
    } else {
      all.support <- reader(file); conv <- FALSE # custom file supplied, don't convert
      if(!is(all.support)[1]=="ChipInfo") {
        stop("alternate.file file contained object not of class ChipInfo")
      }  else {
        options(chip.info=file)
        message("Updated current ChipInfo object to: ",getOption("chip.info"))
      }
    }
#   if(ucsc(all.support)!=build) { conv <- TRUE } # needs converting
    if(is(all.support)[1]=="list") {
      if(length(all.support)>1) { typz <- sapply(lapply(all.support,is),"[",1) } # in case multiple objects in file
    } else {
      typz <- is(all.support)[1]
    }
    if(all(typz=="ChipInfo")) {
      if((build %in% c("hg18","hg38")) & conv) {
        if(build=="hg18") { all.support <- convTo36(all.support) } else { all.support <- convTo38(all.support) }
        tfn <- cat.path(tempdir(),"ImmunoChip_ChipInfo",suf=substr(build,3,4),ext="rda")
        save(all.support,file=tfn) # save this converted object in the tmp directory for this session
        options(chip.info=tfn)
      }
      #old# assign("all.support",value=all.support,envir=globalenv())  # change global environment to namespace of iChip package
    } else {
      stop("object (all.support) in the file",file,
           "should have type ChipInfo, or else object all.support in the global environment has been modified")
    }
  } else {
    file <- getOption("chip.info");
    if(file=="ImmunoChip_37_builtIn") {
      all.support <- humarray::ImmunoChipB37
    } else {
      if(file.exists(file)) { 
        all.support <- reader(file) 
      } else { 
        stop("getOption('chip.info') file did not exist")
      }
      if(!is(all.support)[1]=="ChipInfo") { stop("chip.info file contained object not of ChipInfo class") }
    }
    if(ucsc(all.support)!=build) {
      all.support <- chip.support(build=build,refresh=TRUE,alternate.file=alternate.file,warn.build=warn.build)
      if(ucsc(all.support)!=getOption("ucsc") & warn.build) {
        warning("'build' did not match the current default reference genome. Use options(ucsc=",build,") if you wish to change the default reference genome. otherwise use 'warn.build=FALSE' to hide this warning")
      }
    }
  }
  #old#if(!exists("all.support",envir=globalenv())) { stop("ChipInfo data object 'all.support' not found") }  
  #old#all.support <- get("all.support",envir=globalenv())  
  # change global environment to namespace of iChip package
  if(is(all.support)[1]=="list") {
    nnn <- (names(all.support))
    if(length(nnn)>1) {
      # multiple objects in file, probably different builds
      any37 <- c(grep("37",nnn),grep("hg19",nnn))
      any36 <- c(grep("36",nnn),grep("hg18",nnn))
      if(build=="hg19" & length(any37)>0) { all.support <- all.support[[nnn[any37[1]]]]}
      if(build=="hg38" & length(any37)>0) { all.support <- convTo38(all.support[[nnn[any37[1]]]])}
      if(build=="hg18" & length(any36)>0) { all.support <- all.support[[nnn[any36[1]]]]}
    }
  }
  return(all.support)
}






#' Convert from chip ID labels to dbSNP rs-ids
#' 
#' Most SNPs will have an 'rs-id' from dbSNP/HapMap, and these are often the standard for reporting or
#' annotation lookup. These can differ from the IDs used on the chip. This functions looks at the current
#' snp support (ChipInfo object) and returns rs-ids in place of chip IDs. Currently rs-ids are always 
#' from build37.
#' @param ids character, meant to be a list of chip ids, but if rs-ids are present they will not be altered.
#' @return A character vector of SNP rs-ids, where the input was chip ids, rs-ids or a mixture, any text
#' other than this will result in NA values being returned in the character vector output.
#' @export
#' @seealso \code{\link{rs.to.id}}, \code{\link{GENE.to.ENS}}, \code{\link{ENS.to.GENE}}
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' id.to.rs(c("imm_11_2138800","rs9467354","vh_1_1108138")) # middle one is already a rs.id
id.to.rs <- function(ids) {
  ids <- clean.snp.ids(ids)
  all.support <- chip.support()
  if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
  rsvec <- mcols(all.support)$rs.id[match(ids,rownames(all.support))]
  rsvec2 <- mcols(all.support)$rs.id[match(ids,mcols(all.support)$rs.id)]
  rsvec[is.na(rsvec)] <- rsvec2[is.na(rsvec)]
  return(rsvec)
}


#' Convert from dbSNP rs-ids to chip ID labels
#' 
#' Most SNPs will have an 'rs-id' from dbSNP/HapMap, and these are often the standard for reporting
#' or annotation lookup. These can differ from the IDs used on the chip. This functions looks at 
#' the current snp support (ChipInfo object) and looks up chip IDs based on rs-ids.
#' @param rs.ids character, meant to be a list of rs-ids, but if chip-ids are present they will not be
#' altered.
#' @param manifest logical, if TRUE return the ids as specified by the manifest/official set, or as
#' stored in the column 'chip.id'. If FALSE return the IDs as stored in the rownames of the ChipInfo
#' object, which can differ when the chip.id is an illegal format for an R row/column name.
#' @param multi.list logical, some rs-ids could map to multiple chip ids. It is recommended that if
#' that is the case then a letter should be appended to duplicate rs-ids to make them unique in the
#' ChipInfo object, e.g, rs1234, rs1234b, rs1234c, etc. If multi.list is TRUE, then the id list 
#' will be returned as a list, and any time an rs-id is entered without a letter suffix, all 
#' possible corresponding chip ids will be listed.
#' @param force.flat logical, if 'multi.list' is true, then some rs-ids might map to more than one
#' SNP. If force.flat is TRUE, then multiple SNP listings will be concatenated with a comma. If 
#' FALSE, then a list will be returned, with multiple entries where applicable.
#' @return A character vector of SNP chip-ids, where the input was rs-ids, chip-ids or a mixture, 
#' any text other than this will result in NA values being returned in the character vector output.
#' Or, if multi-list is true, then returns a list instead, which takes more than 1 value where 
#' there are multiple chip-ids with the same rs-id; if there are no such rs-id duplicates the
#' result will still be a list. Currently rs-ids are always from build37.
#' @export
#' @seealso \code{\link{id.to.rs}}, \code{\link{GENE.to.ENS}}, \code{\link{ENS.to.GENE}}
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' rs.to.id(c("rs689","rs9467354","rs61733845"))  # middle one has no chip id
#' \donttest{
#' test.ids <- c("rs61733845","rs2227313","rs11577783","rs3748816",
#'                                     "rs12131065","rs3790567","rs2270614")
#' rs.to.id(test.ids, multi.list=TRUE) # list with duplicates
#' }
rs.to.id <- function(rs.ids,manifest=FALSE,multi.list=TRUE,force.flat=TRUE) {
  rs.ids <- clean.snp.ids(rs.ids)
  all.support <- chip.support()
  if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
  delay.force <- FALSE
  if(manifest) { if(force.flat) { 
    force.flat <- FALSE; #warning("force.flat set to FALSE, cannot be TRUE when manifest=TRUE") 
    multi.list=TRUE; delay.force <- TRUE
  }}
  if(multi.list) {
    X0 <- rmv.trail(rs.ids)
    X1 <- paste0(X0,"b"); X2 <- paste0(X0,"c"); X3 <- paste0(X0,"d"); X4 <- paste0(X0,"a")
    X <- cbind(rs.to.id(X0,multi.list=F),rs.to.id(X1,multi.list=F),
               rs.to.id(X2,multi.list=F),rs.to.id(X3,multi.list=F),rs.to.id(X4,multi.list=F))
    idvec <- apply(X,1,function(x) { unique(narm(x)) })
    if(!is.list(idvec)) { 
      if(!force.flat) {
        idvec <- as.list(idvec) 
      }
    } else {
      if(force.flat) {
        idvec <- sapply(idvec,paste,collapse=",") 
      }
    }
    #warning("'multi.list' option was used, but no duplicate rs-ids found, so returning a vector, not a list") }
  } else {
    #prv(all.support,rs.ids)
    #print(showMethods(rownames))
    #print(rownames(all.support))
    idvec <- rownames(all.support)[match(rs.ids,mcols(all.support)$rs.id)]
    idvec2 <- rownames(all.support)[match(rs.ids,rownames(all.support))]
    #prv(idvec,idvec2)
    idvec[is.na(idvec)] <- idvec2[is.na(idvec)]
  }
  if(manifest) {
    if(!is.null(mcols(all.support)[,"chip.id"])) {
      mano <- function(X) { mcols(all.support)[,"chip.id"][match(X,rownames(all.support))] }
      if(is(idvec)[1]=="list") {
        idvec <- lapply(idvec,mano)
        if(delay.force) {
          idvec <- sapply(idvec,paste,collapse=",") # force.flat if this option was selected
        }
      } else {
        idvec <- mano(idvec)
      }
      return(idvec)
    } else {
      warning("chip.id was NULL, no 'manifest' Ids found, returning rownames"); return(idvec)
    }
  } else {
    return(idvec)
  }
}



#' Convert from chip/rs-ids to manifest chip ID labels
#' 
#' Some SNP-ids aren't legal R row/column names. In order to match datafiles and store annotation
#' as objects, this package converts SNP-names to sanitized versions if necessary, that are
#' legal row/column names. This function converts from such 'legal' versions of the IDs back
#' to the proper names, as per the chip manifest document (or whatever is stored in the chip.id field
#' of the chip.support() object, accessible using chipId()).
#' @param ids character, either a list of rs-ids or chip-ids. chip ids are preferable as they
#' are unique, and rs.ids are not. Using this function is not recommended for rs.id lists that
#' might have entries that map to multiple chip ids, because entries other than the first will be 
#' ignored. For such cases, use 'rs.to.id(manifest=TRUE,multi.list=TRUE,...)'.
#' @return A character vector of SNP chip-ids matching the manifest format.
#' @export
#' @seealso \code{\link{id.to.rs}}, \code{\link{GENE.to.ENS}}, \code{\link{ENS.to.GENE}}
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#'   test.ids <- c("imm_1_898835","rs61733845","rs115005664","rs114582555",
#'       "chr1_20131940","chr1_20133829","rs150992667","rs138231315","rs111577708","rs187104718")
#'   manifest(c("chr1_20131940","ccc_1_67429655_A_G"))
#'   manifest(test.ids)  # even when some are rs-id, still works
#'   data.frame(rs.id=test.ids,legal.id=rs.to.id(test.ids),manifest.id=manifest(test.ids))
manifest <- function(ids) {
  mano <- function(X) { rs.to.id(X,manifest=T,force.flat=FALSE,multi.list=FALSE) }
  if(is(ids)[1]=="list") {
    ids <- lapply(ids,mano)
  } else {
    ids <- mano(ids)
  }
  return(ids)
}



#' Find chromosome for SNP ids, gene name or band
#' 
#' Allows retrieval of the chromosome associated with a SNP-id, HGNC gene label, karyotype band,
#' or vector of such ids. For SNPs the ids can be either chip ids, or rs-ids, but must be contained
#'  in the current annotation. Default behaviour is to assume 'id' are SNP ids, but if none are
#'  found in the SNP annotation, the id's will be passed to functions Pos.gene() and Pos.band() to
#'  see whether a result is found. This latter step will only happen if no SNP ids are retreived in
#'  the first instance, and if snps.only=TRUE, then genes and bands will not be searched and NA's 
#'  returned. If you are repeatedly searching for chromosomes for genes/bands, using the dedicated 
#'  Pos.gene and Pos.band functions would be slightly faster than relying on the fallback behaviour
#'  of the Chr() function.  See documentation for these functions for more information. The build
#'  used will be that in the current ChipInfo object.
#' @param ids character, a vector of rs-ids or chip-ids representing SNPs in the current ChipInfo
#'  annotation, or gene ids, or karyotype bands. Can also be a SnpMatrix object.
#' @param dir character, only relevant when gene or band ids are entered, in this case 'dir' is the location
#' to download gene and cytoband information; if left as NULL, depending on the value of 
#' getOption("save.annot.in.current"), the annotation will either be saved in the working directory to 
#' speed-up subsequent lookups, or deleted after use.
#' @param snps.only logical, if TRUE, only search SNP ids, ignore the possibility of genes/cytobands.
#' @return A character vector of Chromosomes for each ids, with NA values where no result was found.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Pos}}
#' @examples
#' \donttest{
#' setwd(tempdir())
#' Chr(c("rs689","rs9467354","rs61733845"))
#' Chr("CTLA4")
#' Chr("13q21.31")
#' Chr(c("CTLA4","PTPN22"),snps.only=TRUE) # fails as these are genes
#' Chr(c("rs689","PTPN22","13q21.31")) # mixed input, will default to SNPs, as at least 1 was found
#' }
Chr <- function(ids,dir=NULL,snps.only=FALSE) {
  ic.chr <- function(ids) {
    ic.ids <- clean.snp.ids(ids)
    all.support <- chip.support()
    if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found")  }  ## load object: all.support [snp support for whole chip]
    outlist <- chr(all.support)[match(ic.ids,rownames(all.support))]
    return(outlist)
  }
  typ <- is(ids)[1]
  if(typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { ids <- colnames(ids) }
  query <- rs.to.id(ids)
  if(!snps.only & all(is.na(query))) { 
    ## unless the 'snps.only' function is set, then if it looks like we have not been handed
    ## snp ids, then check for band ids or gene ids instead
    chr.ify <- function(X) { tolower(names(X)); return(gsub("chr","",X["chr"],ignore.case=TRUE)) }
    numpqs <- (length(grep("q",ids))+length(grep("p",ids)))
    if(numpqs==length(ids)) { try.band <- T } else { try.band <- F }
    if(try.band) {
      suppressWarnings(test <- Pos.band(ids,dir=dir))
      if(!is.null(test)) { return(chr.ify(test)) }
    }
    suppressWarnings(test <- Pos.gene(ids,dir=dir))
    if(!is.null(test)) { return(chr.ify(test)) } 
  } 
  ic <- ic.chr(query)
  return(ic)
}


#' Find the chromosome position for SNP ids, gene name or band
#' 
#' Allows retrieval of the the chromosome position associated with a SNP-id, HGNC gene label, 
#'  karyotype band, or vector of such ids. For SNPs the ids can be either chip ids, or rs-ids,
#'  but must be contained in the current annotation. Default behaviour is to assume 'id' are 
#'  SNP ids, but if none are found in the SNP annotation, the id's will be passed to functions
#'  Pos.gene() and Pos.band() to see whether a result is found. This latter step will only happen
#'  if no SNP ids are retreived in the first instance, and if snps.only=TRUE, then genes and bands
#'  will not be searched and NA's returned. If you are repeatedly searching for positions for 
#'  genes/bands, using the dedicated Pos.gene() and Pos.band() functions would be slightly faster
#'  than relying on the fallback behaviour of the Pos() function. Note that the position for
#'  genes and bands are not a single point, so the result will be a range with start and end, 
#'  see 'values' below. See documentation for these functions for more information.
#' @param ids character, a vector of rs-ids or chip-ids representing SNPs in the current ChipInfo 
#' annotation,or gene ids, or karyotype bands. Can also be a SnpMatrix object.
#' @param dir character, only relevant when gene or band ids are entered, in this case 'dir' is 
#' the location to download gene and cytoband information; if left as NULL, depending on the 
#' value of getOption("save.annot.in.current"), the annotation will either be saved in the 
#' working directory to speed-up subsequent lookups, or deleted after use.
#' @param snps.only logical, if TRUE, only search SNP ids, ignore the possibility of 
#' genes/cytobands.
#' @return When ids are SNP ids, returns a numeric vector of positions for each id, with NA 
#' values where no result was found. When ids are genes or karyotype bands, will return a 
#' data.frame with columns 'chr' [chromosome], 'start' [starting position of feature], 'end'
#' [end position of feature], and the band without the chromosome prefix, if ids are bands. 
#' Note that this function cannot retrieve multiple ranges for a single gene (e.g, OR2A1 in
#' build 38), which means you'd need to use Pos.gene(). The coordinates used will be of 
#' version getOption(ucsc="hg18"), or ucsc(chip.support()), which should be equivalent.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}
#' @examples
#' \donttest{
#' setwd(tempdir())
#' Pos(c("rs689","rs9467354","rs61733845"))
#' Pos("CTLA4") # returns a range
#' Pos("13q21.31") # returns a range
#' Pos(c("CTLA4","PTPN22"),snps.only=TRUE) # fails as these are genes
#' Pos(c("rs689","PTPN22","13q21.31")) # mixed input, will default to SNPs, as at least 1 was found
#' }
Pos <- function(ids,dir=NULL,snps.only=FALSE) {
  all.support <- chip.support()
  ic.pos <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
    outlist <- start(all.support)[match(ic.ids,rownames(all.support))]
    return(outlist)
  }
  typ <- is(ids)[1]
  if(typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { ids <- colnames(ids) }
  query <- rs.to.id(ids)
  if(!snps.only & all(is.na(query))) { 
    ## unless the 'snps.only' function is set, then if it looks like we have not been handed
    ## snp ids, then check for band ids or gene ids instead
    numpqs <- (length(grep("q",ids))+length(grep("p",ids)))
    if(numpqs==length(ids)) { try.band <- T } else { try.band <- F }
    if(try.band) {
      suppressWarnings(test <- Pos.band(ids,dir=dir))
      if(!is.null(test)) { return(test) }
    }
    suppressWarnings(test <- Pos.gene(ids,dir=dir))
    if(!is.null(test)) { return(test) } 
  } 
  ic <- ic.pos(query)
  return(ic)
}



#' Order rs-ids or ichip ids by chrosome and position
#' 
#' Simple function to sort a character list of SNP ids into genome order.
#' @param ids character, vector of SNP rs-ids or chip-ids, see rs.to.id()
#' @return the same vector 'ids', sorted by genome position
#' @seealso \code{\link{rs.to.id}}, \code{\link{id.to.rs}}, \code{\link{Chr}}, \code{\link{Pos}}
#' @export
#' @examples
#' \donttest{
#' snp.ids <- c("rs3842724","imm_11_2147527","rs689","rs9467354","rs61733845")
#' Chr(snp.ids) # shows each is on a different chromosome
#' Pos(snp.ids)
#' ids.by.pos(snp.ids)
#' Chr(ids.by.pos(snp.ids))
#' Pos(ids.by.pos(snp.ids))
#' }
ids.by.pos <- function(ids) {
  if(!is.character(ids)) { stop("ids must be a character vector") }
  pp <- Pos(ids)
  if(any(is.na(pp))) { stop("invalid id list, 'ids' must all be valid, with position information in the ChipInfo object") }
  ids <- ids[order(pp)]
  cc <- Chr(ids)
  ids <- ids[order.chr(cc)]
  return(ids)
}




#' Find the chromosome, start and end position for gene names
#' 
#' Allows retrieval of the the chromosome position associated with a HGNC gene label, 
#'  or vector of such labels. Note that the position returned for genes is not a 
#'  single point as for SNPs, so the result will be a chromosome, then a position range with
#'  start and end.
#' @param genes character, a vector of gene ids
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, else 
#' a data.frame
#' @param band logical, whether to include band/stripe in returned object
#' @param one.to.one logical, some genes have split ranges, TRUE merges these to give only 1 range 
#' per gene, NB: this is the default behaviour when using the more general Pos() function
#' @param remap.extra logical, if TRUE genes with chromosome annotation 'c6_cox' and 'c6_QBL' will
#'  be mapped to chromosome 6, and 'NT_xxxx' chromosome labels will all be mapped to 'Z_NT', etc
#' @param discard.extra logical, if TRUE then any gene hit with chromosome not in 1:22, X, Y, XY, MT, 
#' will be discarded.
#' @param warnings logical, whether to show warnings when some/all ids are not matched to the 
#' reference
#' @return Returns a data.frame with columns 'chr' [chromosome], 'start' [starting position of the
#'  gene],'end' [end position of the gene], or if bioC=TRUE, then returns a GRanges object with
#'  equivalent information, and if band=TRUE, then an extra column is added with band information
#'  If returning a data.frame, then it will be in the same order as 'genes'. If bioC=TRUE, then
#'  the result will be in genome order, regardless of the order of 'genes'.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}, \code{\link{Pos}}, \code{\link{Pos.band}}, \code{\link{Band}}, 
#' \code{\link{Band.gene}}, \code{\link{Band.pos}}, \code{\link{Gene.pos}}
#' @examples
#' \donttest{
#' setwd(tempdir())
#' Pos.gene(c("CTLA4","PTPN22"))
#' Pos.gene("MYC",build=36)
#' Pos.gene("MYC",build=37)
#' Pos.gene(c("CTLA4","PTPN22"),bioC=TRUE,band=TRUE)
#' Pos.gene(c("CTLA4","OR2A1"),one.to.one=TRUE,build=38) # OR2A1 is split over two ranges
#' Pos.gene(c("CTLA4","OR2A1"),one.to.one=FALSE,build=38)
#' Pos.gene("RNU2-1",one.to.one=FALSE,bioC=TRUE,build=38) # RNU2-1 is split over multiple ranges
#' }
Pos.gene <- function(genes,build=NULL,dir=NULL,bioC=FALSE,band=FALSE,one.to.one=TRUE,
                     remap.extra=FALSE,discard.extra=TRUE,warnings=TRUE) {
  if(!is.character(genes)) { stop("'genes' must be a character vector of gene names") }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  char.lim <- 100
  ga <- suppressWarnings(get.gene.annot(dir=dir,build=build,bioC=bioC,one.to.one=one.to.one,
                       remap.extra=remap.extra,discard.extra=discard.extra,GRanges=TRUE))
  typ <- is(ga)[1]
  if(typ=="GRanges") {  mt <- match(genes,mcols(ga)$gene) } else { mt <- match(genes,ga$gene) }
  failz <- paste(genes[is.na(mt)],collapse=", "); if(nchar(failz)>char.lim) { failz <- paste(substr(failz,1,char.lim),",...",sep="") }
  if(length(mt)<1 | all(is.na(mt))) { 
    if(warnings) { warning("did not find any 'genes' features: ",failz) }; return(NULL) }
  if(any(is.na(mt))) { 
    cnt <- length(which(is.na(mt)))
    if(warnings) { warning("did not find the following ",cnt," 'genes' features: ",failz) }
  }
  if(!one.to.one) { 
    order.important <- FALSE  # not currently implemented
    if(order.important) {
     if(length(genes)>0) {
       mt <- NULL
       for(gg in 1:length(genes)) {
         if(typ=="GRanges") { 
           mt <- c(mt,which(mcols(ga)$gene %in% genes[gg]))
         } else {
           mt <- c(mt,which(ga$gene %in% genes[gg]))
         }
       }
     } else { warning("length of genes entered was zero?") }     
    } else {
      if(typ=="GRanges") { 
        mt <- which(mcols(ga)$gene %in% genes)
      } else {
        mt <- which(ga$gene %in% genes)
      }
    }
  }
  outlist <- ga[(narm(mt)),]
  if(typ=="GRanges") {
    if(!band) { mcols(outlist) <- mcols(outlist)[,-which(colnames(mcols(outlist)) %in% "band")]  }
    cnn <- colnames(mcols(outlist)); rnn <- mcols(outlist)$gene
  } else {
    if(!band) { outlist <- outlist[,-which(colnames(outlist) %in% "band")] }
    cnn <- colnames(outlist); rnn <- outlist$gene
  }
  if(one.to.one & ("gene" %in% cnn) & !anyDuplicated(genes)) {
    rownames(outlist) <- rnn
    if(all(genes %in% rownames(outlist))) {
      outlist <- outlist[match(genes,rnn),]
    }
    # outlist <- outlist[,-which(colnames(outlist) %in% "gene")]
  }
  if(bioC) { outlist <- as(outlist,"GRanges") }
  return(outlist)
}


#' Find the chromosome, start and end position for cytoband names
#' 
#' Allows retrieval of the the chromosome position of a karyotype/cytoband label, 
#'  or vector of such labels. Note that the position returned for bands is not a 
#'  single point as for SNPs, so the result will be a chromosome, then a position range with
#'  start and end, and lastly the band without the chromosome prefix
#' @param bands character, a vector of cytoband labels, chromosome[p/q]xx.xx ; 
#' e.g, 13q21.31, Yq11.221, 6p23, etc
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, else 
#' a data.frame
#' @return Returns a data.frame with columns 'chr' [chromosome], 'start' [starting position of the
#'  gene],'end' [end position of the gene] and 'band' [band without the chromosome prefix],
#'  or if bioC=TRUE, then returns a GRanges object with equivalent information.
#'  If returning a data.frame, then it will be in the same order as 'bands'. If bioC=TRUE, then
#'  the result will be in genome order, regardless of the order of 'bands'.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}, \code{\link{Pos}}, \code{\link{Pos.gene}}, \code{\link{Band}}, 
#' \code{\link{Band.gene}}, \code{\link{Band.pos}}, \code{\link{Gene.pos}}
#' @examples
#' setwd(tempdir())
#' Pos.band("1p13.2")
#' Pos.band("Yq11.221",build=36)
#' Pos.band("Yq11.221",build=37)
#' Pos.band(c("13q21.31","1p13.2","2q33.2","6p23"),bioC=TRUE)
Pos.band <- function(bands,build=NULL,dir=NULL,bioC=FALSE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  char.lim <- 100
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.cyto(build=build,bioC=bioC,dir=dir,GRanges=FALSE)
  mt <- match(bands,rownames(ga))
  failz <- paste(bands[is.na(mt)],collapse=", "); if(nchar(failz)>char.lim) { failz <- paste(substr(failz,1,char.lim),",...",sep="") }
  msg <- ("format for bands is: chromosome[p/q]xx.xx ; e.g, 13q21.31, Yq11.221, 6p23, etc")
  if(length(mt)<1 | all(is.na(mt))) { 
    warning("did not find any 'bands' features: ",failz) ; warning(msg); return(NULL) }
  if(any(is.na(mt))) { 
    cat("format for bands is: chromosome[p/q]xx.xx ; e.g, 13q21.31, Xq27.1, 6p23, etc")
    cnt <- length(which(is.na(mt)))
    warning("did not find the following ",cnt," 'bands' features: ",failz) ; warning(msg)
  }
  outlist <- ga[sort(mt[!is.na(mt)]),]
  if(any(colnames(outlist) %in% "negpos")) { outlist <- outlist[,-which(colnames(outlist) %in% "negpos")] }
  if(all(bands %in% rownames(outlist)) & !bioC) {
    outlist <- outlist[bands,]
  }
  return(outlist)
}


#' Retrieve the cytoband(s) for snp ids, genes or locations
#' 
#' Allows retrieval of the the cytoband/karyotype label, based on multiple
#'  possible input featues, including SNP chip or rs-ids, HGNC gene labels, GRanges or
#'  RangedData object, chromosome and position vectors. The most robust way to use the
#'  function is to use the parameter names to imply the type of input, e.g, use the 'genes'
#'  parameter to input gene labels, the 'snps' parameter to enter SNP ids, etc. However,
#'  if you enter the first argument as a GRanges or RangedData object instead of using the
#'  'ranges' argument, this will be detected and automatically moved to the 'ranges' parameter.
#' @param genes character, an optional vector of gene ids, or RangedData/GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the band from
#' @param ranges optional GRanges or RangedData object describing positions for which we want bands
#' @param snps optional SNP ids, e.g, chip ids or rs-ids, to retrieve the band they fall within
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param ... further arguments to Band.gene if entering gene names, or further arguments to 
#' Band.pos if entering ranges, or chr, pos/start/end
#' @return Returns a vector of bands, if any entries span more than one band, the bands will be
#' concatenated as character type, delimited by semicolons (;)
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}, \code{\link{Pos}}, \code{\link{Pos.gene}}, \code{\link{Band}}, 
#' \code{\link{Band.gene}}, \code{\link{Band.pos}}, \code{\link{Gene.pos}}
#' @examples
#' \donttest{
#' setwd(tempdir())
#' Band(chr=1,pos=1234567) # using chr,pos vectors
#' rd <- RangedData(ranges=IRanges(start=87654321,end=87654321),space=1)
#' gr <- as(rd,"GRanges")
#' Band(rd)    # using RangedData, autodetects this parameter should be 'ranges' not 'genes'
#' Band(ranges=gr) # using GRanges
#' Band("SLC6A4")  # serotonin gene [5-HTT]
#' a.few.snps <- c("rs3842724","imm_11_2147527","rs9467354")
#' Band(a.few.snps) # using SNP ids in the 'genes' parameter (still works!)
#' Band(snps=a.few.snps) # using SNP ids with the dedicated 'snps' parameter is quicker
#' Band(chr="X",pos=8000000)
#' # Band() with longer ranges  #
#' Band(chr=12,start=40000000,end=50000000,build="hg19") # concatenates if range spans multiple bands
#' Band(chr=12,start=40000000,end=50000000,build="hg18") # one extra band in the older annotation
#' }
Band <- function(genes=NULL,chr=NULL,ranges=NULL,snps=NULL,build=NULL,dir=NULL,...) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(all(is.character(genes))) { 
    Band.gene(genes=genes,build=build,dir=dir,...,warnings=FALSE)
  } else {
    if(!is.null(snps)) {
      pp <- Pos(snps); ch <- Chr(snps)
      Band.pos(chr=ch,pos=pp,build=build,dir=dir)
    } else {
      if((is(genes)[1] %in% c("RangedData","GRanges")) & is.null(ranges) ) { ranges <- genes }
      Band.pos(chr=chr,ranges=ranges,build=build,dir=dir,...)
    }
  }
}


#' Retrieve the cytoband(s) for genes labels
#' 
#' Allows retrieval of the the cytoband/karyotype label for HGNC gene labels.
#' @param genes character, an optional vector of gene ids, or RangedData/GRanges object
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param append.chr logical, it is typical that the chromosome character preceeds cytoband labels,
#'  but if this parameter is set to FALSE, it will be left off.
#' @param data.frame logical, if data.frame is true, instead of returning a vector of full cytoband
#' labels, a data.frame will be returned.
#' @param warnings logical, if warnings=FALSE and SNP ids are entered instead of Gene labels,
#' then the function will automatically detect this and return the result of Band(snps='genes')
#' @return Returns a vector of bands, if any entries span more than one band, the bands will be
#' concatenated as character type, delimited by semicolons (;). If data.frame is true, instead of 
#' returning a vector of full cytoband labels, a data.frame will be returned with a 'chr'
#' [chromosome] column, 'band' cytoband label  without the chromosome prefix, and rownames 
#' equal to 'genes'
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}, \code{\link{Pos}}, \code{\link{Pos.gene}}, \code{\link{Band}}, 
#' \code{\link{Band.gene}}, \code{\link{Band.pos}}, \code{\link{Gene.pos}}
#' @examples
#' \donttest{
#' setwd(tempdir())
#' a.few.snps <- c("rs3842724","imm_11_2147527","rs9467354")
#' Band.gene("HLA-C") # using chr,pos vectors
#' Band.gene(a.few.snps)  # fails with warning as these are SNPs, not genes
#' Band.gene(a.few.snps,warnings=FALSE) # with warnings=FALSE this continues with snps entered
#' Band.gene("SLC6A4")  # serotonin gene [5-HTT]
#' Band.gene("SLC6A4",append.chr=FALSE)
#' Band.gene("SLC6A4",data.frame=TRUE)
#' }
Band.gene <- function(genes,build=NULL,dir=getwd(),append.chr=TRUE,data.frame=FALSE,warnings=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  pg <- Pos.gene(genes=genes,build=build,dir=dir,bioC=F,band=TRUE,one.to.one=TRUE,warnings=warnings)
  ## a bit contradictory what i've done here with 'warnings' :
  if(length(pg)<1 & !warnings) { 
    #if(warnings) { warning("perhaps parameter 'genes' should have been 'snps'?") } 
    return(Band(snps=genes)) 
  }
  if(data.frame) {
    if(all(c("start","end") %in% colnames(pg))) { pg <- pg[,-which(colnames(pg) %in% c("start","end"))] }
    if(all(c("gene") %in% colnames(pg))) { rownames(pg) <- pg[["gene"]] ; pg <- pg[,-which(colnames(pg) %in% c("gene"))] }
    out <- pg
  } else {
    if(append.chr) {
      out <- paste(pg[["chr"]],pg[["band"]],sep="")
    } else {
      out <- pg[["band"]]
    }
  }
  return(out)    
}



#' Find the gene(s) overlapping a chromosome location
#' 
#' Allows retrieval of genes intersected by a chromosome and position, which can be entered
#' using chr, pos/start/end vectors, or a RangedData or GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the possible overlapping gene(s)
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' start or end if this is entered, and vice-versa
#' @param start integer, an optional vector of start points for chromosome ranges
#' @param end integer, an optional vector of end points for chromosome ranges
#' @param ranges optional GRanges or RangedData object describing positions for which we want genes,
#' removing the need to enter chr, pos, start or end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, or
#' RangedData if 'ranges' is RangedData, else a data.frame
#' @param one.to.one logical, whether to concatenate multiple hits for the same range into one result,
#' or spread the result over multiple lines, one for each gene overlapped
#' @return Returns a set of genes separated by semicolons (if more than one) for each range entered.
#' If bioC=TRUE, returns the equivalent as a GRanges object, unless a RangedData object was used
#' for the ranges parameter, in which case a RangedData object would be returned. If one.to.one is
#' FALSE, then instead of concatenating multiple genes into one line per range, each is listed 
#' separately as a new row, with an index added to correspond to the original input order of ranges,
#' if bioC=TRUE; or just adds additional elements to the resulting vector if bioC=FALSE.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}, \code{\link{Pos}}, \code{\link{Pos.band}}, \code{\link{Band}}, 
#' \code{\link{Band.gene}}, \code{\link{Band.pos}}, \code{\link{Gene.pos}}
#' @examples
#' \donttest{
#' setwd(tempdir())
#' Gene.pos(chr=6, start=31459636, end=31462760)
#' Gene.pos(chr=22, pos=3452345) # no gene here
#' Gene.pos(Chr("rs689"),Pos("rs689")) # combine with Chr() and Pos() to find gene(s) for SNP rs689
#' Gene.pos(chr=1,start=114000000,end=115000000,build="hg19") # multiple genes in range
#' Gene.pos(chr=1,start=114000000,end=115000000,one.to.one=FALSE) # list separately
#' ii <- Pos.gene(c("CTLA4","PTPN22"))
#' Gene.pos(ii$chr,ii$start,ii$end,bioC=FALSE) # returns same genes inputted on line above
#' }
Gene.pos <- function(chr=NA,pos=NA,start=NA,end=NA,ranges=NULL,build=NULL,dir=NULL,bioC=FALSE,one.to.one=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is(chr)[1] %in% c("RangedData","GRanges"))  { ranges <- chr } # in case first parameter is used
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) {
    if(any(!is.na(pos))) { start <- pos; end <- pos }
    if(length(chr)==1 & length(start)>1) { chr <- rep(chr,times=length(start)) }
    if(length(chr)!=length(start)) { stop("chr vector must have same length as pos or start/end") }
    if(any(is.na(chr))) { stop("cannot have chr=NA") }
    Pos <- matrix(ncol=2,nrow=length(start))
    for (cc in 1:length(start)) {
      Pos[cc,] <- force.chr.pos(Pos=c(start[cc],end[cc]),Chr=chr[cc],dir=dir,build=build)
    }
    #if(any(tolower(substr(chr,1,3))!="chr")) { chr <- gsub("chr",chr,sep="") }
    #chr <- gsub("chrchr","chr",chr)
    testData <- RangedData(ranges=IRanges(start=Pos[,1],end=Pos[,2]),space=chr,index=1:length(chr),universe=build[1])
    testData <- toGenomeOrder2(testData,strict=T)
  } else {
    testData <- ranges # set.chr.to.char(ranges)
    if(typ=="GRanges") { testData <- as(ranges,"RangedData") }
    #chr <- chr(testData)
    if(!bioC) { warning("bioC was set false, but ranges argument in use so overriding") ; bioC <- T }
    if("index" %in% colnames(testData)) { 
      warning("'index' is a reserved column name for ranges objects passed to this function so will be replaced. Consider renaming this column if this is undesired") }
    testData[["index"]] <- 1:nrow(testData)
  }
  testData <- set.chr.to.numeric(testData,keep=T)
  #return(testData)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- suppressWarnings(get.gene.annot(build=build,dir=dir,GRanges=FALSE))
  #ga <- set.chr.to.numeric(ga,keep=F)
  newDataList <- vector("list",nrow(testData))
  overlaps <- findOverlaps(testData,ga)
  genez <- ga$gene[subjectHits(overlaps)]
  indexz <- queryHits(overlaps)
  #prv(indexz,genez)
  if(length(indexz)<1) { return(NA) }
  if(!one.to.one) {
    newData <- testData[queryHits(overlaps),]
    newData[["gnm.index"]] <- indexz
    newData[["gene"]] <- genez
  } else {
    out <- tapply(genez,factor(indexz),c,simplify=FALSE)
    out <- sapply(out,function(X) { X <- narm(unique(X)); X <- X[X!=""] ; paste(X,collapse=";") })
    newData <- testData
    newData[["gene"]] <- rep("intergenic",nrow(newData))
    if(!is.null(names(out))) {
      newData[["gene"]][as.numeric(names(out))] <- out
    } else {
      newData[["gene"]] <- out
    }
  }
  if(bioC) {
    if(typ!="RangedData") { newData <- as(newData,"GRanges") }
    return(newData)
  } else {
    OO <- newData[["gene"]][order(newData[["index"]])]
    OO <- narm(unique(OO)); OO <- OO[OO!=""]
    return(OO)
  }
  #if(!all(chr %in% chr2(ga))) { stop("invalid chromosome(s) entered") } # redundant i think
}





#' Find the cytoband(s) overlapping a chromosome location
#' 
#' Allows retrieval of cytobands/karyotypes intersected by a chromosome and position, which can be 
#' entered using chr, pos/start/end vectors, or a RangedData or GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the possible overlapping cytoband(s)
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' start or end if this is entered, and vice-versa
#' @param start integer, an optional vector of start points for chromosome ranges
#' @param end integer, an optional vector of end points for chromosome ranges
#' @param ranges optional GRanges or RangedData object describing positions for which we want bands,
#' removing the need to enter chr, pos, start or end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, or
#' RangedData if 'ranges' is RangedData, else a data.frame
#' @param one.to.one logical, whether to concatenate multiple hits for the same range into one result,
#' or spread the result over multiple lines, one for each cytoband overlapped
#' @return Returns a set of cytobands separated by semicolons (if more than one) for each range entered.
#' If bioC=TRUE, returns the equivalent as a GRanges object, unless a RangedData object was used
#' for the ranges parameter, in which case a RangedData object would be returned. If one.to.one is
#' FALSE, then instead of concatenating multiple cytobands into one line per range, each is listed 
#' separately as a new row, with an index added to correspond to the original input order of ranges,
#' if bioC=TRUE; or just adds additional elements to the resulting vector if bioC=FALSE.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}, \code{\link{Pos}}, \code{\link{Pos.band}}, \code{\link{Band}}, 
#' \code{\link{Band.gene}}, \code{\link{Band.pos}}, \code{\link{Gene.pos}}
#' @examples
#' \donttest{
#' setwd(tempdir())
#' Band.pos(chr=6, start=31459636, end=31462760)
#' Band.pos(chr=22, pos=3452345) 
#' Band.pos(Chr("rs689"),Pos("rs689")) # combine Chr(), Pos() to find the cytoband for SNP rs689
#' Band.pos(chr=1,start=110000000,end=120000000,build="hg19") # multiple cytobands in range
#' Band.pos(chr=1,start=110000000,end=120000000,one.to.one=FALSE) # list separately
#' Band.pos(Pos.band(c("13q21.31","1p13.2"),bioC=TRUE)) # use ranges object returned by Pos.band()
#' # note that 3 ranges are returned for each entry as the start/end overlap the adjacent ranges
#' }
Band.pos <- function(chr=NA,pos=NA,start=NA,end=NA,ranges=NULL,build=NULL,dir=NULL,bioC=FALSE,one.to.one=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is(chr)[1] %in% c("RangedData","GRanges"))  { ranges <- chr } # in case first parameter is used
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) {
    if(any(!is.na(pos))) { start <- pos; end <- pos }
    if(length(chr)==1 & length(start)>1) { chr <- rep(chr,times=length(start)) }
    if(length(chr)!=length(start)) { stop("chr vector must have same length as pos or start/end") }
    if(any(is.na(chr))) { stop("cannot have chr=NA") }
    Pos <- matrix(ncol=2,nrow=length(start))
    for (cc in 1:length(start)) {
      Pos[cc,] <- force.chr.pos(Pos=c(start[cc],end[cc]),Chr=chr[cc],dir=dir,build=build)
    }
    #if(any(tolower(substr(chr,1,3))!="chr")) { chr <- gsub("chr",chr,sep="") }
    #chr <- gsub("chrchr","chr",chr)
    testData <- RangedData(ranges=IRanges(start=Pos[,1],end=Pos[,2]),space=chr,index=1:length(chr),universe=build[1])
    testData <- toGenomeOrder2(testData,strict=T)
  } else {
    if(typ=="GRanges") { ranges <- as(ranges,"RangedData") }
    testData <- ranges # set.chr.to.char(ranges)
    if("index" %in% colnames(testData)) { warning("'index' is a reserved column name for ranges objects passed to this function so will be replaced. Consider renaming this column if this is undesired") }
    testData[["index"]] <- 1:nrow(testData)
  }
  testData <- set.chr.to.numeric(testData)
  #return(testData)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  cyto <- get.cyto(build=build,bioC=TRUE,dir=dir,GRanges=FALSE)
  cyto <- set.chr.to.numeric(cyto,keep=F)
  #if(any(colnames(cyto) %in% "negpos")) { cyto <- cyto[,-which(colnames(cyto) %in% "negpos")] }
  newDataList <- vector("list",nrow(testData))
  overlaps <- findOverlaps(testData,cyto)
  bandz <- rownames(cyto)[subjectHits(overlaps)]
  indexz <- queryHits(overlaps)
  if(length(indexz)<1) { return(NA) }
  if(!one.to.one) {
    newData <- testData[queryHits(overlaps),]
    newData[["gnm.index"]] <- indexz
    newData[["band"]] <- bandz
  } else {
    out <- tapply(bandz,factor(indexz),c,simplify=FALSE)
    if(one.to.one) { 
      out <- sapply(out,function(X) { paste(X,collapse=";") }) 
      newData <- testData
      #newData[["band"]] <- out 
      newData[["band"]] <- rep("",nrow(newData))
      if(!is.null(names(out))) {
        newData[["band"]][as.numeric(names(out))] <- out
      } else {
        newData[["band"]] <- out
      }
    }
  }
  if(bioC) {
    if(typ!="RangedData") { newData <- as(newData,"GRanges") }
    return(newData)
  } else {
    return(newData[["band"]][order(newData[["index"]])])
  }
  #if(!all(chr %in% chr2(ga))) { stop("invalid chromosome(s) entered") } # redundant i think
}



#' Returns the A and B allele for SNP ids
#' 
#' For a set of chip ids or rs ids, returns a two column matrix containing the A and B allele. 
#' For snpStats objects the default is that A,B are coded in alphabetical order, so A,C; A,T; 
#' C,T; C,G are possible A,B pairs. Allele codes are specific to each dataset, so you should
#' upload your allele codes into the current ChipInfo object to make the alleles produced by
#' this function meaningful.
#' @param ids character, a list of chip ids or rs-ids as contained in the current ChipInfo object
#' @return Returns a two column matrix containing the A and B allele.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso \code{\link{Chr}}, \code{\link{Pos}}, \code{\link{Pos.band}}, \code{\link{Band}}, 
#' \code{\link{Band.gene}}, \code{\link{Band.pos}}, \code{\link{Gene.pos}}
#' @examples
#' \donttest{
#' snp.ids <- c("rs3842724","rs9729550","rs1815606","rs114582555","rs1240708","rs6603785")
#' AB(snp.ids) 
#' }
AB <- function(ids) {
  all.support <- chip.support()
  ic.ab <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { all.support <- chip.support() }  ## load object: all.support [snp support for whole chip]
    outlist <- cbind(A1(all.support)[match(ic.ids,rownames(all.support))],A2(all.support)[match(ic.ids,rownames(all.support))])
    return(outlist)
  }
  out <- matrix(nrow=length(ids),ncol=2)
  ic.ab.id <- ic.ab(rs.to.id(ids))
  return(ic.ab.id)
}


#' Retrieve SNP ids or positions in specified range
#' 
#' This function will always use the build in getOption('ucsc'), so use options() if it needs to 
#' change.
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' Alternatively chr can be a RangedData or GRanges object in which case SNP lists will be returned
#' in a list for each row of the ranges object.
#' @param start integer, genomic position to define the start of the range to look for SNPs,
#'  should be between 1 and the length of the chromosome 'chr'
#' @param end integer, genomic position to define the end of the range to look for SNPs,
#'  should be between 1 and the length of the chromosome 'chr', and >= start
#' @param ids logical, if TRUE will return snp ids (chip ids, for rs-ids, use id.to.rs on the output), 
#' or if FALSE will return the chromosome positions of the SNPs.
#' @export
#' @return Set of SNP ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr', that
#' fall within the genomic range described by the chr, start, and end parameters. Alternatively, if
#' chr is a RangedData or GRanges object then multiple SNP lists will be returned
#' in a list for each row of the ranges object.
#' @examples
#' snps.in.range(1,9000000,10000000)
#' snps.in.range(10,19000000,20000000,ids=TRUE)
#' snps.in.range(10,19000000,20000000,ids=FALSE) # return positions instead of rs-ids
snps.in.range <- function(chr, start=NA, end=start, ids=TRUE) { 
  # ids - whether to return ichip SNP ids or positions
  if(is(chr)[1]=="RangedData" | is(chr)[1]=="GRanges") {
    chrz <- chr2(chr); stz <- start(chr); enz <- end(chr)
    output <- vector("list",nrow(chr))
    for(cc in 1:nrow(chr)) {
      output[[cc]] <- snps.in.range(chrz[cc],stz[cc],enz[cc],ids=ids)
    }
    names(output) <- rownames(chr)
    return(output)
  }
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(start)>1) { warning("start should be length 1, using only first entry"); start <- start[1] }
  if(length(end)>1) { warning("end should be length 1, using only first entry"); end <- end[1] }
  if(start>end) { warning("start was higher than end, so switching") }
  the.range <- sort(c(start,end))
  all.support <- chip.support()
  #if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
 # if(!exists("all.support")) { all.support <- chip.support() }  ## load object: all.support [snp support for whole chip]
  all.chr <- chr(all.support)
  all.pos <- start(all.support)[all.chr %in% chr]
  if(length(all.pos)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  validz <- which(all.pos>=the.range[1] & all.pos<=the.range[2])
  if(ids) {
    out <- rownames(all.support)[all.chr %in% chr][validz]
  } else {
    out <- start(all.support)[(all.chr %in% chr)][validz]
  }
  return(out)
}


#' Retrieve the 'n' closest SNP ids or positions near specified locus
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' @param pos integer, genomic position, should be between 1 and the length of the chromosome 'chr'
#' @param n integer, the number of nearest SNPs to seek, if there aren't enough in the annotation
#' then NAs will fill the gaps to force the return value length to equal 'n'
#' @param side character, can be 'either', 'left' or 'right' and specifies which side of the 'pos'
#' to look for nearest snps (where left is decreasing genomic position and right is increasing)
#' @param ids logical, if TRUE will return snp ids (chip ids, for rs-ids, use id.to.rs on the output), 
#' or if FALSE will return the chromosome positions of the SNPs.
#' @param limit integer, a limit on the maximum distance from the position 'pos' can be specified
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @export
#' @seealso \code{\link{expand.nsnp}}, \code{\link{nearest.gene}}
#' @return Set of SNP ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr'.
#' If the number of SNPs on the chromosome or the bounds of the 'side' and 'limit' parameters
#' restrict the number returned to less than 'n' then the return value will be padded with NAs.
#' @examples
#' nearest.snp(1,159000000,n=10) # return ids
#' nearest.snp(1,159000000,n=10,build=37)
#' nearest.snp(1,159000000,n=10,build=36,ids=FALSE) # return positions
#' \donttest{
#' nearest.snp(1,159000000,n=10,build=37,ids=FALSE)
#' nearest.snp(6,25000000,n=10,build=37,ids=FALSE,side="left")  # only SNPs to the left of the locus
#' nearest.snp(6,25000000,n=10,build=37,ids=FALSE,side="right") # only SNPs to the right of the locus
#' }
nearest.snp <- function(chr, pos, n=1, side=c("either","left","right"),ids=TRUE,limit=NULL,build=NULL) { 
  # ids - whether to return ichip SNP ids or positions
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(pos)>1) { warning("pos should be length 1, using only first entry"); pos <- pos[1] }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  all.support <- chip.support(build=build,warn.build=FALSE)
  if(!exists("all.support")) { all.support <- chip.support() }  ## load object: all.support [snp support for whole chip]
  side <- tolower(side[1]); 
  if(!side %in% c("either","left","right")) {
    side <- "either"; warning("invalid side argument, defaulting to 'either'") }
  if(!is.null(limit)) { if(!is.numeric(limit)) { limit <- NULL; warning("invalid limit argument, defaulting to NULL") } }
  all.chr <- chr(all.support)
  all.pos <- start(all.support)[all.chr %in% chr]
  if(length(all.pos)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  difz <- pos-all.pos
  all.true <- difz==difz
  if(is.null(limit)) { lfilt <- all.true } else { lfilt <- abs(difz)<=limit }
  if(side=="left") { filt <- difz>0 & lfilt }
  if(side=="right") { filt <- difz<0 & lfilt }
  if(side=="either") { filt <- all.true & lfilt }
  Difz <- abs(difz[filt])
  if(length(Difz)<n)  { warning("fewer than ",n," positions found for 'chr' specified (within 'limit'), NAs returned") }
  indx <- order(Difz)[1:n]
  if(ids) {
    out <- rownames(all.support)[all.chr %in% chr][filt][indx]
  } else {
    out <- start(all.support)[(all.chr %in% chr)][filt][indx]
  }
  return(out)
}



#' Obtain nearby SNP-lists within a recombination window
#' 
#' For a snp.id (or list), extend a window around that chromosome location in recombination
#' units (centimorgans) and return the list of SNPs from the current ChipInfo object that
#' lie in this window. This is a way of extracting SNPs in linkage disequilibrium with an 
#' index SNP, that could also be plausible causal candidates. Runs fastest for build 36,
#' otherwise internal conversion takes place (runs using build based on getOptions('ucsc')).
#' @param snpid.list character, list of snp-ids (e.g, rs-id or chip id) to obtain lists for.
#' SNPs must all be from the same chromosome - if ranges for SNPs spanning multiple ranges
#' are desired, you must use multiple calls. A warning will be given if SNPs from the same
#' karyotype band are entered as index SNPs, as in a typical GWAS analysis only one SNP would
#' be used like this from each region, ignore the warning if this is not the case for your
#' application.
#' @param cM numeric, the number of centimorgans to extend the window either side of each SNP
#' @param bp.ext numeric, optional number of base-pairs to extend the window by in addition
#' to the centimorgan extension
#' @param excl.snps character, a list of rs-id or chip-ids of SNPs to exclude from the list
#' returned, as, for instance, they may have failed quality control such as call-rate.
#' @param name.by.bands give labels to each sublist returned by the karotype/cytoband name, 
#' but faster not to do this
#' @return Returns a list of vectors of snp-ids falling within the window(s) specified and not
#' in 'excl.snps'. Each snp in 'snpid.list' will correspond to an element in the list returned.
#' If name.by.bands is TRUE, then these list elements will each be named using the local
#' karyotype/cytoband location
#' @export
#' @seealso \code{\link{snps.in.range}}, \code{\link{get.recombination.map}}, 
#' \code{\link{recomWindow}}, \code{\link{conv.37.36}}, \code{\link{conv.36.37}}, \code{\link{expand.nsnp}}
#' @examples
#' # examples not run as too slow
#' \donttest{
#' result <- get.nearby.snp.lists("rs900569")
#' # trick below to extract SNPs within 0.1-0.2cM
#' get.nearby.snp.lists("rs900569",cM=0.2,excl.snps=result[[1]]) 
#' # note that the same query can return a different set with build 36 versus 37
#' get.nearby.snp.lists(c("rs689","rs4909944"),cM=0.001,name.by.bands=FALSE) 
#' }
get.nearby.snp.lists <- function(snpid.list,cM=0.1,bp.ext=0,excl.snps=NULL,name.by.bands=TRUE) {
  #if(!exists("all.support")) { print(load("all.support.RData")) }
  all.support <- chip.support()
#  if(is.null(build)) { build <- getOption("ucsc") }
#  build <- ucsc.sanitizer(build)
  build <- ucsc.sanitizer(getOption("ucsc"))
  if(!build %in% c("hg19","hg18","hg38")) { stop("only builds 36,37,38 are supported") }
  snpic.list <- rs.to.id(snpid.list)
  if(length(snpic.list)<1) { stop("snpic.list must contain at least 1 id")}
  cyto <- get.cyto(dir=getwd(),GRanges=FALSE); cyto[["gene"]] <- rownames(cyto)
  #which.snps <- match(snpid.list,mcols(all.support)$rs.id)
  #if(any(is.na(which.snps))) { stop(paste("NAs in dbSNP match:",paste(snpid.list[is.na(which.snps)],collapse=","))) }
  snps.locs <- Pos(snpid.list)
  snps.chrs <- Chr(snpid.list)
  if(build=="hg18") {
    snps.locs36 <- snps.locs;  
    if(length(snps.locs)>1) {
      snps.locs37 <- conv.36.37(chr=snps.chrs,pos=snps.locs)[,"start"] 
    } else {
      snps.locs37 <- conv.36.37(chr=snps.chrs,pos=snps.locs)["start"] 
    }
  } else { 
    if(build=="hg38") {
      snps.locs38 <- snps.locs; 
      if(length(snps.locs)>1) {
        snps.locs37 <- conv.38.37(chr=snps.chrs,pos=snps.locs)[,"start"]
        snps.locs36 <- conv.37.36(chr=snps.chrs,pos=snps.locs)[,"start"]
      } else {
        snps.locs37 <- conv.38.37(chr=snps.chrs,pos=snps.locs)["start"]
        snps.locs36 <- conv.37.36(chr=snps.chrs,pos=snps.locs)["start"]
      }
    } else {
      # assume hg19
      snps.locs37 <- snps.locs; 
      if(length(snps.locs)>1) {
        snps.locs36 <- conv.37.36(chr=snps.chrs,pos=snps.locs)[,"start"]
      } else {
        snps.locs36 <- conv.37.36(chr=snps.chrs,pos=snps.locs)["start"]
      }
    }
  }
  next.chr <- unique(Chr(snpid.list)); if(length(next.chr)>1) { stop("enter snpids from only 1 chromosome at a time!") }
  if(any(snps.locs!=sort(snps.locs))) { 
    warning("snp-ids not in position order, rearrangement is preferred but will attempt to continue")
    sort.back <- match(snps.locs,sort(snps.locs))
  } else { sort.back <- 1:length(snps.locs) }
  ddz <- snpic.list[duplicated(snpic.list)]
  if(length(ddz)>0) { warning("dup SNPs:",ddz,"\n") }
  #prv(snps.locs,snps.locs36)
  snp.rd <- RangedData(ranges=IRanges(start=snps.locs,end=snps.locs,names=snpic.list),
                       space=rep(next.chr,length(snps.locs)))
  snp.rd <- toGenomeOrder2(snp.rd,strict=T) # think it autosorts anyway, but just in case
  if(name.by.bands) {
    #snp.rd <- annot.cnv(snp.rd,gs=cyto,quiet=TRUE); colnames(snp.rd) <- "band"
    bands <- Band.pos(ranges=snp.rd, build=build) #   snp.rd$band
  }
  #prv(next.chr,cM,bp.ext)
  ## recomWindow uses build36 only, so convert back afterwards
  nxt.window <- lapply(snps.locs36, function(X,...) { recomWindow(start=as.numeric(X),...) },
                       chr=next.chr,window=cM,bp.ext=bp.ext,info=FALSE)
  if(build=="hg18") {
    st.window <- sapply(nxt.window, "[",1)
    en.window <- sapply(nxt.window, "[",2)
  } else {
    #prv(next.chr,nxt.window)
    nncc <- rep(next.chr,length(nxt.window))
    if(length(snps.locs)>1) {
      st.window <- conv.36.37(chr=nncc,pos=sapply(nxt.window, "[",1))[,"start"]
      en.window <- conv.36.37(chr=nncc,pos=sapply(nxt.window, "[",2))[,"start"]
    } else {
      st.window <- conv.36.37(chr=nncc,pos=sapply(nxt.window, "[",1))["start"]
      en.window <- conv.36.37(chr=nncc,pos=sapply(nxt.window, "[",2))["start"]
    }
    if(build=="hg38") {
      st.window <- conv.37.38(chr=nncc,pos=st.window)
      en.window <- conv.37.38(chr=nncc,pos=en.window)
    }
  }
  pozz <- start(all.support)
  n.snps <- vector("list",length(st.window))
  for(cc in 1:length(st.window)) {
    n.snps[[cc]] <- which(chr(all.support)==next.chr &
                            pozz>=st.window[cc] & 
                            pozz<=en.window[cc] &
                            (!rownames(all.support) %in% excl.snps) &
                            (!mcols(all.support)$rs.id %in% excl.snps) 
    )
  }
  grp.labs <- lapply(n.snps,function(X) { rownames(all.support)[X] })
  if(name.by.bands) {
    if(length(unique(bands))!=length(bands)) { warning("these bands are not unique ==> ",paste(bands[duplicated(bands)],collapse=",")) }
    grpz <- 1:length(bands)
    names(grp.labs) <- paste(grpz,bands,sep=":")
  }
  grp.labs <- grp.labs[sort.back]
  return(grp.labs)
}




################## end support ##########################



