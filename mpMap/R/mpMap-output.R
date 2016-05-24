#' Output mpcross objects to other file formats
#'
#' Outputs the genotype information from an 'mpcross' object to files which can be read in to either R/qtl cross format or R/happy.hbrem format
#' @rdname mpMap-output
#' @aliases write2cross write2happy write.mpcross
#' @export
#' @usage write2cross(object, filestem, chr, ...)
#' write2happy(object, filestem, chr, ...)
#' write.mpcross(object, filestem="mp", chr, format=c("qtl", "happy"), ...)
#' @param object Object of class \code{mpcross}
#' @param filestem Filestem for all files which are output; may include directory
#' @param chr Subset of chromosomes to be output; default is all
#' @param format Output format - R/qtl or happy.hbrem
#' @param ... Additional arguments
#' @return 
#' For R/qtl cross format, two files are output: filestem.ril.csv and filestem.founder.csv which can be then read into the R/qtl package using 
#' \code{\link[qtl]{readMWril}}. 

#' For R/happy.hbrem format, two files are output for each chromosome: e.g. filestem.chr1.alleles and filestem.chr1.data. These can be directly input to \code{happy}.
#' @seealso \code{\link[qtl]{readMWril}}, \code{happy}

write.mpcross <- 
function(object, filestem="mp", chr, format=c("qtl", "happy"), ...)
{
  # creates output files which can be input to R/qtl
  if (missing(object))
	stop("Missing a required argument for this function")
  if (missing(filestem))
	stop("Missing a required argument for this function")

  if (missing(format)) 
	stop("Must give a format for output files")
  
  if (missing(chr))
	chr <- c(1:length(object$map))

  ### 
  if (format=="qtl")
	write2cross(object, filestem, chr, ...)

  if (format=="happy")
	write2happy(object, filestem, chr, ...)
}

write2cross <- 
function(object, filestem, chr, ...)
{
  # creates output files which can be input to R/qtl
  if (missing(object))
	stop("Missing a required argument for this function")
  if (missing(filestem))
	stop("Missing a required argument for this function")

  if (is.null(object$pheno)) {
    cat("Must have phenotype to output Rqtl files - inserting random values")
    object$pheno$pheno <- rnorm(nrow(object$finals))
    object$pheno <- as.data.frame(object$pheno) 
  }
  if (is.null(object$map))
    stop("Must have map to output Rqtl files")

  if (missing(chr))
	chr <- c(1:length(object$map))

  if (is.character(chr)) chr <- match(chr, names(object$map))

  obj <- subset(object, chr=chr)

  n.founders <- nrow(obj$founders)
  n.pheno <- ncol(as.matrix(obj$pheno))

  strains <- LETTERS[1:n.founders]

  rilfile <- paste(filestem, ".ril.csv", sep="")
  founderfile <- paste(filestem, ".founder.csv", sep="")

  # rilfile needs to include a phenotype cross with funnels
  # construct cross phenotype
  cross <- apply(CR_cross(obj), 1, function(x) return(paste(strains[x], collapse="")))

#sapply(object$id, function(x) return(paste(strains[as.numeric(pedtofun(object$pedigree,x))],collapse="")))
  chrnam <- rep(names(obj$map), unlist(lapply(obj$map, length))) 
 
  ril <- data.frame(obj$pheno, cross, obj$finals)
  ril[,n.pheno+1] <- as.character(ril[,n.pheno+1])
  names(ril)[(n.pheno+2):ncol(ril)] <- colnames(obj$finals) 

  vec <- c(rep("", n.pheno+1), as.character(chrnam))
  names(vec) <- colnames(ril)
  write.csv(t(vec), rilfile, quote=FALSE, row.names=FALSE)
  write.table(t(c(rep("", n.pheno+1), unlist(obj$map))), rilfile, col.names=FALSE, quote=FALSE, row.names=FALSE, append=TRUE, sep=",")

#  write.csv(ril, rilfile, row.names=FALSE, quote=FALSE)
  write.table(ril, rilfile, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")

  write.csv(cbind(strains,obj$founders), founderfile, row.names=FALSE, quote=FALSE)
}

write2happy <- 
function(object, filestem, chr, ...)
{
  # creates output files which can be input to HAPPY
  if (missing(object))
	stop("Missing a required argument for this function")
  if (missing(filestem))
	stop("Missing a required argument for this function")

  if (is.null(object$pheno))
	stop("Must have phenotype to output HAPPY files")
  if (is.null(object$map))
	stop("Must have map to output HAPPY files")

  if (missing(chr))
	chr <- c(1:length(object$map))

  n.founders <- nrow(object$founders)

  datafile <- paste(filestem, ".data", sep="")
  allelesfile <- paste(filestem, ".alleles", sep="")

  write(paste("markers", length(unlist(object$map)), "strains", n.founders), allelesfile)
  write(paste("strain_names", paste(object$fid, collapse=" ")), allelesfile, append=TRUE)

    mrkchr <- 1:ncol(object$finals)
    elem <- list(id=object$id, phe=round(object$pheno,4), gen=object$finals[,mrkchr[as.numeric(gl(length(mrkchr),2))]])
    datamat <- as.data.frame(do.call(cbind, elem))
    write.table(datamat, datafile, row.names=FALSE, col.names=FALSE)

  for (j in chr)
  {
    mrkchr <- match(names(object$map[[j]]), colnames(object$founders))
    mrkchr <- mrkchr[!is.na(mrkchr)]

    for (i in 1:length(mrkchr))
    {
 	allfreq <- names(table(object$founders[,mrkchr[i]]))
	n.alleles <- length(allfreq)+1
    	write(paste("marker", colnames(object$founders)[mrkchr[i]], n.alleles, j, object$map[[j]][i]), allelesfile, append=TRUE)
	# deal with missing data
	write(paste("allele", "NA", paste(rep(1/n.founders, n.founders), collapse=" ")), allelesfile, append=TRUE)
	for (k in 1:(n.alleles-1))
	  write(paste("allele", allfreq[k], paste(as.numeric(object$founders[,mrkchr[i]]==as.numeric(allfreq[k]))/sum(object$founders[,mrkchr[i]]==allfreq[k]), collapse=" ")), allelesfile, append=TRUE)
    }
  }

}
