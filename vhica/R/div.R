div <- function
(file=NULL, sequence=NULL, sqs=NULL, method="LWL85", pairwise=TRUE, max.lim=3)
{
	stopifnot(
		!(is.null(file) && is.null(sequence)),
		method[1] %in% c("LWL85"),
		requireNamespace("gtools", quietly=TRUE))

	if (!is.null(file)) {
		if (!requireNamespace("seqinr", quietly=TRUE)) {
			stop("Reading FASTA files require package seqinr")
		}		
		sequence <- seqinr::read.fasta(file)
	}
	sequence <- .checkseq(sequence, gene.name=if (is.null(file)) "" else file)
	if (is.null(sqs)) {
		sqs <- names(sequence)
	}
	combn <- gtools::combinations(n=length(sqs), r=2, v=sqs)
	if (method[1]=="LWL85") {
		return(data.frame(div=.LWL85(sequence, combn[,1], combn[,2], pairwise=pairwise, max.lim=max.lim), sp1=combn[,1], sp2=combn[,2]))
	} 
	stop("Method ", method, " unknown.") # This should never happen
}
