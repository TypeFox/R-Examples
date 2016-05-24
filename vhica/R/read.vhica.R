read.vhica <-
function (gene.fasta=NULL, target.fasta=NULL, cb.filename=NULL, div.filename=NULL, 
	reference = "Gene", divergence = "dS", CUB.method="ENC", div.method="LWL85", div.pairwise=TRUE, div.max.lim=3, species.sep="_", gene.sep=".", family.sep=".", ...) 
{
	stopifnot( 
		!(is.null(gene.fasta) && is.null(target.fasta)) || 
		!(is.null(cb.filename) && is.null(div.filename)))
    vhica.obj <- list()
    if (!is.null(gene.fasta)) {
		vhica.obj$cbias <- 
			.seq.codon.bias(gene.fasta=gene.fasta, target.fasta=target.fasta, method=CUB.method, species.sep=species.sep, family.sep=family.sep)
		vhica.obj$div <- 
			.seq.divergence(sequence.fasta=c(gene.fasta, target.fasta), method=div.method, pairwise=div.pairwise, max.lim=div.max.lim, species.sep=species.sep, family.sep=family.sep)
		if (!is.null(cb.filename))
			write.table(vhica.obj$cbias, file=cb.filename, sep="\t", quote=FALSE, row.names=TRUE)
		if (!is.null(div.filename)) 
			write.table(vhica.obj$div, file=div.filename, sep="\t", quote=FALSE, row.names=FALSE)
    } else {
		vhica.obj$cbias <- .read.codon.bias(file = cb.filename, reference = reference)
		vhica.obj$div <- .read.divergence(file = div.filename, divergence = divergence)    
    }
    vhica.obj$reg <- .reference.regression(vhica.obj$cbias, vhica.obj$div, 
        reference = reference, divergence = divergence, family.sep=family.sep, ...)
    vhica.obj$reference <- reference
    tmp.target <- levels(vhica.obj$cbias[, "Type"])
    vhica.obj$target <- tmp.target[tmp.target != reference][1]
    vhica.obj$divergence <- divergence
    vhica.obj$family.sep=family.sep
    class(vhica.obj) <- c("vhica", class(vhica.obj))
    return(vhica.obj)
}
