# The bedr package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

download.datasets <- function(datasets = "all", data.dir = paste0(Sys.getenv("HOME"),"/bedr/data")) {

	if (datasets == "all") {
		datasets <- c("refgene","hg19","b37","hugo", "cosmic","clinvar", "agilent", "nimblegen","gap");
		}

	old.dir <- getwd();
	setwd(data.dir);

	# get the fasta files from the broads ftp.  thank you!
	if ("hg19" %in% datasets && !file.exists(paste0(data.dir,"/ucsc.hg19.fasta.gz"))) {
		download.file("ftp://ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.gz",destfile="ucsc.hg19.fasta.gz", extra="--user gsapubftp-anonymous", method = "wget")
		download.file("ftp://ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.fai.gz",destfile="ucsc.hg19.fasta.fai.gz", extra="--user gsapubftp-anonymous", method = "wget")
		}

	if ("b37" %in% datasets && !file.exists(paste0(data.dir, "/human_g1k_v37.fasta.gz"))) {
		download.file("ftp://ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta.gz",destfile="human_g1k_v37.fasta.gz", extra="--user gsapubftp-anonymous", method = "wget")
		download.file("ftp://ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta.fai.gz",destfile="ucsc.hg19.fasta.fai.gz", extra="--user gsapubftp-anonymous", method = "wget")
		}

	# get the gene data from the 
	if ("refgene" %in% datasets && !file.exists(paste0(data.dir, "/refGene.txt.gz"))) {
		download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz", destfile="refGene.txt.gz");
		download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.sql", destfile="refGene.sql");
		}

	# parse refgene
	if ("refgene" %in% datasets && !file.exists(paste0(data.dir, "/gene_regions.txt"))) {
		refgene    <- query.ucsc(paste0(data.dir, "/refGene.txt.gz"));
		# only select pr coding genes
		refgene.nm <- refgene[grepl("^NM",refgene$name),c("chrom","txStart","txEnd","name2","strand")]
		# only include canonical chr
		chr <- paste0("chr", c(1:22,"X","Y"));
		refgene.nm <- refgene.nm[refgene.nm[,1] %in% chr,]
		# sort and merge
		refgene.nm <- snm(refgene.nm,check.chr = FALSE);
		# remove genes with multiple positions
		duplicated.gene <- duplicated(refgene.nm[,4]) | duplicated(rev(refgene.nm[,4]));
		refgene.nm <- refgene.nm[!duplicated.gene,];
		write.table(refgene.nm, paste0(data.dir, "/gene_regions.txt"), sep = "\t", quote = FALSE, row.names = FALSE);
		}

	# centromeres/telomeres and problem regions	
	if ("gap" %in% datasets && !file.exists(paste0(data.dir, "/gap.txt.gz"))) {
		download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz", destfile="gap.txt.gz");
		download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.sql", destfile="gap.sql");
		}

	# repeatMasker data
	if ("repeatmasker" %in% datasets && !file.exists(paste0(data.dir, "/rmsk.txt.gz"))) {
		download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz", destfile="rmsk.txt.gz");
		download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.sql", destfile="rmsk.sql");
		}

	# hugo gene names with alias's
	if ("hugo" %in% datasets && !file.exists(paste0(data.dir, "/hugo_protein_coding_genes.txt.gz"))) {
		download.file("ftp://ftp.ebi.ac.uk/pub/databases/genenames/locus_groups/protein-coding_gene.txt.gz", destfile = "hugo_protein_coding_genes.txt.gz");
		}
	
	# cosmic
	if ("cosmic" %in% datasets && !file.exists(paste0(data.dir, "/CosmicCodingMuts_v68.vcf.gz"))) {
		download.file("ftp://ngs.sanger.ac.uk/production/cosmic/CosmicCodingMuts_v68.vcf.gz", destfile ="CosmicCodingMuts_v68.vcf.gz");
		download.file("ftp://ngs.sanger.ac.uk/production/cosmic/CosmicNonCodingVariants_v68.vcf.gz", destfile = "CosmicNonCodingVariants_v68.vcf.gz");
		}

	# clinvar
	if ("clinvar" %in% datasets && !file.exists(paste0(data.dir, "/clinvar.vcf.gz"))) {
		download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_00-latest.vcf.gz", destfile = "clinvar.vcf.gz");
		download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_00-latest.vcf.gz.tbi", destfile = "clinvar.vcf.gz.tbi");
		}

	setwd(old.dir)

	return(0);
	}
