library(testthat)
library(Lahman)
library(RSQLite)

test.db.1 <- function()
{
    samples <- data.frame(sample_id=as.integer(sapply(1:100, function(x) rep(x,2))), wave=as.integer(sapply(1:100, function(x) 1:2)), did_collect=sample(c("Y", "N"), size=100, replace=T))
    
    dna <- samples[samples$did_collect == "Y",c("sample_id", "wave")]
    dna$lab_id <- paste("dna", paste(dna$sample_id, dna$wave, sep="_"), sep="_")
    dna$lab_id[sample.int(nrow(dna), size=10)] <- NA
    dna$ng.ul <- abs(rnorm(n=nrow(dna), mean=146, sd=98))
    dna$ng.ul[is.na(dna$lab_id)] <- NA
        
    clinical <- data.frame(sample_id=1:100, sex=sample(c("M", "F"), size=100, replace=T), age=sample(1:20, 100, replace=T), status=sample(c(0L,1L), size=100, replace=T), var_wave_1=rnorm(n=100), var_wave_2=rnorm(n=100))
    clinical <- clinical[-sample.int(nrow(clinical), 12),]
    
    return(list(samples=samples, dna=dna, clinical=clinical))
}

test.db.2 <- eval(structure(list(allele = structure(list(allele_id = 1:10, alleles = c(".", 
"G", "A", ".", "A", "T", ".", "C", "A", "."), allele_num = c(-1L, 
0L, 1L, -1L, 0L, 1L, -1L, 0L, 1L, -1L), ref_id = c(1L, 1L, 1L, 
2L, 2L, 2L, 3L, 3L, 3L, 4L)), .Names = c("allele_id", "alleles", 
"allele_num", "ref_id"), row.names = c(NA, -10L), class = "data.frame"), 
    genotype = structure(list(geno_id = c(1L, 2L, 3L, 4L, 5L, 
    6L, 7L, 8L, 9L, 11L), geno_chr = c(1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L), allele_num = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L), strain = c("NODShiLtJ", "NODShiLtJ", "NODShiLtJ", 
    "NODShiLtJ", "NODShiLtJ", "NODShiLtJ", "NODShiLtJ", "NODShiLtJ", 
    "NODShiLtJ", "NODShiLtJ"), ref_id = c(1L, 2L, 3L, 4L, 5L, 
    6L, 7L, 8L, 9L, 11L)), .Names = c("geno_id", "geno_chr", 
    "allele_num", "strain", "ref_id"), row.names = c(NA, -10L
    ), class = "data.frame"), probe_align = structure(list(probe_align_id = 1:10, 
        probe_chr = c("1", "1", "1", "1", "1", "1", "1", "1", 
        "1", "1"), probe_start = c(3214633L, 3214895L, 3215576L, 
        3216284L, 3216483L, 3216537L, 3216591L, 3421710L, 3421740L, 
        3421808L), probe_end = c(3214657L, 3214919L, 3215600L, 
        3216308L, 3216507L, 3216561L, 3216615L, 3421734L, 3421764L, 
        3421832L), probe_ind = 23944:23953), .Names = c("probe_align_id", 
    "probe_chr", "probe_start", "probe_end", "probe_ind"), row.names = c(NA, 
    -10L), class = "data.frame"), probe_info = structure(list(
        probe_ind = 1:10, fasta_name = c("chr1:3044331-3044355;+;GGGCTCAAATTCACCTGTGTTAACT", 
        "chr1:3044334-3044358;+;AGGGGGCTCAAATTCACCTGTGTTA", "chr1:3044335-3044359;+;GAGGGGGCTCAAATTCACCTGTGTT", 
        "chr1:3044336-3044360;+;GGAGGGGGCTCAAATTCACCTGTGT", "chr1:3044337-3044361;+;GGGAGGGGGCTCAAATTCACCTGTG", 
        "chr1:3044504-3044528;+;AAGGCTTGGGCTATGCGCTCTCCAG", "chr1:3044557-3044581;+;AAGGTCCTCGCTCCTCTTTCATATA", 
        "chr1:3044558-3044582;+;GAAGGTCCTCGCTCCTCTTTCATAT", "chr1:3044559-3044583;+;GGAAGGTCCTCGCTCCTCTTTCATA", 
        "chr1:3044560-3044584;+;AGGAAGGTCCTCGCTCCTCTTTCAT"), 
        probe_id = c(271822L, 866088L, 240008L, 34772L, 396933L, 
        613322L, 984143L, 826914L, 1075207L, 160847L), align_status = c("MultiMapped", 
        "MultiMapped", "MultiMapped", "MultiMapped", "MultiMapped", 
        "MultiMapped", "MultiMapped", "MultiMapped", "MultiMapped", 
        "MultiMapped")), .Names = c("probe_ind", "fasta_name", 
    "probe_id", "align_status"), row.names = c(NA, -10L), class = "data.frame"), 
    probe_to_snp = structure(list(probe_snp_id = 1:10, ref_id = c(1L, 
    2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 9L), probe_align_id = c(137L, 
    137L, 252L, 263L, 561L, 601L, 606L, 622L, 635L, 636L)), .Names = c("probe_snp_id", 
    "ref_id", "probe_align_id"), row.names = c(NA, -10L), class = "data.frame"), 
    reference = structure(list(ref_id = c(1L, 2L, 3L, 4L, 5L, 
    6L, 7L, 8L, 9L, 11L), seqnames = c("1", "1", "1", "1", "1", 
    "1", "1", "1", "1", "1"), start = c(4785683L, 4785703L, 8595537L, 
    8682000L, 10719780L, 12990759L, 13114284L, 13122474L, 13125617L, 
    13150850L), end = c(4785683L, 4785703L, 8595537L, 8682000L, 
    10719780L, 12990759L, 13114284L, 13122474L, 13125617L, 13150850L
    ), filter = c("TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", 
    "TRUE", "TRUE", "TRUE", "TRUE"), vcf_annot_id = c(1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L)), .Names = c("ref_id", "seqnames", 
    "start", "end", "filter", "vcf_annot_id"), row.names = c(NA, 
    -10L), class = "data.frame"), vcf_annot = structure(list(
        vcf_annot_id = 1L, vcf_name = "/Users/bottomly/Desktop/resources/vcfs/mgp.v3.snps.rsIDdbSNPv137.vcf.gz", 
        type = "SNV"), .Names = c("vcf_annot_id", "vcf_name", 
    "type"), row.names = c(NA, -1L), class = "data.frame")), .Names = c("allele", 
"genotype", "probe_align", "probe_info", "probe_to_snp", "reference", 
"vcf_annot")))

test.schema.2 <- list(probe_info=list(db.cols=c("probe_ind", "fasta_name", "probe_id", "align_status"),
                                db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "INTEGER", "TEXT"),
                                db.constr="CONSTRAINT probe_idx UNIQUE (fasta_name)",
                                dta.func=function(x) x, should.ignore=FALSE, foreign.keys=NULL),
                probe_align=list(db.cols=c("probe_align_id", "probe_chr", "probe_start", "probe_end", "probe_ind"),
                                db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "INTEGER", "INTEGER", "INTEGER"),
                                db.constr="CONSTRAINT geno_idx UNIQUE (probe_chr, probe_start, probe_end, probe_ind)",
                                dta.func=function(x) x, should.ignore=FALSE, foreign.keys=list(probe_info=list(local.keys="probe_ind", ext.keys="fasta_name"))),
                vcf_annot=list(db.cols=c("vcf_annot_id", "vcf_name", "type"),
                               db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "TEXT"),
                               db.constr="CONSTRAINT probe_idx UNIQUE (vcf_name, type)",
                               dta.func=function(x) x,
                               should.ignore=TRUE, foreign.keys=NULL),
                reference=list(db.cols=c("ref_id", "seqnames", "start", "end", "filter", "vcf_annot_id"),
                               db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "INTEGER", "INTEGER", "TEXT", "INTEGER"),
                               db.constr="CONSTRAINT ref_idx UNIQUE (seqnames, start, end, vcf_annot_id)",
                               dta.func=function(x) x, should.ignore=TRUE, foreign.keys=list(vcf_annot=list(local.keys="vcf_annot_id", ext.keys=c("vcf_name", "type")))),
                allele=list(db.cols=c("allele_id", "alleles", "allele_num", "ref_id"),
                            db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "INTEGER", "INTEGER"),
                            db.constr="CONSTRAINT alelle_idx UNIQUE (alleles, allele_num, ref_id)",
                            dta.func=function(x) x, should.ignore=TRUE, foreign.keys=list(vcf_annot=list(local.keys="vcf_annot_id", ext.keys=c("vcf_name", "type")),
                                                                                            reference=list(local.keys="ref_id", ext.keys=c("seqnames", "start", "end", "vcf_annot_id")))),
                genotype=list(db.cols=c("geno_id", "geno_chr", "allele_num","strain", "ref_id"),
                              db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "INTEGER", "INTEGER", "TEXT", "INTEGER"),
                              db.constr="CONSTRAINT geno_idx UNIQUE (ref_id, strain, geno_chr, allele_num)",
                              dta.func=function(x) x, should.ignore=TRUE, foreign.keys=list(vcf_annot=list(local.keys="vcf_annot_id", ext.keys=c("vcf_name", "type")),
                                                                                                reference=list(local.keys="ref_id", ext.keys=c("seqnames", "start", "end", "vcf_annot_id")))),
                probe_to_snp=list(db.cols=c("probe_snp_id", "ref_id", "probe_align_id"),
                                  db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "INTEGER", "INTEGER"),
                                  db.constr="CONSTRAINT p_s_idx UNIQUE (ref_id, probe_align_id)",
                                  dta.func=function(x) x, should.ignore=TRUE, foreign.keys=list(vcf_annot=list(local.keys="vcf_annot_id", ext.keys=c("vcf_name", "type")),
                                                                                                    reference=list(local.keys="ref_id", ext.keys=c("seqnames", "start", "end", "vcf_annot_id")),
                                                                                                    probe_align=list(local.keys="probe_align_id", ext.keys=c("probe_chr", "probe_start", "probe_end")))))





test_check("poplite")
