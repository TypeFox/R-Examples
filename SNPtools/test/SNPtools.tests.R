################################################################################
# SNPtools tests.
# Daniel Gatti
# Dan.Gatti@jax.org
# Mar. 25, 2013
################################################################################
library(SNPtools)


snp.file = "/Users/dgatti/Documents/SNP/Sanger/Build38/SNPIndels/cc.snps.NCBI38.txt.gz"
indel.file = "/Users/dgatti/Documents/SNP/Sanger/Build38/SNPIndels/cc.indels.NCBI38.txt.gz"
sv.file = "/Users/dgatti/Documents/SNP/Sanger/Build38/SV/cc.svs.NCBI38.txt.gz"
mgi.file = "/Users/dgatti/Documents/MGI/MGI.20130305.sorted.txt.gz"

# Get strains for SNPs.
as = get.strains(file = snp.file)

# Get strains for Indels.
as = get.strains(file = indel.file)

# Get strains for SVs.
as = get.strains(file = sv.file)

# Get SNPs.
snps = get.variants(file = snp.file, chr = 12, start = 104, end = 106)

# Get Indels.
indels = get.variants(file = indel.file, chr = 12, start = 104, end = 106,
         type = "indel")

# Get SVs.
sv = get.variants(file = sv.file, chr = 12, start = 104, end = 106,
     type = "sv")

# Catagorize variants.
varcat = categorize.variants(snps, mgi.file = mgi.file)

# Catagorize variants.
varcat = categorize.variants(indels, mgi.file = mgi.file)

# Catagorize variants.
varcat = categorize.variants(sv, mgi.file = mgi.file)

# Variant plot for SNPs.
variant.plot(var.file = snp.file, mgi.file = mgi.file,
             chr = 12, start = 104, end = 106, type = "snp")

variant.plot(var.file = snp.file, mgi.file = mgi.file,
             chr = 12, start = 104, end = 106, type = "snp", 
             pattern = c("A/J", "C57BL/6J", "129S1/SvImJ", "NZO/HlLtJ",
             "PWK/PhJ"))

variant.plot(var.file = snp.file, mgi.file = mgi.file,
             chr = 12, start = 104, end = 106, type = "snp", qtl = qtl)

variant.plot(var.file = snp.file, mgi.file = mgi.file,
             chr = 12, start = 104, end = 106, type = "snp", pattern =
             c("A/J", "C57BL/6J", "129S1/SvImJ", "NZO/HlLtJ",
             "PWK/PhJ"), qtl = qtl)



# Gene Plot
mgi = get.mgi.features(file = mgi.file, chr = 12, start = 104, end = 106,
      source = "MGI", type = "gene")
gene.plot(mgi = mgi, highlight = c("Serpina3a", "Gsc", "Dicer1"))