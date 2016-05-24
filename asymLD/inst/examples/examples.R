library(asymLD)

# An example using haplotype frequencies from Wilson(2010)
data(hla.freqs)
hla.a_b <- hla.freqs[hla.freqs$locus1=="A" & hla.freqs$locus2=="B",]
compute.ALD(hla.a_b)
hla.freqs$locus <- paste(hla.freqs$locus1, hla.freqs$locus2, sep="-")
compute.ALD(hla.freqs[hla.freqs$locus=="C-B",])

# An example using genotype data from the haplo.stats package
require(haplo.stats)
data(hla.demo)
geno <- hla.demo[,5:8]  #DPB-DPA 
label <- unique(gsub(".a(1|2)", "", colnames(geno)))
label <- paste("HLA*",label,sep="")
keep <- !apply(is.na(geno) | geno==0, 1, any)
em.keep  <- haplo.em(geno=geno[keep,], locus.label=label)
hapfreqs.df <- cbind(em.keep$haplotype, em.keep$hap.prob) 
#format dataframe for ALD function
names(hapfreqs.df)[dim(hapfreqs.df)[2]] <- "haplo.freq"
names(hapfreqs.df)[1] <- "allele1"
names(hapfreqs.df)[2] <- "allele2"
hapfreqs.df$locus1 <- label[1]
hapfreqs.df$locus2 <- label[2]
head(hapfreqs.df)
compute.ALD(hapfreqs.df)
# Note that there is substantially less variablity (higher ALD) for HLA*DPA1 
# conditional on HLA*DPB1 than for HLA*DPB1 conditional on HLA*DPA1, indicating 
# that the overall variation for DPA1 is relatively low given specific DPB1 alleles

# An example using SNP data where results are symmetric and equal to the ordinary correlation measure (r)
data(snp.freqs)
snps <- c("rs1548306", "rs6923504", "rs4434496", "rs7766854")
compute.ALD(snp.freqs[snp.freqs$locus1==snps[2] & snp.freqs$locus2==snps[3],])

snp.freqs$locus <- paste(snp.freqs$locus1, snp.freqs$locus2, sep="-")
by(snp.freqs,list(snp.freqs$locus),compute.ALD)

# SNP1 & SNP2 : the r correlation & ALD measures are equivalent due to symmetry for bi-allelic SNPs
p.AB <- snp.freqs$haplo.freq[1]
p.Ab <- snp.freqs$haplo.freq[2]
p.aB <- snp.freqs$haplo.freq[3]
p.ab <- snp.freqs$haplo.freq[4]
p.A <- p.AB + p.Ab
p.B <- p.AB + p.aB 
r.squared <- (p.AB - p.A*p.B)^2 / (p.A*(1-p.A)*p.B*(1-p.B))
sqrt(r.squared) #the r correlation measure
compute.ALD(snp.freqs[snp.freqs$locus1==snps[1] & snp.freqs$locus2==snps[2],])

