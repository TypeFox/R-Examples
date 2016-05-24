setClass("region.data", representation(

populations      = "list",
populations2     = "list",
popmissing       = "list", 
outgroup         = "list",
outgroup2        = "list",

reading.frame    = "list",
rev.strand       = "list",
Coding.matrix    = "list",
Coding.matrix2   = "list",
UTR.matrix       = "list",
Intron.matrix    = "list",
Exon.matrix      = "list",
Gene.matrix      = "list",

#Coding.info      = "list",
#UTR.info         = "list",
#Intron.info      = "list",
#Exon.info        = "list",
#Gene.info        = "list",


CodingSNPS       = "list",
UTRSNPS          = "list",
IntronSNPS       = "list",
ExonSNPS         = "list",
GeneSNPS         = "list",

transitions      = "list",  # matrix_sv  transition war eine 1
biallelic.matrix = "list",  # matrix_pol
biallelic.sites  = "list",  # matrix_pos
biallelic.sites2 = "list",  # important for SNP DATA
reference        = "list",
matrix_codonpos  = "list", # codonpos. of biallelics
synonymous       = "list", # synnonsyn
matrix_freq      = "list",
n.singletons     = "list", # unic
trans.transv.ratio     = "list",
n.valid.sites        = "list", # algsites
n.biallelic.sites       = "list",
polyallelic.sites = "list", # mhitbp
n.nucleotides    = "list", # sum_sam
biallelic.compositions  = "list", # TCGA
biallelic.substitutions = "list", # subst
minor.alleles    = "list", # mutations
codons           = "list",
sites.with.gaps  = "list", # gaps
sites.with.unknowns = "list"
 
))

#### SHOW ######
setMethod("show", "region.data",
 function(object){
 # print(summary(object))
 cat("-----\n")
 cat("SLOTS:\n")
 cat("-----\n")
 out <- data.frame(Slots=c("populations","populations2","outgroup","transitions","biallelic.matrix","n.singletons",
                   "biallelic.sites","reference","n.nucleotides","biallelic.compositions","synonymous","biallelic.substitutions","polyallelic.sites","sites.with.gaps","sites.with.unknowns",
"minor.alleles","codons","IntronSNPS","UTRSNPS","CodingSNPS","ExonSNPS","GeneSNPS"),         
               

                Description=c("Samples of each population (rows)","Samples of each population (names)","Samples of outgroup",
               "Biallelic site transitions","Biallelic matrix","Number of singletons",
               "Position of biallelic sites","SNP reference","Number of nucleotides per sequence","Nucleotides per sequence (biallelic)",
                "Synonymous biallelic sites","Biallelic substitutions","Sites with >2 nucleotides","Sites with gap positions","Sites with unknown positions","Minor alleles","Codons of biallelic substitutions","SNPs in intron region","SNPs in UTR region","SNPs in coding region","SNPs in exon region","SNPs in gene region") )
 print(out)
 cat("\n---------------\n")
 cat("These are the Slots (class region.data) \n")
  
})

#setMethod("summary","DATA",function(object){
#cat("Some basic informations: \n")
#cat("---------------------------- \n")
#cat("Path of Gene: ",object@genename,"\n")
#cat("unic:       number of single mutations:",object@unic,"\n")
#cat("totalmis:   number of missing positions",object@totalmis,"\n")
#cat("bial_sites: number of biallelic positions",object@bial_sites,"\n")
#})
