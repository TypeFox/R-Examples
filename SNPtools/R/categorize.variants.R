categorize.variants <-
function (variants, mgi.file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz") 
{
  if(missing(variants) || nrow(variants) == 0) {
    return(NULL)
  } # if(missing(variants) || nrow(variants) == 0)

  type = attr(variants, "type")

  hdr = variants[,1:5]
  chr = hdr$CHR[1]
  start = 0
  end = 0

  if(type %in% c("snp", "indel")) {
    start = min(hdr$POS)
    end   = max(hdr$POS)
  } else if(type == "sv") {
    start = min(hdr$START)
    end   = max(hdr$END)
  } # else if(type == "sv")

  mgi = get.mgi.features(mgi.file, chr = chr, start = start, 
        end = end, source = "all", type = c("transcript", "mRNA", 
        "gene", "exon", "three_prime_UTR", "five_prime_UTR"))
  output = cbind(hdr, symbol = rep(NA, nrow(hdr)), id = rep(NA, 
           nrow(hdr)), type = rep(NA, nrow(hdr)))
  output$type = "intergenic"

  # First classify SNPs in genes.
  gene.rows = which(mgi$type == "gene")
#  This consumed too much memory.
#  inter = outer(hdr$POS, mgi$start[gene.rows], ">=") & outer(hdr$POS, 
#          mgi$stop[gene.rows], "<=")
#  colnames(inter) = gene.rows
#  inter = apply(inter, 2, which)
#  inter = inter[sapply(inter, length) > 0]

  if(type %in% c("snp", "indel")) {
    for(i in gene.rows) {
      inter = which(hdr$POS >= mgi$start[i] & hdr$POS <= mgi$stop[i])
      if(length(inter) > 0) {
        name = get.gene.name(mgi$ID[i], mgi)
        output$symbol[inter] = name
        output$id[inter] = mgi$ID[gene.rows[i]]
        output$type[inter] = "intron"
      } # if(length(inter) > 0)
    } # for(i)
  } else if(type == "sv") {
    for(i in gene.rows) {
      inter = which(hdr$START >= mgi$start[i] & hdr$END <= mgi$stop[i])
      if(length(inter) > 0) {
        name = get.gene.name(mgi$ID[i], mgi)
        output$symbol[inter] = name
        output$id[inter] = mgi$ID[gene.rows[i]]
        output$type[inter] = "intron"
      } # if(length(inter) > 0)
    } # for(i)
  } # else if(type == "sv")

  # Exons.
  exon.rows = which(mgi$type == "exon")
  if(type %in% c("snp", "indel")) {
    for(i in exon.rows) {
      inter = which(hdr$START >= mgi$start[i] & hdr$END <= mgi$stop[i])
      if(length(inter) > 0) {
       output$type[inter] = "exon"
      } # if(length(inter) > 0)
    } # for(i)
  } else if(type == "sv") {
    for(i in exon.rows) {
      inter = which(hdr$START >= mgi$start[i] & hdr$END <= mgi$stop[i])
      if(length(inter) > 0) {
       output$type[inter] = "exon"
      } # if(length(inter) > 0)
    } # for(i)
  } # else if(type == "sv")

  # 3' UTR
  utr.rows = which(mgi$type == "three_prime_UTR")
  if(type %in% c("snp", "indel")) {
    for(i in utr.rows) {
      inter = which(hdr$POS >= mgi$start[i] & hdr$POS <= mgi$stop[i])
      if(length(inter) > 0) {
        output$type[inter] = "3UTR"
      } # if(length(inter) > 0)
    } # for(i)
  } else if(type == "sv") {
    for(i in utr.rows) {
      inter = which(hdr$START >= mgi$start[i] & hdr$END <= mgi$stop[i])
      if(length(inter) > 0) {
        output$type[inter] = "3UTR"
      } # if(length(inter) > 0)
    } # for(i)
  } # else if(type == "sv")

  # 5' UTR
  utr.rows = which(mgi$type == "five_prime_UTR")
  if(type %in% c("snp", "indel")) {
    for(i in utr.rows) {
      inter = which(hdr$POS >= mgi$start[i] & hdr$POS <= mgi$stop[i])
      if(length(inter) > 0) {
        output$type[inter] = "5UTR"
      } # if(length(inter) > 0)
    } # for(i)
  } else if(type == "sv") {
    for(i in utr.rows) {
      inter = which(hdr$START >= mgi$start[i] & hdr$END <= mgi$stop[i])
      if(length(inter) > 0) {
        output$type[inter] = "5UTR"
      } # if(length(inter) > 0)
    } # for(i)
  } # else if(type == "sv")

  return(output)
} # categorize.variants()


get.gene.name = function(value, mgi) {
  first3 = substr(value, 1, 3)
  name = NA
  if(first3 == "MGI") {
    name = mgi$Name[mgi$ID == value]
  } else {
    m = intersect(which(mgi$source == "MGI"), grep(value, mgi$Dbxref))
    if(length(m) > 0) {
      name = mgi$Name[m]
    } else {
      name = value
    } # else
  } # else

  return(name)
} # get.gene.name()
