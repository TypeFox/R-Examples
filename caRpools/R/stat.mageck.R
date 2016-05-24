stat.mageck=function(untreated.list, treated.list, namecolumn=1, fullmatchcolumn=2, norm.fun="median", extractpattern=expression("^(.+?)_.+"), mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr", filename=NULL, fdr.pval=0.05){
  
  # mageckfolder: where to store analysis files
  # by default: getwd

  # Dataset will be passed on to Mageck Python Script

  # check for identity
  # take NAME element of each dataset and check whether it is identical in the others
  
  non=lapply(
    c(untreated.list,treated.list),
    function(x) stopifnot(
      identical(
        x[,namecolumn],
        treated.list[[1]][,namecolumn]
      )
    )
  )
  # extract gene names
  gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
  design.names=treated.list[[1]][,namecolumn]
  
  
  
  # provide the necessary dataset combination used for mageck
  
  # sgRNA name  sample.treated.reads (replicate)  sample.untreated.reads (replicate)
  
  untreated.list<-do.call(
      "cbind",
      lapply(
        untreated.list,
        FUN= function (x) return(x[,fullmatchcolumn])
      )
    )
  treated.list<-do.call(
    "cbind",
    lapply(
      treated.list,
      FUN= function (x) return(x[,fullmatchcolumn])
    )
  )
  
  # get the number of replicates
  ncol.untreated = ncol(untreated.list)
  ncol.treated = ncol(treated.list)
  
  #print(design.names)
  designs = as.character(design.names)
  # create data.frame for use
  dataset.combined = data.frame( 
    designs = designs,
    genes = gene.names,
    stringsAsFactors=FALSE)
  
  for(i in 1:ncol.treated)
  {
    dataset.combined[,2+i] = as.numeric(treated.list[,i])
  }
  for(i in 1:ncol.untreated)
  {
    dataset.combined[,2+ncol.treated+i] = as.numeric(untreated.list[,i]) 
  }
  
  
  # Write file to pass on to Mageck
  if(is.null(filename)) { filename="mageckanalysisfile"}
  if(is.null(mageckfolder)) {dirstore = getwd()} else {dirstore = mageckfolder}
  
  dataset.combined.file=paste(dirstore,"/",filename, "_MAGeCK_sgRNA.tab", sep="")
  write.table(dataset.combined, file=dataset.combined.file, row.names=FALSE,quote=FALSE)
  
  #print(dataset.combined$designs)
  
  # pass on the dataset to Mageck via rPython
#   
#   usage: mageck test [-h] -k COUNT_TABLE -t TREATMENT_ID [-c CONTROL_ID]
#   [-n OUTPUT_PREFIX] [--norm-method {none,median,total}]
#   [--normcounts-to-file]
#   [--gene-test-fdr-threshold GENE_TEST_FDR_THRESHOLD]
#   [--adjust-method {fdr,holm}] [--variance-from-all-samples]
#   [--sort-criteria {neg,pos}] [--keep-tmp]
  
  
# Treated samples start in 3rd row
treated.samples = 0
for(i in seq(from = 1, to = (0 + ncol.treated-1), by = 1))
{
  treated.samples = paste(treated.samples, i, sep=",")
}

# Untreated samples start at 3rd + ncol.treated +1
untreated.samples = 0+ncol.treated
for(i in seq(from = (0+ncol.treated+1), to = (0+ncol.treated+ncol.untreated-1), by = 1))
{
  untreated.samples = paste(untreated.samples, i, sep=",")
}


# # Create String
#  if(ncol.treated== 2)
#  {
#    treated.samples = paste(0, 1, sep=",")
#    if(ncol.untreated== 2)
#    {
#      
#      untreated.samples = paste(ncol.treated, ncol.treated+1, sep=",")
#    }
#    else {
#      untreated.samples = ncol.treated
#    }
#  }
#  else {
#    treated.samples = 0 
#    
#    if(ncol.untreated== 2)
#    {
#      
#      untreated.samples = paste(1, 2, sep=",")
#    }
#    else {
#      untreated.samples = 1
#    }
#    
#  }



# str(dataset.combined.file)
# str(treated.samples)
# str(untreated.samples)
# str(norm.fun)
# str(sort.criteria)
# str(adjust.method)
# str(filename)

mageckstring = paste("mageck test", "--count-table", dataset.combined.file, "--treatment-id", treated.samples, "--control-id", untreated.samples, "--norm-method", norm.fun, "--sort-criteria", sort.criteria, "--gene-test-fdr-threshold ", fdr.pval, "--adjust-method", adjust.method, "--output-prefix", filename, sep=" ")
# Pass on to Mageck
system(mageckstring)

# load files created by Mageck
# Filenames:
# sample1.gene_summary.txt
# sample1.sgrna_summary.txt


data.mageck.genes = load.file(paste(filename, "gene_summary.txt", sep="." ))
data.mageck.sgrna = load.file(paste(filename, "sgrna_summary.txt", sep="." ))


# Layout gene file:
#id  num	neg|score	neg|p-value	neg|fdr	neg|rank	neg|goodsgrna	pos|score	pos|p-value	pos|fdr	pos|rank	pos|goodsgrna
#ENSG00000171552	30	1.0203e-05	0.00020577	0.084158	1	18	0.33353	0.70983	0.998632	254	2

# Layout sgRNA file:
#sgrna   Gene   control_count   treatment_count control_mean    treat_mean      control_var     adj_var score   p.low   p.high  p.twosided      FDR     high_in_treatment
#INO80B_m74682554   INO80B        0.0/0.0 1220.1598778/1476.14096301      0.810860655738  1348.15042041   0.0     19.0767988005   308.478081895   1.0     1.11022302463e-16       2.22044604925e-16       1.57651669497e-14       True

dataset.return = data.frame( 
  row.names = data.mageck.genes$id,
  genes = data.mageck.genes$id,
  pos = data.mageck.genes[,paste("pos.", adjust.method, sep="")],
  rank.pos = data.mageck.genes[,"pos.rank"],
  neg = data.mageck.genes[,paste("neg.", adjust.method, sep="")],
  rank.neg = as.numeric(data.mageck.genes[,"neg.rank"]),
  #sgrna.neg = as.numeric(sgRNA.sig.neg$x),
 # sgrna.pos = as.numeric(sgRNA.sig.pos$x),
  sgrna.neg.good = as.numeric(data.mageck.genes[,"neg.goodsgrna"]),
  sgrna.pos.good = as.numeric(data.mageck.genes[,"pos.goodsgrna"]),
  stringsAsFactors=FALSE)

# now add the correct number of significant sgRNAs
#dataset.return$sgrna.neg = apply(dataset.return, 1, function(u) return(sgRNA.sig.neg[sgRNA.sig.neg$Group.1==u["genes"],"x"]))
#dataset.return$sgrna.pos = apply(dataset.return, 1, function(u) return(sgRNA.sig.pos[sgRNA.sig.pos$Group.1==u["genes"],"x"]))

  
return(list(genes = dataset.return, sgRNA = data.mageck.sgrna))
  
  
}