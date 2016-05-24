


GSA.read.gmt=function(filename){
#
## Read in and parse a gmt file (gene set file) from the  Broad institute
# this is tricky, because each lines (geneset) has a variable length
#  I read the file twice, first to pick up the geneset name and description
# in the   first two  columns, then I read it all in as a long string

# The beginning and end of each gene set in the string
# is determined by matching
# BOTH on  geneset name and description (since geneset names sometimes
# occur as genenames elsewhere in the file)

a=scan(filename,what=list("",""),sep="\t", quote=NULL, fill=T, flush=T,multi.line=F)
geneset.names=a[1][[1]]

geneset.descriptions=a[2][[1]]

dd=scan(filename,what="",sep="\t", quote=NULL)


nn=length(geneset.names)
n=length(dd)
ox=rep(NA,nn)

ii=1
for(i in 1:nn){
cat(i)
 while((dd[ii]!=geneset.names[i]) | (dd[ii+1]!=geneset.descriptions[i]) ){
  ii=ii+1
  }
 ox[i]=ii
 ii=ii+1
}

genesets=vector("list",nn)

for(i in 1:(nn-1)){
cat(i,fill=T)
i1=ox[i]+2
i2=ox[i+1]-1
geneset.descriptions[i]=dd[ox[i]+1]
 genesets[[i]]=dd[i1:i2]
}

geneset.descriptions[nn]=dd[ox[nn]+1]
 genesets[[nn]]=dd[(ox[nn]+2):n]
out=list(genesets=genesets,geneset.names=geneset.names, geneset.descriptions=geneset.descriptions)
class(out)="GSA.genesets"
return(out)
}
