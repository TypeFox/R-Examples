jacklociArpRANDOM <-
function(datafile, n.retain)
{
input=read.table(file=datafile, sep='\t', colClasses='character', quote="'")
noext= gsub("[.].*$","",datafile)
input[10,1]='MissingData=\'?\''

nloci=ncol(input)-2

subset=sample(1:nloci, n.retain, replace=F)

subset2=c(1,2,(subset+2))
new.input=input[,subset2]

filename=paste(noext,'_',n.retain,'loci.arp',sep='')
write.table(new.input, file=filename, quote=F, row.names=F, col.names=F, sep='\t')

}

