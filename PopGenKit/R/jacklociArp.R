jacklociArp <-
function(datafile, rm.loci)
{
input=read.table(file=datafile, sep='\t', colClasses='character', quote="'")
noext= gsub("[.].*$","",datafile)
input[10,1]='MissingData=\'?\''

nloci=ncol(input)-2

new.input=input[,-(rm.loci+2)]

new.nloci=nloci-length(rm.loci)

filename=paste(noext,'_jackloci.arp',sep='')
write.table(new.input, file=filename, quote=F, row.names=F, col.names=F, sep='\t')

}

