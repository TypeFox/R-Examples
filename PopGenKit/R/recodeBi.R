recodeBi <-
function(datafile, major.alleles)
{
input=read.table(file=datafile, sep='\t', colClasses='character', quote="'")
noext= gsub("[.].*$","",datafile)
input[10,1]='MissingData=\'?\''

inputBi=input

nloci=ncol(input)-2

for (a in 1:nloci)
{
 for (b in 1:nrow(input))
{
if (input[b,(a+2)]!='') {
if (input[b,(a+2)]!='?') {
if (input[b,(a+2)]!=major.alleles[a]) inputBi[b,(a+2)]='300' }}
}}

filename=paste(noext,'_Bi.arp',sep='')
write.table(inputBi, file=filename, quote=F, row.names=F, col.names=F, sep='\t')

}

