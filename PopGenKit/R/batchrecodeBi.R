batchrecodeBi <-
function(major.alleles)
{

filelist=list.files(pattern='[.arp]$')
nbfiles=length(filelist)
for (t in 1 :nbfiles)
{
recodeBi(filelist[t],major.alleles)
}
#créer le fichier batch pour Arlequin
arlfilelist=gsub(pattern='[.]arp', replacement='_Bi.arp',filelist)
write.table(arlfilelist, file='batchBi.arb', quote=F, row.names=F, col.names=F, sep='\r')
}

