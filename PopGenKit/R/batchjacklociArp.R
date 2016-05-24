batchjacklociArp <-
function(rm.loci)
{

filelist=list.files(pattern='[.arp]$')
nbfiles=length(filelist)
for (t in 1 :nbfiles)
{
jacklociArp(filelist[t],rm.loci)
}
#créer le fichier batch pour Arlequin
arlfilelist=gsub(pattern='[.]arp', replacement='_jackloci.arp',filelist)
write.table(arlfilelist, file='batchjackloci.arb', quote=F, row.names=F, col.names=F, sep='\r')
}

