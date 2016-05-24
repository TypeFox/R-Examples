batchconvertjack <-
function(ndigit=3, rm.loci)
{

filelist=list.files(pattern='[.gen]$')
nbfiles=length(filelist)
for (t in 1 :nbfiles)
{
convertjack(filelist[t],ndigit,rm.loci)
}
#créer le fichier batch pour Arlequin
arlfilelist=gsub(pattern='[.]gen', replacement='_jack.arp',filelist)
write.table(arlfilelist, file='batchjack.arb', quote=F, row.names=F, col.names=F, sep='\r')
}

