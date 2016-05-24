batchconvert <-
function(ndigit=3)
{

filelist=list.files(pattern='[.gen]$')
nbfiles=length(filelist)
for (t in 1 :nbfiles)
{
convert(filelist[t],ndigit)
}
#créer le fichier batch pour Arlequin
arlfilelist=gsub(pattern='[.]gen', replacement='.arp',filelist)
write.table(arlfilelist, file='batch.arb', quote=F, row.names=F, col.names=F, sep='\r')
}

