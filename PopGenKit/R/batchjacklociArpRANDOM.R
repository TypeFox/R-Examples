batchjacklociArpRANDOM <-
function(n.retain)
{

filelist=list.files(pattern='[.arp]$')
nbfiles=length(filelist)
for (t in 1 :nbfiles)
{
jacklociArpRANDOM(filelist[t],n.retain)
}
#créer le fichier batch pour Arlequin
filename=paste('_',n.retain,'loci.arp',sep='')
arlfilelist=gsub(pattern='[.]arp', replacement=filename,filelist)
write.table(arlfilelist, file='batchjackloci.arb', quote=F, row.names=F, col.names=F, sep='\r')
}

