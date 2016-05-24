sldrop<-function(statslist, varname, type = 1) 
{
    newstatsl <- statslist
    for (i in 1:length(statslist)) {
       vnms.ind<-dimnames(statslist[[i]][[type]])[[3]]
       vdrop<-which(vnms.ind%in%varname)
       newstatsl[[i]][[type]]<-statslist[[i]][[type]][,,-vdrop]
    }
    newstatsl
}
