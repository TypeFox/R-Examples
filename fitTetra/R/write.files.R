write.files <-
function(mrkresult, modelfile, scorefile, diploscorefile, logfile, allmodelsfile) {
  if (class(mrkresult$scores)!="logical") { # test for NA
    if (file.exists(scorefile))
      write.table(mrkresult$scores,file=scorefile,sep="\t",na="",row.names=F,col.names=F,quote=F,append=T)
    else write.table(mrkresult$scores,file=scorefile,sep="\t",na="",row.names=F,col.names=T,quote=F)
  }
  if (diploscorefile!="" && class(mrkresult$diploscores)!="logical") { # test for NA)
    if (file.exists(diploscorefile))
      write.table(mrkresult$diploscores,file=diploscorefile,sep="\t",na="",row.names=F,col.names=F,quote=F,append=T)
    else write.table(mrkresult$diploscores,file=diploscorefile,sep="\t",na="",row.names=F,col.names=T,quote=F)
  }  
  if (allmodelsfile!="" && class(mrkresult$allmodeldata)!="logical") { # test for NA
    if (file.exists(allmodelsfile))
      write.table(mrkresult$allmodeldata,file=allmodelsfile,sep="\t",na="",row.names=F,col.names=F,quote=F,append=T)
    else write.table(mrkresult$allmodeldata,file=allmodelsfile,sep="\t",na="",row.names=F,col.names=T,quote=F)
  }
  if (file.exists(modelfile))
         write.table(mrkresult$modeldata,file=modelfile,sep="\t",na="",row.names=F,col.names=F,quote=F,append=T)
  else write.table(mrkresult$modeldata,file=modelfile,sep="\t",na="",row.names=F,col.names=T,quote=F)
  if (logfile!="") write(mrkresult$log, ncolumns=1, file=logfile, append=T)
}
