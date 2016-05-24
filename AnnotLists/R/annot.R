#annot_file=""; data_dir=""; res_path=""; display=TRUE
annot <-
function(annot_file="", data_dir="", res_path="", display=TRUE)
{  
  annotlist<-function(annot, list_id, res_path, display)
  {
    liste = read.table(file=list_id, header=TRUE, sep="\t")
    res = annot[duplicated(c(liste[,1], annot[,1]))[(nrow(liste)+1):(nrow(annot)+nrow(liste))],]
    if(ncol(liste)>2) res = cbind(res[order(res[,1]),], liste[order(liste[,1]),])

    if(nrow(res)!=nrow(liste))
    {
      reste = liste[!duplicated(c(res[,1], liste[,1]))[(nrow(res)+1):(nrow(res)+nrow(liste))],]
      write.table(res, row.names=FALSE, file = paste(res_path, "/Rest_", substr(basename(list_id), 0, (nchar(basename(list_id))-4)), ".txt", sep=""), sep="\t")
    }
    
    if(nrow(res)!=0)
    {
      write.table(res, row.names=FALSE, file = paste(res_path, "/Annot_", substr(basename(list_id), 0, (nchar(basename(list_id))-4)), ".txt", sep=""), sep="\t")
    }
    if(display)  write(paste("\t\t", substr(basename(list_id), 0, (nchar(basename(list_id))-4)), ": Annotations=", (round(nrow(res)/nrow(liste), digits=2)*100), "%", sep=""), file="")
  }
  
  if(display)  
  {
    write("\t########################################################", file="")
    write("\t#                                                      #", file="")
    write("\t#                     Annot (v1.2)                     #", file="")
    write("\t#                                                      #", file="")
    write("\t########################################################", file="")
    write("\t[Run man.annot() for quick help]", file="")
    write("\tPrototype: annot(annot_file=\"\", data_dir=\"\")\n", file="")
  }
  if(res_path=="")
  {
    res_path = paste(getwd(), "/AnnotLists_", format(Sys.time(), "(%H-%M-%S)_%a_%d_%b_%Y"), sep="")
    dir.create(res_path)
  }
  
  if(display)  write("\tAnnotations file", file="")
  flush.console()
  if(annot_file!="")
  {
    annot = read.table(file=annot_file, header=TRUE, sep="\t")
  }else{
    if(Sys.info()["sysname"]=="Windows")
    {
      annot = read.table(file=file.choose(), header=TRUE, sep="\t")
    }else{
      require(tcltk)
      annot = read.table(file=tk_choose.files(), header=TRUE, sep="\t")
    }
  }

  if(display)  write("\tLists folder\n", file="")
  flush.console()
  if(data_dir!="")
  {
    listes = as.matrix(list.files(path = data_dir, full.names = TRUE)) 
  }else{
    if(Sys.info()["sysname"]=="Windows")
    {
      listes = as.matrix(list.files(path = choose.dir(), full.names = TRUE))
    }else{
      listes = as.matrix(list.files(path = tk_choose.dir(), full.names = TRUE))
    }
  }
  apply(listes, 1, function(x) annotlist(annot, x, res_path, display))

  if(display)  write(paste("\n\tResulting files have been placed @ ", res_path, sep=""), file="")
}

