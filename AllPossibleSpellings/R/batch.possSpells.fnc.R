batch.possSpells.fnc <-
function(fn=list.files(pattern=".*\\.rda")){
  for(file in fn){
    sink(file=paste(gsub("(.*)\\.rda","\\1",file),"_log.txt",sep=""),split=TRUE)
    cat("loading file",file,"\n")
    load(file)
    possible.spellings=possSpells.fnc(words=words)
    write(possible.spellings,file=paste(gsub("(.*)\\.rda","\\1",file),".txt",sep=""))
    sink(file=NULL)
  }
}

