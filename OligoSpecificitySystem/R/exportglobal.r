"exportglobal"<-function(n){
fileName<-tclvalue(tkgetSaveFile())
    filename<-paste(fileName,".txt",sep="")
    if (filename == "")
    return()
    if(n==1) write(txt1,file=filename,sep=";")
    if(n==2) write(txt2,file=filename,sep=";")
    if(n==3) write(txt3,file=filename,sep=";")
    if(n==4) write(newdb,file=filename,sep=";")
}