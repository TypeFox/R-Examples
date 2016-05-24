#' @export
dupmodel<-function(path, 
                   path.new,
                   modsp.file="model.spec.csv",
                   model.file, 
                   cols.file, 
                   data,
                   bat.file,
                   model.file.new="test.mdl", 
                   cols.file.new="cols1.txt", 
                   data.new="data1.txt"){
  rootwd=getwd()
  modspfl=paste(path, modsp.file, sep="/")
  if(!file.exists(path.new)) dir.create(path.new)
  
  if(file.exists(modspfl)){
    mod.spec=read.csv(modspfl)
    model.file=mod.spec$model.file[1]
    cols.file=mod.spec$cols.file[1]
    data=mod.spec$data[1]
    bat.file=mod.spec$bat.file}
  
    ctlfl<-paste(path,model.file,sep="/")  
    colsfl<-paste(path,cols.file,sep="/")
    datafl<-paste(path,data,sep="/")

    file.copy(ctlfl,paste(path.new,model.file.new,sep="/"))  
    file.copy(colsfl,paste(path.new,cols.file.new,sep="/"))
    file.copy(datafl,paste(path.new,data.new,sep="/")) 
    file.copy(paste(path,bat.file,sep="/"),paste(path.new,bat.file,sep="/")) 
  
}