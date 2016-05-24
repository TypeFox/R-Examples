#' @export
"phxvpc.sim"<-function(path, 
                         vpcpath=NULL, 
                         ivar="t",
                         nsim=200,
                         pstrat=NULL,
                         setseed=NULL,
                         pred.corr=NULL,
                         var.corr=FALSE,
                         pi=c(0.025, 0.5, 0.975),
                         pi.ci=c(0.025, 0.975),
                         bin.option=NULL,
                         bin.bound=NULL,
                         bin.center=NULL,
                         modsp.file="model.spec.csv",
                         out.file="out0001.txt",
                         clean=FALSE,
                         hold=FALSE){
  rootwd=getwd()
  modspfl=paste(path, modsp.file, sep="/")
  outfl=paste(path, out.file, sep="/")
  if(all(file.exists(c(modspfl,outfl)))){
    ######checking model.file, cols.file, and data
    mod.spec=read.csv(modspfl)
    model.file=mod.spec$model.file[1]
    cols.file=mod.spec$cols.file[1]
    data=mod.spec$data[1]
    est.method=mod.spec$method
    est.iterlimit=mod.spec$iterlimit
    bat.file=mod.spec$bat.file
    ctlfl<-paste(path,model.file,sep="/")  
    colsfl<-paste(path,cols.file,sep="/")
    datafl<-paste(path,data,sep="/")
    
    if(all(file.exists(c(ctlfl,colsfl,datafl)))){
      
      if(is.null(vpcpath)) vpcpath=paste(path,"vpc_1",sep="/")
      
      if(file.exists(vpcpath)){
        write(paste("Warning: VPC path",vpcpath,"already exists."),"")

        while(file.exists(vpcpath)){
          dir.name=dirname(vpcpath)
          folder.name=basename(vpcpath)
          last.elm=strsplit(folder.name, split="_")[[1]][2]
          other.elm=strsplit(folder.name, split="_")[[1]][1]
          if(last.elm %in% as.character(c(1:999))) { 
            last.elm=as.numeric(last.elm)+1
            folder.name=paste(other.elm, last.elm, sep="_")
          }else folder.name=paste(folder.name,"new", sep=".")
          vpcpath=paste(dir.name, folder.name, sep="/")
        }
        
        write(paste("New VPC path", vpcpath,"automatically assigned"),"")
      }

      dir.create(vpcpath)
      
      model.text<-(appmodest(path,model.file,out.file=out.file))$ctltxt
      newctlfl<-paste(vpcpath,"test.mdl",sep="/")
      writeLines(model.text,con=newctlfl,sep="\n",useBytes=FALSE)  
      file.copy(colsfl,paste(vpcpath,"cols1.txt",sep="/"))
      file.copy(datafl,paste(vpcpath,"data1.txt",sep="/")) 
      file.copy(paste(path,bat.file,sep="/"),paste(vpcpath,bat.file,sep="/")) 
      
      phxnlme(path=vpcpath,model.file="test.mdl", cols.file="cols1.txt", data="data1.txt",
                   method=est.method, iterlimit=est.iterlimit)
      setwd(rootwd)
      simmodel(vpcpath=vpcpath, nsim=nsim, pstrat=pstrat, setseed=setseed, pred.corr=pred.corr,
                    var.corr=var.corr, pi=pi, pi.ci=pi.ci, bin.option=bin.option, 
                    bin.bound=bin.bound, bin.center=bin.center, ivar=ivar, clean=clean, hold=hold
      )
      
    }else{
      write("Error: Required file(s) missing. Check: path, modsp.file, model.file, cols.file, data","")
      setwd(rootwd)
      return(NULL)}
    
  }else{ write("Error: Required file(s) model.spec.csv or out0001.txt missing. Execute phxnlme.R first","")
         setwd(rootwd)
         return(NULL)}
}