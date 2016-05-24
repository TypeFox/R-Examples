utils::globalVariables(c("dmp.txt"))

phxdata <- function(){
  
  # Read model specification for fit
  if(all(file.exists(c("model.spec.csv","dmp.txt","out0001.txt","predtable.dat")))){
    mod.spec = read.csv("model.spec.csv")
    data=as.character(mod.spec$data)
    
    # read predtable.dat and out.txt
    predtable = read.table("predtable.dat",header=TRUE)
    
    source("dmp.txt")

    # count number of records
    rcount = nrow(predtable)
    
    out.txt = read.table("out0001.txt",skip=(length(readLines("out0001.txt"))-rcount-1),header=TRUE)
    
    # read etacov.txt
    ideta = read.csv("IdEta.txt",header = FALSE)
    ideta = ideta[,5:ncol(ideta)]
    names(ideta) = c("ID","param","value")
    ideta = reshape(ideta, v.names = "value", idvar = "ID",
                    timevar = "param", direction = "wide")
    names(ideta) = sub("value. ","",names(ideta))
    
    ## Read in data and add posthoc etas  
    if(length(grep(".csv",data))!=0){
      phxnlme.data = read.csv(data,header=TRUE)
    } else{
      phxnlme.data = read.table(data,header=FALSE,sep=",")
      col.names = read.delim(data,header = FALSE,nrow=1)
      col.names = col.names[!is.na(col.names)]
      
      if(length(col.names)==length(phxnlme.data)) {
        names(phxnlme.data)= col.names
      } else {
        #names(phxnlme.data)= col.names[-1]
        col.names = read.delim(data,header = FALSE,nrow=1,sep=",")
        names(phxnlme.data)= col.names
      }      
    }
        
    names(phxnlme.data)[1] = "ID" #rename to ID for merging with output
    
    # read parmtable and add to data together with typical values
    param = read.csv("parmtable.csv",header = TRUE)
    names(param)[[2]] = "ID"
    phxnlme.data = merge(phxnlme.data,ideta)
    phxnlme.data = merge(phxnlme.data,param)
    phxnlme.data = data.frame(phxnlme.data,t(as.matrix(dmp.txt$coefficients$fixed)))
    
    # create eta table
    etatable = as.data.frame(dmp.txt$coefficients$random$Subject)
    
    # One per ID for parameter vs non-varying covariate
    phxnlme.data.first = phxnlme.data[!duplicated(phxnlme.data$ID),]
    
    # Create phxnlme output object
    phxd = list(etatable=etatable,out.txt=out.txt,phxnlme.data=phxnlme.data,
                phxnlme.data.first=phxnlme.data.first,param=param,dmp.txt=dmp.txt)
      
    return(phxd)
            
  }else{
      write("Error: Required files are missing.")
  }
}

