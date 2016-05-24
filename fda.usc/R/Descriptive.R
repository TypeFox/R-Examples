#  for functional data
func.mean<-function(fdataobj){
if (!is.fdata(fdataobj)) fdataobj<-fdata(fdataobj)
fdataobj[["data"]]<-matrix(apply(fdataobj[["data"]],2,mean,na.rm=TRUE),nrow=1)
fdataobj$names$main<-"mean"
fdataobj
}


func.var<-function(fdataobj){
if (!is.fdata(fdataobj)) fdataobj<-fdata(fdataobj)
n<-dim(fdataobj)[1]
fdataobj[["data"]]<-(n-1)*apply(fdataobj[["data"]],2,var)/n
fdataobj[["data"]]<-matrix(fdataobj[["data"]],nrow=1)
fdataobj$names$main<-"var"
fdataobj
}

func.trim.FM=function(fdataobj,...){depth.FM(fdataobj,...)$mtrim}
func.trim.mode=function(fdataobj,...){depth.mode(fdataobj,...)$mtrim}
func.trim.RP=function(fdataobj,...){depth.RP(fdataobj,...)$mtrim} 
func.trim.RT=function(fdataobj,...){depth.RT(fdataobj,...)$mtrim} 
func.trim.RPD=function(fdataobj,...){depth.RPD(fdataobj,...)$mtrim}

func.med.FM=function(fdataobj,...){depth.FM(fdataobj,...)$median} 
func.med.mode=function(fdataobj,...){depth.mode(fdataobj,...)$median}
func.med.RP=function(fdataobj,...){ depth.RP(fdataobj,...)$median}
func.med.RT=function(fdataobj,...){ depth.RT(fdataobj,...)$median}
func.med.RPD=function(fdataobj,...){ depth.RPD(fdataobj,...)$median}

func.trimvar.FM=function(fdataobj,...){
  lista=depth.FM(fdataobj,...)$ltrim
  func.var(fdataobj[lista,])
}


func.trimvar.mode=function(fdataobj,...){
  lista=depth.mode(fdataobj,...)$ltrim
  func.var(fdataobj[lista,])
  }

func.trimvar.RP=function(fdataobj,...){
 lista=depth.RP(fdataobj,...)$ltrim
 func.var(fdataobj[lista,])}

func.trimvar.RPD=function(fdataobj,...){
 lista=depth.RPD(fdataobj,...)$ltrim
 func.var(fdataobj[lista,])}


func.trim.RT=function(fdataobj,...){depth.RT(fdataobj,...)$mtrim}
func.med.RT=function(fdataobj,...){ depth.RT(fdataobj,...)$median}
func.trimvar.RT=function(fdataobj,...){
 lista=depth.RT(fdataobj,...)$ltrim
 func.var(fdataobj[lista,])}

#  for multivarate data
# func.trim.SD=function(fdataobj,...){depth.SD(fdataobj,...)$mtrim}
# func.trim.PD=function(fdataobj,...){depth.PD(fdataobj,...)$mtrim}
# func.trim.HD=function(fdataobj,...){depth.HD(fdataobj,...)$mtrim} 
# func.trim.MhD=function(fdataobj,...){depth.MhD(fdataobj,...)$mtrim}

# func.med.SD=function(fdataobj,...){depth.SD(fdataobj,...)$median} 
# func.med.PD=function(fdataobj,...){depth.PD(fdataobj,...)$median}
# func.med.HD=function(fdataobj,...){ depth.HD(fdataobj,...)$median}
# func.med.MhD=function(fdataobj,...){ depth.MhD(fdataobj,...)$median}


