# 13-4-2006, MRC Epid

twinan90 <- function(mzdat,dzdat,vname='mzdz',xlamb=1,const=0,vmiss=-9,path=1,
           ped=0,nvar=1,form='((1x,a1,5x,F6.2))')
{
  nMZ <- dim(mzdat)[1]
  nDZ <- dim(dzdat)[1]
  max.twin.pairs=max(nMZ,nDZ)
  nc=nchar(vname)
  nf=nchar(form)
  logfile <- paste(vname,".log",sep="")
  outfile <- paste(vname,".out",sep="")
  pedfile <- paste(vname,".ped",sep="")
  unlink(logfile)
  unlink(outfile)
  unlink(pedfile)
  names <- paste(vname,logfile,outfile,pedfile,form,sep="")
  z<-.Fortran("newtw5",d1=as.double(mzdat),d2=as.double(dzdat),
              nmz=as.integer(nMZ),ndz=as.integer(nDZ),mp=as.integer(max.twin.pairs),
              xlamb=as.double(xlamb),const=as.double(const),vmiss=as.double(vmiss),
              path=as.double(path),ped=as.double(ped),
              names=as.character(names),nchar=as.integer(nc),
              nvar=as.integer(nvar),nf=as.integer(nf),PACKAGE="gap")
  cat("\n")
  cat("Output and log files are ",outfile," and ",logfile,"\n")
  if(ped==1) cat("FISHER data file is ",pedfile,"\n")
  cat("\n")
  x <- replace(mzdat,mzdat==vmiss,NA)
  y <- replace(dzdat,dzdat==vmiss,NA)
  all.MZ <- apply(is.na(mzdat),1,any)
  nMZ <- dim(mzdat[!all.MZ,])[1]
  all.DZ <- apply(is.na(dzdat),1,any)
  nDZ <- dim(dzdat[!all.DZ,])[1]
  rMZ <- cor(x[,1],x[,2],use="complete.obs")
  rDZ <- cor(y[,1],y[,2],use="complete.obs")
  h2 <- 2*(rMZ-rDZ)
  var.rMZ <- (1-rMZ^2)^2/nMZ
  var.rDZ <- (1-rMZ^2)^2/nDZ
  seh2 <- sqrt(var.rMZ+var.rDZ)
  covMZ <- cov(x,use="complete.obs")
  covDZ <- cov(y,use="complete.obs")
  invisible(list(h2=h2,seh2=seh2,nMZ=nMZ,nDZ=nDZ,rMZ=rMZ,rDZ=rDZ,covMZ=covMZ,covDZ=covDZ))
}
