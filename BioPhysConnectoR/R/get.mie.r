#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  keul(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################



get.mie<-function(aln,method="ORMI",gapchar=NULL,nullmod=NULL,logMI=FALSE){

  logMI<-ifelse(logMI,1,0)
  
  transcode<-function(a,gc){
    aa<-unique(as.vector(a))
    gc<-gc[gc%in%aa]
    if(length(gc)>0){
      for(i in 1:length(gc)){
	ind<-which(aa==gc[i])
	aa<-c(aa[-ind],aa[ind])
      }
    }
    aaa<-a
    for(i in 1:length(aa)){
      aaa[which(a==aa[i])]<-i-1
    }
    return(list(seq=aaa,gc=length(gc)))
  }

  tc<-transcode(a=aln,gc=gapchar)
  seqs<-tc$seq
  gap_chars<-tc$gc
  alphabetlength<-length(unique(as.vector(seqs)))
  n<-dim(aln)[2] # Anzahl der Aminosäre-positionen pro Sequenz
  s<-dim(aln)[1] # Anzahl der Sequenzen
  # Initialisierung der an C zu übergebenen Vektoren und Matrizen
  MI_avg<-matrix(0,n,n) 
  if(method=="ORMI"){
    out<-.C("computeORMI",seqs=as.integer(seqs),n=as.integer(n),s=as.integer(s),MI_avg=as.double(MI_avg),alphabetlength=as.integer(alphabetlength),logMI=as.integer(logMI),PACKAGE="BioPhysConnectoR")
  }
  if(method=="DEMI"){
    out<-.C("computeDEMI",seqs=as.integer(seqs),n=as.integer(n),s=as.integer(s),MI_avg=as.double(MI_avg),gap_chars=as.integer(gap_chars),alphabetlength=as.integer(alphabetlength),logMI=as.integer(logMI),PACKAGE="BioPhysConnectoR")
  }
  if(method=="SUMI"){
    out<-.C("computeSUMI",seqs=as.integer(seqs),n=as.integer(n),s=as.integer(s),MI_avg=as.double(MI_avg),gap_chars=as.integer(gap_chars),alphabetlength=as.integer(alphabetlength),logMI=as.integer(logMI),PACKAGE="BioPhysConnectoR")
  }
  if(method=="ESMI"){
    out<-.C("computeESMI",seqs=as.integer(seqs),n=as.integer(n),s=as.integer(s),MI_avg=as.double(MI_avg),gap_chars=as.integer(gap_chars),alphabetlength=as.integer(alphabetlength),logMI=as.integer(logMI),PACKAGE="BioPhysConnectoR")
  }

  out<-matrix(out$MI_avg,ncol=n)
  out[which(is.nan(out))]<-0
  MI<-out

  if(!is.null(nullmod)){
    if(!is.numeric(nullmod)){
      cat("The argument \"nullmod\" has to be nummeric!\n")
      stop(75)
    }
    MI_avg<-mi<-mi2<-var<-matrix(0,n,n)
    for(i in 1:nullmod){
      seqs<-matrix(seqs,s,n)
      for(p in 1:n){
	seqs[,p]<-sample(seqs[,p])
      }
      seqs<-as.vector(seqs)
      if(method=="ORMI"){
	out<-.C("computeORMI",seqs=as.integer(seqs),n=as.integer(n),s=as.integer(s),MI_avg=as.double(MI_avg),alphabetlength=as.integer(alphabetlength),logMI=as.integer(logMI),PACKAGE="BioPhysConnectoR")
	mi<-mi+matrix(out$MI_avg,ncol=n)
	mi2<-mi2+(matrix(out$MI_avg,ncol=n)^2) 
      }
      if(method=="DEMI"){
	out<-.C("computeDEMI",seqs=as.integer(seqs),n=as.integer(n),s=as.integer(s),MI_avg=as.double(MI_avg),gap_chars=as.integer(gap_chars),alphabetlength=as.integer(alphabetlength),logMI=as.integer(logMI),PACKAGE="BioPhysConnectoR")
	mi<-mi+matrix(out$MI_avg,ncol=n)
	mi2<-mi2+(matrix(out$MI_avg,ncol=n)^2)
      }
      if(method=="SUMI"){
	out<-.C("computeSUMI",seqs=as.integer(seqs),n=as.integer(n),s=as.integer(s),MI_avg=as.double(MI_avg),gap_chars=as.integer(gap_chars),alphabetlength=as.integer(alphabetlength),logMI=as.integer(logMI),PACKAGE="BioPhysConnectoR")
	mi<-mi+matrix(out$MI_avg,ncol=n)
	mi2<-mi2+(matrix(out$MI_avg,ncol=n)^2)
      }
      if(method=="ESMI"){
	out<-.C("computeESMI",seqs=as.integer(seqs),n=as.integer(n),s=as.integer(s),MI_avg=as.double(MI_avg),gap_chars=as.integer(gap_chars),alphabetlength=as.integer(alphabetlength),logMI=as.integer(logMI),PACKAGE="BioPhysConnectoR")
	mi<-mi+matrix(out$MI_avg,ncol=n)
	mi2<-mi2+(matrix(out$MI_avg,ncol=n)^2)
      }
    }
    mi<-mi/nullmod
    mi2<-mi2/nullmod
    var<-(mi2)-(mi^2)
    out<-list(mi=MI,nullmodel=mi,nullsquare=mi2,nullvar=var)
  }
  return(out)
}
