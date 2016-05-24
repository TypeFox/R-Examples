mia<-function(hapfile="hap.out",assfile="assign.out",miafile="mia.out",so=0,ns=0,mi=0,allsnps=0,sas=0)
{
  # to call up mi.inference here

  z<-.C("mia",hapfile=as.character(hapfile),assfile=as.character(assfile),
        miafile=as.character(miafile),so=as.integer(so),ns=as.integer(ns),mi=as.integer(mi),
        allsnps=as.integer(allsnps),sas=as.integer(sas),PACKAGE="gap")
}
