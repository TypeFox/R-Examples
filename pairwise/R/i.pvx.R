#-----------internal function ------------------
pvx<-function(theta,thres,xm=NULL){
  # func. by joerg-henrik heine jhheine(at)googlemail.com
  # nichts geändert aber zum merken theta: einzelne zahl; thres: thurstonian threshold eines items
  # korrigierte formel aus markus buch seite 330 siehe auch s. 522 2006
  s<-0:length(thres)
  thres<-c(0,thres)
  oben<- exp((s*theta)-cumsum(thres))
  unten<- sum(exp((s*theta)-cumsum(thres)))
  px<-oben / unten
  names(px)<-0:(length(thres)-1)
  P<-sum(oben / unten)
  #if(length(xm)==0){return(list(px=px,P=P))}
  if(length(xm)==0){return(px)}
  if(length(xm)!=0){px[xm]}
}     
#---------------------------------------------------  