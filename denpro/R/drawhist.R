drawhist<-function(dendat,binlkm,epsi=0,plkm){
#piirtaa 2-ulotteisessa tapauksessa histogramma estimaattorin kuvaajan 
#
#plkm on kuvaajan hilan pisteiden lkm
#
#dendat<-matrix(rnorm(20),10) 
#koe<-drawhist(dendat,binlkm=3,plk=30)
#persp(koe$x,koe$y,koe$z,phi=30,theta=60)
#
hi<-histo(dendat,binlkm,epsi)
recs<-hi$recs
values<-hi$values   
#
ans<-drawgene(values,recs,plkm)
return(list(x=ans$x,y=ans$y,z=ans$z))
}                                      








