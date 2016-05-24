histo<-function(dendat,binlkm,epsi=0)
{
# Constructs a histogram estimate: result is given by giving level
# sets of the estimate

supp<-support(dendat,epsi)
regdat<-den2reg(dendat,binlkm,supp)
palvak<-makehis(regdat)
values<-palvak$values
recs<-palvak$recs

integ<-0
recnum<-length(values)
for (i in 1:recnum){
   integ<-integ+values[i]*massone(recs[i,])
}

values<-values/integ

return(list(values=values,recs=recs))
}




