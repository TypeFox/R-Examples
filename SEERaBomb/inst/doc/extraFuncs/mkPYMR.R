mkPYMR <- function(PY,yearEnd) {
  py = PY[,1]
  age=PY[,2] 
  year=PY[,3]
  print(range(py))
  print(range(age))
  print(range(year))
  nr = length(py);
  yrs=1973:yearEnd
  nyears=length(yrs)
  ages=0.5:125.5
  #   ages=0.5:135.5
  nages=length(ages)
  PYM=matrix(0,ncol=nyears,nrow=nages)
  #    PYM=numeric(nyears*nages)
  if (nr>=1) for (i in 1:nr) {
    pystrip=py[i]
    rem = pystrip-floor(pystrip)
    quo = floor(pystrip)
    strtAgeIndx=round(age[i]+.5)
    strtYrIndx=round(year[i]-1972)
    for (j in 0:quo)   
      if ( (strtAgeIndx+j<=nages) & (strtYrIndx+j<=nyears))
        if (j==quo) 
          PYM[strtAgeIndx+j,strtYrIndx+j]=PYM[strtAgeIndx+j,strtYrIndx+j]+rem   else
            PYM[strtAgeIndx+j,strtYrIndx+j]=PYM[strtAgeIndx+j,strtYrIndx+j]+1
  }
#   print(nyears)
#   print(nages)
#   print(nages*nyears)
#   print(head(PYM))
#   print(tail(PYM))
#   print(dim(PYM))
  # length(PYM)<-nages*nyears
  #   dim(PYM)=c(nages,nyears)
  colnames(PYM)=yrs
  rownames(PYM)=ages
  PYM
}


