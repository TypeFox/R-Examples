findobsR<-function(dendat, ordobs, leftrecend, coordi)
{
  maara<-length(ordobs)
  j<-1
  ind<-ordobs[j]
  havakoor<-dendat[ind,coordi]

  lcount<-0
  while ((havakoor<=leftrecend) && (j<maara)){
          lcount<-lcount+1
          j<-j+1
          ind<-ordobs[j]
          havakoor<-dendat[ind,coordi]
  }
  ind<-ordobs[maara]
  havakoor<-dendat[ind,coordi] 
  if (havakoor<=leftrecend) lcount<-lcount+1

  leftend<-lcount  #/* could be zero */

return(leftend)
}

