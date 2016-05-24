"con.comp" <-
function(comat){
  nc <- ncol(comat)
#  print(nc)
  ccn <- rep(0, times=nc) # con.comp number for each point
  fhist <- rep(FALSE, times=nc) # indicator if point had been under consideration  
  stn <- 0                  # current cc number
  pn <- 1                  # point nr. to which similar objects are looked for
  while(pn>0){
    stn <- stn+1
    repeat{
      sm <- 0              # smallest new point nr.
      ccn[pn] <- stn
      fhist[pn] <- TRUE
      if(nc>1)
      {
  	for(i in 2:nc)
        {
# print(comat[i,])
# cat(i, pn, ccn[i], comat[i,pn],"\n")
  	  if((ccn[i]==0) & (comat[i,pn]))
            ccn[i] <- stn
  	  if ((sm==0) & (ccn[i]==stn) & (fhist[i]==FALSE))
            sm <- i
  	} # for i
# cat("stn=", stn, "sm=", sm, "\n")
      } # if nc>1
      if (sm>0)
        pn <- sm
      else
        break
    } # repeat
#    print("repeat terminated")
    pn <- 0
    i <- 2
    while(i<=nc){
      if(ccn[i]==0){
        pn <- i
        i <- nc
      } # if
      i <- i+1
    } # while i
  } # while pn>0 (stn-loop)
  ccn
}
