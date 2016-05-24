genOrd <-
function(nObs, ordPmat, binObj){
ep0 = generate.binary( nObs, binObj$pvec, binObj$del.next)
Mydata= BinToOrd(binObj$pvec, ordPmat, binObj$Mlocation, ep0)
return(Mydata)
}
