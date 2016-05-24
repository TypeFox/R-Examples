sumcol <-
function(rowdat){

  output<-c(sum(rowdat[3:10]), sum(rowdat[11:12]), sum(rowdat[13:16]), sum(rowdat[17]))
  return(t(output))
}

