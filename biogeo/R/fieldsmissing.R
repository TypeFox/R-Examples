fieldsmissing <-
function (dat, fields) 
{
  ck <- checkdatastr(dat)
  if (any(ck$Present == FALSE)==TRUE) {
    fa <- ck[ck$Present == FALSE,1]
    fd<-as.character(fa)
    msg<-paste(fd,collapse=", ")
    stop(paste("Fields missing: ", msg))
}
}
