`convert2strings` <-
function(pattmat){
################################################################
# initializes resulting string
# uses p2string
################################################################
  strvec<-NULL
  for (i in 1:nrow(pattmat))
      strvec[i]<-p2string(pattmat[i,])
  strvec
}
