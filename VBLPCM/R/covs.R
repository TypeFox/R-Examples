vblpcmcovs<-function(N, model, Y, edgecovs=NULL, sendcovs=NULL, receivecovs=NULL, socialcovs=NULL)
  {
  XX_e<-matrix(rep(1,N^2),ncol=1) # all get the intercept term
  if (!is.null(sendcovs)) # model sendcovs as edgecovs
    {
    if (!is.matrix(sendcovs))
      {
      print("Please specify the sender covariates as a matrix with N rows")
      print("running model without sender covariates")
      } else
    for (i in 1:ncol(sendcovs))
      XX_e<-cbind(XX_e,c(t(matrix(rep(sendcovs[,i],N),N))))
    }
  if (!is.null(receivecovs)) # model receivecovs as edgecovs
    {
    if (!is.matrix(receivecovs))
      {
      print("Please specify the receiver covariates as a matrix with N rows")
      print("running model without receiver covariates")
      } else
    for (i in 1:ncol(receivecovs))
      XX_e<-cbind(XX_e,rep(receivecovs[,i],N))
    }
  if (!is.null(socialcovs)) # model socialcovs as edgecovs
    {
    if (!is.matrix(socialcovs))
      {
      print("Please specify the social covariates as a matrix with N rows")
      print("running model without social covariates")
      } else
    for (i in 1:ncol(socialcovs))
      XX_e<-cbind(XX_e,c(t(matrix(rep(socialcovs[,i],N),N)))+rep(socialcovs[,i],N))
    }
  if (!is.null(edgecovs))
    XX_e<-cbind(XX_e,edgecovs)
  return(list("XX_e"=XX_e))
  }
  
  
