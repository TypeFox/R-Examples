prepdata <-
function(rawdata) {
  indata = list(SNP=as.character(rawdata$Name))
  indata$SampleID = gsub(".Theta","",names(rawdata)[grep(".Theta",names(rawdata))])
  indata$P = length(indata$SNP)
  indata$N = length(indata$SampleID)
  indata$Theta = as.matrix(rawdata[,grep("Theta",names(rawdata))])
  indata$R = as.matrix(rawdata[,grep("Theta",names(rawdata))+1])
  indata$GType = as.matrix(rawdata[,setdiff(grep("GType",names(rawdata)),grep("Custom.GType",names(rawdata)))])
  indata$Score = as.matrix(rawdata[,setdiff(grep("Score",names(rawdata)),grep("GenTrain.Score",names(rawdata)))])
  indata$X.Raw = as.matrix(rawdata[,grep("X.Raw",names(rawdata))])
  indata$Y.Raw = as.matrix(rawdata[,grep("Y.Raw",names(rawdata))])
  indata$X = as.matrix(rawdata[,grep(".X",names(rawdata))[!is.element(grep(".X",
    names(rawdata)),grep(".Raw",names(rawdata)))]])
  indata$Y = as.matrix(rawdata[,grep(".Y",names(rawdata))[!is.element(grep(".Y",
    names(rawdata)),grep(".Raw",names(rawdata)))]])
  indata$logratio = indata$Theta
  for (j in 1:indata$N) {
    if (sum(indata$X[,j]>0,na.rm=TRUE) > 1) {

      tmpmod <- lm(indata$X[,j][indata$X[,j]>0] ~ 
        indata$X.Raw[,j][indata$X[,j]>0])
      xtmp = indata$X[,j] - tmpmod$coefficients[1]
      bool = indata$X[,j]==0
      bool[is.na(bool)] = FALSE
      xtmp[bool] = indata$X.Raw[,j][bool]*tmpmod$coefficients[2]
      tmpmod = lm(indata$Y[,j][indata$Y[,j]>0] ~ 
        indata$Y.Raw[,j][indata$Y[,j]>0])
      ytmp <- indata$Y[,j] - tmpmod$coefficients[1]
      bool <- indata$Y[,j]==0
      bool[is.na(bool)] = FALSE
      ytmp[bool] <- indata$Y.Raw[,j][bool]*tmpmod$coefficients[2]
      indata$logratio[,j] = log(ytmp/xtmp)
      rm(bool,tmpmod,xtmp,ytmp)
    } else {
      indata$logratio[,j] = rep(NA,indata$P)
    }
  }
  indata$R.trans = indata$R
  for (p in 1:indata$P) {
    if (sum(indata$R[p,]>.05,na.rm=TRUE) > 1) {
      tranpar = boxcox(indata$R[p,][indata$R[p,]>.05] ~ indata$logratio[p,][indata$R[p,]>.05],plotit=FALSE)
      lambda = tranpar$x[tranpar$y==max(tranpar$y)]
      if (is.na(lambda)) lambda = 0
      if (lambda==0) {indata$R.trans[p,] = log(indata$R[p,])} else {
        indata$R.trans[p,] = (indata$R[p,]^lambda-1)/lambda
      }
      rm(lambda,tranpar)
    }
  }
  return(indata)
}
