library(gains)

data(ciaScores)
g1<-gains(actual=ciaScores$CellPhonesPP[ciaScores$train==1],predicted=ciaScores$PredOLS[ciaScores$train==1],
      optimal=TRUE)
g2<-gains(actual=ciaScores$CellPhonesPP[ciaScores$train==0],predicted=ciaScores$PredOLS[ciaScores$train==0],
      conf="t")
print(g1)
print(g2)
