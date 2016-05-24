# plot ROC with ggplot2
library(ggplot2)
ggroc <- function(data=data,bin=0.01,roccol="green",sp=19,output="roc.pdf")
{
  pn <- nrow(subset(data,data[,2]== 1))
  fn <- nrow(data) - pn
  diag = data.frame(x = seq(0,1,by=bin),y = seq(0,1,by=bin))
  cutoffs <- seq(0,1,by=bin)
  x = 0
  y = 0
  for (i in cutoffs)
  {
    tpn <- nrow(subset(data, data[,1] >= i & data[,2] == 1))
    fpn <- nrow(subset(data, data[,1] >= i & data[,2] == 0))
    fnn <- pn - tpn
    tnn <- fn - fpn
    tpr <- tpn/(tpn + fnn)
    fpr <- fpn/(fpn + tnn)
    x <-c(x,fpr)
    y <- c(y,tpr)
  }
  FPR = ''
  TPR = ''
  rocdata <- data.frame("FPR" = x, "TPR" = y)
  p <- ggplot(data=rocdata, aes(x = FPR, y = TPR)) + geom_point(color=roccol) + geom_line(color=roccol) + geom_line(data=diag, aes(x=x,y=y),color="red") 
  f <- p + geom_point(data=diag, aes(x=x,y=y),color="red",shape=sp) + theme(axis.text = element_text(size=16), title=element_text(size=18)) + labs(x="False Positive Rate", y="True Positive Rate", title="ROC curve")
  ggsave(f, filename = output, width=8,height=6,units=c("in"))
}

