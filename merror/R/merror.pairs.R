"merror.pairs" <- function(df,labels=names(df))
{
  pairs(df,xlim=range(df,na.rm=TRUE),ylim=range(df,na.rm=TRUE),
    upper.panel=panel.merror,lower.panel=NULL,labels=labels)
}
