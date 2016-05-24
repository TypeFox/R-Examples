plotTabCont0<-function(Table){
M<-max(Table)
barplot(Table,beside=TRUE,ylim=c(0,M*5/4),legend.text=TRUE,args.legend = list(x = "top",cex=0.75),col=rev(brewer.pal(3,name="PuRd")))
devAskNewPage(ask = TRUE)
assoc(Table,gp = shading_hcl)
devAskNewPage(ask = TRUE)
mosaicplot(Table,col=rev(brewer.pal(3,name="PuRd")),las=1,main="")
  if (min(dim(Table)) > 2) {
devAskNewPage(ask = TRUE)
        biplot(corresp(Table, nf = 2))
    }
}
