boxMANA <-
function(comp,type="perc",ymax=NULL,ymin=NULL,yunit=NULL,leg="topright",   #bootstrap or jackknife object
  cex_ylab=1,cex_yaxis=1,cex_names=1) {
  mater<- grep("maternal", colnames(comp))
  add<- grep("additive", colnames(comp))
  nonadd<- grep("nonadd", colnames(comp))
if (type == "perc") {
  dat<- stack(100*comp[,c(add,nonadd,mater)]/ comp$Total)
  colnames(dat)<- c("variance","component")
  label<- "phenotypic variance (%)" }
if (type == "raw") {
  dat<- stack(comp[,c(add,nonadd,mater)])
  colnames(dat)<- c("variance","component")
  label<- "phenotypic variance" }
dat$component<- factor(dat$component, levels(dat$component)[c(1,3,2)])
 if (is.null(ymax)) { ymax<- max(dat$variance) }
 if (!is.null(ymax)) { ymax<- ymax }
 if (is.null(ymin)) { ymin<- 0 }
 if (!is.null(ymin)) { ymin<- ymin }
if (is.null(comp$trait)) { num<-1; name_lab<- ""
 box_plot<- boxplot(variance~ component, dat, pch=20, las=2, xaxt='n', yaxt='n',
   ylim=c(ymin,ymax),ylab=label,cex.lab=cex_ylab,col=c("gray55","gray","gray95"),at=1:3)
 legend(paste(leg),c("additive","non-additive","maternal"),fill=c("gray55","gray","gray95")) }
if (!is.null(comp$trait)) { num <- length(levels(comp$trait)); name_lab<- levels(comp$trait)
  dat$trait<- rep(levels(comp$trait),tapply(comp[,1],comp$trait,length))
  full<- 1:(num*4);rem<- seq(4,num*4,4); loc<- full[-rem]
box_plot<- boxplot(variance~ component + trait, dat, pch=20, las=2, xaxt='n', yaxt='n',
   ylim=c(ymin,ymax),ylab=label,cex.lab=cex_ylab,col=c("gray55","gray","gray95"),at=loc)
  legend(paste(leg),c("additive","non-additive","maternal"),fill=c("gray55","gray","gray95"))
  axis(1, at=seq(2,num*4,4),labels=name_lab,cex.axis=cex_names)  }
  if (is.null(yunit)) { yunit<- (ymax-ymin)/5 }
  if (!is.null(yunit)) { yunit<- yunit }
  axis(2, at=seq(ymin,ymax,yunit),labels=seq(ymin,ymax,yunit),las=1,cex.axis=cex_yaxis)
}
