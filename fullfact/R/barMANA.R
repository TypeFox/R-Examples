barMANA <-
function(ci_dat,type="perc",bar_len=0.1,ymax=NULL,ymin=NULL,yunit=NULL,leg="topright",     #ci object
  cex_ylab=1,cex_yaxis=1,cex_names=1) {
  error.bar <- function(x, y, upper, lower=upper, length=bar_len,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)}
if (type == "perc") { dat<- ci_dat$percentage; label<- "phenotypic variance (%)" }
if (type == "raw") { dat<- ci_dat$raw; label<- "phenotypic variance" }
if (is.null(dat$trait)) { num <- 1; name_lab<- ""
  ord<- matrix(0,ncol=1,nrow=3)
  ord[,1][1]<- which(dat$component=="additive")
  ord[,1][2]<- which(dat$component=="nonadd")
  ord[,1][3]<- which(dat$component=="maternal") }
if (!is.null(dat$trait)) { num <- length(levels(dat$trait)); name_lab<- levels(dat$trait)
  ord<- matrix(0,ncol=1,nrow=3*num)
  ord[,1][seq(1,3*num,3)]<- which(dat$component=="additive")
  ord[,1][seq(2,3*num,3)]<- which(dat$component=="nonadd")
  ord[,1][seq(3,3*num,3)]<- which(dat$component=="maternal") }
 dat_ci<- matrix(dat[,3][ord],ncol=num,nrow=3) #median/mean
 lwr_ci<- matrix(dat$lower[ord],ncol=num,nrow=3)
 upp_ci<- matrix(dat$upper[ord],ncol=num,nrow=3)
 if (is.null(ymax)) { ymax<- max(dat$upper[ord]) }
 if (!is.null(ymax)) { ymax<- ymax }
 if (is.null(ymin)) { ymin<- 0 }
 if (!is.null(ymin)) { ymin<- ymin }
ci_plot<- barplot(dat_ci,beside=T,ylab=label,col=c("gray55","gray","gray95"),
  names.arg=name_lab,cex.names=cex_names,
  yaxt='n',ylim=c(ymin,ymax),cex.lab=cex_ylab)
  error.bar(ci_plot,dat_ci,upper=upp_ci,lower=lwr_ci)
  legend(paste(leg),c("additive","non-additive","maternal"),fill=c("gray55","gray","gray95"))
  if (is.null(yunit)) { yunit<- (ymax-ymin)/5 }
  if (!is.null(yunit)) { yunit<- yunit }
  axis(1, at=c(0,4*num),labels=FALSE)
  axis(2, at=seq(ymin,ymax,yunit),labels=seq(ymin,ymax,yunit),las=1,cex.axis=cex_yaxis)
}
