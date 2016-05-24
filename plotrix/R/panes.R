panes<-function(mat=NULL,widths=rep(1,ncol(mat)),heights=rep(1,nrow(mat)),
 nrow=2,ncol=2,mar=c(0,0,1.6,0),oma=c(2.5,1,1,1)) {

 oldpar<-par("mar","mfrow","oma")
 if(is.null(mat)) par(mfrow=c(nrow,ncol),mar=mar,oma=oma)
 else layout(mat,widths=widths,heights=heights)
 return(oldpar)
}
