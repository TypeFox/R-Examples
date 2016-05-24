plot_bg<-function(col="lightgray") {
 plim<-par("usr")
 rect(plim[1],plim[3],plim[2],plim[4],col=col)
}
