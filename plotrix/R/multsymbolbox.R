multsymbolbox<-function(x1,y1,x2,y2,tot,relw=0.8,fg=par("fg"),bg=par("bg"),
 box=TRUE,debug=FALSE,...) {

 x1 <- rep(x1, length(tot))
 y1 <- rep(y1, length(tot))
 x2 <- rep(x2, length(tot))
 y2 <- rep(y2, length(tot))
 fg <- rep(fg, length(tot))
 bg <- rep(bg, length(tot))
 for (i in 1:length(tot)) {
  if (tot[i] > 0)
   symbolbox(x1[i],y1[i],x2[i],y2[i],tot[i],relw=relw,
    fg=fg[i],bg=bg[i],box=box,debug=debug,...)
 }
}
