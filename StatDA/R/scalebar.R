scalebar <-
function(Xlowerleft,Ylowerleft,Xupperright,Yupperright,shifttext,shiftkm,sizetext)
{
# draw unit bar on map
#
rect(Xlowerleft,Ylowerleft,Xupperright,Yupperright)
rect(Xlowerleft+(Xupperright-Xlowerleft)/2,Ylowerleft,Xupperright,Yupperright,col=1)
mtext("0",side=1,at=Xlowerleft,line=shifttext, cex=sizetext)
mtext("50",side=1,at=Xlowerleft+(Xupperright-Xlowerleft)/2,line=shifttext, cex=sizetext)
mtext("100",side=1,at=Xupperright,line=shifttext, cex=sizetext)
mtext("km",side=1,at=Xupperright+shiftkm,line=shifttext, cex=sizetext)
invisible()
}

