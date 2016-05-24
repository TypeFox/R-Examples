Northarrow <-
function(Xbottom,Ybottom,Xtop,Ytop,Xtext,Ytext,Alength,Aangle,Alwd,Tcex)
{
# Draw North arrow on a map
#
arrows(Xbottom,Ybottom,Xtop,Ytop,length=Alength,angle=Aangle,lwd=Alwd)
text(Xtext,Ytext,"N",cex=Tcex)
invisible()
}

