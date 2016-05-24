
# Draws A as a heatmap, with rows and columns reordered as bicluster rows and
# columns
drawHeatmap=function(x, bicResult=NULL, number=NA, local=TRUE, beamercolor=FALSE,paleta,...)
  {
  n=dim(x)[1]
  m=dim(x)[2]
  
  #Color palette
  numColores=255*2
  gvect=c(array(255:0),array(0,dim=255))
  rvect=c(array(0,dim=255),array(0:255))
  bvect=array(0,dim=numColores)
  if(beamercolor)
  {paleta=paleta}
  else
  {paleta=rgb(rvect, gvect, bvect, 255, maxColorValue=255)}
  
  oldpar=par(c("mai", "mar", "mgp", "xpd"))
  	  	  
  if(is.null(bicResult))
    {#draw just the matrix
	par(mar=c(0,0,0,0), mai=c(0,0,0,0))
	image(1:m,1:n, t(x), col=paleta, axes=FALSE)
    }
  else
    {
    if(is.na(number) || number>bicResult@Number || number<=0)
	    {
	    stop("Error: the bicluster does not exist in the result set")
	    break
	    }
		
	bicRows = which(bicResult@RowxNumber[,number])
	bicCols = which(bicResult@NumberxCol[number, ])
	
	if(local)#Only the heatmap of the bicluster
		{
		x2=x[bicRows,bicCols]
		x3=length(paleta)*(x2-min(x))/(max(x)-min(x))
		paleta2=paleta[min(x3):max(x3)]
		par(mar=c(1,6,7,0)+0.1, mgp=c(0,1,0))
		#par(xpd=T)
		
		image(1:length(bicCols),1:length(bicRows), t(x2), col=paleta2, axes=F, 
				xlab=paste("Bicluster ",number, " (size ", length(bicRows), "x", length(bicCols),")"), 
				ylab="")
		axis(side=3, 1L:length(bicCols), labels = colnames(x2), las = 2, line = -0.5, tick = 0)
		axis(side=2,  1L:length(bicRows), labels = rownames(x2), las = 2, line = -0.5, tick = 0)		
		#mtext(colnames(x2), side=3, 1L:length(bicCols),  las = 2, line = -0.5, srt=45)
		#mtext(rownames(x2), side=2,  1L:length(bicRows), las = 2, line = -0.5, srt=45)
		#text(x=1L:length(bicCols),labels=colnames(x2), srt=45)		
		} 	 
	else #the whole matrix with the bicluster reordered (could be very small sometimes)
	   {
		par(mar=c(0,0,0,0), mai=c(0,0,0,0))
		image(1:m,1:n,
	         t(x[c(setdiff(c(1:n),bicRows), bicRows),
	         c(bicCols,setdiff(c(1:m),bicCols))]),
	         col=paleta, axes=FALSE)
	  
	   desp=(n-length(bicRows))/n
	   grid.lines(x=unit(c(0,1),"npc"),y=unit(c(desp,desp),"npc"), gp=gpar(col="yellow"))
	   desp=length(bicCols)/m
	   grid.lines(y=unit(c(0,1),"npc"),x=unit(c(desp,desp),"npc"), gp=gpar(col="yellow"))
	   }
    }
  par(oldpar)
  }
