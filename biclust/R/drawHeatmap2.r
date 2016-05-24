
# Draws A as a heatmap, with rows and columns reordered as bicluster rows and
# columns
drawHeatmap2=function(x, bicResult=NULL, number=NA, plotAll=FALSE)
  {
  if(is.null(bicResult))
    {#draw just the matrix
    n=dim(x)[1]
    m=dim(x)[2]

    #Color palette
    numColores=255*2
    gvect=c(array(255:0),array(0,dim=255))
    rvect=c(array(0,dim=255),array(0:255))
    bvect=array(0,dim=numColores)
    paleta=rgb(rvect, gvect, bvect, 255, maxColorValue=255)

     oldmai=par("mai")
     oldmar=par("mar")
     par(mai=c(0,0,0,0),mar=c(0,0,0,0))

    image(1:m,1:n, t(x), col=paleta, axes=FALSE)
    
    par(mai=oldmai, mar=oldmar)
    }
  else
    {
    
    if(plotAll)
      {
            n = dim(x)[1]
            m = dim(x)[2]
            numColores = 255 * 2
            gvect = c(array(255:0), array(0, dim = 255))
            rvect = c(array(0, dim = 255), array(0:255))
            bvect = array(0, dim = numColores)
            paleta = rgb(rvect, gvect, bvect, 255, maxColorValue = 255)
            oldmai = par("mai")
            oldmar = par("mar")
            par(mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
            bicRows = row(matrix(bicResult@RowxNumber[, 1]))[bicResult@RowxNumber[,1] == T]
            bicCols = row(matrix(bicResult@NumberxCol[1,]))[bicResult@NumberxCol[1, ] == T]
            rowlength<-1:bicResult@Number
            rowlength[1]<-length(bicRows)
            collength<-1:bicResult@Number
            collength[1]<-length(bicCols)
            if(bicResult@Number>=2) 
            {
            for (i in 2:bicResult@Number)
            {
            bicRows = c(bicRows,row(matrix(bicResult@RowxNumber[,i]))[bicResult@RowxNumber[,i] == T])
            rowlength[i]<-length(bicRows)
            bicCols = c(bicCols,row(matrix(bicResult@NumberxCol[i,]))[bicResult@NumberxCol[i, ] == T])
            collength[i]<-length(bicCols)
            }
            }
            image(1:m, 1:n, t(x[c(setdiff(c(1:n), bicRows), bicRows),
                c(bicCols, setdiff(c(1:m), bicCols))]), col = paleta,
                axes = FALSE)


            desp = (n-rowlength[1])/n
            grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(desp, 
                desp), "npc"), gp = gpar(col = "yellow"))
            desp = (collength[1])/m
            grid.lines(y = unit(c(0, 1), "npc"), x = unit(c(desp, 
                desp), "npc"), gp = gpar(col = "yellow"))
                
            for (i in 2:bicResult@Number)
            {
            desp = (n- rowlength[i])/n
            grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(desp,
                desp), "npc"), gp = gpar(col = "yellow"))
            desp = (collength[i]-collength[i-1])/m
            grid.lines(y = unit(c(0, 1), "npc"), x = unit(c(desp,
                desp), "npc"), gp = gpar(col = "yellow"))

            }

            par(mai = oldmai, mar = oldmar)
      }
    else
      {
      if(is.na(number) || number>bicResult@Number || number<=0)
      {
      stop("Error: the bicluster does not exist in the result set",call.=FALSE)
      }
      n=dim(x)[1]
      m=dim(x)[2]
  
      #Color palette
      numColores=255*2
      gvect=c(array(255:0),array(0,dim=255))
      rvect=c(array(0,dim=255),array(0:255))
      bvect=array(0,dim=numColores)
      paleta=rgb(rvect, gvect, bvect, 255, maxColorValue=255)
  
      oldmai=par("mai")
      oldmar=par("mar")
      par(mai=c(0,0,0,0),mar=c(0,0,0,0))
  
      bicRows=row(matrix(bicResult@RowxNumber[,number]))[bicResult@RowxNumber[,number]==T]
      bicCols=row(matrix(bicResult@NumberxCol[number,]))[bicResult@NumberxCol[number,]==T]
      image(1:m,1:n,
             t(x[c(setdiff(c(1:n),bicRows), bicRows),
                c(bicCols,setdiff(c(1:m),bicCols))]),
            col=paleta, axes=FALSE)
  
  
  
      desp=(n-length(bicRows))/n
      grid.lines(x=unit(c(0,1),"npc"),y=unit(c(desp,desp),"npc"), gp=gpar(col="yellow"))
      desp=length(bicCols)/m
      grid.lines(y=unit(c(0,1),"npc"),x=unit(c(desp,desp),"npc"), gp=gpar(col="yellow"))
      
      par(mai=oldmai, mar=oldmar)
      }
    }
  }