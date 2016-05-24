draw.brachiopod<-function(BK=matrix(), x=0, y=0,  col="black", fill=NULL, ... )
  {
    if(missing(x)) x = 0
    if(missing(y)) y=0
    if(missing(col)) col="black"
    if(missing(fill)) fill=NULL
    if(missing(BK))
      {
        BK =  get.brachiopod()
        ##############   default brachiopod looks like this
####   bx =  c(-.5,  -0.4648232373807,   -0.3765825678890,-0.2144105266610,   -0.0593931343107,
####     0.0593931343107,0.2144105266610,0.3765825678890,0.4648232373807, 0.5)
####    by =    c(0.0 , 0.139203811333, 0.248963900427, 0.372735915788,0.419442336678,
####     0.419442336678,0.372735915788,0.248963900427,0.139203811333,  0.0)
####   BK = cbind(bx, by)     
      }

    if(!is.null(fill))
      {
        
        polygon(BK[,1], BK[,2], col=fill, border=NA)


      }


    polygon(BK[,1], BK[,2], col=NA, border=col  )

    cenx = mean(c(BK[1,1], BK[10,1]))
    ceny =  mean(c(BK[1,2], BK[10,2]))

    segments(cenx     , ceny   ,
             BK[,1], BK[,2] , col=col, ... )
    

  }

