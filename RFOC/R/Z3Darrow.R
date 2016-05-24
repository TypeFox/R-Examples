`Z3Darrow` <-
function(len = 1 , basethick = 0.1, headlen = .6, headlip=.1 )
  {

    if(missing(headlen)) { headlen = .6  }
    if(missing(len)) {  len = 1 }
    if(missing(basethick)) {  basethick = 0.1 }
    if(missing(headlip)) { headlip=.1  }

    headthick = basethick+headlip
    
    aglyph = list()
    anorm = list()

###########  base rect part
###########  Bottom 
    i = 1
    aglyph[[i]] = matrix(c(-basethick, -basethick, 0,
            -basethick, basethick,0,
            basethick, basethick,0,
            basethick, -basethick, 0) ,  ncol=3, nrow=4, byrow=TRUE)
    anorm[[i]] = c(0,0,-1)

    ##  upper
    i = 1+i
    aglyph[[i]] = matrix(c(   -basethick, basethick,0,
            -basethick, basethick,len,
            basethick, basethick, len,
            basethick,  basethick, 0) ,ncol=3, nrow=4, byrow=TRUE)
    anorm[[i]] = c(0,1,0)

    ###  left
    i = 1+i
    aglyph[[i]] = matrix(c(  -basethick, -basethick,0,
            -basethick, -basethick,len,
            -basethick,  basethick,len,
            -basethick, basethick, 0) ,ncol=3, nrow=4, byrow=TRUE)
    anorm[[i]] = c(-1,0,0)

    ### lower
    i = 1+i
    aglyph[[i]] = matrix(c( basethick, -basethick,0,
            basethick, -basethick,len,
            -basethick,  -basethick,len,
            -basethick,  -basethick, 0 ) ,ncol=3, nrow=4, byrow=TRUE)
    anorm[[i]] = c(0,-1,0)


    
    i = 1+i
    aglyph[[i]] = matrix(c( basethick, -basethick, 0,
            basethick, basethick, 0,
            basethick, basethick,len,
            basethick, -basethick, len ) ,ncol=3, nrow=4, byrow=TRUE)
  anorm[[i]] = c(1,0,0)

    
###########  top
    i = 1+i
    aglyph[[i]] = matrix(c(-basethick,-basethick, len,
                           basethick , -basethick,len,
                           basethick , basethick,len,
                           -basethick, basethick,len) ,  ncol=3, nrow=4, byrow=TRUE)
     anorm[[i]] = c(0,0,1)
 

###########  head part (5 pieces)

   ## if(FALSE){
    ###  bottom of pyramid
    i = 1+i
    aglyph[[i]] = matrix(c(   -headthick, -headthick, len,
            -headthick,  headthick, len,
            headthick,  headthick, len,
            headthick, -headthick, len) ,  ncol=3, nrow=4, byrow=TRUE)
    anorm[[i]] = c(0,0,-1)
    ##  right panel
    i = 1+i
    aglyph[[i]] = matrix(c(  0,  0, len+headlen,
            headthick,-headthick, len,
            headthick,headthick, len
            ) ,  ncol=3, nrow=3, byrow=TRUE)
    anorm[[i]]  = RSEIS::xprod( aglyph[[i]][2,]-aglyph[[i]][1,], aglyph[[i]][3,]-aglyph[[i]][1,])
    ##  upper
    i = 1+i
    aglyph[[i]] = matrix(c(    0,         0, len+headlen,
            headthick, headthick, len,
            -headthick,  headthick, len) ,  ncol=3, nrow=3, byrow=TRUE)
   anorm[[i]]  = RSEIS::xprod( aglyph[[i]][2,]-aglyph[[i]][1,], aglyph[[i]][3,]-aglyph[[i]][1,])
  
    ###  left
    i = 1+i
    aglyph[[i]] = matrix(c(     0, 0,len+headlen,
            -headthick , headthick, len,
            -headthick,  -headthick,len) ,  ncol=3, nrow=3, byrow=TRUE)
   anorm[[i]]  = RSEIS::xprod( aglyph[[i]][2,]-aglyph[[i]][1,], aglyph[[i]][3,]-aglyph[[i]][1,])
  
    ##  lower
    i = 1+i
    aglyph[[i]] = matrix(c( 0,0, len+headlen,
            -headthick ,-headthick, len,
            headthick,  -headthick,  len) ,  ncol=3, nrow=3, byrow=TRUE)
   anorm[[i]]  = RSEIS::xprod( aglyph[[i]][2,]-aglyph[[i]][1,], aglyph[[i]][3,]-aglyph[[i]][1,])
  
 ## }

    centroid = rep(0, length=3)
    np = 0
    for(i in 1:length(aglyph))
      {
        centroid[1] = centroid[1]+sum(aglyph[[i]][,1])
        centroid[2] = centroid[2]+sum(aglyph[[i]][,2])
        centroid[3] = centroid[3]+sum(aglyph[[i]][,3])
        np = np+length(aglyph[[i]][,1])
      }

    centroid = centroid/np
    
    attr(aglyph, "centroid") = centroid
    attr(aglyph, "np") = np
 #####   if(TRUE)
   #####   {
#########  test to make sure all the vectors are point in correct direction
     #####   for(i in 1:length(aglyph))
       #####   {
        #####    XX = RSEIS::xprod( aglyph[[i]][2,]-aglyph[[i]][1,], aglyph[[i]][3,]-aglyph[[i]][2,])
            #####print(paste(sep=' ', c(i, XX) ))
        #####  }
   #####   }

    return(list(aglyph=aglyph, anorm=anorm ) )

  }

