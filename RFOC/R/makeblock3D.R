`makeblock3D` <-
function(block1)
  {

    gblock= list()
     
    anorm = list()

    i = 1
    v = c(1,2,3,4)
    gblock[[i]] =  matrix(c(block1[v[1],1], block1[v[1],2], block1[v[1],3],
            block1[v[2],1], block1[v[2],2], block1[v[2],3],
            block1[v[3],1], block1[v[3],2], block1[v[3],3],
            block1[v[4],1], block1[v[4],2], block1[v[4],3] ) ,  ncol=3, nrow=4, byrow=TRUE)
     anorm[[i]]  = RSEIS::xprod( gblock[[i]][2,]-gblock[[i]][1,], gblock[[i]][3,]-gblock[[i]][2,])
     anorm[[i]] = anorm[[i]]/sqrt(sum(anorm[[i]]^2))
  
    i = i+1
    v = c(1,5,6,2)
    gblock[[i]] =  matrix(c(block1[v[1],1], block1[v[1],2], block1[v[1],3],
            block1[v[2],1], block1[v[2],2], block1[v[2],3],
            block1[v[3],1], block1[v[3],2], block1[v[3],3],
            block1[v[4],1], block1[v[4],2], block1[v[4],3] ) ,  ncol=3, nrow=4, byrow=TRUE)
        anorm[[i]]  = RSEIS::xprod( gblock[[i]][2,]-gblock[[i]][1,], gblock[[i]][3,]-gblock[[i]][2,])
    anorm[[i]] = anorm[[i]]/sqrt(sum(anorm[[i]]^2))

    i = i+1
    v = c(1,4,8,5)
    gblock[[i]] =  matrix(c(block1[v[1],1], block1[v[1],2], block1[v[1],3],
            block1[v[2],1], block1[v[2],2], block1[v[2],3],
            block1[v[3],1], block1[v[3],2], block1[v[3],3],
            block1[v[4],1], block1[v[4],2], block1[v[4],3] ) ,  ncol=3, nrow=4, byrow=TRUE)
        anorm[[i]]  = RSEIS::xprod( gblock[[i]][2,]-gblock[[i]][1,], gblock[[i]][3,]-gblock[[i]][2,])
    anorm[[i]] = anorm[[i]]/sqrt(sum(anorm[[i]]^2))
    i = i+1
    v = c(5,8,7,6)
    gblock[[i]] =  matrix(c(block1[v[1],1], block1[v[1],2], block1[v[1],3],
            block1[v[2],1], block1[v[2],2], block1[v[2],3],
            block1[v[3],1], block1[v[3],2], block1[v[3],3],
            block1[v[4],1], block1[v[4],2], block1[v[4],3] ) ,  ncol=3, nrow=4, byrow=TRUE)
        anorm[[i]]  = RSEIS::xprod( gblock[[i]][2,]-gblock[[i]][1,], gblock[[i]][3,]-gblock[[i]][2,])
    anorm[[i]] = anorm[[i]]/sqrt(sum(anorm[[i]]^2))
    
    i = i+1
    v = c(4,3,7,8)
    gblock[[i]] =  matrix(c(block1[v[1],1], block1[v[1],2], block1[v[1],3],
            block1[v[2],1], block1[v[2],2], block1[v[2],3],
            block1[v[3],1], block1[v[3],2], block1[v[3],3],
            block1[v[4],1], block1[v[4],2], block1[v[4],3] ) ,  ncol=3, nrow=4, byrow=TRUE)
        anorm[[i]]  = RSEIS::xprod( gblock[[i]][2,]-gblock[[i]][1,], gblock[[i]][3,]-gblock[[i]][2,])
    anorm[[i]] = anorm[[i]]/sqrt(sum(anorm[[i]]^2))


    i = i+1
    v = c(2,6,7,3)
    gblock[[i]] =  matrix(c(block1[v[1],1], block1[v[1],2], block1[v[1],3],
            block1[v[2],1], block1[v[2],2], block1[v[2],3],
            block1[v[3],1], block1[v[3],2], block1[v[3],3],
            block1[v[4],1], block1[v[4],2], block1[v[4],3] ) ,  ncol=3, nrow=4, byrow=TRUE)
        anorm[[i]]  = RSEIS::xprod( gblock[[i]][2,]-gblock[[i]][1,], gblock[[i]][3,]-gblock[[i]][2,])
    anorm[[i]] = anorm[[i]]/sqrt(sum(anorm[[i]]^2))
    return(list(aglyph=gblock, anorm=anorm ))

  }

