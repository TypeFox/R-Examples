NoOverlap<-function(x, y, focsiz , SEL=0, OLDx=0, OLDy=0, cenx=0, ceny=0)
  {

######   
    if(missing(cenx)) cenx = mean(x)
    if(missing(ceny)) ceny = mean(y)
    if(missing(SEL)) SEL = 1:length(x)

    if(missing(OLDx)) OLDx=rep(cenx, length(x))
    if(missing(OLDy)) OLDy=rep(ceny, length(y))

    newx = x
    newy = y


    if( any(is.na(newx)) ) {
      print("BAD INPUT")
      ix = which(is.na(newx))
      print(ix)
       return(list(x=newx, y=newy))
      
    }
   if( any(is.na(newy)) ) {
      print("BAD INPUT")
      ix = which(is.na(newy))
      print(ix)
       return(list(x=newx, y=newy))
      
    }


    if( any(is.na(OLDx)) ) {
      print("BAD INPUT")
      ix = which(is.na(OLDx))
      print(ix)
       return(list(x=newx, y=newy))
      
    }
   if( any(is.na(OLDy)) ) {
      print("BAD INPUT")
      ix = which(is.na(OLDy))
      print(ix)
       return(list(x=newx, y=newy))
      
    }
    
###     print(paste('LEN OLDX=', length(OLDx),'LEN OLDY=', length(OLDy)   ))
    
###     print(SEL)
    
###  points(newx[SEL],newy[SEL],pch=6,col=grey(.8))
    
    for(j in 1:length(SEL))
      {
        i = SEL[j]
        points(newx[i],newy[i],pch=6,col="green")
###     print(paste(sep=" ", "check", j, i, newx[i] , OLDx[i] , newy[i],  OLDy[i]))
        
        jx = newx[-i]
        jy = newy[-i]

        if( newx[i]==OLDx[i] & newy[i] == OLDy[i])
          {
            vec =   c( (newx[i]-mean(OLDx) ) , (newy[i]-mean(OLDy) ))
            vec = vec/sqrt(vec[1]^2 + vec[2]^2)
          }
        else
          {
            vec =   c( (newx[i]-OLDx[i]) , (newy[i]-OLDy[i]))
            vec = vec/sqrt(vec[1]^2 + vec[2]^2)
          }
        
        
###     print(paste("##################     working on ", j, i, vec[1], vec[2]  ))
###    cat(cbind(jx, jy), sep="\n")
        
###  for each point keep looping until it is nudged away from all other points
        while(TRUE){

          px = newx[i]
          py = newy[i]
          ##   draw.circ(newx[i], newy[i], focsiz, col='green')
          disi = sqrt( (jx-px)^2 + (jy-py)^2)
          kdis = disi
          
          if(any(is.na(kdis)))
            {
              rdis = c(j, i, px, py, length(kdis) , length(which( is.na(kdis) ))  , range(kdis, na.rm=TRUE))
              
              cat(paste(rdis, collapse=" "), sep="\n")

              ## return(list(x=newx, y=newy))
              break
              
            }
          
          
          if(any(kdis<2*focsiz)   )
            {
### nudge the location of newx newy away from cenx ceny


              
              newx[i]= newx[i]+focsiz*vec[1]
              newy[i]= newy[i]+focsiz*vec[2]



              
            }
          else
            {
              break
            }
        }
        
      }

    return(list(x=newx, y=newy))
    
  }



