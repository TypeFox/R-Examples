mohrleg <-
function(ES)
  {

    u = par('usr')


    for(i in 1:length(ES$values))
      {

        tex1 = substitute(sigma[y]==x , list(x=ES$values[i], y=i) )
        
        mtext(tex1, side = 3, line = -i, at=u[2], adj=1 )

      }

  }

