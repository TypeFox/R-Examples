loadMAT <-
function(fn, pos=1)
  {
    load(fn)


####  what is in the data set:
    NG = names(G)

    GIMPUT = G
    NG = names(GIMPUT)

    G = GIMPUT$G

    for(i in 1:length(NG))
      {
        assign(NG[i],GIMPUT[[i]],  pos=pos )
      }



  }
