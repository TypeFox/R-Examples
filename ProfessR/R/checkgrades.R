checkgrades<-function(D1, id=NULL, names=NULL)
  {
    if(missing(id)) { id = NULL }
    if(missing(names)) { names= NULL }
   
    o = order(D1$grades)
    DF = data.frame(grades=D1$grades[o], lett=D1$lett[o], scor=D1$scor[o])
    
    if(!is.null(id))
      {
        DF$id = id[o]
      }
    
    if(!is.null(names))
      {
        DF$name = names[o]
      }

    
    print(DF )

    
  }
