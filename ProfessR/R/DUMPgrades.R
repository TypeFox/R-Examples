DUMPgrades<-function(D1, file= NULL, id=NULL, names=NULL)
  {
    if(missing(id)) { id = NULL }
    if(missing(names)) { names= NULL }
    if(missing(file)) { file = NULL }

    DF = data.frame(grades=D1$grades, lett=D1$lett, scor=D1$scor)
    
    if(!is.null(id))
      {
        DF$id = id
      }
    
    if(!is.null(names))
      {
        DF$name = names
      }

    
    write.table(file=paste(sep=".", file, "csv"),  DF , sep=",", row.names =FALSE )

    
  }
