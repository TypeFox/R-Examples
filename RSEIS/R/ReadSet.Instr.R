`ReadSet.Instr` <-
function(file)
  {
    ##  cmd = paste(sep=" ","cat ", file)
    ##    INSTF = system(cmd, intern=TRUE)

    if(length(file)==1 & file.exists(file[1]) )
      {
        ACON = file(description=file, open = "r")
        INSTF = readLines(con = ACON, n = -1, ok = TRUE, warn = TRUE)
        close(ACON)
      }
    else
      {
        INSTF =file

      }

    poles = NULL
    zeros = NULL
    Knorm = NULL
    Sense = NULL
    nz=0
    np=0

    
    a = unlist(strsplit(INSTF[1],split=' '))
    nz = as.numeric(a[2])
    if(nz>0)
      {
        zeros = vector(length=nz, mode="complex")
        for( i in 1:nz)
          {
            a = unlist(strsplit(INSTF[1+i],split='\\ '))
            a = a[a!=""]
            zeros[i] = complex(real=as.numeric(a[1]), imaginary=as.numeric(a[2]))
          }
      }
    ip= 1+nz+1

    
    a = unlist(strsplit(INSTF[ip],split=' '))
    np = as.numeric(a[2])
    if(np>0)
      {
        poles = vector(length=np, mode="complex")
        for( i in 1:np)
          {
            a = unlist(strsplit(INSTF[ip+i],split='\\ '))
            a = a[a!=""]
            poles[i] = complex(real=as.numeric(a[1]), imaginary=as.numeric(a[2]))
          }
      }
    ip= np+nz+2+1
    a = unlist(strsplit(INSTF[ip],split=' '))
    Knorm =  as.numeric(a[2])
    ip= ip+1
    a = unlist(strsplit(INSTF[ip],split=' '))
    Sense =  as.numeric(a[2])
    return(list(np=np, poles=poles, nz=nz, zeros=zeros, Knorm=Knorm, Sense=Sense))
    
  }

