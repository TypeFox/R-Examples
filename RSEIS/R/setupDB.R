setupDB<-function(DB, token=TRUE, split="\\.")
  {
    if(missing(token)) { token=TRUE }
    if(missing(split)) {   split="\\."  }

    
    ruler(DB$fn[1])

    
    if(token)
      {
        SP1 = strsplit(DB$fn[1], split=split)

        print(SP1)

        readline(prompt = "First HIT return, then answer question:")
        
        stacolmns = readline(prompt = "Type in the token location of the station and channel: ")

        SC1 = as.numeric(unlist( strsplit(stacolmns, split=" ")  ))


        SPLIfn = strsplit(DB$fn, split="\\.")

        das1 = as.vector(unlist(lapply(SPLIfn, getmem, SC1[1])))


        comp1 = as.vector(unlist(lapply(SPLIfn, getmem, SC1[2])))

      }
    else
      {

        readline(prompt = "HIT return")
        
        stacolmns = readline(prompt = "Type in the columns for the station name:")
        
        SC1 = as.numeric(unlist( strsplit(stacolmns, split=" ")  ))
        
        
        das1  = substr(DB$fn, SC1[1], SC1[2])
        
        compcolmns = readline(prompt = "Type in the columns for the component name:")
        
        SC2 = as.numeric(unlist( strsplit(compcolmns, split=" ")  ))
        
        comp1 = substr(DB$fn, SC2[1], SC2[2])
################################### 
      }

    
################################### 
    
#############  set up beginning and ending times for each file:
   origyr = min(DB$yr)
    attr(DB, "origyr")<- origyr

    
    eday = EPOCHday(DB$yr, jd =  DB$jd, origyr=origyr)
    DB$t1 = eday$jday + DB$hr/24 + DB$mi/(24*60) + DB$sec/(24*3600)
    DB$t2 = DB$t1 + DB$dur/(24*3600)
    
    DB$sta = das1
    DB$comp  = comp1
    return(DB)

  }

