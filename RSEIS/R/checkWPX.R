checkWPX<-function(wpx)
  {
    ###    returns 0 for no problems; 1,2,3,4 for problems
    ####  return error code
    SUCCESScode = 0
    
    zpx = cleanWPX()
    nam1 = names(zpx)
    nam2 = names(wpx)
###  check to make sure all the elements of a basic wpx list are there	
    m1 = match(nam1, nam2)
    
    if(any(is.na(m1)))
      {
        print("Error WPX list is incomplete")
        ww = which(is.na(m1))
        p1 = paste(collapse=" ", nam1[ww])
        print(paste("MISSING:",p1))
        SUCCESScode = 1
      
      }
    
### check station names and components

    if(any(is.na(wpx$name)))
       {
         print("Error WPX: station names")
         SUCCESScode = 2
       }
    if(any(is.na(wpx$comp)))
       {
         print("Error WPX: component names")
         SUCCESScode = 3
       }

    ###  check the dates and times
    if(any(is.na(c(wpx$yr, wpx$jd, wpx$hr, wpx$mi, wpx$sec)         ) ))
       {
         print("Error WPX: incomplete dates")
         SUCCESScode = 4
       }


    
  ##   WPX = data.frame(WPX, stringsAsFactors = FALSE)
    invisible(SUCCESScode)
  }
