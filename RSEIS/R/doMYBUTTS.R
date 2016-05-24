`doMYBUTTS` <-
function(butt="",clicks=NULL, x=NULL)
  {

    print(paste(sep=' ', "Hi there you dummy", butt))
    print(clicks)
    
    if(identical("MED", butt) )
          {
            print(median(x))
          }
    if(identical("AVE", butt) )
          {
            print(mean(x))
          }
    if(identical("MIN", butt) )
          {
            print(min(x))
          }

    
  }

