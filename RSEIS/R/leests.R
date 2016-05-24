`leests` <-
function(a, dt=0.008)
  {
   if(missing(dt)) { dt=1 }

   myts = list()
   
    if(is.ts(a))
      {
        y = a 
        dt = deltat(a)

        myts = list(y=y, dt=dt)

        return(myts)
      }

   if(is.list(a) )
     {
       
       if( any(names(a)=="y") )
         {
           y = a$y
         }
       
       dt = 1
       if( !any(names(a)=="dt") )
         {  
           if( !any(names(a)=="deltat") )
             {
               dt=1
             }else{
               dt = a$deltat
             }
        
         }else{
           dt = a$dt
         }
       
       
     }
   else
     {

       y = a
     }
   
   myts = list(y=y, dt=dt)

   return(myts)

  }

