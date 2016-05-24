  dms<-function(d1)
          {
            left1 =  d1*(3600);
            
            degs = trunc( left1 / (3600));
            left2 =  left1 - degs*(3600);
            mins = trunc(left2/60.0);
            osec = left2 - mins*60;
            
            return(list(d=degs, m=mins, s=osec))
          }
