getmagsize<-function( mag, minsize=1,   slope=1, minmag=0,  maxmag=8, style=1)
  {
    if(missing(minsize)) minsize=1
    mdiff = (maxmag - minmag);

    hsize = minsize;
    if(style==0)
      {
        
        hsize = minsize;
      }


    if(style==1)
      {
        
        hsize = minsize + 0.008 * exp(mag);
        hsize[hsize > 50.0] =   37.5 + .25 * hsize[hsize > 50.0] 
      }

    if(style==2)
      {
        
        if((maxmag - minmag) == 0)
          { hsize = minsize;}
        else
          {
            hsize = minsize + slope * (mag - minmag) / (maxmag - minmag)
          }
      }

    return(hsize);
}
