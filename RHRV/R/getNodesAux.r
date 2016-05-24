
getNodesAux=function(f,sampling,error,type,deltaOriginalInterval,relative)
{   # auxliary variables
    found=FALSE;
    i=1;
    n_interval_in=0;
    # Iterate until the node is found
    while(found==FALSE)
    {
          bj=n_interval_in
          # width of the interval
          delta = sampling/(2^(i+1));
          # searching
          for (j in (bj*2):(bj*2+1))
          {
            interval= c(j*delta,(j+1)*delta);
            if (f %rhrv_in% interval)
                 n_interval_in=j
            found =  getError(f,interval,type,deltaOriginalInterval,relative)<error;
            if (found)
               break       
           

          }
          if (found==FALSE)
              i=i+1;

    }
    return(c(i,j))
}