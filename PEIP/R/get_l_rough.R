get_l_rough <-
function(n,deg)
  {
   ##  require(Matrix)
    
    if(deg < 0 | floor(deg) != deg)
      {
        print('degree must be a non-negative integer');
        return(NULL)
      }


    if(deg==0)
      {
        L=diag(n)
        return(L)
      }

    ## % let df approximate the first derivative

    df=c(-1,1,rep(0,length=deg-1))   


    for(i in 2:deg)
      {
        ## % take the difference of the lower order derivative and itself shifted left 
        ## % to get a derivative one order higher
        df=c(0,df[1:deg])-c(df[1:deg],0);
      }

    dn=n-deg;
    L=Matrix::sparseMatrix(i = n-deg, j= n, x=0);

    ## % add the ith element of df to L i-1 elements to the right of the diagonal
    for(i in 1:(deg+1))
      {

        
        L=L+Matrix::sparseMatrix(i=1:(n-deg),
          j= c(1:dn)+i-1,
          x = df[i]*rep(1,length=dn),
          dims=c(dn,n)   );
      }

    
    ## L=unpack(L);

    return(L)

  }
