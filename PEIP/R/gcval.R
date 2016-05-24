gcval <-
function(U,s,b,npoints)
  {
### %  [reg_min,g,alpha] = gcval
### % Smallest regularization parameter.
    
    smin_ratio = 16*.Machine$double.eps;

### % get matrix sizes
    m = dim(U);
    m = m[1]
    p = dim(s); 
    p = p[1]

### % project the data into a more useful space
    beta = t(U) %*% b;; 

### % sort s in the opposite order and divide the first element in each row by
### % the second to evaluate gamma


    gamma = s[seq(from=p, by=-1, to=1),1] /  s[  seq(from=p, by=-1, to=1) ,2]


    beta = beta[seq(from=m, by=-1, to=1)]

    
### % Vector of regularization parameters. 
    alpha = rep(0, npoints );
### %initiatize g
    g = rep(0, npoints );
    gamma2 = gamma^2; 

### % the last alpha is the larger of the smallest ratio, or the largest ratio 
### % times a very small number


    alpha[npoints]  = max(c( gamma[p],gamma[1]*smin_ratio) );


    r = (gamma[1]/alpha[npoints])^(1/(npoints-1)); 
### % logarithmically distribute the alphas to be used
    for(i in seq(from=npoints-1, by=-1, to=1))
      {
        alpha[i] = r*alpha[i+1];
      }
    
### % Evaluate GCV function values. 
    for(i in 1:npoints )
      {
        g[i] = gcv_function(alpha[i],gamma2,beta); 
      }   

### % find the minimal value of g and save it into a output variable
    mingi = which.min(g);
    reg_min=alpha[mingi];


    return(list( reg_min=reg_min, g=g,  alpha=alpha ) ) 

  }
