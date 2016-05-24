#Function to calculate the chi square statistic to the maximum likelihood estimates of the generalized Poisson and the Poisson model.
calc_chisq_statistic <- function(x,lambda,theta)
{
  n = length(x);
  x_bar = mean(x);
  
  lx = max(x);

  obs = rep(0,lx+1);

  for(i in 1:n)
  {
    obs[x[i]+1] = obs[x[i]+1] + 1;
  }

  mark1 = 0;
  mark2 = 0;
  flag1 = 1; # I don't know what to do with these things
  flag2 = 1;

  chisq1 = 0;
  df1 = 0;
 
  chisq2 = 0;
  df2 = 0;

  expected_gp = rep(0,lx+1);
  expected_p = rep(0,lx+1);

  dont_combine_gp = rep(0,lx+1);
  dont_combine_p = rep(0,lx+1);


#Generalized Poisson Model

  zs = rep(0,lx+1);
  p = rep(0,lx+1);

  zs = theta + (0:lx)*lambda;
  
  zsnew=zs[zs>0];
  len =length(zsnew)-1;

  p = -1*log(gamma(1:(len+1))) + log(theta) +(-1:(len-1))*log(zs[1:(len+1)])-zs[1:(len+1)];

  expected_gp = exp(p)*n;

  dont_combine_gp = expected_gp >= 5;
  sum_cells = sum(dont_combine_gp);
  
  if(sum_cells > 2 && flag1 == 1)
  {
    mark1 = 1;
    df1 = sum_cells-2;

    temp_new = 0;
    temp_new = sum(((obs[dont_combine_gp] - expected_gp[dont_combine_gp])^2)/expected_gp[dont_combine_gp]);
    so = sum(obs[dont_combine_gp]);
    se = sum(expected_gp[dont_combine_gp]);
  
    temp_new = temp_new + ((se-so)^2)/(n-se);
    
    if((n-se) == 0)
    {
      mark1 = 0;
    }
    chisq1 = temp_new;
  }
      
#Poisson Model

  p = log(dpois(0:lx,x_bar));

  expected_p = exp(p)*n;
   
  dont_combine_p = expected_p >=5;

  sum_cells = sum(dont_combine_p);

  if(sum_cells > 1 && flag2 == 1)
  {
    mark2 = 1;
    df2 = sum_cells-1;
    
    temp_new = 0;
    temp_new = sum(((obs[dont_combine_p] - expected_p[dont_combine_p])^2)/expected_p[dont_combine_p]);
    so = sum(obs[dont_combine_p]);
    se = sum(expected_p[dont_combine_p]);

   
    if((n-se) == 0)
    {
      mark2 = 0;
    }
    else
    {
      temp_new = temp_new + ((se-so)^2)/(n-se);
    }
    
    chisq2 = temp_new; 
  }
  
  return(list(mark1=mark1,mark2=mark2,df1=df1,df2=df2,chisq1=chisq1,chisq2=chisq2));
}

likelihood_ratio_generalized_poisson_exon_gene <- function(z,theta1,lambda1,x,theta2,lambda2,v,theta3,lambda3,y,theta4,lambda4)
{
  m = length(z); # = length(v)
  n = length(x); # = length(y)

  lz = max(z);
  obsz = rep(0,lz+1);

  lv = max(v);
  obsv = rep(0,lv+1);

  for(i in 1:m)
  {
    obsz[z[i]+1] = obsz[z[i]+1] + 1;
    obsv[v[i]+1] = obsv[v[i]+1] + 1;
  }

  lx = max(x);
  obsx = rep(0,lx+1);
  
  ly = max(y);
  obsy = rep(0,ly+1);

  for(i in 1:n)
  {
    obsx[x[i]+1] = obsx[x[i]+1] + 1;
    obsy[y[i]+1] = obsy[y[i]+1] + 1;
  }

#Initialization
  s1 = theta2;
  s2 = theta4;
  b = (theta1/theta2 + theta3/theta4)/2;

  mark = 0;
  flag = 1;
#Newton Raphson to get the MLE
  for(its in 1:10000)
  {
    zs1 = b*s1 +(0:lz)*lambda1;
    zs2 = s1 +(0:lx)*lambda2;
    zs3 = b*s2 + (0:lv)*lambda3;
    zs4 = s2 + (0:ly)*lambda4;

    if(sum(zs1==0) > 0 || 
       sum(zs2==0) > 0 ||
       sum(zs3==0) > 0 ||
       sum(zs4==0) > 0 )
    {
      flag = 0;
      break;
    }

    f1 = (m+n)/s1 -m*b -n + sum(obsz[1:(lz+1)]*(-1:(lz-1))*b/zs1) + sum(obsx[1:(lx+1)]*(-1:(lx-1))/zs2);
    f2 = (m+n)/s2 - m*b -n + sum(obsv[1:(lv+1)]*(-1:(lv-1))*b/zs3) + sum(obsy[1:(ly+1)]*(-1:(ly-1))/zs4);
    f3 = (2*m/b) - m*(s1+s2) + sum(obsz[1:(lz+1)]*(-1:(lz-1))*s1/zs1) + sum(obsv[1:(lv+1)]*(-1:(lv-1))*s2/zs3);

    a = -1*(m+n)/(s1^2) + sum(obsz[1:(lz+1)]*(-1:(lz-1))*(-1)*(b^2)/(zs1^2)) + sum(obsx[1:(lx+1)]*(-1:(lx-1))*(-1)/(zs2^2));
    p = -1*m + sum(obsz[1:(lz+1)]*(-1:(lz-1))*(0:lz)*lambda1/(zs1^2));
    c = -1*(m+n)/(s2^2) + sum(obsv[1:(lv+1)]*(-1:(lv-1))*(-1)*(b^2)/(zs3^2)) + sum(obsy[1:(ly+1)]*(-1:(ly-1))*(-1)/(zs4^2));
    d = -1*m + sum(obsv[1:(lv+1)]*(-1:(lv-1))*(0:lv)*lambda3/(zs3^2));
    e = -2*m/(b^2) + sum(obsz[1:(lz+1)]*(-1:(lz-1))*(-1)*(s1^2)/(zs1^2)) + sum(obsv[1:(lv+1)]*(-1:(lv-1))*(-1)*(s2^2)/(zs3^2));

    zs = a*c*e - a*d*d - p*p*c;
    if(is.na(zs))
    {
      flag = 0;
      break;
    }
    if( flag == 1 && zs != 0)
    {
      ds1 = ((-1)*f1*(c*e-d^2)-f2*p*d+f3*p*c)/zs;
      ds2 = ((-1)*f1*p*d - f2*(a*e-p*p)+f3*a*d)/zs;
      db = (f1*p*c + f2*a*d - f3*a*c)/zs;
      err = 0;
    
      zs1 = max(abs(s1),0.01);
      zs2 = abs(ds1)/zs1;
      err = max(err,zs2);
      s1 = s1+ds1;

      zs1 = max(abs(s2),0.01);
      zs2 = abs(ds2)/zs1;
      err = max(err,zs2);
      s2 = s2+ds2;

      zs1 = max(abs(b),0.01);
      zs2 = abs(db)/zs1;
      err = max(err,zs2);
      b = b+db;

      if(err<=1.0e-6 && s1 > 0 && s2 > 0 && b > 0)
      {
        mark = 1;
        break;
      }
    }
    else
    {
      break;
    }
  }

  if(lambda1 < 0)
  {
    tm = floor(-1*s1*b/lambda1);
    if(tm < 4)
    {
      mark = 0;
    }
    zs = -1 *s1*b/tm;
    if(zs < -1 && lambda1 < -1)
    {
      mark = 0;
    }
  }

  if(lambda2 < 0)
  {
    tm = floor(-1*s1/lambda2);
    if(tm < 4)
    {
      mark = 0;
    }
    zs = -1*s1/tm;
    if(zs < -1 && lambda2 < -1)
    {
      mark = 0;
    }
  }

  if(lambda3 < 0)
  {
    tm = floor(-1*s2*b/lambda3);
    if(tm < 4)
    {
      mark = 0;
    }
    zs = -1 *s2*b/tm;
    if(zs < -1 && lambda3 < -1)
    {
      mark = 0;
    }
  }
 
  if(lambda4 < 0)
  {
    tm = floor(-1*s2/lambda4);
    if(tm < 4)
    {
      mark = 0;
    }
    zs = -1*s2/tm;
    if(zs < -1 && lambda4 < -1)
    {
      mark = 0;
    }
  }

#The likelihood ratio test
  if(mark==1)
  {
    zs3_1 = b*s1+(0:lz)*lambda1;
    zs3_2 = theta1 +(0:lz)*lambda1;
    zs3_3 = s1+(0:lx)*lambda2;
    zs3_4 = theta2 + (0:lx)*lambda2;
    zs3_5 = b*s2 + (0:lv)*lambda3;
    zs3_6 = theta3 +(0:lv)*lambda3;
    zs3_7 = s2 + (0:ly)*lambda4;
    zs3_8 = theta4 + (0:ly)*lambda4;

    if(min(zs3_1,zs3_2,zs3_3,zs3_4,zs3_5,zs3_6,zs3_7,zs3_8) <= 0)
    {
      mark = 0;
      break;
    }

    LR1 = (n+m)*log(s1) -(n+b*m)*s1 + 2*m*log(b) + (n+m)*log(s2) - (n+b*m)*s2;
    LR2 = m*log(theta1) + n*log(theta2) + m*log(theta3) + n*log(theta4) - m*theta1 -n*theta2 - m*theta3 - n*theta4;

    LR1 = LR1 + sum(obsz[1:(lz+1)]*(-1:(lz-1))*log(zs3_1)) + sum(obsx[1:(lx+1)]*(-1:(lx-1))*log(zs3_3)) + sum(obsv[1:(lv+1)]*(-1:(lv-1))*log(zs3_5)) + sum(obsy[1:(ly+1)]*(-1:(ly-1))*log(zs3_7));
    LR2 = LR2 + sum(obsz[1:(lz+1)]*(-1:(lz-1))*log(zs3_2)) + sum(obsx[1:(lx+1)]*(-1:(lx-1))*log(zs3_4)) + sum(obsv[1:(lv+1)]*(-1:(lv-1))*log(zs3_6)) + sum(obsy[1:(ly+1)]*(-1:(ly-1))*log(zs3_8));

    test = -2*(LR1-LR2);
    return(list(mark=mark,Gptest=test));
  }
  else
  {
    return(list(mark=0,Gptest=-1));
  }
}

    
        
    
    
likelihood_ratio_poisson_exon_gene<-function(z,x,v,y)
{
  m = length(z);
  n = length(x);

  l5 = mean(z);
  l6 = mean(x);
  l7 = mean(v);
  l8 = mean(y);

  sumz = sum(z);
  sumv = sum(v);
  sumx = sum(x);
  sumy = sum(y);

  b = (sumz+sumv)*n/(m*(sumx+sumy));
  
  l1 = (sumz+sumx)/(n+m*b);
  l2 = (sumv+sumy)/(n+m*b);

  LR1 = sumz*log(l1) + sumx*log(l1) - (n+m*b)*l1 + sumz*log(b) + sumv*log(l2) + sumy*log(l2) - (n+m*b)*l2 + sumv*log(b);
  LR2 = sumz*log(l5) + sumx*log(l6) - m*l5 - n*l6 + sumv*log(l7) + sumy*log(l8) -m*l7 - n*l8;

  test = -2*(LR1-LR2);
  return(list(Ptest=test));
}

likelihood_ratio_tissue_generalized_poisson <- function(x,lambda1,theta1,y,lambda2,theta2,w)
{
  n = length(x);
  m = length(y);

  lx = max(x);
  ly = max(y);

  obsx = rep(0,lx+1);
  obsy = rep(0,ly+1);

#Populating the counts for x
  for(i in 1:n)
  {
    obsx[x[i]+1] = obsx[x[i]+1] + 1;
  }

#Populating the counts for y
  for(i in 1:m)
  {
    obsy[y[i]+1] = obsy[y[i]+1] + 1;
  }

#Calculating parameters for the NULL:theta1 = theta2

  theta = (theta1*n + theta2*m)/(n+m);

  mark = 0;
  flag = 1;

#Newton Raphson Method to Estimate theta
  for(its in 1:10000)
  {
    zs1 = rep(0,lx+1);
    zs2 = rep(0,ly+1);

    zs1 = theta + (0:lx)*lambda1;
    zs2 = w*theta + (0:ly)*lambda2;

    if( sum(zs1 == 0) > 0 || sum(zs2 == 0) > 0)
    {
      flag = 0;
      break;
    }

    hs = (n+m)/theta - n - m*w + sum(obsx[1:(lx+1)]*((-1:(lx-1))/zs1[1:(lx+1)])) + sum(obsy[1:(ly+1)]*((-1:(ly-1))*w/zs2[1:(ly+1)]));
    
    is = -1*(n+m)/(theta^2) + sum(obsx[1:(lx+1)]*(-1)*(-1:(lx-1))/(zs1[1:(lx+1)]^2)) + sum(obsy[1:(ly+1)]*(-1)*(-1:(ly-1))*(w^2)/(zs2[1:(ly+1)]^2));

    if(flag == 1 && is != 0)
    {
      dtheta = -1*hs/is;
      zs1 = max(abs(theta),0.01);
      zs2 = abs(dtheta)/zs1;
      
      theta = dtheta + theta;
      if(zs2 <= 1.0e-6 && theta > 0)
      {
        #Converged
        mark = 1;
        break;
      }
    }
    else
    {
      break;
    }
  }
  

  if(lambda1 < 0)
  {
    tm = floor((-1)*theta/lambda1);
    if(tm < 4)
    {
      mark = 0;
    }
    zs1 = -1*theta/tm;
    if(zs1 < -1 && lambda1 < -1)
    {
      mark = 0;
    }
  }

  if(lambda2 < 0)
  {
    tm = floor((-1)*theta*w/lambda2);
    if(tm < 4)
    {
      mark = 0;
    }
    zs1 = -1 * theta * w/tm;
    if(zs1 < -1 && lambda2 < -1)
    {
      mark = 0;
    }
  }

  if(mark == 1)
  {
  # Perform the likelihood ratio test

    zs3_1 = theta+(0:lx)*lambda1;
    zs3_2 = theta1 + (0:lx)*lambda1;
    zs4_1 = w*theta + (0:ly)*lambda2;
    zs4_2 = theta2 + (0:ly)*lambda2;

    if(sum(zs3_1 == 0) > 0 || sum(zs3_2 == 0) > 0 || sum(zs4_1 == 0) > 0|| sum(zs4_2 == 0) > 0)
    {
      mark = 0;
      LR_0 = 0;
      LR_1 = 0;
    }
    else
    {
      LR_0 = (n+m)*log(theta) - (n+w*m)*theta + m*log(w) + sum(obsx[1:(lx+1)]*(-1:(lx-1))*log(zs3_1[1:(lx+1)])) + sum(obsy[1:(ly+1)]*(-1:(ly-1))*log(zs4_1[1:(ly+1)])) ;
      
      LR_1 = n*log(theta1) + m*log(theta2) - n*theta1 - m*theta2 +sum(obsx[1:(lx+1)]*(-1:(lx-1))*log(zs3_2[1:(lx+1)])) + sum(obsy[1:(ly+1)]*(-1:(ly-1))*log(zs4_2[1:(ly+1)]));
    } 
    
    test = -2*(LR_0 - LR_1);
    return(list(mark=mark,Gptest=test));
  }
  else
  {
    return(list(mark=0,Gptest=0));
  }
}


likelihood_ratio_tissue_poisson <- function(x,lambda1,y,lambda2,w)
{
  n = length(x);
  m = length(y);

  sumx = sum(x);
  sumy = sum(y);

  lambda = (sumx+sumy)/(n+w*m);

  LR0 = sumx*log(lambda) + sumy*log(lambda) - (n+w*m)*lambda + sumy*log(w);
 
  LR1 = sumx*log(lambda1) + sumy*log(lambda2) - n*lambda1 - m*lambda2;

  test = -2*(LR0 - LR1);

  return(list(Ptest=test));
}

generalized_poisson_likelihood <- function(y)
{
  n = length(y);
  y_bar = mean(y);
 
  var_y = var(y);

  ly = max(y);

  if(var_y == 0 || ly < 2)
  {
    return(list(mark=0,theta=-1,lambda=-1,y_bar=y_bar,length=n));
  }
  else
  {
    l = 1-sqrt(y_bar/var_y);

    obs = rep(0,ly+1);
  
    for(i in 1:n)
    {
      obs[y[i]+1] = obs[y[i]+1] + 1;
    }
  }
  
  mark = 0;
  flag = 1;

# Running the Newton Raphson Algorithm to caculate the Maximum Likelihood Estimate
  for(its in 1:10000)
  {
    hl_i = rep(0,ly+1);
    il_i = rep(0,ly+1);
    tmp = rep(0,ly+1);

    tmp = y_bar + ((1:ly)-y_bar)*l;

    if(sum(tmp[2:ly]==0) == 0)
    {
      hl_i = (1:(ly-1))*(2:ly)*obs[3:(ly+1)]/tmp[2:ly];
      il_i = (1:(ly-1))*(2:ly)*obs[3:(ly+1)]*((2:ly)-y_bar)/(tmp[2:ly]^2);
    }
    else
    {
      flag = 0;
      break;
    }
    hl = sum(hl_i);
    il = sum(il_i);

    hl = hl-n*y_bar;
    if(flag == 1 && il != 0)
    {
      delta_l = hl/il;
      z1 = max(abs(l),0.01);
      z2 = abs(delta_l)/z1;
      l = l + delta_l;
      if(z2 <= 1.0e-6 && l < 1.0)
      {
        mark = 1;
        break;
      }
    }
    else
    {
      break;
    } 
  }

  theta = y_bar*(1-l);
  if(l < 0)
  {
    tm = floor(-1*theta/l);
    if(tm < 4)
    {
      mark = 0;
    }
    z1 = -1*theta/tm;
    if(z1 < -1 && l < -1)
    {
      mark = 0;
    }
  }
  output = c(mark,theta,l,y_bar);
  return(list(mark=mark,theta=theta,lambda=l,y_bar=y_bar,length=n));
}

#perform_permutation_de <- function(x,y,num_permute)
#{
#  gp_permute = matrix(0,num_permute,2);
#
#  for(i in 1:num_permute)
#  {
#    sum = x+y;
#    xnew = sapply(sum,rbinom,n=1,prob=0.5);
#    ynew = sum-xnew;
#    par1 = generalized_poisson_likelihood(xnew);
#    par2 = generalized_poisson_likelihood(ynew);
#
#    if(par1$mark == 1 && par2$mark == 1)
#    {
#      out = likelihood_ratio_tissue_generalized_poisson(xnew,par1$l,par1$theta,ynew,par2$l,par2$theta,1);
#      gp_permute[i,1] = out$mark;
#      gp_permute[i,2] = out$Gptest;
#    }
#  }
#  
#  #Estimate the gamma parameters
#  subset = gp_permute[,1] == 1;
#  m = mean(gp_permute[subset,2]);
#  v = var(gp_permute[subset,2]);
#  scale = v/m;
#  shape = m/scale;
#
#  return(list(scale = scale,shape=shape));
#}

perform_permutation_de <- function(x,y,num_permute)
{
  sum = x+y;

  xnew = sapply(sum,rbinom,n=num_permute,p=0.5)
  ynew = -1*apply(xnew,1,'-',sum)
  ynew = t(ynew)
  parx = apply(xnew,1,generalized_poisson_likelihood)
  pary = apply(ynew,1,generalized_poisson_likelihood)
  markx = sapply(parx,function(x) x[[1]])
  thetax = sapply(parx,function(x) x[[2]])
  lambdax = sapply(parx,function(x) x[[3]])
  marky = sapply(pary,function(x) x[[1]])
  thetay = sapply(pary,function(x) x[[2]])
  lambday = sapply(pary,function(x) x[[3]])

  subset = markx==1 & marky==1
  xnew = xnew[subset,]
  ynew = ynew[subset,]
  out = mapply(likelihood_ratio_tissue_generalized_poisson,split(xnew,1:nrow(xnew)),lambdax[subset],thetax[subset],split(ynew,1:nrow(ynew)),lambday[subset],thetay[subset],1)

  subset2 = out[1,] == 1
  y = sapply(out[2,subset2],mean);
  m = mean(y);
  v = var(y);
  scale = v/m;
  shape = m/scale;

  return(list(scale = scale,shape = shape));
}

perform_permutation_ds <- function(x,y,exon_start,exon_end,n,num_permute)
{
  sum = x+y;
  xnew = sapply(sum,rbinom,n=num_permute,p=0.5)
  ynew = -1*apply(xnew,1,'-',sum)
  ynew = t(ynew)
 
  znew = xnew[,exon_start:exon_end];
  vnew = ynew[,exon_start:exon_end];
  
  add_end = n-length(znew[1,]);
  mat = matrix(0,num_permute,add_end);
  znew = cbind(znew,mat);
  vnew = cbind(vnew,mat);

  parx = try(apply(xnew,1,generalized_poisson_likelihood));
  pary = try(apply(ynew,1,generalized_poisson_likelihood));
  parz = try(apply(znew,1,generalized_poisson_likelihood));
  parv = try(apply(vnew,1,generalized_poisson_likelihood));
   
  markx = sapply(parx,function(x) x[[1]])
  thetax = sapply(parx,function(x) x[[2]])
  lambdax = sapply(parx,function(x) x[[3]])
  marky = sapply(pary,function(x) x[[1]])
  thetay = sapply(pary,function(x) x[[2]])
  lambday = sapply(pary,function(x) x[[3]])
  markz = sapply(parz,function(x) x[[1]])
  thetaz = sapply(parz,function(x) x[[2]])
  lambdaz = sapply(parz,function(x) x[[3]])
  markv = sapply(parv,function(x) x[[1]])
  thetav = sapply(parv,function(x) x[[2]])
  lambdav = sapply(parv,function(x) x[[3]])

  subset = markx == 1 & marky == 1 & markz == 1 & markv == 1
  xnew = xnew[subset,];
  ynew = ynew[subset,]
  vnew = vnew[subset,]
  znew = znew[subset,]

  out = mapply(likelihood_ratio_generalized_poisson_exon_gene,split(znew,1:nrow(znew)),thetaz[subset],lambdaz[subset],split(xnew,1:nrow(xnew)),thetax[subset],lambdax[subset],split(vnew,1:nrow(vnew)),thetav[subset],lambdav[subset],split(ynew,1:nrow(ynew)),thetay[subset],lambday[subset]);

  subset2 = out[1,] == 1;
  y = sapply(out[2,subset2],mean);
  m = mean(y);
  v = var(y);
  scale = v/m;
  shape = m/scale;

  return(list(scale = scale,shape = shape));
}
 
#perform_permutation_ds <- function(x,y,exon_start,exon_end,n,num_permute)
#{
#  gp_permute = matrix(0,num_permute,2);
#
#  for(i in 1:num_permute)
#  {
#    sum = x+y;
#    xnew = sapply(sum,rbinom,n=1,prob=0.5);
#    ynew = sum-xnew;
# 
#    znew = xnew[exon_start:exon_end];
#    vnew = ynew[exon_start:exon_end];
#    
#    add_end = n-length(znew);
#    znew = c(znew,rep(0,add_end));
#    vnew = c(vnew,rep(0,add_end));
#
#    parx = generalized_poisson_likelihood(xnew);
#    pary = generalized_poisson_likelihood(ynew);
#    parz = generalized_poisson_likelihood(znew);
#    parv = generalized_poisson_likelihood(vnew);
#
#    if(parx$mark == 1 && pary$mark == 1 & parz$mark == 1 && parv$mark == 1)
#    {
#      out = likelihood_ratio_generalized_poisson_exon_gene(znew,parz$theta,parz$l,xnew,parx$theta,parx$l,vnew,parv$theta,parv$l,ynew,pary$theta,pary$l);
#      gp_permute[i,1] = out$mark;
#      gp_permute[i,2] = out$Gptest;
#    }
#  }
#  #Estimate the gamma parameters
#  subset = gp_permute[,1] == 1;
#  m = mean(gp_permute[subset,2]);
#  v = var(gp_permute[subset,2]);
#  scale = v/m;
#  shape = m/scale;
#
#  return(list(scale = scale,shape=shape));
#}
