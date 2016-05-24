evaluateLogLik <- function(v,theta_x,theta_y,alpha_x,alpha_y,J,Sx,Sy,Nx,Ny,KDA_X,KDA_Y)
{
  nx = v;
  ny = 1-v;
  I_X = alpha_x * nx * (J-1) / (1 - alpha_x * nx - alpha_y * ny);
  I_Y = alpha_y * ny * (J-1) / (1 - alpha_x * nx - alpha_y * ny);

  a = lgamma(J+1);
  b = rep(0,length(I_X));

  for(cnt in 1:length(I_X)) {
    b[cnt] = lgamma(I_X[cnt] + I_Y[cnt] + J) - lgamma(I_X[cnt] + I_Y[cnt]);
  }
  c = a - b;


  l = ((theta_x/2) - 1) * log(nx) + ((theta_y/2) - 1) * log(ny);

  s = calcSumKDA2(Sx,Nx,I_X,(theta_x/2),KDA_X);
  z = calcSumKDA2(Sy,Ny,I_Y,(theta_y/2),KDA_Y);
  
  result = c(c + l + s + z);
  output <- rep(-Inf,length(v));
  for(i in 1:length(v)) output[i] = result[i];
  return(output);
}