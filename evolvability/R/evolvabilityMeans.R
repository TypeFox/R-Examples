#' @export
evolvabilityMeans = function(G, means=1)  {
  if(means[1]==1) means=rep(1, nrow(G))
  G = G/(means%*%t(means))
  e_mean = mean(eigen(G)$values)
  e_max = max(eigen(G)$values)
  e_min = min(eigen(G)$values)
  Heig = 1/mean(1/eigen(G)$values)
  Ieig = var(eigen(G)$values)/(mean(eigen(G)$values))^2
  Ieig2 =  var(eigen(G)$values^2)/(mean(eigen(G)$values^2))^2
  Iinveig = var(1/eigen(G)$values)/(mean(1/eigen(G)$values))^2
  k = nrow(G)
  c_mean = Heig*(1+(2*Iinveig)/(k+2))
  r_mean = sqrt(mean(eigen(G)$values^2))*(1-(Ieig2/(4*(k+2))))
  a_mean = (Heig/e_mean)*(1+2*(Ieig+Iinveig-1+(Heig/e_mean)+2*Ieig*Iinveig/(k+2))/(k+2))
  i_mean = 1-a_mean
  c(e_mean = e_mean, e_min = e_min, e_max = e_max, c_mean = c_mean, 
    r_mean = r_mean, a_mean = a_mean, i_mean = i_mean)
}

