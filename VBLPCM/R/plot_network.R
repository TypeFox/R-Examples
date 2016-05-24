vblpcmdrawpie <- function(center,radius,probs,n=50,colours=1:length(probs))
{
  x <- c(0,cumsum(probs)/sum(probs))
  dx <- diff(x)
  np <- length(probs)
  for (i in 1:np)
  {
    t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
    xc <- center[1] + c(cos(t2p), 0) * radius
    yc <- center[2] + c(sin(t2p), 0) * radius
    polygon(xc, yc, border = FALSE, col = colours[i])
  }
}  

plot.vblpcm<-function(x, ..., R2=0.2, main="Variational-Bayes Positions", alpha=0.5, colours=1:x$G, RET=FALSE)
  {
  model<-x$model
  d<-x$d
  N<-x$N
  P_n<-x$P_n
  G<-x$G
  Y<-x$Y
  V_xi_n<-x$V_xi_n
  V_z<-x$V_z
  V_eta<-x$V_eta
  V_lambda<-x$V_lambda
  V_alpha<-x$V_alpha
  inv_sigma02<-x$inv_sigma02
  omega2<-x$omega2
  mu_nought<-0
  if (d<=2) 
    {
    V_z_2D<-V_z
    V_eta_2D<-V_eta
    }
  if (d>2) 
    {
    cat("Dimension of latent space is greater than 2: using cmdscale of positions for 2D plot\n")
    V_z_2D<-cmdscale(dist(V_z))
    V_eta_2D<-matrix(NaN,G,2)
    for (g in 1:G)
      {
       for (dd in 1:2)
        {
        tmpsum1 = 0
        tmpsum2 = 0
        for (i in 1:N)
          {
          tmpsum1 = tmpsum1 + 0.5*V_lambda[g,i]*inv_sigma02*V_alpha[g]*V_z_2D[i,dd]
          tmpsum2 = tmpsum2 + V_lambda[g,i]*0.5*inv_sigma02*V_alpha[g]
          }
        V_eta_2D[g,dd] = (tmpsum1 + 0.5*mu_nought/ omega2)/(tmpsum2+0.5/ omega2)
        }
      }
    }
  XLIM=c(min(V_z_2D[,1]),max(V_z_2D[,1]))
  YLIM=c(min(V_z_2D[,2]),max(V_z_2D[,2]))
  XLIM=XLIM*1.1
  YLIM=YLIM*1.1
  pad<-formals(plot.network.default)$pad
  #if (model=="plain") 
  object.scale = formals(plot.network.default)$object.scale
  piesize <- R2 
  if (model=="rreceiver" | model=="rsender")
    {
    #piesize <- R2*(exp(scale(V_xi_n)))
    tmp<-scale(V_xi_n)
    piesize <- R2+R2*tmp
    }
  if (model=="rsocial") 
    {
    #piesize <- R2*(exp(scale(V_xi_n[seq(1,2*N,2)]+V_xi_n[seq(2,2*N,2)])))
    #piesize <- R2*(exp(scale(V_xi_n[,1]+V_xi_n[,2])))
    tmp<-scale(V_xi_n[,1]+V_xi_n[,2])
    piesize <- R2+R2*tmp
    }
  pie.order <- order(piesize, decreasing = TRUE)
  vertex.cex=0 
  if (model!="plain") vertex.cex=piesize*3.5 else vertex.cex=(piesize*diff(XLIM))/10
  plot.network(x$net, coord = V_z_2D, main = main,
               xlim = XLIM, ylim = YLIM, 
  	       suppress.axes = 0, edge.col = rgb(t(rep(190/255,3)),alpha=alpha), vertex.cex=vertex.cex, vertex.col=colours, ...)
  # add unknown links to plot
  for (i in 1:N)
    for (j in 1:N)
      if (is.na(Y[i,j])) #
        points(V_z_2D[c(i,j),1], V_z_2D[c(i,j),2], col=8, lty=2, t='l')
  plot_matrix<-cbind(V_z_2D,vertex.cex/2,t(V_lambda))
  plot_func<-function(x) vblpcmdrawpie(center=x[1:2], radius=x[2+1], probs=x[(2+2):(2+1+G)], n = 20, colours=colours)
  apply(plot_matrix, 1, plot_func)
  points(V_eta_2D, col=colours, pch=3)
  symbols(V_eta_2D, circles=2*sqrt(1/(x$inv_sigma02*x$V_alpha)), add=1, fg=colours, inches=FALSE)
  #symbols(V_eta_2D, circles=sqrt(x$V_lambda%*%x$V_sigma2*x$sigma02), add=1, fg=colours, inches=FALSE)
  #points(0,0,pch=3)
  if (RET) 
    return(list(positions=V_z_2D,clusters=V_eta_2D))
  }

