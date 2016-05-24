Youden3Grp.Normal <-
function(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus)
  {
    ###This function calculates the resulting Youden Index when the true parameters of the three groups distr. are unknown and have to be estimated by MLE (sample mean/var)
         
    Se <-1-pnorm(t.plus, mean=mu.plus, sd=s.plus, lower.tail = TRUE)
    Sp <- pnorm(t.minus, mean=mu.minus, sd=s.minus, lower.tail = TRUE)
    Sm <- pnorm(t.plus, mean=mu0, sd=s0, lower.tail = TRUE)-pnorm(t.minus, mean=mu0, sd=s0, lower.tail = TRUE)
    
    youden <- 0.5*(Se+Sp+Sm-1)
    
    return(data.frame(t.minus=t.minus,t.plus=t.plus,Se=Se,Sp=Sp,Sm=Sm,youden=youden))
  }

