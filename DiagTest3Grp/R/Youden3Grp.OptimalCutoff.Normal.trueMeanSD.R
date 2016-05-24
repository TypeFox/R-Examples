Youden3Grp.OptimalCutoff.Normal.trueMeanSD <-
function(mu.minus,mu0,mu.plus,s.minus,s0,s.plus)
  {
    #######Under normal assumptions, give the optimal cutoff points t.minus, t.plus which maximize Youden index for 3 categorical groups
    
    ######If want to compute 2 groups, leave mu.plus, s.plus as NULL or missing and if the two groups have equal variance, then the cutoff
    #####      is located as the mid-point of the two means
    ######mu.minus,mu0,mu.plus=NULL,s.minus,s0,s.plus: if true para in the group normal distr are known, use them, otherwise use sample mean and sd

    Youden3Grp.Normal.trueMeanSD <- function(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus)
      {
        ###give Youden Index when the true parameters of the three groups distr. are known
        Se <-1-pnorm(t.plus, mean=mu.plus, sd=s.plus, lower.tail = TRUE)
        Sp <- pnorm(t.minus, mean=mu.minus, sd=s.minus, lower.tail = TRUE)
        Sm <- pnorm(t.plus, mean=mu0, sd=s0, lower.tail = TRUE)-pnorm(t.minus, mean=mu0, sd=s0, lower.tail = TRUE)
        youden <- 0.5*(Se+Sp+Sm-1)
        return(youden)
      }

    var.minus <- s.minus^2
    var0 <- s0^2
    
    if(s.minus==s0) t.minus <- 0.5*(mu.minus+mu0)
    else
      {
        t.minus <- (mu0*var.minus-mu.minus*var0)-s.minus*s0*sqrt((mu.minus-mu0)^2+(var.minus-var0)*log(var.minus/var0))
        t.minus <- t.minus/(var.minus-var0)
      }
   
    if(s0==s.plus) t.plus <- 0.5*(mu0+mu.plus)
    else
      {
        
        var.plus <- s.plus^2
        t.plus <- (mu.plus*var0-mu0*var.plus)-s0*s.plus*sqrt((mu0-mu.plus)^2+(var0-var.plus)*log(var0/var.plus))
        t.plus <- t.plus/(var0-var.plus)
      }
    youden <- Youden3Grp.Normal.trueMeanSD(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,t.minus,t.plus)
    
    return(data.frame(t.minus=t.minus,t.plus=t.plus,youden=youden))
    
  }

