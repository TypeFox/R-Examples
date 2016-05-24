Sp.Sm.Se <- function(x,y,z,t.minus,t.plus)
  {
    ###This function calculates the empirical correct classification probabilities, once identified the optimal cut-points from VUS or Youden3Grp analyses.
    #Specificity=Pr(x<=t.minus),Sensitivity=Pr(z>=t.plus) and Sm=Pr(t.minus<=y<=t.plus)
    x <- na.exclude(x)
    y <- na.exclude(y)
    z <- na.exclude(z)
    
    Sp <- function(x,t.minus){sum(x<=t.minus)/length(x)}##specificity Pr(x<=t-|D-)
    Se <- function(z,t.plus){sum(z>=t.plus)/length(z)}##sensitivity Pr(z>=t+|D+)
    Sm <- function(y,t.minus,t.plus){sum(y>t.minus & y<t.plus)/length(y)}## Pr(t-<=y<=t+|D0)

    sp <- Sp(x,t.minus)
    se <- Se(z,t.plus)
    sm <- Sm(y,t.minus,t.plus)
    return(c(Sp=sp,Sm=sm,Se=se))
  }
