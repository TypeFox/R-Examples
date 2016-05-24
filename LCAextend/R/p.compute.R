# Modification version 1.1: ajout du terme de décalage pour les covariables
# comme argument supplémentaire optionnel
p.compute <-
function(alpha,decal=0)
{
    p <- NULL
    if((length(alpha)==1)&&(is.na(alpha))) p <- 1
    else
    {
        for(i in 1:length(alpha))
        {
            if(alpha[i]==-Inf) p[i] <- 0
            else if(alpha[i]==Inf) p[i] <- 1-ifelse(i==1,0,sum(p[1:(i-1)]))
            else 
            {
                cum.p <- ifelse(i==1,0,sum(p[1:(i-1)]))
                cum.alpha <- sum(alpha[1:i][!is.infinite(alpha[1:i])])+decal
                p[i] <- exp(cum.alpha)/(1+exp(cum.alpha))-cum.p
            }
        }
        p[length(alpha)+1] <- 1-sum(p)
    }
    p
}

