ldknn.odds <-
function(type, probability, reference){
    lr<-summary(glm(type==reference~-1+probability,family=binomial))
    odds<-exp(coefficients(lr)['probability','Estimate'])
    z.value<-coefficients(lr)['probability','z value']
    p.value<-coefficients(lr)['probability','Pr(>|z|)']
    return(list(odds=odds, 'z.value'=z.value, 'Pr(>|z|)'=p.value))
}

