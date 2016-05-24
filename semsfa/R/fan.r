fan<-function(lambda_fan,resp,Ey,ineffD){
    n=length(resp)
    sigma_fan=sqrt(mean((resp-Ey)^2)/(1-(2*lambda_fan^2)/(pi*(1+lambda_fan^2))))
    mu=sqrt(2/pi)*sigma_fan*lambda_fan/sqrt(1+lambda_fan^2)

    if(ineffD)
    { eps=resp-Ey-mu
    rr=+(n/2)*log(2/pi)-n*log(sigma_fan) + sum(pnorm(-eps*lambda_fan/sigma_fan,log.p=TRUE))
          -(1/(2*sigma_fan^2))*sum(eps^2)}
    else
    { eps=resp-Ey+mu
    rr=+(n/2)*log(2/pi)-n*log(sigma_fan) + sum(pnorm(eps*lambda_fan/sigma_fan,log.p=TRUE))
          -(1/(2*sigma_fan^2))*sum(eps^2)
    }
    rr
}