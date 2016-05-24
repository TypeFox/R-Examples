
dt.lap=function(x, df, ncp=0, log=FALSE, normalize=c('central', 'integral','none'),...)
{
    x0=(sqrt(4*df*df+4*x*x*df+x*x*ncp*ncp)+x*ncp)/2/sqrt(x*x+df)
    logC=df/2*log(df/2)-df*ncp*ncp/2/(x*x+df)-.5*log(pi/2)-lgamma(df/2)-(df+1)/2*log(x*x+df)
    logNCT=logC+df*log(x0)-.5*(x0-ncp*x/sqrt(x*x+df))^2+.5*log(2*pi*x0*x0/(x0*x0+df))+
            pnorm(0,x0,sqrt(1/(1+df/x0/x0)),lower.tail=FALSE,log.p=TRUE)

    normalization=match.arg(normalize)
    if (normalization=='integral'){
          warning('integral normalization not implemented yet; set normalize="central"')
         normalization='central'
    }
    logNormalizer=if(normalization=='central') lgamma(df/2+1/2)+df/2-.5*log(2*pi)-df/2*log(df/2)-
                                               pnorm(0,sqrt(df),1/sqrt(2),lower.tail=FALSE,log.p=TRUE)
                  else if (normalization=='none') 0 

    logAns=logNormalizer+logNCT
    if(log)return(logAns) else return(exp(logAns))

}

