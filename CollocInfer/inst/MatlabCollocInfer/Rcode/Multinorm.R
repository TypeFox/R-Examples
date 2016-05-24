#### A general multivariate normal function; some of this 
# should be coded as C
make.multinorm <- function()
{

#########################################################################

multinorm = list()
multinorm$fn <- function(y,times,x,pars,more)
{
    #  computes negative log density for multivariate normal distribution
    F = more$fn(times,x,pars,more$f.more)    
    S = more$var.fn(times,x,pars,more$v.more)   
    d = rep(0, dim(F)[1])
    if( dim(S)[1] == 1) {
      #  if variance is constant, inverse of S is computed once
      Sinv = solve(S[1,,])     
      lognormconst = -0.5*log(det(Sinv)) 
    }                                                     
    for(i in 1:dim(F)[1]){
        if( dim(S)[1] > 1) { 
          Sinv = solve(S[i,,]) 
          lognormconst = -0.5*log(det(Sinv)) 
        }
        resi = y[i,]-F[i,]
        d[i] = 0.5*t(resi) %*% Sinv %*% resi + lognormconst
    }
    return(d)
}

#########################################################################

multinorm$dfdy <- function(y,times,x,pars,more)
{
     F = more$fn(times,x,pars,more$f.more)
     S = more$var.fn(times,x,pars,more$v.more)

     d = array(0,dim(y))

     if( dim(S)[1] == 1){ SS = solve(S[1,,]) }
     
     for(i in 1:dim(y)[1]){
        if( dim(S)[1] > 1){ SS = solve(S[i,,]) }
        d[i,] = SS%*%(y[i,]-F[i,])
     }
     return(d)
}

#########################################################################

multinorm$dfdx <- function(y,times,x,pars,more)
{
     F = more$fn(times,x,pars,more$f.more)
     dF = more$dfdx(times,x,pars,more$f.more)
     
     S = more$var.fn(times,x,pars,more$v.more)
     dS = more$var.dfdx(times,x,pars,more$v.more)

     d = array(0,dim(x))
     
     
     if( dim(S)[1] == 1){ 
        SS = solve(S[1,,]) 
        dSS = array(dS[1,,,],dim(dS)[2:4])
    }

     for(i in 1:dim(x)[1]){
        if( dim(S)[1] > 1){
            SS = solve(S[i,,])
            dSS = array(dS[i,,,],dim(dS)[2:4])
        }
        wdifs = SS%*%(y[i,]-F[i,])

        d[i,] = -t(t(array(dF[i,,],dim(dF)[2:3]))%*%wdifs) 
            for(j in 1:dim(x)[2]){

              d[i,j] = d[i,j] - 0.5*t(wdifs)%*%dSS[,,j]%*%wdifs +
                    0.5*sum(diag(SS%*%dSS[,,j])) 
            }
     }
     return(d)
}

#########################################################################

multinorm$dfdp <- function(y,times,x,pars,more)
{
     F = more$fn(times,x,pars,more$f.more)
     dF = more$dfdp(times,x,pars,more$f.more)
     
     S = more$var.fn(times,x,pars,more$v.more)
     dS = more$var.dfdp(times,x,pars,more$v.more)

     d = matrix(0,dim(x)[1],length(pars))

     if( dim(S)[1] == 1){ 
        SS = solve(S[1,,]) 
        dSS = array(dS[1,,,],dim(dS)[2:4])
    }

     for(i in 1:dim(x)[1]){
        if( dim(S)[1] > 1){
            SS = solve(S[i,,])
            dSS = array(dS[i,,,],dim(dS)[2:4])
        }
        wdifs = SS%*%(y[i,]-F[i,])
        d[i,] = -t(array(dF[i,,],dim(dF)[2:3]))%*%wdifs 
           for(j in 1:length(pars)){
                d[i,j] = d[i,j] - 0.5*t(wdifs)%*%dSS[,,j]%*%wdifs +
                    0.5*sum(diag(SS%*%dSS[,,j])) 
           }    
    }

    return(d)

}

#########################################################################

multinorm$d2fdx2 <- function(y,times,x,pars,more)
{
     F = more$fn(times,x,pars,more$f.more)
     dF = more$dfdx(times,x,pars,more$f.more)
     d2F = more$d2fdx2(times,x,pars,more$f.more)
     
     S = more$var.fn(times,x,pars,more$v.more)
     dS = more$var.dfdx(times,x,pars,more$v.more)
     d2S = more$var.d2fdx2(times,x,pars,more$v.more)

     d = array(0,c(dim(x)[1],dim(x)[2],dim(x)[2]))

     if( dim(S)[1] == 1){
        SS = solve(S[1,,])
        dSS = array(dS[1,,,],dim(dS)[2:4])
        d2SS = array(d2S[1,,,,],dim(d2S)[2:5])
     }

     for(i in 1:dim(x)[1]){
        if( dim(S)[1] > 1){
            SS = solve(S[i,,])
            dSS = array(dS[i,,,],dim(dS)[2:4])
            d2SS = array(d2S[i,,,,],dim(d2S)[2:5])
        }
        wdifs = SS%*%(y[i,]-F[i,])
        
        for(j in 1:dim(x)[2]){     
            for(k in j:dim(x)[2]){
                d[i,j,k] = -(d2F[i,,j,k]-dF[i,,j]%*%SS%*%dSS[,,k]-dF[i,,k]%*%SS%*%dSS[,,j])%*%wdifs +
                   0.5*t(wdifs)%*%(dSS[,,j]%*%SS%*%dSS[,,k]+dSS[,,k]%*%SS%*%dSS[,,j]-d2SS[,,j,k])%*%wdifs +
                   t(dF[i,,j])%*%SS%*%dF[i,,k] - 0.5*sum(diag(SS%*%(dSS[,,j]%*%SS%*%dSS[,,k] - d2SS[,,j,k])))
                d[i,k,j] = d[i,j,k]
            }     
        }
     }
     return(d)
}

#########################################################################

multinorm$d2fdy2 <- function(y,times,x,pars,more)
{
    S = more$var.fn(times,x,pars,more$v.more)
    
    d = array(0,c(dim(y),dim(y)[2]))
    
    if( dim(S)[1] == 1){ SS = solve(S[1,,]) }  
    
    for(i in 1:dim(x)[1]){ 
        if( dim(S)[1] > 1){ SS = solve(S[i,,]) }
        d[i,,] = 0.5*(t(SS) + SS) 
    }
    
    return(d)
    
}

#########################################################################

multinorm$d2fdxdy <- function(y,times,x,pars,more)
{
     F = more$fn(times,x,pars,more$f.more)
     dF = more$dfdx(times,x,pars,more$f.more)
 
     S = more$var.fn(times,x,pars,more$v.more)
     dS = more$var.dfdx(times,x,pars,more$v.more)
     d = array(0,c(dim(x)[1],dim(x)[2],dim(y)[2]))

     if( dim(S)[1] == 1){
        SS = solve(S[1,,])
        dSS = array(dS[1,,,],dim(dS)[2:4])
     }
     for(i in 1:dim(x)[1]){
        if( dim(S)[1] > 1){
            SS = solve(S[i,,])
            dSS = array(dS[i,,,],dim(dS)[2:5])
        }
        for(j in 1:dim(x)[2]){
            d[i,j,] = -SS%*%dSS[,,j]%*%SS%*%(y[i,]-F[i,]) - SS%*%dF[i,,j]        
        }     
     }

    return(d)
}

#########################################################################

multinorm$d2fdxdp <- function(y,times,x,pars,more)
{
     F = more$fn(times,x,pars,more$f.more)
     dFx = more$dfdx(times,x,pars,more$f.more)
     dFp = more$dfdp(times,x,pars,more$f.more)
     d2F = more$d2fdxdp(times,x,pars,more$f.more)
     
     S = more$var.fn(times,x,pars,more$v.more)
     dSx = more$var.dfdx(times,x,pars,more$v.more)
     dSp = more$var.dfdp(times,x,pars,more$v.more)
     d2S = more$var.d2fdxdp(times,x,pars,more$v.more)

     d = array(0,c(dim(x)[1],dim(x)[2],length(pars)))

     if( dim(S)[1] == 1){
       SS = solve(S[1,,])
       dSSx = array(dSx[1,,,],dim(dSx)[2:4])
       dSSp = array(dSp[1,,,],dim(dSp)[2:4])
       d2SS = array( d2S[1,,,,],dim(d2S)[2:5])
     }

     for(i in 1:dim(x)[1]){
        if( dim(S)[1] > 1){
            SS = solve(S[i,,])
            dSSx = array(dSx[i,,,],dim(dSx)[2:4])
            dSSp = array(dSp[i,,,],dim(dSp)[2:4])
            d2SS = array(d2S[i,,,,],dim(d2S)[2:5])
        }
        wdifs = SS%*%(y[i,]-F[i,])
        for(j in 1:dim(x)[2]){     
            for(k in 1:length(pars)){
                d[i,j,k] = -(d2F[i,,j,k]-dFx[i,,j]%*%SS%*%dSSp[,,k]-dFp[i,,k]%*%SS%*%dSSx[,,j])%*%wdifs +
                   0.5*t(wdifs)%*%(dSSx[,,j]%*%SS%*%dSSp[,,k]+dSSp[,,k]%*%SS%*%dSSx[,,j]-d2SS[,,j,k])%*%wdifs + 
                   t(dFx[i,,j])%*%SS%*%dFp[i,,k] - 0.5*sum(diag(SS%*%(dSSx[,,j]%*%SS%*%dSSp[,,k] - d2SS[,,j,k])))
            }     
        }
     }
     return(d)
}

#########################################################################

multinorm$d2fdydp <- function(y,times,x,pars,more)
{
     F = more$fn(times,x,pars,more$f.more)
     dF = more$dfdp(times,x,pars,more$f.more)
   
     S = more$var.fn(times,x,pars,more$v.more)
     dS = more$var.dfdp(times,x,pars,more$v.more)

     d = array(0,c(dim(y)[1],dim(y)[2],length(pars)))

     if( dim(S)[1] == 1){ 
        SS = solve(S[1,,]) 
        dSS = array(dS[1,,,],dim(dS)[2:4])
     }     
     for(i in 1:dim(y)[1]){
        if( dim(S)[1] > 1){ 
            SS = solve(S[i,,]) 
            dSS = array(dS[i,,,],dim(dS)[2:4])
        }
        for(j in 1:length(pars)){
            d[i,,j] = -SS%*%(dSS[,,j]%*%SS%*%(y[i,]-F[i,]) + dF[i,,j])
        }
     }
     return(d)
}

  return(multinorm)
    
}
