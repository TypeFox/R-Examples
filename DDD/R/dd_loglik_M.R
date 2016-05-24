lambdamu = function(n,pars,ddep)
{
    lnn = length(n)
    zeros = rep(0,lnn)
    ones = rep(1,lnn)
    la = pars[1]
    mu = pars[2]
    K = pars[3]
    r = pars[4]
    lavec = la * ones
    muvec = mu * ones
    n0 = (ddep == 2 | ddep == 4)
    if(ddep == 1)
    {
       Kprime = la / (la - mu) * K
       lavec = pmax(zeros,la * (1 - n/Kprime))
    } else if(ddep == 1.3)
    {
       lavec = pmax(zeros,la * (1 - n/K))
    } else if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
    {
       y = -(log(la/mu)/log(K+n0))^(ddep != 2.2)
       lavec = pmax(zeros,la * (n + n0)^y)
    } else if(ddep == 2.3)
    {
       y = K
       lavec = pmax(zeros,la * (n + n0)^y)
    } else if(ddep == 3)
    {
       lavec = la * ones
       muvec = mu + (la - mu) * n/K
    } else if(ddep == 4 | ddep == 4.1 | ddep == 4.2)
    {
       y = (log(la/mu)/log(K+n0))^(ddep != 4.2)
       muvec = mu * (n + n0)^y
    } else if(ddep == 5)
    { 
       lavec = pmax(zeros,la - 1/(r + 1)*(la - mu)/K * n)
       muvec = muvec = mu + r/(r + 1)*(la - mu)/K * n
    }
    return(list(lavec,muvec))
}

dd_loglik_M_aux = function(pars,lx,k,ddep)
{
    nvec = 0:(lx - 1);
    lambdamu_nk = lambdamu(nvec + k,pars,ddep)
    lambda_nk = lambdamu_nk[[1]]
    mu_nk = lambdamu_nk[[2]];
    i1 = 2:lx; j1 = 1:(lx - 1); x1 = lambda_nk[1:(lx - 1)] * (nvec[1:(lx - 1)] + 2 * k)
    i2 = 1:(lx - 1); j2 = 2:lx; x2 = mu_nk[2:lx] * nvec[2:lx];
    i3 = 1:lx; j3 = i3; x3 = -(lambda_nk + mu_nk) * (nvec + k)
    MM = Matrix::sparseMatrix(i = c(i1,i2,i3),j = c(j1,j2,j3),x = c(x1,x2,x3), dims = c(lx,lx));
    MM[lx,lx] = -mu_nk[lx] * (nvec[lx] + k)
    return(MM)
}

changepars = function(pars)
{
    if(length(pars) <= 3)
    {
       sel = 1:2
    } else {
       sel = c(1:2,4:5)
    }
    if(sum(pars[sel] == Inf) > 0)
    {
       pars[which(pars[sel] == Inf)] = 1E+10
    }
    if(sum(pars[sel] == 0) > 0)
    {
       pars[which(pars[sel] == 0)] = 1E-14
    }
    if(sum(pars[sel] == 1) > 0)
    {
       pars[which(pars[sel] == 1)] = pars[which(pars[sel] == 1)] - 1E-14
    }
    return(pars)
}

dd_loglik_M = function(pars,lx,k,ddep,tt,p)
{
    pars = changepars(pars)
    MM = dd_loglik_M_aux(pars,lx,k,ddep)
    p = expoRkit::expv(x = MM,v = p,t = tt,m = 50L)
    return(p)
}

dd_loglik_M_bw_aux = function(pars,lx,k,ddep)
{
    nvec = 0:(lx - 1);
    lambdamu_nk = lambdamu(nvec + k,pars,ddep)
    lambda_nk = lambdamu_nk[[1]]
    mu_nk = lambdamu_nk[[2]];
    i1 = 2:lx; j1 = 1:(lx - 1); x1 = mu_nk[1:(lx - 1) + 1] * nvec[1:(lx - 1) + 1];
    i2 = 1:(lx - 1); j2 = 2:lx; x2 = lambda_nk[2:lx - 1] * (nvec[2:lx] + 2 * k - 1);
    i3 = 1:lx; j3 = i3; x3 = -(lambda_nk + mu_nk) * (nvec + k)
    MM = Matrix::sparseMatrix(i = c(i1,i2,i3),j = c(j1,j2,j3),x = c(x1,x2,x3), dims = c(lx,lx));
    MM[lx,lx] = -mu_nk[lx] * (nvec[lx] + k)
    return(MM)
}

dd_loglik_M_bw = function(pars,lx,k,ddep,tt,p)
{
    pars = changepars(pars)
    MM = dd_loglik_M_bw_aux(pars,lx,k,ddep)
    p = expoRkit::expv(x = MM,v = p,t = tt,m = 50L)
    return(p)
}

lambdamu2 = function(n,pars,ddep)
{
    lnn = length(n)
    laM = pars[1]
    muM = pars[2]
    KM = pars[3]
    n0 = (ddep == 2 | ddep == 4)
    nx1 = rep(0:(lnn - 1),lnn)
    dim(nx1) = c(lnn,lnn) # row index = number of species in first group 
    nx2 = t(nx1) # column index = number of species in second group
    nxt = nx1 + nx2
    if(ddep == 1) 
    { 
        lavec = pmax(matrix(0,lnn,lnn),laM - (laM-muM)/KM * nxt)
        muvec = muM * matrix(1,lnn,lnn)
    } else if(ddep == 1.3) 
    { 
        lavec = pmax(matrix(0,lnn,lnn),laM * (1 - nxt/KM))
        muvec = muM * matrix(1,lnn,lnn)
    } else if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
    { 
        x = -(log(laM/muM)/log(KM+n0))^(ddep != 2.2)
        lavec = pmax(matrix(0,lnn,lnn),laM * (nxt + n0)^x)
        muvec = muM * matrix(1,lnn,lnn)
    } else if(ddep == 2.3)
    { 
        x = KM
        lavec = pmax(matrix(0,lnn,lnn),laM * (nxt + n0)^x)
        muvec = muM * matrix(1,lnn,lnn)
    } else if(ddep == 3)
    {
        lavec = laM * matrix(1,lnn,lnn)
        muvec = muM + (laM - muM)/KM * nxt
    } else if(ddep == 4 | ddep == 4.1 | ddep == 4.2)
    {
        lavec = laM * matrix(1,lnn,lnn)
        x = (log(laM/muM)/log(KM+n0))^(ddep != 4.2)
        muvec = (nxt + n0)^x
    }
    return(list(lavec,muvec))
}

dd_loglik_M2_aux = function(pars,lx,ddep)
{
    nvec = 0:(lx - 1);
    lambdamu_n = lambdamu2(nvec,pars,ddep)
    lambda_n = lambdamu_n[[1]]
    mu_n = lambdamu_n[[2]];
    ly = lx^2
    # set elements
    auxM1 = rep(0:(lx - 1),times = lx) + rep(0:(lx - 1),each = lx)
    auxM2 = auxM1; auxM2[seq(lx,ly,by = lx)] = 0;
    auxM3 = rep(1,ly); auxM3[seq(lx,ly,by = lx)] = 0;
    # subdiagonal
    i1 = 2:ly; j1 = 1:(ly - 1); x1 = lambda_n[1 + auxM1[1:(ly - 1)]] * (auxM2[1:(ly - 1)] > 0) * rep(0:(lx - 1),lx)[-ly];
    # sublxdiagonal
    i2 = (lx + 1):ly; j2 = 1:(ly - lx); x2 = lambda_n[1 + auxM1[1:(ly - lx)]] * (auxM2[1:(ly - lx)] > 0) * rep(0:(lx - 2),each = lx);
    # superdiagonal
    i3 = 1:(ly - 1); j3 = 2:ly; x3 = mu_n[1 + auxM1[2:ly]] * rep(0:(lx - 1), times = lx)[-1]; 
    # superlxdiagonal
    i4 = 1:(ly - lx); j4 = (lx + 1):ly; x4 = (mu_n[1 + auxM1[1:ly]] * rep(0:(lx - 1),each = lx))[-c(1:lx)];
    i5 = 1:ly; j5 = i5; x5 = -(lambda_n[1 + auxM1[1:ly]] + mu_n[1 + auxM1[1:ly]]) * auxM1[1:ly];
    MM = Matrix::sparseMatrix(i = c(i1,i2,i3,i4,i5),j = c(j1,j2,j3,j4,j5),x = c(x1,x2,x3,x4,x5));
    return(MM)
}

dd_loglik_M2 = function(pars,lx,ddep,tt,p)
{
    pars = changepars(pars)
    MM = dd_loglik_M2_aux(pars,lx,ddep)
    p = expoRkit::expv(x = MM,v = p,t = tt,m = 50L)
    return(p)
}

lambdamu3 = function(n,pars,ddep,kM,kS)
{
    lnn = length(n) + kM + kS
    laM = pars[1]
    muM = pars[2]
    K = pars[3]
    laS = pars[4]
    muS = pars[5]
    n0 = (ddep == 2 | ddep == 4)
    nx1 = rep(0:(lnn - 1),lnn)
    dim(nx1) = c(lnn,lnn) # row index = number of species in first group 
    nx2 = t(nx1) # column index = number of species in second group
    nxt = nx1 + nx2 + kM + kS
    muvecM = muM * matrix(1,lnn,lnn)
    muvecS = muS * matrix(1,lnn,lnn)    
    if(ddep == 1.3) 
    { 
        lavecM = pmax(matrix(0,lnn,lnn),laM * (1 - nxt/K))
        lavecS = pmax(matrix(0,lnn,lnn),laS * (1 - nxt/K))
    } else if(ddep == 2.3)
    { 
        x = K
        lavecM = pmax(matrix(0,lnn,lnn),laM * (nxt + n0)^x)
        lavecS = pmax(matrix(0,lnn,lnn),laS * (nxt + n0)^x)
    }     
    return(list(lavecM,muvecM,lavecS,muvecS))
}

dd_loglik_M3_aux = function(pars,lx,ddep,kM,kS)
{
    nvec = 0:(lx - 1);
    lambdamu_n = lambdamu3(nvec,pars,ddep,kM,kS)
    lambdaM_n = lambdamu_n[[1]]
    muM_n = lambdamu_n[[2]];
    lambdaS_n = lambdamu_n[[3]]
    muS_n = lambdamu_n[[4]];
    ly = lx^2
    # set elements
    auxM1 = rep(0:(lx - 1),times = lx) + rep(0:(lx - 1),each = lx)
    auxM2 = auxM1; auxM2[seq(lx,ly,by = lx)] = 0;
    auxM3 = rep(1,ly); auxM3[seq(lx,ly,by = lx)] = 0;
    # subdiagonal
    i1 = 2:ly; j1 = 1:(ly - 1); x1 = lambdaS_n[1 + auxM1[1:(ly - 1)]] * (kS + auxM2[1:(ly - 1)] > 0) * (rep(0:(lx - 1),lx)[-ly] + 2 * kS);
    # sublxdiagonal
    i2 = (lx + 1):ly; j2 = 1:(ly - lx); x2 = lambdaM_n[1 + auxM1[1:(ly - lx)]] * (kM + auxM2[1:(ly - lx)] > 0) * (rep(0:(lx - 2),each = lx) + 2 * kM);
    # superdiagonal
    i3 = 1:(ly - 1); j3 = 2:ly; x3 = muS_n[1 + auxM1[2:ly]] * rep(0:(lx - 1), times = lx)[-1];
    # superlxdiagonal
    i4 = 1:(ly - lx); j4 = (lx + 1):ly; x4 = (muM_n[1 + auxM1[1:ly]] * rep(0:(lx - 1),each = lx))[-c(1:lx)];
    # diagonal
    i5 = 1:ly; j5 = i5; x5 = -(lambdaM_n[1 + auxM1[1:ly]] + muM_n[1 + auxM1[1:ly]]) * (rep(0:(lx - 1),each = lx) + kM) - (lambdaS_n[1 + auxM1[1:ly]] + muS_n[1 + auxM1[1:ly]]) * (rep(0:(lx - 1),times = lx) + kS)
    
    MM = Matrix::sparseMatrix(i = c(i1,i2,i3,i4,i5),j = c(j1,j2,j3,j4,j5),x = c(x1,x2,x3,x4,x5));
    return(MM)
}

dd_loglik_M3 = function(pars,lx,ddep,tt,p,kM,kS)
{
    pars = changepars(pars)
    MM = dd_loglik_M3_aux(pars,lx,ddep,kM,kS)
    p = expoRkit::expv(x = MM,v = p,t = tt,m = 50L)
    return(p)
}
