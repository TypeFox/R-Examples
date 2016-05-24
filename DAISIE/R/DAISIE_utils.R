roundn = function(x, digits = 0)
{
    fac = 10^digits
    return(trunc(fac * x + 0.5)/fac)
}

countspecies = function(datalistelement)
{
    N = length(datalistelement$branching_times) - 1 + datalistelement$missing_species
}

counttype1 = function(datalistelement)
{
    N1 = 0
    if(length(datalistelement$type1or2) > 0)
    {
        N1 = (datalistelement$type1or2 == 1)
    }
}

countspeciestype1 = function(datalistelement)
{
    N1 = 0
    if(length(datalistelement$type1or2) > 0)
    {
        if(datalistelement$type1or2 == 1)
        {
           N1 = length(datalistelement$branching_times) - 1 + datalistelement$missing_species
        }
    }
}

countimmi = function(datalistelement)
{
    datalistelement$stac != 2
}

fconstr13 = function(x,pars1,x_E,age)
{
    lac = pars1[1]
    laa = pars1[5]
    ga = pars1[4]
    A = x - lac
    C = ga + laa + 2 * lac
    ff = (1 + A/C * (1 - exp(-C * age))) * exp(-A * age) - (1 - x_E)
    return(ff)
}

fconstr15 = function(x,pars1,x_E,x_I,age)
{
    lac = pars1[1]
    laa = pars1[5]
    A = x - lac
    B_c = -1/age * log(1 - x_I)
    ga = B_c - x - laa - lac
    C = ga + laa + 2 * lac
    ff = (1 + A/C * (1 - exp(-C * age))) * exp(-A * age) - (1 - x_E)
    return(ff)
}

calcMN = function(datalist,pars1)
{
    N = sum(unlist(lapply(datalist,countspecies)))
    if(is.null(datalist[[1]]$not_present))
    {
        M = datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
        if(!is.na(pars1[6]))
        {   
           if(is.na(pars1[11]))
           {
              M = datalist[[1]]$not_present_type1 + sum(unlist(lapply(datalist,counttype1)))
           } else {
              M = M - max(0,roundn(pars1[11] * M)) 
           }
           N = sum(unlist(lapply(datalist,countspeciestype1)))      
        }
    } else {
        M = datalist[[1]]$not_present + length(datalist) - 1
    }
    return(c(M,N))
}

DAISIE_eq = function(datalist,pars1,pars2)
{
    eqmodel = pars2[5]
    ddep = pars2[2]
    MN = calcMN(datalist,pars1)
    M = MN[1]
    N = MN[2]
    I = sum(unlist(lapply(datalist,countimmi)))
    rNM = N/M
    rIM = I/(M - I)
    rIN = I/(N - I)
    clado = pars1[1] * ((1 - N/pars1[3])^(ddep == 1 || ddep == 11)) * (exp(-N/pars1[3]))^(ddep == 2 || ddep == 21)    
    ana = pars1[5]
    # Equilibrium based on deterministic model in terms of N
    if(eqmodel == 1)
    {
        immi = pars1[4] * ((1 - N/pars1[3])^(ddep == 11)) * (exp(-N/pars1[3]))^(ddep == 21)
        ext = clado + immi * (1/rNM - 1)       
        pars1[2] = ext
    }
    # Equilibrium model based on deterministic model in terms of E and I
    if(eqmodel == 2) # Only eq for N
    {
        ext = pars1[2]
        immitot = 1/(1/rNM * 1/(ext - clado) - 1/(ana + clado + ext))
        immi = immitot / ((1 - N/pars1[3])^(ddep == 11) * (exp(-N/pars1[3]))^(ddep == 21))
        pars1[4] = immi
    }
    if(eqmodel == 3) # Only eq for E            
    {
        immi = pars1[4] * ((1 - N/pars1[3])^(ddep == 11)) * (exp(-N/pars1[3]))^(ddep == 21)
        ext = clado + (ana + 2 * clado) * rIN        
        pars1[2] = ext
    }
    if(eqmodel == 4) # Only eq for I            
    {
        ext = pars1[2]
        immitot = (ext + ana + clado) * rIM
        immi = immitot / ((1 - N/pars1[3])^(ddep == 11) * (exp(-N/pars1[3]))^(ddep == 21))
        pars1[4] = immi
    }
    if(eqmodel == 5) # Eq for E and I            
    {
        ext = clado + (ana + 2 * clado) * rIN        
        immitot = (ext + ana + clado) * rIM
        immi = immitot / ((1 - N/pars1[3])^(ddep == 11) * (exp(-N/pars1[3]))^(ddep == 21))
        pars1[2] = ext
        pars1[4] = immi
    }             
    if(eqmodel == 13) # Within x_E of equilibrium for E - diversity-dependence not implemented
    {        
        x_E = pars2[10]
        x_I = pars2[11]
        age = datalist[[1]]$island_age
        pars1[2] = uniroot(f = fconstr13,interval = c(pars1[1] + 1E-6, pars1[1] + 10),pars1 = pars1,x_E = x_E, age = age)$root
        ga_c = -1/age * log(1 - x_I) - pars1[1] - pars1[2] - pars1[5]
        if(pars1[4] < ga_c)
        {
            cat("The non-endemics do not satisfy the equilibrium criterion for these parameters.\n")
        } 
    }
    if(eqmodel == 15) # Within x_E and x_I of equilibrium for both E and I - diversity-dependence not implemented
    {        
        x_E = pars2[10]
        x_I = pars2[11]
        age = datalist[[1]]$island_age
        pars1[2] = uniroot(f = fconstr15,interval = c(pars1[1] + 1E-6, pars1[1] + 10),pars1 = pars1,x_E = x_E, x_I = x_I, age = age)$root 
        pars1[4] = -1/age * log(1 - x_I) - pars1[1] - pars1[2] - pars1[5]
    }                                                                               
    return(pars1)
}