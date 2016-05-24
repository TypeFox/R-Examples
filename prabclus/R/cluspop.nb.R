"cluspop.nb" <-
function(neighbors,p.nb=0.5,n.species,clus.specs,reg.group,
                       grouppf=10, n.regions=length(neighbors),
                       vector.species=rep(1,n.species),
                       pdf.regions=rep(1/n.regions,n.regions),count=TRUE,
                       pdfnb=FALSE){
# print(vector.species)
  out <- matrix(0,ncol=n.species,nrow=n.regions)
  if (pdfnb)
  {
    for (i in 1:n.regions)
      pdf.regions[i] <- pdf.regions[i]/max(1,length(neighbors[[i]]))
    pdf.regions <- pdf.regions/sum(pdf.regions)
  }
  pdf.group <- pdf.complement <- pdf.regions
  pdf.groupstart <- pdf.cstart <- rep(0,n.regions)
  sp <- spc <- spr <- 0
  prob.group <- sum(pdf.group[reg.group])
  pdf.complement[reg.group] <- pdf.regions[reg.group]/grouppf
  spc <- sum(pdf.complement[reg.group])
  sp <- spc*grouppf
  reg.c <- (1:n.regions)[-reg.group]
  pdf.group[reg.c] <- pdf.regions[reg.c]/grouppf
  pdf.cstart[reg.c] <- pdf.regions[reg.c]/(1-sp)
  pdf.complement[reg.c] <- pdf.cstart[reg.c]*(1-spc)
  spr <- sum(pdf.group[reg.c])
  pdf.groupstart[reg.group] <- pdf.regions[reg.group]/sp
  pdf.group[reg.group] <- pdf.groupstart[reg.group]*(1-spr)
  cdf.local <- cdf.regions <- cdf.groupstart <- cdf.cstart <- c()
  for (i in 1:n.regions){
    cdf.regions[i] <- sum(pdf.regions[1:i])
    cdf.groupstart[i] <- sum(pdf.groupstart[1:i])
    cdf.cstart[i] <- sum(pdf.cstart[1:i])
  }
# print(pdf.groupstart)    
# print(cdf.groupstart)    
# regular species
  for (i in 1:(n.species-clus.specs))
  {
    if(count)
      cat("Species ",i,"\n")
    spec.regind <- spec.neighb <- rep(FALSE,n.regions)
    nsize <- vector.species[1+floor(length(vector.species)*runif(1))]
# print(nsize)
    r1 <- runif(1)
    reg <- 1+sum(r1>cdf.regions)
# print(reg)
    spec.regind[reg] <- TRUE
    for (k in neighbors[[reg]])
      spec.neighb[k] <- TRUE
    out[reg,i] <- 1
    if(nsize>1)
      for (j in 2:nsize)
        if ((sum(spec.neighb)==0) | (sum(pdf.regions[spec.neighb])<1e-8) |
             (sum(spec.neighb | spec.regind)==n.regions))
# no further neighbors or only neighbors, i.e., next region is drawn from all
# remaining
        {
          nreg <- sum(!spec.regind)
          pdf.local <- pdf.regions[!spec.regind]
          pdf.local <- pdf.local/sum(pdf.local)
          for (l in 1:nreg)
            cdf.local[l] <- sum(pdf.local[1:l])
# cat(nreg, "\n")
          r1 <- runif(1)
          zz <- 1+sum(r1>cdf.local[1:nreg])
# cat(zz,"\n")
          reg <- (1:n.regions)[!spec.regind][zz]
# cat("reg, all ",reg,"\n")  
          spec.regind[reg] <- TRUE
          spec.neighb[reg] <- FALSE
          for (k in neighbors[[reg]])
            spec.neighb[k] <- !(spec.regind[k])
          out[reg,i] <- 1
        }
        else
          if (runif(1)<p.nb)
# next region is drawn from non-neighbors (jump)
          { 
            regs <- !(spec.regind | spec.neighb)
            nreg <- sum(regs)
            pdf.local <- pdf.regions[regs]
            pdf.local <- pdf.local/sum(pdf.local)
            for (l in 1:nreg)
              cdf.local[l] <- sum(pdf.local[1:l])
            r1 <- runif(1)
            zz <- 1+sum(r1>cdf.local[1:nreg])
# cat(nreg," ",zz,"\n")
            reg <- (1:n.regions)[regs][zz]  
# cat("reg, jump ",reg,"\n")  
            spec.regind[reg] <- TRUE
            for (k in neighbors[[reg]])
              spec.neighb[k] <- !(spec.regind[k])
            out[reg,i] <- 1
# if (sum(out[,i])!=sum(spec.regind))
#   cat("error: sum= ",sum(out[,i])," ind=",sum(spec.regind),"\n")
          }
          else
          {
# next region is drawn from neighbors
            nreg <- sum(spec.neighb)
            pdf.local <- pdf.regions[spec.neighb]
# print(pdf.local)
            pdf.local <- pdf.local/sum(pdf.local)
            for (l in 1:nreg)
              cdf.local[l] <- sum(pdf.local[1:l])
# print(cdf.local)
            r1 <- runif(1)
            zz <- 1+sum(r1>cdf.local[1:nreg])
# cat("nreg= ",nreg," zz =",zz,"\n")
            reg <- (1:n.regions)[spec.neighb][zz]  
# cat("reg, neighbor ",reg,"\n")  
            spec.regind[reg] <- TRUE
            spec.neighb[reg] <- FALSE
            for (k in neighbors[[reg]])
              spec.neighb[k] <- !(spec.regind[k])
            out[reg,i] <- 1
# if (sum(out[,i])!=sum(spec.regind))
#   cat("error: sum= ",sum(out[,i])," ind=",sum(spec.regind),"\n")
          }
# end if nsize>1 for j 
# cat("out=",sum(out[,i]),"  ind=",sum(spec.regind),"  nb=",sum(spec.neighb),
#    "  nni=",sum(!(spec.regind | spec.neighb)),"\n")
  } # for i - regular species
# species from reg.group
  for (i in 1:clus.specs)
  {
    ind <-i+n.species-clus.specs 
    if(count)
      cat("Clustered species ",ind,"\n")
    groupind <- runif(1)<prob.group
    if (groupind){
      spec.regind <- spec.neighb <- rep(FALSE,n.regions)
      nsize <- vector.species[1+floor(length(vector.species)*runif(1))]
  # print(nsize)
      r1 <- runif(1)
      reg <- 1+sum(r1>cdf.groupstart)
  # print(reg)
      spec.regind[reg] <- TRUE
      for (k in neighbors[[reg]])
        spec.neighb[k] <- TRUE
      out[reg,ind] <- 1
      if(nsize>1)
        for (j in 2:nsize)
          if ((sum(spec.neighb)==0) | (sum(pdf.group[spec.neighb])<1e-8) |
               (sum(spec.neighb | spec.regind)==n.regions))
  # no further neighbors or only neighbors, i.e., next region is drawn from all
  # remaining
          {
            nreg <- sum(!spec.regind)
            pdf.local <- pdf.group[!spec.regind]
            pdf.local <- pdf.local/sum(pdf.local)
            for (l in 1:nreg)
              cdf.local[l] <- sum(pdf.local[1:l])
  # cat(nreg, "\n")
            r1 <- runif(1)
            zz <- 1+sum(r1>cdf.local[1:nreg])
  # cat(zz,"\n")
            reg <- (1:n.regions)[!spec.regind][zz]
  # cat("reg, all ",reg,"\n")  
            spec.regind[reg] <- TRUE
            spec.neighb[reg] <- FALSE
            for (k in neighbors[[reg]])
              spec.neighb[k] <- !(spec.regind[k])
            out[reg,ind] <- 1
          }
          else
            if (runif(1)<p.nb)
  # next region is drawn from non-neighbors (jump)
            { 
              regs <- !(spec.regind | spec.neighb)
              nreg <- sum(regs)
              pdf.local <- pdf.group[regs]
              pdf.local <- pdf.local/sum(pdf.local)
              for (l in 1:nreg)
                cdf.local[l] <- sum(pdf.local[1:l])
              r1 <- runif(1)
              zz <- 1+sum(r1>cdf.local[1:nreg])
  # cat(nreg," ",zz,"\n")
              reg <- (1:n.regions)[regs][zz]  
  # cat("reg, jump ",reg,"\n")  
              spec.regind[reg] <- TRUE
              for (k in neighbors[[reg]])
                spec.neighb[k] <- !(spec.regind[k])
              out[reg,ind] <- 1
  # if (sum(out[,i])!=sum(spec.regind))
  #   cat("error: sum= ",sum(out[,i])," ind=",sum(spec.regind),"\n")
            }
            else
            {
  # next region is drawn from neighbors
              nreg <- sum(spec.neighb)
              pdf.local <- pdf.group[spec.neighb]
  # print(pdf.local)
              pdf.local <- pdf.local/sum(pdf.local)
              for (l in 1:nreg)
                cdf.local[l] <- sum(pdf.local[1:l])
  # print(cdf.local)
              r1 <- runif(1)
              zz <- 1+sum(r1>cdf.local[1:nreg])
  # cat("nreg= ",nreg," zz =",zz,"\n")
              reg <- (1:n.regions)[spec.neighb][zz]  
  # cat("reg, neighbor ",reg,"\n")  
              spec.regind[reg] <- TRUE
              spec.neighb[reg] <- FALSE
              for (k in neighbors[[reg]])
                spec.neighb[k] <- !(spec.regind[k])
              out[reg,ind] <- 1
  # if (sum(out[,i])!=sum(spec.regind))
  #   cat("error: sum= ",sum(out[,i])," ind=",sum(spec.regind),"\n")
            }
  # end if nsize>1 for j 
  # cat("out=",sum(out[,i]),"  ind=",sum(spec.regind),"  nb=",sum(spec.neighb),
  #    "  nni=",sum(!(spec.regind | spec.neighb)),"\n")
    } # if groupind           
# species from complement
    else{
      spec.regind <- spec.neighb <- rep(FALSE,n.regions)
      nsize <- vector.species[1+floor(length(vector.species)*runif(1))]
  # print(nsize)
      r1 <- runif(1)
      reg <- 1+sum(r1>cdf.cstart)
  # print(reg)
      spec.regind[reg] <- TRUE
      for (k in neighbors[[reg]])
        spec.neighb[k] <- TRUE
      out[reg,ind] <- 1
      if(nsize>1)
        for (j in 2:nsize)
          if ((sum(spec.neighb)==0) | (sum(pdf.complement[spec.neighb])<1e-8) |
               (sum(spec.neighb | spec.regind)==n.regions))
  # no further neighbors or only neighbors, i.e., next region is drawn from all
  # remaining
          {
            nreg <- sum(!spec.regind)
            pdf.local <- pdf.complement[!spec.regind]
            pdf.local <- pdf.local/sum(pdf.local)
            for (l in 1:nreg)
              cdf.local[l] <- sum(pdf.local[1:l])
  # cat(nreg, "\n")
            r1 <- runif(1)
            zz <- 1+sum(r1>cdf.local[1:nreg])
  # cat(zz,"\n")
            reg <- (1:n.regions)[!spec.regind][zz]
  # cat("reg, all ",reg,"\n")  
            spec.regind[reg] <- TRUE
            spec.neighb[reg] <- FALSE
            for (k in neighbors[[reg]])
              spec.neighb[k] <- !(spec.regind[k])
            out[reg,ind] <- 1
          }
          else
            if (runif(1)<p.nb)
  # next region is drawn from non-neighbors (jump)
            { 
              regs <- !(spec.regind | spec.neighb)
              nreg <- sum(regs)
              pdf.local <- pdf.complement[regs]
              pdf.local <- pdf.local/sum(pdf.local)
              for (l in 1:nreg)
                cdf.local[l] <- sum(pdf.local[1:l])
              r1 <- runif(1)
              zz <- 1+sum(r1>cdf.local[1:nreg])
  # cat(nreg," ",zz,"\n")
              reg <- (1:n.regions)[regs][zz]  
  # cat("reg, jump ",reg,"\n")  
              spec.regind[reg] <- TRUE
              for (k in neighbors[[reg]])
                spec.neighb[k] <- !(spec.regind[k])
              out[reg,ind] <- 1
  # if (sum(out[,i])!=sum(spec.regind))
  #   cat("error: sum= ",sum(out[,i])," ind=",sum(spec.regind),"\n")
            }
            else
            {
  # next region is drawn from neighbors
              nreg <- sum(spec.neighb)
              pdf.local <- pdf.complement[spec.neighb]
  # print(pdf.local)
              pdf.local <- pdf.local/sum(pdf.local)
              for (l in 1:nreg)
                cdf.local[l] <- sum(pdf.local[1:l])
  # print(cdf.local)
              r1 <- runif(1)
              zz <- 1+sum(r1>cdf.local[1:nreg])
  # cat("nreg= ",nreg," zz =",zz,"\n")
              reg <- (1:n.regions)[spec.neighb][zz]  
  # cat("reg, neighbor ",reg,"\n")  
              spec.regind[reg] <- TRUE
              spec.neighb[reg] <- FALSE
              for (k in neighbors[[reg]])
                spec.neighb[k] <- !(spec.regind[k])
              out[reg,ind] <- 1
  # if (sum(out[,i])!=sum(spec.regind))
  #   cat("error: sum= ",sum(out[,i])," ind=",sum(spec.regind),"\n")
            }
  # end if nsize>1 for j 
  # cat("out=",sum(out[,i]),"  ind=",sum(spec.regind),"  nb=",sum(spec.neighb),
  #    "  nni=",sum(!(spec.regind | spec.neighb)),"\n")
    } # else (complement)
  } # for i
  out
}
