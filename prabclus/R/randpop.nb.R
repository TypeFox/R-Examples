"randpop.nb" <-
function(neighbors,p.nb=0.5,n.species,
                       n.regions=length(neighbors),
                       vector.species=rep(1,n.species), species.fixed=FALSE,
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
  cdf.local <- cdf.regions <- c()
  for (i in 1:n.regions)
    cdf.regions[i] <- sum(pdf.regions[1:i])
  for (i in 1:n.species)
  {
    if(count)
      cat("Species ",i,"\n")
    spec.regind <- spec.neighb <- rep(FALSE,n.regions)
    nsize <- ifelse(species.fixed, vector.species[i],
                    vector.species[1+floor(length(vector.species)*runif(1))])
#    print(vector.species)
#    print(nsize)
    r1 <- runif(1)
    reg <- 1+sum(r1>cdf.regions)
# print(reg)
#    print(cdf.regions)
    spec.regind[reg] <- TRUE
    for (k in neighbors[[reg]])
      spec.neighb[k] <- TRUE
    out[reg,i] <- 1
#    print(out)
#    print(nsize)
    if(nsize>1)
      for (j in 2:nsize)
        if (all(!spec.neighb) | all(pdf.regions[spec.neighb]<1e-8) |
             all(spec.neighb | spec.regind) |
            all(pdf.regions[!(spec.regind | spec.neighb)]<1e-8))
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
#          if (spec.regind[reg]) cat("reg, all ",reg,"\n")  
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
#            if (is.na(spec.regind[reg]))  cat("reg, jump ",r1," ",cdf.local,"\n")  
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
#            if (spec.regind[reg]) cat("reg, neighbor ",reg,"\n")  
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
#   "  nni=",sum(!(spec.regind | spec.neighb)),"\n")
  } # for i           
  out
}
