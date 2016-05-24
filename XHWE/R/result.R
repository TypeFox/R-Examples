result <-
function(ped,loci,dv,start.rho,simuno,status_missing,allele_missing,header,itertime,SNP_name) {

  n2f <- sum(ped[,5]==2 & ped[,7]==1 & ped[,8]==1) 

  n1f <- sum((ped[,5]==2 & ped[,7]==1 & ped[,8]==2) | (ped[,5]==2 & ped[,7]==2 & ped[,8]==1))   
   
  n0f <- sum(ped[,5]==2 & ped[,7]==2 & ped[,8]==2)
   
  n1m <- sum(ped[,5]==1 & ped[,7]==1 )  
             
  n0m <- sum(ped[,5]==1 & ped[,7]==2 )  

  nf <- n2f + n1f + n0f

  nm <- n1m + n0m

  if(nm > 0 & nf < 0.1) {

  print(cbind(SNP_name,"Warning: all the female founders are missing!")) 

  pm.male <- n1m / nm

  test.testa <- data.frame(LRT.test = 0, LRT.test2 = 0, LRT.test1 = 0, z0=0, z1=0, z2=0)

  test.pvalue <- data.frame(Pvalue = 1, Pvalue.boot = 1, Pvalue2 = 1, Pvalue1 = 1, Pvalue1.boot = 1, Pvalue.z0 = 1, Pvalue.z1 = 1, Pvalue.z2 = 1)

  test.para <- data.frame(p.0=NA,p.01=NA, rho.01=NA, pf.z.01=NA, rho.z.01=NA, pm.02=pm.male, pf.02=NA, pm=pm.male, pf1=NA, rho1=NA, pf=NA, rho.z=NA)

  }

  else if(nm < 0.1 & nf > 0) {

    print(cbind(SNP_name,"Warning: all the male founders are missing!")) 

    #############################################
    # Z2-stat with df being 1

    nfT2 <- 2 * nf
  
    pf <- (n2f * 2 + n1f) / nfT2
  
    pf2 <- (pf)^2

    delta <- (n2f / nf) - pf2
  
    pq <- pf * (1 - pf)
  
    Z2 <- sqrt(nf) * (delta + pq / nfT2) / pq
  
    z2 <- (Z2)^2
  
    #############################################

    rho.z <- delta/pq
    
    emf1 <- emf(0, 0, n2f, n1f, n0f, 0, nf, start.rho, dv, itertime)

    pf1.last <- emf1$pf.last
  
    rho1.last <- emf1$rho.last

    pf1qf1 <- pf1.last * (1 - pf1.last)

    tt1 <- rho1.last * pf1qf1
  
    p2.e <- pf1.last^2 + tt1
  
    p1.e <- 2 * (pf1qf1 - tt1)
  
    #############################################
    # Chi-square test with df being 1
  
    p.all <- pf
  
    p2.all <- p.all^2
  
    p1.all <- 2 * p.all * (1 - p.all)
  
    LRT.num <- Likelihoodfun(0, p2.e, p1.e, 0, 0, n2f, n1f, n0f)
  
    LRT.den <- Likelihoodfun(0, p2.all, p1.all, 0, 0, n2f, n1f, n0f)
  
    LRT.test1 <- 2 * (LRT.num - LRT.den)

  
    ###########################################################
    #Z2 with df being 1
    Pvalue.z2 <- pchisq(z2, df = 1,lower.tail=F)

    # Chi-square test with df being 1 for rho=0
    Pvalue1 <- pchisq(LRT.test1, df = 1,lower.tail=F)

    test.testa <- data.frame(LRT.test = 0, LRT.test2 = 0, LRT.test1 = LRT.test1, z0=0, z1=0, z2=z2)

    test.pvalue <- data.frame(Pvalue = 1, Pvalue.boot = 1, Pvalue2 = 1, Pvalue1 = Pvalue1, Pvalue1.boot = 1, Pvalue.z0 = 1, Pvalue.z1 = 1, Pvalue.z2 = Pvalue.z2)

    test.para <- data.frame(p.0=NA, p.01=NA, rho.01=NA, pf.z.01=NA, rho.z.01=rho.z, pm.02=NA, pf.02=pf, pm=NA, pf1=pf1.last, rho1=rho1.last, pf=pf, rho.z=rho.z)

 }

  
  else if (nf > 0 & nm > 0) {

    LRT.array4 <- array(0, dim=c(simuno, 1))
  
    LRT.array6 <- array(0, dim=c(simuno, 1))

    #############################################
    # Z1-stat with df being 1

    nfT2 <- 2 * nf
  
    pf <- (n2f * 2 + n1f) / nfT2
  
    pm <- n1m / nm
  
    pf2 <- (pf)^2
  
    var1 <- pm * (1-pm) / nm 
  
    var2 <- (pf - 2 * pf2 + n2f / nf) / nfT2
  
    Z1 <- (pm - pf) / (sqrt(var1 + var2))
  
    z1 <- (Z1)^2
  
    #############################################
    # Z2-stat with df being 1
  
    delta <- (n2f / nf) - pf2
  
    pq <- pf * (1 - pf)
  
    Z2 <- sqrt(nf) * (delta + pq / nfT2) / pq
  
    z2 <- (Z2)^2
  
    #############################################
    # Z0-stat with df being 2
  
    z0 <- z1 + z2
  
    rho.z <- delta/pq
    
    emf1 <- emf(n1m, n0m, n2f, n1f, n0f, nm, nf, start.rho, dv, itertime)
  
    pm1.last <- emf1$pm.last
  
    pf1.last <- emf1$pf.last
  
    rho1.last <- emf1$rho.last

    pf1qf1 <- pf1.last * (1 - pf1.last)

    tt1 <- rho1.last * pf1qf1
  
    p2.e <- pf1.last^2 + tt1
  
    p1.e <- 2 * (pf1qf1 - tt1)
  
    #############################################
    # Chi-square test with df being 2

    sum1 <- nm + nfT2
  
    p.all <- (n1m + n2f * 2 + n1f)/sum1
  
    p2.all <- p.all^2
  
    p1.all <- 2 * p.all * (1 - p.all)
  
    LRT.num <- Likelihoodfun(pm1.last, p2.e, p1.e, n1m, n0m, n2f, n1f, n0f)
  
    LRT.den <- Likelihoodfun(p.all, p2.all, p1.all, n1m, n0m, n2f, n1f, n0f)
  
    LRT.test <- 2 * (LRT.num - LRT.den)
 
  
    ###########################################################
    # Chi-square test with df being 1 for rho = 0
  
    pf1 <- 2 * pq
  
    LRT.den1 <- Likelihoodfun(pm, pf2, pf1, n1m, n0m, n2f, n1f, n0f)
  
    LRT.test1 <- 2 * (LRT.num - LRT.den1)
  
  
    ###########################################################
    # Chi-square test with df being 1 for pm = pf
    em2 <- emc(n1m, n0m, n2f, n1f, n0f, nm, nf, start.rho, dv, itertime)
  
    pc2.last <- em2$pc.last
  
    rho2.last <- em2$rho.last

    pcqc <- pc2.last * (1 - pc2.last)

    tt2 <- rho2.last * pcqc
  
    p2f2 <- pc2.last^2 + tt2
  
    p1f2 <- 2 * (pcqc - tt2)
  
    LRT.den2 <- Likelihoodfun(pc2.last, p2f2, p1f2, n1m, n0m, n2f, n1f, n0f)
  
    LRT.test2 <- 2 * (LRT.num - LRT.den2)
  
    ###########################################################
    #BOOTSTRAP

    #rho=0 and pm=pf    
    p0.all <- 1 - p2.all - p1.all
    
    pf0 <- 1 - pf2 - pf1 

    nm.all <- array(0,dim=c(simuno,2))

    nm.a0 <- rbinom(simuno,nm,p.all)

    nm.all[,1] <- nm.a0

    nm.all[,2] <- nm - nm.a0 
    
    nf.all <- rmultinom(simuno, size = nf, prob = c(p2.all,p1.all,p0.all))

    n.all <- cbind(nm.all,t(nf.all))

    nf.all02 <- rmultinom(simuno, size = nf, prob = c(pf0,pf1,pf2))

    nf.all2 <- t(nf.all02)
     
    for(j in 1:simuno){

      n1ma<-n.all[j,1]

      n0ma<-n.all[j,2]

      n2f.a<-n.all[j,3]

      n1f.a<-n.all[j,4]

      n0f.a<-n.all[j,5]
    
      pfa <- (n1ma + n2f.a * 2 + n1f.a)/sum1
    
      pf2a <- (pfa)^2
    
      pf1a <- 2 * pfa * (1 - pfa)
    
      emfa <- emf(n1ma, n0ma, n2f.a, n1f.a, n0f.a, nm, nf, start.rho, dv, itertime)
    
      pma.last <- emfa$pm.last
    
      pfa.last <- emfa$pf.last
    
      rhoa.last <- emfa$rho.last

      pfaqfa <- pfa.last * (1 - pfa.last)

      tt3 <- rhoa.last * pfaqfa
    
      p2.aa <- pfa.last^2 + tt3
    
      p1.aa <- 2 * (pfaqfa - tt3)
    
      LRT.dena <- Likelihoodfun(pfa, pf2a, pf1a, n1ma, n0ma, n2f.a, n1f.a, n0f.a)
    
      LRT.numa <- Likelihoodfun(pma.last, p2.aa, p1.aa, n1ma, n0ma, n2f.a, n1f.a, n0f.a)
    
      LRT.testa <- 2 * (LRT.numa - LRT.dena)
    
      LRT.array6[j,1] <- LRT.testa
    
    
      #rho=0 
    
      n0fb <-nf.all2[j,1]

      n1fb <-nf.all2[j,2]

      n2fb <-nf.all2[j,3]

      pfb <- (n2fb * 2 + n1fb) / nfT2
    
      pf2b <- (pfb)^2
    
      pf1b <- 2 * pfb * (1 - pfb)
    
      emfb <- emf(n1m, n0m, n2fb, n1fb, n0fb, nm, nf, start.rho, dv, itertime)
    
      pmb.last <- emfb$pm.last
    
      pfb.last <- emfb$pf.last
    
      rhob.last <- emfb$rho.last

      pfbqfb <- pfb.last * (1 - pfb.last)
   
      tt4 <- rhob.last * pfbqfb 
    
      p2.b <- pfb.last^2 + tt4
    
      p1.b <- 2 * (pfbqfb - tt4)
    
      LRT.numb <- Likelihoodfun(pmb.last, p2.b, p1.b, n1m, n0m, n2fb, n1fb, n0fb)
    
      LRT.denb <- Likelihoodfun(pm, pf2b, pf1b, n1m, n0m, n2fb, n1fb, n0fb)
    
      LRT.testb <- 2 * (LRT.numb - LRT.denb)
    
      LRT.array4[j,1] <- LRT.testb
    
    }
  
  
    #Z1 with df being 1
    Pvalue.z1 <- pchisq(z1, df = 1,lower.tail=F)
  
    #Z2 with df being 1
    Pvalue.z2 <- pchisq(z2, df = 1,lower.tail=F)
  
    #Z0 with df being 2
    Pvalue.z0 <-  pchisq(z0, df = 2,lower.tail=F)
    
    # Chi-square test with df being 2
    Pvalue <- pchisq(LRT.test, df = 2, lower.tail=F)

    Pvalue.boot <- sum(LRT.array6[,1] >= LRT.test) / simuno
    
    # Chi-square test with df being 1 for rho=0
    Pvalue1 <- pchisq(LRT.test1, df = 1,lower.tail=F)

    Pvalue1.boot <- sum(LRT.array4[,1] >= LRT.test1) / simuno
    
    # Chi-square test with df being 1 for pm = pf
    Pvalue2 <- pchisq(LRT.test2, df = 1,lower.tail=F)
  
    test.testa <- data.frame(LRT.test = LRT.test, LRT.test2 = LRT.test2, LRT.test1 = LRT.test1, z0=z0, z1=z1, z2=z2)

    test.pvalue <- data.frame(Pvalue = Pvalue, Pvalue.boot = Pvalue.boot, Pvalue2 = Pvalue2, Pvalue1 = Pvalue1, Pvalue1.boot = Pvalue1.boot, Pvalue.z0 = Pvalue.z0, Pvalue.z1 = Pvalue.z1, Pvalue.z2 = Pvalue.z2)

    test.para <- data.frame(p.0=p.all, p.01=pc2.last, rho.01=rho2.last, pf.z.01=p.all, rho.z.01=rho.z, pm.02=pm, pf.02=pf, pm=pm1.last, pf1=pf1.last, rho1=rho1.last, pf=pf, rho.z=rho.z)
  
  }


else {

    test.testa <- data.frame(LRT.test = 0, LRT.test2 = 0, LRT.test1 = 0, z0=0, z1=0, z2=0)

    test.pvalue <- data.frame(Pvalue = 1, Pvalue.boot = 1, Pvalue2 = 1, Pvalue1 = 1, Pvalue1.boot = 1, Pvalue.z0 = 1, Pvalue.z1 = 1, Pvalue.z2 = 1)

    test.para <- data.frame(p.0=NA, p.01=NA, rho.01=NA, pf.z.01=NA, rho.z.01=NA, pm.02=NA, pf.02=NA, pm=NA, pf1=NA, rho1=NA, pf=NA, rho.z=NA)


    }

list(test.testa, test.pvalue, test.para)
  
}
