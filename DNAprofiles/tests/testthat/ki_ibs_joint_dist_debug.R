Zdebug.ki.ibs.dist<-function(a,b,hyp.1,hyp.2="UN",hyp.true="UN",f.ki,f.true=f.ki,theta.ki=0,theta.true=theta.ki){
  # function computes the conditional joint dist of the ibs and ki (hyp.1 vs hyp.2) under hyp.true
  # a,b is genotype of profile @ locus
  # this function is not exported; only the function treating any number of loci is exported
  
  # look up allele freqs
  f.a.ki  <- as.vector(f.ki)[a];  f.b.ki  <- as.vector(f.ki)[b]
  f.a.hyp.true <- as.vector(f.true)[a]; f.b.hyp.true <- as.vector(f.true)[b]
  
  # look up ibd probs
  k.hyp.1 <- ibdprobs(hyp.1)
  k.hyp.2 <- ibdprobs(hyp.2)
  k.hyp.true <- ibdprobs(hyp.true)
  
  # compute dist
  # homozygous case gives three possibilities, heterozygous six
  if (a==b){ # x has a/a
    f.z.hyp.true <- (1-f.a.hyp.true)
    f.z.ki <- 1-f.a.ki
    f0.ki <- c(f.a.ki,f.z.ki)
    f0.true <- c(f.a.hyp.true,f.z.hyp.true)
    
    ## separate code for theta==0 and theta>0 --> useful for debugging
    # lr (hyp 1 vs unr)
    if (theta.ki==0){ 
      x1  <- c(k.hyp.1[1],                                              # z/z
               k.hyp.1[1]+k.hyp.1[2]*(1/2)*(1/f.a.ki),                   # a/z
               k.hyp.1[1]+k.hyp.1[2]*(1/f.a.ki)+k.hyp.1[3]*(1/f.a.ki^2))  # a/a
      # lr (v unr.) under hyp.2
      x2  <- c(k.hyp.2[1],                                              # z/z
               k.hyp.2[1]+k.hyp.2[2]*(1/2)*(1/f.a.ki),                   # a/z
               k.hyp.2[1]+k.hyp.2[2]*(1/f.a.ki)+k.hyp.2[3]*(1/f.a.ki^2))  # a/a      
    }else{
      x1 <- c(k.hyp.1[1], # z/z
              k.hyp.1[1]+k.hyp.1[2]/2*1/(pr.next.allele(i=1,seen=t(c(2,1,1)),f=f0.ki,theta=theta.ki)), #a/z
              k.hyp.1[1]+k.hyp.1[2]*1/(pr.next.allele(i=1,seen=t(c(1,1,1)),f=f0.ki,theta=theta.ki))+ #a/a
                k.hyp.1[3]*1/(pr.next.alleles(ij=t(c(1,1)),seen=t(c(1,1)),f=f0.ki,theta=theta.ki)))    
      x2 <- c(k.hyp.2[1], # z/z
              k.hyp.2[1]+k.hyp.2[2]/2*1/(pr.next.allele(i=1,seen=t(c(2,1,1)),f=f0.ki,theta=theta.ki)), #a/z
              k.hyp.2[1]+k.hyp.2[2]*1/(pr.next.allele(i=1,seen=t(c(1,1,1)),f=f0.ki,theta=theta.ki))+ #a/a
                k.hyp.2[3]*1/(pr.next.alleles(ij=t(c(1,1)),seen=t(c(1,1)),f=f0.ki,theta=theta.ki)))              
    }
    
    x <- x1 # possible ki's
    if (!identical(k.hyp.2,ibdprobs("UN"))) x<- x/x2
    
    ibs <- c(0,1,2)
    if (theta.true==0){
      p.x <- c(k.hyp.true[1]*f.z.hyp.true^2,                            # z/z
               k.hyp.true[1]*2*f.a.hyp.true*f.z.hyp.true+k.hyp.true[2]*f.z.hyp.true,   # a/z
               k.hyp.true[1]*f.a.hyp.true^2+k.hyp.true[2]*f.a.hyp.true+k.hyp.true[3])  # a/a      
    }else{
      p.x <- c(k.hyp.true[1]*pr.next.alleles(ij=t(c(2,2)),seen=t(c(1,1)),f=f0.true,theta=theta.true), #z/z
               k.hyp.true[1]*2*pr.next.alleles(ij=t(c(1,2)),seen=t(c(1,1)),f=f0.true,theta=theta.true)+#a/z
                 k.hyp.true[2]*pr.next.alleles(ij=t(c(2)),seen=t(c(1,1)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*pr.next.alleles(ij=t(c(1,1)),seen=t(c(1,1)),f=f0.true,theta=theta.true)+#a/a
                 k.hyp.true[2]*pr.next.alleles(ij=t(c(1)),seen=t(c(1,1)),f=f0.true,theta=theta.true)+
                 k.hyp.true[3])      
    }
    
  }else{ # x has a/b
    f.z.hyp.true <- (1-f.a.hyp.true-f.b.hyp.true)
    f.z.ki <- 1-f.a.ki-f.b.ki
    f0.ki <- c(f.a.ki,f.b.ki,f.z.ki)
    f0.true <- c(f.a.hyp.true,f.b.hyp.true,f.z.hyp.true)
    # lr (hyp.1 vs unr.)
    if (theta.ki==0){
      x1 <- c(k.hyp.1[1],                             # z/z
              k.hyp.1[1]+k.hyp.1[2]*(1/4)*(1/f.a.ki),  # a/z
              k.hyp.1[1]+k.hyp.1[2]*(1/2)*(1/f.a.ki),  # a/a
              k.hyp.1[1]+k.hyp.1[2]*(1/4)*(1/f.b.ki),  # b/z
              k.hyp.1[1]+k.hyp.1[2]*(1/2)*(1/f.b.ki),  # b/b
              k.hyp.1[1]+k.hyp.1[2]*(1/4)*(1/f.b.ki+1/f.a.ki)+k.hyp.1[3]*(1/2)*(1/(f.a.ki*f.b.ki))) # a/b      
      x2 <- c(k.hyp.2[1],                             # z/z
              k.hyp.2[1]+k.hyp.2[2]*(1/4)*(1/f.a.ki),  # a/z
              k.hyp.2[1]+k.hyp.2[2]*(1/2)*(1/f.a.ki),  # a/a
              k.hyp.2[1]+k.hyp.2[2]*(1/4)*(1/f.b.ki),  # b/z
              k.hyp.2[1]+k.hyp.2[2]*(1/2)*(1/f.b.ki),  # b/b
              k.hyp.2[1]+k.hyp.2[2]*(1/4)*(1/f.b.ki+1/f.a.ki)+k.hyp.2[3]*(1/2)*(1/(f.a.ki*f.b.ki))) # a/b
    }else{
      x1 <- c(k.hyp.1[1], # z/z
              k.hyp.1[1]+k.hyp.1[2]/4/pr.next.alleles(ij=t(c(1)),seen=t(c(3,1,2)),f=f0.ki,theta=theta.ki), # a/z
              k.hyp.1[1]+k.hyp.1[2]/2/pr.next.alleles(ij=t(c(1)),seen=t(c(1,1,2)),f=f0.ki,theta=theta.ki), # a/a
              k.hyp.1[1]+k.hyp.1[2]/4/pr.next.alleles(ij=t(c(2)),seen=t(c(3,1,2)),f=f0.ki,theta=theta.ki), # b/z
              k.hyp.1[1]+k.hyp.1[2]/2/pr.next.alleles(ij=t(c(2)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki), # b/b
              k.hyp.1[1]+k.hyp.1[2]/4/pr.next.alleles(ij=t(c(1)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki)+ # a/b
                k.hyp.1[2]/4/pr.next.alleles(ij=t(c(2)),seen=t(c(2,1,1)),f=f0.ki,theta=theta.ki)+
                k.hyp.1[3]/2/pr.next.alleles(ij=t(c(1,2)),seen=t(c(1,2)),f=f0.ki,theta=theta.ki))      
      x2 <- c(k.hyp.2[1], # z/z
              k.hyp.2[1]+k.hyp.2[2]/4/pr.next.alleles(ij=t(c(1)),seen=t(c(3,1,2)),f=f0.ki,theta=theta.ki), # a/z
              k.hyp.2[1]+k.hyp.2[2]/2/pr.next.alleles(ij=t(c(1)),seen=t(c(1,1,2)),f=f0.ki,theta=theta.ki), # a/a
              k.hyp.2[1]+k.hyp.2[2]/4/pr.next.alleles(ij=t(c(2)),seen=t(c(3,1,2)),f=f0.ki,theta=theta.ki), # b/z
              k.hyp.2[1]+k.hyp.2[2]/2/pr.next.alleles(ij=t(c(2)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki), # b/b
              k.hyp.2[1]+k.hyp.2[2]/4/pr.next.alleles(ij=t(c(1)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki)+ # a/b
                k.hyp.2[2]/4/pr.next.alleles(ij=t(c(2)),seen=t(c(2,1,1)),f=f0.ki,theta=theta.ki)+
                k.hyp.2[3]/2/pr.next.alleles(ij=t(c(1,2)),seen=t(c(1,2)),f=f0.ki,theta=theta.ki))
    }
    
    x <- x1 # possible ki's
    if (!identical(k.hyp.2,ibdprobs("UN"))) x<- x/x2
    
    ibs <- c(0,1,1,1,1,2)
    
    if (theta.true==0){
      p.x <- c(k.hyp.true[1]*f.z.hyp.true^2,                                # z/z
               k.hyp.true[1]*2*f.a.hyp.true*f.z.hyp.true+k.hyp.true[2]*(1/2)*f.z.hyp.true, # a/z
               k.hyp.true[1]*f.a.hyp.true^2+k.hyp.true[2]*(1/2)*f.a.hyp.true,         # a/a
               k.hyp.true[1]*2*f.b.hyp.true*f.z.hyp.true+k.hyp.true[2]*(1/2)*f.z.hyp.true, # b/z
               k.hyp.true[1]*f.b.hyp.true^2+k.hyp.true[2]*(1/2)*f.b.hyp.true,         # b/b
               k.hyp.true[1]*2*f.a.hyp.true*f.b.hyp.true+k.hyp.true[2]*((1/2)*(f.a.hyp.true+f.b.hyp.true))+k.hyp.true[3]) # a/b      
    }else{
      p.x <- c(k.hyp.true[1]*pr.next.alleles(ij=t(c(3,3)),seen=t(c(1,2)),f=f0.true,theta=theta.true),# z/z 
               k.hyp.true[1]*2*pr.next.alleles(ij=t(c(1,3)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # a/z
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(3)),seen=t(c(1,2)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*1*pr.next.alleles(ij=t(c(1,1)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # a/a
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(1)),seen=t(c(1,2)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*2*pr.next.alleles(ij=t(c(2,3)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # b/z
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(3)),seen=t(c(1,2)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*1*pr.next.alleles(ij=t(c(2,2)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # b/b
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(2)),seen=t(c(1,2)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*2*pr.next.alleles(ij=t(c(1,2)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # a/b
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(2)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(1)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+
                 k.hyp.true[3] )      
    }
    
  }
  
  if (any(is.na(x))) stop("NA or NaNs encountered")
  
  p.nonzero <- p.x>0 # only retain the events with non-zero probability
  
  cbind(p.x[p.nonzero],x[p.nonzero],ibs[p.nonzero])
}
