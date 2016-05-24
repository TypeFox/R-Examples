## Computes the sum: sum{i<j<k} p_i^I p_j^J p_k^K, where c(I,J,K) = pow
Pijk <- function(p,pow=c(1,1,1)){
  m <- length(p)
  p.inner <- unlist(lapply(2:(m-1),function(x) sum((p^pow[2])[x:(m-1)]*(sum(p^pow[3])-cumsum(p^pow[3])[x:(m-1)]))))
  sum((p^pow[1])[1:(m-2)]*p.inner)
}

## Computes the sum: sum{i<j<k<l} p_i^I p_j^J p_k^K p_l^L, where c(I,J,K,L) = pow
Pijkl <- function(p,pow=c(1,1,1,1)){
  m <- length(p)
  sum((p^pow[1])[1:(m-3)]*unlist(lapply(1:(m-3),function(x) Pijk(p[-(1:x)],pow=pow[-1]))))
}

## Computes the P_{0/0}, P_{0/1}, P_{1/0} for a given locus
RPs <- function(p,t=0,r=0,k=rep(0,3),browse=FALSE){
  ## Modified versions of S_n to account for no summation over R-alleles:
  s1 <- sum(p); s2 <- sum(p^2); s3 <- sum(p^3); s4 <- sum(p^4)
  ## ASSUMES THETA=0
  ## Profile 1: HOM = AiAi
  pAiAi.AiAi <- s4 ## AiAi,AiAi
  pAiAi.AjAj <- s2^2-s4 ## AiAi,AjAj
  pAiAi.AiAj <- s3-s4 ## AiAi,AiAj
  pAiAi.AjAk_ijk <- Pijk(p,pow=c(2,1,1))
  pAiAi.AjAk_jik <- Pijk(p,pow=c(1,2,1))
  pAiAi.AjAk_jki <- Pijk(p,pow=c(1,1,2))
  pAiAi.AjAk <- 2*(pAiAi.AjAk_ijk + pAiAi.AjAk_jik + pAiAi.AjAk_jki) ## AiAi,AjAk  ## in the sum i<j<k ; j<i<k; j<k<i each permutating of j and k gives 2 outcomes
  ## Profile 1: HET = AiAj
  pAiAj.AiAj <- s2^2-s4 ## AiAj,AiAj
  pAiAj.AiAk_ijk <- Pijk(p,pow=c(2,1,1)) 
  pAiAj.AiAk_jik <- Pijk(p,pow=c(1,2,1))
  pAiAj.AiAk_jki <- Pijk(p,pow=c(1,1,2))
  pAiAj.AiAk <- 2*(pAiAj.AiAk_ijk + pAiAj.AiAk_jik + pAiAj.AiAk_jki) ## AiAj,AiAk ## Equals sum_{i,j,k}^{!=} pi^2 pj pk
  pAiAj.AkAl <- Pijkl(p,pow=c(1,1,1,1)) ## AiAj,AkAl ## Equals sum_{i,j,k,l}^{!=} pi pj pk pl ## There are 4 ways to choose 1st index, 3 for 2nd, 2 for 3rd = 4*3*2*1 = 24
  ##
  weir2007 <- c("AiAi.AiAi" = pAiAi.AiAi, ## Done
                "AiAi.AjAj" = pAiAi.AjAj, ## Done
                "AiAi.AiAj" = 2*2*pAiAi.AiAj, ## Done
                "AiAi.AjAk" = 2*pAiAi.AjAk, ## 
                "AiAj.AiAj" = 2*pAiAj.AiAj, ## Done
                "AiAj.AiAk" = 4*pAiAj.AiAk, ## 
                "AiAj.AkAl" = 24*pAiAj.AkAl) ## Done
  ## if(abs(sum(weir2007)-1)>10^(-10)) cat("Don't sum to one")
  ## Compute the probabilities of shared alleles
  p0 <- 0 + ## AiAi.AiAi
    (1-r)^2*pAiAi.AjAj + ## AiAi.AjAj ## If no R allele - no match
      0 + ## 4 AiAi.AiAj
        (1-r)^3*2*pAiAi.AjAk + (1-r)^2*r*2*pAiAi.AjAk + ## 2 AiAi.AjAk ## If no R allele - no match ## 
          0 + ## 2 AiAj.AiAj
            0 + ## 4 AiAj.AiAk
              (1-r)^4*24*pAiAj.AkAl + (4*(1-r)^2*r^2 + 4*(1-r)^3*r)*4*pAiAj.AkAl ## 24 AiAj.AkAl
  p1 <- 0 + ## AiAi,AiAi
    0 + ## AiAi.AjAj
      (1-r)^2*4*pAiAi.AiAj + ## 4 AiAi.AiAj ## if no R or if Aj=R - partial match
        (1-r)^2*r*2*pAiAi.AjAk + ## 2 AiAi.AjAk
          0 + ## 2 AiAj.AiAj
            (1-r)^3*4*pAiAj.AiAk + r*(1-r)^2*4*pAiAj.AiAk + (4*r*(1-r)^2 + 2*(1-r)*r^2)*4*pAiAj.AiAk_jik +  ## 4 AiAj.AiAk ## If no R - partial macth
              (4*(1-r)^2*r^2 + 4*(1-r)^3*r)*20*pAiAj.AkAl ## 24 AiAj.AkAl
  p2 <- pAiAi.AiAi + ## AiAi.AiAi
    (2*r*(1-r) + r^2)*pAiAi.AjAj + ## AiAi.AjAj ## If one or both profs have R - match
      (r^2 + 2*(1-r)*r)*4*pAiAi.AiAj + ## 4 AiAi.AiAj ## If both Ai=R and Aj=R or only Ai=R - match
        (3*(1-r)*r^2 + (1-r)^2*r + r^3)*2*pAiAi.AjAk + ## 2 AiAi.AjAk
          2*pAiAj.AiAj + ## 2 AiAj.AiAj
            (4*r*(1-r)^2 + 2*(1-r)*r^2)*4*(pAiAj.AiAk_ijk + pAiAj.AiAk_jki) + (2*(1-r)*r^2 + r^3)*4*pAiAj.AiAk + ## 4 AiAj.AiAk // 
              (2*(1-r)^2*r^2 + 4*(1-r)*r^3 + r^4)*24*pAiAj.AkAl ## 24 AiAj.AkAl
  if(browse) {
    comb <- c("AiAi.AiAi" = 0 + 0 + pAiAi.AiAi, ## AiAi.AiAi 
              "AiAi.AjAj" = (1-r)^2*pAiAi.AjAj + 0 + (2*r*(1-r) + r^2)*pAiAi.AjAj, ## AiAi.AjAj 
              "AiAi.AiAj" =0 + (1-r)^2*4*pAiAi.AiAj + (r^2 + 2*(1-r)*r)*4*pAiAi.AiAj, ## 4 AiAi.AiAj 
              "AiAi.AjAk" = (1-r)^3*2*pAiAi.AjAk + (1-r)^2*r*2*pAiAi.AjAk + (1-r)^2*r*2*pAiAi.AjAk + (3*(1-r)*r^2 + (1-r)^2*r + r^3)*2*pAiAi.AjAk, ## 2 AiAi.AjAk
              "AiAj.AiAj" = 0 + 0 + 2*pAiAj.AiAj, ## 2 AiAj.AiAj
              "AiAj.AiAk" = 0 + (1-r)^3*4*pAiAj.AiAk + r*(1-r)^2*4*pAiAj.AiAk + (4*r*(1-r)^2 + 2*(1-r)*r^2)*4*pAiAj.AiAk_jik +
              (4*r*(1-r)^2 + 2*(1-r)*r^2)*4*(pAiAj.AiAk_ijk + pAiAj.AiAk_jki) + (2*(1-r)*r^2 + r^3)*4*pAiAj.AiAk, ## 4 AiAj.AiAk 
              "AiAj.AkAl" = (1-r)^4*24*pAiAj.AkAl + (4*(1-r)^2*r^2 + 4*(1-r)^3*r)*4*pAiAj.AkAl + (4*(1-r)^2*r^2 + 4*(1-r)^3*r)*20*pAiAj.AkAl +
              (2*(1-r)^2*r^2 + 4*(1-r)*r^3 + r^4)*24*pAiAj.AkAl) ## 24 AiAj.AkAl 
    return(list(weir=weir2007,comb=comb))
  }
  c(p0,p1,p2)
}
