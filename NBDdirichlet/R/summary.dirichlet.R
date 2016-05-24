"summary.dirichlet" <-
function(object,      # object of "dirichlet" class 
                              t=1,      # multiple of time period
                              type=c("buy","freq","heavy","dup"),
                              digits=2, # decimal digits to keep for output
                              ## following parameters are specific to each summary type
                              freq.cutoff=5, # cutoff point for freq dist
                              heavy.limit=1:6, # freq range defining buyers with such purchase freq
                              dup.brand=1, # which focal brand to study its "duplication" statistics
                              ...
                              ) {
  obj <- object
  obj$period.set(t)  # change the time period

  result <- list()

  for (tt in type) {
    if (tt == "buy") {
      r <- data.frame (
                       ## Observed Market Share (Per Capita Purchase Rate)
                       ## Theoreticl Brand Penetration is:                                        
                       pen.brand=round(sapply(1:obj$nbrand, obj$brand.pen), digits),
                       ## Purchase Rate of the Brand
                       pur.brand=round(sapply(1:obj$nbrand, obj$brand.buyrate), digits),
                       ## [1] 1.8 1.8 1.7 1.7 1.7 1.7 1.6 1.6  (close to Observed Brand Purchase Rate)
                       ## Purchase Rate of the Category by the Brand Buyers
                       pur.cat=round(sapply(1:obj$nbrand, obj$wp), digits),
                       ## 3.2 3.2 3.3 3.3 3.3 3.3 3.4 3.4
                       row.names=obj$brand.name
                       )
      result[[tt]] <- round(r,digits)
    }
    else if (tt =="freq") {
      ## Freq Dist of Purchases of the Individual Brands (Table 4)
      ## p(r) = \sum_{n \geq r}^{n^{*}} p(r, n)

      prob.r <- function(r,j) {
        sum(sapply(r:obj$nstar, function(n,r,j) {obj$Pn(n) * obj$p.rj.n(r,n,j)}, r=r,j=j))
      }

      r <- matrix(0,obj$nbrand,freq.cutoff+2)

      for (j in 1:obj$nbrand)
        r[j,] <- c(sapply(0:freq.cutoff, prob.r,j=j), sum(sapply((freq.cutoff+1):obj$nstar, prob.r,j=j)))
                    
      dimnames(r) <- list(obj$brand.name, c(0:freq.cutoff, paste(freq.cutoff+1, "+",sep="")))
      result[[tt]] <- round(r,digits)
      }
    else if (tt == "heavy"){
      ## The penetration (b) and average purchase frequency (w) among frequent or infrequent buyers of
      ## the Category (Table 5)
      ## "heavy.limit" : the Category Buyers whose Buying Frequency falls in "limit"
      
      r <- matrix(0,obj$nbrand,2)
      Pn.sum <- sum(sapply(heavy.limit, obj$Pn))       #Prob of Category Purchase Freq in Set "limit"
      for (j in 1:obj$nbrand) {

        ## P(buy B at least once | buy Category "limit" (say 1-6) times)
        ## = [P(buy Cat 1-6 times) - P(buy cat 1-6 times, NO B)] / P(buy Cat 1-6 times)
        ## = 1 - \sum_{n=1}^6 P(n) P(0|n) / \sum_{n=1}^6 P(n)

        if (obj$check == T) {cat("compute penetration\n")}

        p0 <- 1 - obj$brand.pen(j, limit=heavy.limit)

        r[j,1] <- 1 - p0 / Pn.sum

        if (obj$check == T) {cat("compute purchase rate\n")}
        ## Given a "heavy" buyer (buying category 1-6 times), average brand purchase freq.
        ## = \sum_{n=1}^6 {P(n) [\sum_{r=1}^n r p(r|n)]} / \sum_{n=1}^6 { P(n) [1 - p(0|n)] }
        r[j,2] <- obj$brand.buyrate(j, heavy.limit) * obj$brand.pen(j) / (Pn.sum - p0)
      }
      dimnames(r) <- list(obj$brand.name, c("Penetration", "Avg Purchase Freq"))
      result[[tt]] <- round(r,digits)
    }
    else if (tt == "dup"){
      k <- dup.brand
      r <- rep(0,obj$nbrand)            # store result for Brand Duplication
      r[k] <- 1                         # Brand Duplication with Itself is 100%
      
      b.k <- obj$brand.pen(k)           # penetration for Brand k (dup.brand)
      others <- c(1:obj$nbrand)[-k] # list of Other Brands that the Focal Brand Buyer may buy

      for (j in others){
        ## The composite brand is [dup.brand=k, j]. Its penetration is:
        ## 1 - \sum_{n=0}^{-\Infinity} P_n * p_k(0|n) * p_j(0|n)

        ## p(0): Overall Prob. of Not Buying Brand k and j                                        
        p0 <- sum( sapply(0:obj$nstar,
                          function(i,x) {obj$Pn(i) * obj$p.rj.n(0,i,x)},
                          x=c(k,j)))
        
        b.j.k <- 1 - p0                               # Composite Brand [k,j]  Penetration
        b.j <- obj$brand.pen(j)
        b.jk <- b.j + b.k - b.j.k # proportion buying BOTH brand j and k.
        b.j.given.k <- b.jk / b.k

        r[j] <- b.j.given.k
      }
      names(r) <- obj$brand.name
      result[[tt]] <- round(r,digits)
    }
  }
  result
}

