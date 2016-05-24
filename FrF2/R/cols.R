cols <- function(k, gen){
    main <- c(2^(0:(k-1)),gen)
    hilf <- combn(k+length(gen),2)
    fi2s <- rep(NA,ncol(hilf))
    for (i in 1:length(fi2s)){
          hilf2 <- which(sapply(Yates[1:(2^k-1)],function(obj2){ 
                                     obj <- mult.gen(Yates[main[hilf[,i]]])
                                     hh <- FALSE 
                                     if (length(obj2)==length(obj)) 
                                          if (all(obj2==obj)) hh <- TRUE
                                     hh}))
          if (length(hilf2)>0)
          fi2s[i]<-hilf2
          }
   # fi2s <- sort(table(fi2s))
    list(main=main,fi2s=fi2s,freecols=setdiff(1:(2^k-1),c(main,fi2s)))
}

