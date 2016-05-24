#(size n.three way nested. model 4 c)
# Section 3.4.2.3. test factor C in B in A

# Three way nested classification. Model IV
# Factor A and C fixed, B random. Determining n
# Testing hypothesis about factor C in B in A
size_n.three_way_nested.model_4_c <-
  function(alpha, beta, delta, a, b, c, cases) ## , maxrekur=20)
{
  ##rekur <- 0
  if(is.na(b)){
    b <- 2
    n <- size_n.three_way_nested.model_4_c(alpha, beta, delta, a, b, c, cases)
    ##cat(paste("debug:",n*b,"\n"))
    min <- b*n
    b.min <- b
    n.min <- n
    while(!is.na(n)){## & rekur<maxrekur){
      ##rekur <- rekur+1
      ##if(rekur>=25) browser()
      b <- b+1
      n <- size_n.three_way_nested.model_4_c(alpha, beta, delta, a, b, c, cases)
      if(is.na(n) | b*n>4*10^4)break
      ##cat(paste("debug:",n*b,"\n"))
      if(n*b<min){
        min <- n*b
        n.min <- n
        b.min <- b
      }
    }
    cat(paste("\nminimum value of b*n is ",b.min*n.min,"\n\n"))
    list(b=b.min,n=n.min)

  } else {
    n <- 2
    dfn <- a*b*(c-1)
    dfd <- a*b*c*(n-1)
    if (cases == "maximin")
      {
        lambda <- 0.5*n*delta*delta
      }
    else if (cases == "minimin")
      {
        lambda <- 0.25*a*b*c*n*delta*delta
      }
    beta.calculated <- Beta(alpha, dfn, dfd, lambda)
    if (is.nan(beta.calculated) || beta.calculated < beta )
      {
        warning(paste("Given parameter will result in too high power.",
                      "To continue either increase the precision or ",
                      "decrease the level of factors."))
        return(NA)
      }
    else
      {
        n <- 5    
        n.new <- 1000
        while (abs(n - n.new)>1e-4)
          ## on windows accuracy only allows ">1e-4", on other systems ">1e-6" works !
          {
            n <- n.new
            dfn <- a*b*(c-1)
            dfd <- a*b*c*(n-1)
            lambda <- ncp(dfn,dfd,alpha,beta)
            if (cases == "maximin")
              {
                n.new <- 2*lambda/(delta*delta)
              }
            else if (cases == "minimin")
              {
                n.new <- 4*lambda/(a*b*c*delta*delta)
              }
          }  
        return(ceiling(n.new))
      }
  }
}


# example
# size.3_4_2_3.test_factor_CinBinA(0.05, 0.1, 0.5, 6, 5, 4, "maximin")
# size.3_4_2_3.test_factor_CinBinA(0.05, 0.1, 0.5, 6, 5, 4, "minimin")


optim_size_n.three_way_nested.model_4_c <- function(alpha, beta, delta, a, b, c, cases){
  n <- size_n.three_way_nested.model_4_c(alpha, beta, delta, a, b, c, cases)
  #print(b*n)
  min <- b*n
  b.min <- b
  n.min <- n
  while(!is.na(n)){
    b <- b+1
    n <- size_n.three_way_nested.model_4_c(alpha, beta, delta, a, b, c, cases)
    if(is.na(n) | b*n>4*10^4)break
    #print(n*b)
    if(n*b<min){
      min <- n*b
      n.min <- n
      b.min <- b
    }
  }
  cat(paste("\nminimum value of b*n is ",b.min*n.min,"\n\n"))
  list(b=b.min,n=n.min)
}
