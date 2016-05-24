GT.decision <-
function(TestStatistic, alpha=0.05, eta=alpha)
  {
    G <- length(TestStatistic)
    group.fdr.star.g <- array(0, G)
    Rgs <- array(0, G)
    ## Go in to each group, calculate the potential number of rejections and mark the corresponding hypothesis, calculate group.fdr.star.g
    for(g in 1:G)
      {
        mg <- TestStatistic[[g]]$mg
        order.fdr.j.g <- sort( TestStatistic[[g]]$fdr.j.g, decreasing=FALSE )
        R_g <- max ( ( cumsum( order.fdr.j.g)/c(1:mg) <= eta )* c(1:mg) )
        if(R_g>0){
          TestStatistic[[g]]$within.group.rej <- ( TestStatistic[[g]]$fdr.j.g <= order.fdr.j.g[R_g] )*1
          TestStatistic[[g]]$eta.g <- sum( order.fdr.j.g[1:R_g] )/R_g
          Rgs[g] <- R_g
        }else{
          TestStatistic[[g]]$within.group.rej <- ( TestStatistic[[g]]$fdr.j.g )*0
          TestStatistic[[g]]$eta.g <- 0
          Rgs[g] <- 0
        }
        TestStatistic[[g]]$fdr.star.g <- 1-( 1- TestStatistic[[g]]$eta.g)*( 1-TestStatistic[[g]]$fdr.g )
        group.fdr.star.g[g] <- TestStatistic[[g]]$fdr.star.g
      }

    ## Step 2. between-group decision
    ORD.g <- order( group.fdr.star.g, decreasing=FALSE )
    number.rej.group <- max ( ( cumsum( Rgs[ORD.g]*group.fdr.star.g[ORD.g] ) <= alpha*cumsum( Rgs[ORD.g] ) ) * c(1:G)  )
    if( number.rej.group>0 )
      {
        ind.rej.group <- ( group.fdr.star.g <= group.fdr.star.g[ ORD.g[ number.rej.group ] ] )*1
        for(g in 1:G){
          if( sum(TestStatistic[[g]]$within.group.rej)!=0){
            TestStatistic[[g]]$between.group.rej <- ind.rej.group[g]
            TestStatistic[[g]]$within.group.rej <- TestStatistic[[g]]$within.group.rej * ind.rej.group[g]
          }else{
            TestStatistic[[g]]$between.group.rej <- 0
            TestStatistic[[g]]$within.group.rej <- TestStatistic[[g]]$within.group.rej * 0
          }
        }
      }else{
        for(g in 1:G){
          TestStatistic[[g]]$between.group.rej <- 0
          TestStatistic[[g]]$within.group.rej <- 0*TestStatistic[[g]]$within.group.rej
        }
      }
    TestStatistic
  }
