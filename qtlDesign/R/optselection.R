"optselection" <- function(cost,d=0,G=NULL,cross)
  {
    if(cross=="bc")
      optselection.bc(cost,d,G)
    else if(cross=="f2")
      optselection.f2(cost,d,G)
    else
      stop("Unknown cross ", cross, ".")
  }

"optselection.bc" <-
function(cost,d=0,G=NULL)
  {
    optimize(f=info2cost.bc,interval=c(0.0001,0.9999),maximum=TRUE,
             G=G,d=d,cost=cost)$maximum
  }

"optselection.f2" <-
function(cost,d=0,G=NULL)
  {
    optimize(f=info2cost.f2,interval=c(0.0001,0.9999),maximum=TRUE,
             G=G,d=d,cost=cost)$maximum
  }

