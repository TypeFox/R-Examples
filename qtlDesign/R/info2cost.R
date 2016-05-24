"info2cost" <- function(sel.frac,cost,d=0,G=NULL,cross)
  {
    if(cross=="bc")
      info2cost.bc(sel.frac,cost,d,G)
    else if(cross=="f2")
      info2cost.f2(sel.frac,cost,d,G)
    else
      stop("Unknown cross ", cross, ".")
  }


"info2cost.bc" <-
function(sel.frac,cost,d=0,G=NULL)
  {
    if((d==0) & is.null(G))
      {
        ans <- info.bc(sel.frac,theta=0)/(1+cost*sel.frac)
      }
    else
      {
        if((d==0)|is.null(G)|(G<=0))
          {
            stop("Cannot compute with given d and G.")
          }
        else
          {
            theta <- recomb(d/100)
            ans <- info.bc(sel.frac,theta)/(1+cost*sel.frac*G/d)
          }
      }
    ans
  }


"info2cost.f2" <-
function(sel.frac,cost,d=0,G=NULL)
  {
    if((d==0) & is.null(G))
      {
        ans <- info.f2(sel.frac,theta=0)/(1+cost*sel.frac)
      }
    else
      {
        if((d==0)|is.null(G)|(G<=0))
          {
            stop("Cannot compute with given d and G.")
          }
        else
          {
            theta <- recomb(d/100)
            ans <- info.f2(sel.frac,theta)$add/(1+cost*sel.frac*G/d)
          }
      }
    ans
  }
