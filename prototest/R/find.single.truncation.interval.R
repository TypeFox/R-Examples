#### function for finding the intervals on x > 0 such that
#### q.sqrt(x) + r.sqrt(1+x) + s <= 0
####  input:
####      - q, r, s = Coefficients for the truncation region function
find.single.truncation.interval <-
function(q, r, s, verbose=FALSE){
  if (verbose){
    print (paste('q = ', q, sep=''), quote=FALSE)
    print (paste('r = ', r, sep=''), quote=FALSE)
    print (paste('s = ', s, sep=''), quote=FALSE)
  }
  
  ## check a whole bunch of cases
  if (sign(q) == sign(r)){ # derivative never becomes 0
    #print ('sign(q) = sign(r) -- monotone function')
    if (q > 0){ # derivative always positive, so increasing
      return (interval.if.derivative.always.positive(q, r, s, start=0, verbose=verbose))
    }else{ # derivatve always negative, so decreasing
      return (interval.if.derivative.always.negative(q, r, s, start=0, verbose=verbose))
    }
  }else{ # derivative becomes 0 once
    if (abs(r) <= abs(q)){ # does not turn on the positive half line
      if (q < r/sqrt(2)){ # derivative always negative - testing it at x = 1 (same for all x > 0)
        return (interval.if.derivative.always.negative(q, r, s, start=0, verbose=verbose))
      }else{ # derivative always positive
        return (interval.if.derivative.always.positive(q, r, s, start=0, verbose=verbose))
      }
    }else{ # have to check two intervals for roots
      #grad.at.0.sign = sign (gradient.truncation.region.function(0, q, r))
      func.at.0 = truncation.region.function(0, q, r, s)
      func.at.0.sign = sign (func.at.0)
      func.at.turning.point = truncation.region.function(q^2/(r^2-q^2), q, r, s)
      func.at.turning.point.sign = sign (func.at.turning.point)
      
      #print (paste('At 0: ', func.at.0, sep=''))
      #print (paste('At turning point (',q^2/(r^2-q^2), '): ', func.at.turning.point, sep=''))
      
      if (func.at.turning.point.sign < 0 & func.at.0.sign < 0 & func.at.turning.point >= func.at.0){ # turns below zero and will keep going negative
        if (verbose){print ("Turns below zero; will continue with negative gradient. Returning half line.", quote=FALSE)}
        return (Intervals(c(0, Inf)))
      }
      if (func.at.turning.point.sign > 0 & func.at.0.sign > 0 & func.at.turning.point == func.at.0){
        if (verbose){print ("Turns at zero; will continue with negative gradient.", quote=FALSE)}
        return (interval.if.derivative.always.negative(q, r, s, start=0, verbose=verbose))
      }
      if (func.at.turning.point.sign > 0 & func.at.0.sign > 0 & func.at.turning.point < func.at.0){
        if (verbose){print ("Turns above zero; will continue with positive gradient. Returning empty set.", quote=FALSE)}
        return (Intervals(c(0, 0)))
      }
      if(func.at.turning.point.sign < 0 & func.at.0.sign < 0 & func.at.turning.point < func.at.0){ # will eventually go positive
        return (interval.if.derivative.always.positive(q, r, s, start=q^2/(r^2-q^2), verbose=verbose))
      }
      if(func.at.turning.point.sign > 0 & func.at.0.sign > 0 & func.at.turning.point > func.at.0){ # will eventually go negative
        return (interval.if.derivative.always.negative(q, r, s, start=q^2/(r^2-q^2), verbose=verbose))
      }
      if(func.at.turning.point.sign > 0 & func.at.0.sign <= 0){ # negative over two intervals
        if (func.at.0.sign == 0){
          left.star = 0
        }else{
          left.star = uniroot (f=function(x){truncation.region.function(x, q, r, s)}, lower=0, upper=q^2/(r^2-q^2))$root
        }
        if (verbose){print(paste('Found root left of turning point: ', left.star, sep=''), quote=FALSE)}
        right.star = find.root (q, r, s, start=q^2/(r^2-q^2), verbose=verbose)
        return (Intervals(rbind(c(0, left.star), c(right.star, Inf))))
      }
      if(func.at.turning.point.sign < 0 & func.at.0.sign >= 0){ # just one negative interval, straddling turning point
        if (func.at.0.sign == 0){
          left.star = 0
        }else{
          left.star = uniroot (f=function(x){truncation.region.function(x, q, r, s)}, lower=0, upper=q^2/(r^2-q^2))$root
        }
        if (verbose){print(paste('Found root left of turning point: ', left.star, sep=''), quote=FALSE)}
        right.star = find.root (q, r, s, start=q^2/(r^2-q^2), verbose=verbose)
        return (Intervals(c(left.star, right.star)))
      }
    }
  }
}
