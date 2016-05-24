CounterClock <- "Cntclck"; NoneClock <- "noneclck"; Clockwise <- "clckws"

IsCounterCl2 <- function( U, V, ref ) {
  return  ((U == V) | ((U < V) == (abs( U - V ) <= ref/2)))
} ## end  IsCounterCl2

IsCounterCl3 <- function( U, V, W, ref ) {
  if ((U == V) | (V == W)) IsCounterCl2( U, W, ref ) else {
    if (U == W) IsCounterCl2( U, V, ref ) else {
      if (U < V) ((V < W) | (W < U)) else  ## U<V<W OR V<W<U OR W<U<V *)
        ((V < W) & (W < U))
    }
  }	## end  ## If *)
} ## end  IsCounterCl3

ClockSense2 <- function( U, V, ref ) {
  if (U == V) NoneClock else {
    if ((U < V) == (abs( U - V ) <= ref*0.5)) CounterClock else
      Clockwise
  }
} ## end  ClockSense2

ClockSense3 <- function( U, V, W, ref ) {
  if((U == V) | (V == W)) ClockSense2( U, W, ref ) else {
    if (U == W) NoneClock else {
      if (U < V) {           ## U<V<W OR V<W<U OR W<U<V *)
        if ((V < W) | (W < U)) CounterClock else
          Clockwise
      } else {  ##if*)
      if ((V < W) & (W < U)) CounterClock else
          Clockwise
      }
    }
  }
} ## end  ClockSense3

