DO.HALFSPACE <-
function()
  {
    k = 1e-6
    T0 = 1000
    x = seq(0,20, length=1000)
    y =  heat.sol(x, T0, k, seq(from=0,to=1000*24*3600, by=100*24*3600))
  }

