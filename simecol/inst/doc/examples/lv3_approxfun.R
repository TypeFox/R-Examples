
lv3 <- new("odeModel",
  main = function(time, init, parms, inputs) {
    s.in <- inputs(time)
    with(as.list(c(init, parms)),{
      ds <- s.in  - b*s*p + g*k
      dp <- c*s*p - d*k*p
      dk <- e*p*k - f*k
      list(c(ds, dp, dk), s.in = s.in)
    })
  },
  parms = c(b = 0.1, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0),
  times  = c(from = 0, to = 200, by = .1),
  inputs = approxfun(
    c(0,   99, 100,  101, 200), # time
    c(0.1, 0.1, 0.5, 0.1, 0.1), # s.in
    rule = 2
  ),
  init = c(s = 1, p = 1, k = 1), # substrate, producer, consumer
  solver = "lsoda"
)
