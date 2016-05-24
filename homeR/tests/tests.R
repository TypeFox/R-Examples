pmv.test <- function() {
  vote <- pmv(clo=1.0,
              met=1.2,
              air.temp=19,
              saturation=40)
  stopifnot(!is.null(vote))
  stopifnot(vote > -0.6 && vote < -0.5)
}

pmv.test.vectors <- function() {
  votes <- pmv(clo=1.0,
               met=1.2,
               air.temp=c(19,30),
               sat=40)
  stopifnot(length(votes)==2)
  stopifnot(all.equal(votes,
                      c(pmv(clo=1.0,
                            met=1.2,
                            air.temp=19,
                            sat=40),
                        pmv(clo=1.0,
                            met=1.2,
                            air.temp=30,
                            sat=40)),
                      tolerance=1e-4))
                        
}

run.tests <- function() {
  pmv.test()
  pmv.test.vectors()
}

library(homeR)
run.tests()
