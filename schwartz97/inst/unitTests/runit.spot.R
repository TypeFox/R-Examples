test.r.sim <- function(){
  obj <-  schwartz2f(s0 = 60, delta0 = 0.2, alpha = 0.1)

  traj <- sapply(1:2000, function(dummy, obj)simstate(5, 0.3, obj)[5,],
                 obj = obj)
  end.points <- rstate(2000, 0.3, obj)

  abs.diff <- abs(rowMeans(traj) - colMeans(rstate(2000, 0.3, obj)))
  checkTrue(abs.diff[1] < 1,  "spot: simstate must converge to rstate")
  checkTrue(abs.diff[2] < .05,  "convenience yield: simstate must converge to rstate")
}

test.mean <- function(){
##  library(schwartz97)
##  library(RUnit)

  obj <- schwartz2f(mu = 0, alpha = 0, delta0 = 0, s0 = 1)
  checkTrue(mean(obj, 1.5)[1] < 1, "Mean must be < 1!")
  checkTrue(mean(obj, 1.5)[2] == 0, "Mean must be == 0!")

  obj <- schwartz2f(mu = 0, alpha = -1, delta0 = 0, s0 = 1)
  checkTrue(mean(obj, 0.5)[1] > 1, "Mean must be > 1!")
  checkTrue(mean(obj, 0.5)[2] < 0, "Mean must be < 0!")

  checkTrue(mean(schwartz2f(alpha = 0.2))[1] <
            mean(schwartz2f(alpha = 0.1))[1],
            "First mean must be smaller than second!")

  checkTrue(mean(schwartz2f(sigmaE = 0.2))[1] <
            mean(schwartz2f(sigmaE = 0.1))[1],
            "First mean must be smaller than second!")

  checkTrue(mean(schwartz2f(mu = 0.1))[1] <
            mean(schwartz2f(mu = 0.2))[1],
            "First mean must be smaller than second!")

}

test.vcov <- function(){

  obj.1 <- schwartz2f(kappa = 1)
  obj.2 <- schwartz2f(kappa = 2)

  checkTrue(vcov(obj.1)[2,2] > vcov(obj.2)[2,2],
            "Faster decay rate of conv. yield should result in less volatility of the conv.yield!")


  obj.1 <- schwartz2f(kappa = 1, rho = 0)
  obj.2 <- schwartz2f(kappa = 2, rho = 0)

  checkTrue(vcov(obj.1)[1,1] > vcov(obj.2)[1,1],
            "Faster decay rate of conv. yield should result in less volatility in the spot!")

}




