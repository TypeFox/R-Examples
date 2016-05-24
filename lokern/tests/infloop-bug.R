### One of the examples that gave an infinite loop:

require(lokern)

ep <- 2e-6
x <- c(rep(118 + ep, 2), 150, 207, 240, 240, 240,
       263-ep, 270, 275-ep, 290-ep)
##x <- c(rep(-32+ep, 2), 0, 57, 90, 90, 90, 113-ep, 120, 125-ep, 140-ep)
y <- c(19, 27, 23, 27, 4, 15, 21, 33, 5, 15, 0)

plot(x,y)

g3  <- glkerns(x, y, x.inOut=FALSE) # used to "freeze up" , i.e. Infinite loop --- now fine
lok <- lokerns(x, y, x.inOut=FALSE) # ditto.
## 'x.inOut=FALSE' - remain back-compatible

with(g3,  lines(est ~ x.out, col = 2))
with(lok, lines(est ~ x.out, col = 3, lwd=2))

g.e <- c(26.370, 26.088, 25.807, 25.526, 25.244, 24.963, 24.682, 24.400, 24.119,
         23.838, 23.556, 23.275, 22.994, 22.712, 22.431, 22.109, 21.552, 20.801,
         20.050, 19.299, 18.548, 17.797, 17.046, 16.296, 15.545, 14.794, 14.043,
         13.292, 12.541, 11.790)
l.e <- c(26.407, 26.120, 25.831, 25.542, 25.253, 24.966, 24.681, 24.399, 24.122,
         23.849, 23.580, 23.315, 23.054, 22.796, 22.537, 22.246, 21.997, 21.669,
         21.204, 20.523, 19.614, 18.656, 17.659, 16.631, 15.579, 14.508, 13.421,
         12.324, 11.221, 10.115)
ii <- 1+10*(0:29)
stopifnot(all.equal( g3$est[ii], g.e, tol = 4e-5),
          all.equal(lok$est[ii], l.e, tol = 4e-5))

if(!interactive()) warnings()

cat('Time elapsed: ', proc.time(),'\n') # "stats"

