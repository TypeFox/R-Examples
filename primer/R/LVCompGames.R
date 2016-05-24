LVCompGames <- function(Alpha=NULL, r1=.1, r2=.1, num=20, time=10, step=1){

### INPUT THE COMPETITION COEFFICIENTS a, and r
#
r1 <- r1; r2 <- r2

### INPUT THE NUMBER OF INITIAL POINTS and TIME, and intervals
if(is.null(Alpha)) a <- matrix(10^runif(4, min=-3, max=-1), ncol=2) else a <- Alpha

#####################################################################
# set up graph with isocline 1

n2iso <- expression(1/a[2,2] - (a[2,1]/a[2,2])*x2)
x2 <- seq(0,1/a[2,1], length=2)

plot(x2, eval(n2iso), type='n', axes=FALSE,
      ylim=c(0, max(1/a[2,2], 1/a[1,2])), xlim=c(0, max(1/a[1,1], 1/a[2,1])),
      ylab=expression("N"[2]), xlab=expression("N"[1]))

# add pretty axes and box

axis(1, c(0, 1/a[1,1], 1/a[2,1]),
          c(0, expression(1/alpha[11]), NA) )
          
axis(1, 1/a[2,1], expression(1/alpha[21]),
     line=1, tick=FALSE)
         
axis(2, c(0, 1/a[2,2], 1/a[1,2]),
          c(0, expression(1/alpha[22]), NA)
          )
axis(2, 1/a[1,2],expression(1/alpha[12]), tick=FALSE, line=1 )
#legend('right', c("N2","N1"), lty=1:2, lwd=2, bty="n")
box()

# add isoclines

lines(x2, eval(n2iso), lty=1)
  
n1iso <- expression(1/a[1,2] - (a[1,1]/a[1,2])*x1)
x1 <- seq(0, 1/a[1,1], length=2)
lines(x1, eval(n1iso), lty=2)

###
### Which isocline is for which species?
### What is the maximum size of each population?
###############################################################

# make the initial (random) points, and a vector of parameters.
maxN2 <- 1/min(a[,2])
maxN1 <- 1/min(a[,1])

init.N1 <- runif(num, min=0, max=maxN1)
init.N2 <- runif(num, min=0, max=maxN2)

paramsST <- c(r1=r1, r2=r2,
            a11 = a[1,1], a12 = a[1,2],
            a21 = a[2,1], a22 = a[2,2])
times <- seq(1, time, by=step)
# for each initial point -- add the points to the graph
for( i in 1:num ) {points(init.N1[i], init.N2[i], pch=19)}

### Predict where Each population will go.
###########################################################################
# on.exit(devAskNewPage(oask), TRUE)
# Now run the solver and add the trajectories
# add isoclines and points and lines

#ask(message="Is there a two-species equilibrium, and if so, is it stable? Hit Return when you have decided.")

readline("Press <Enter> to continue")
for( i in 1:num ) {Ns <- ode(c(init.N1[i], init.N2[i]), times,
                                      lvcomp2, paramsST )
                   lines(Ns[,2],Ns[,3]) }

}
