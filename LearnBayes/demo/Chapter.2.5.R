#####################################
# Section 2.5 Using a Histogram Prior
#####################################

 library(LearnBayes)

 midpt = seq(0.05, 0.95, by = 0.1)
 prior = c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
 prior = prior/sum(prior)
 
 curve(histprior(x,midpt,prior), from=0, to=1,
   ylab="Prior density",ylim=c(0,.3))

 s = 11
 f = 16

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 curve(histprior(x,midpt,prior) * dbeta(x,s+1,f+1), 
   from=0, to=1, ylab="Posterior density")

S=readline(prompt="Type  <Return>   to continue : ")

 p = seq(0, 1, length=500)
 post = histprior(p, midpt, prior) *
        dbeta(p, s+1, f+1)
 post = post/sum(post)
 ps = sample(p, replace = TRUE, prob = post)

 windows()
 hist(ps, xlab="p", main="")

