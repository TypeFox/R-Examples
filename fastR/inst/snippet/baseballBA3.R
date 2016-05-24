# suppress warnings from log(0) 
oldopt <- options(warn=-1)          
summary(nlmax(loglik,p=c(1,1),x=ba))
# get just the mle
nlmax(loglik,p=c(1,1),x=ba)$estimate
options(oldopt)                         # reset options
