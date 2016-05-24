`ss.aipe.sem.path` <-
function(model, Sigma, desired.width, which.path, conf.level=.95, assurance=NULL, ...){

if(!requireNamespace("sem", quietly = TRUE)) stop("The package 'sem' is needed; please install the package and try again.")
  
Sigma[upper.tri(Sigma)]<-0 
N<-1000000
result.sem <- sem::sem(model, Sigma, N, ...)
p <- result.sem$n
H <- result.sem$vcov * N

j <- which(which.path==names(result.sem$coeff))
if(j==0) stop("The path of interest is not included in the model; make sure the parameter names are specified correctly")
h.jj <- H[j,j]
  
omega<- desired.width
alpha<- 1- conf.level
z <- qnorm(1-alpha/2)

if(is.null(assurance)){
    N <- 4*z^2* h.jj /omega^2
    N <- ceiling(N)
    } 
    
if(!is.null(assurance)){
    N0 <- ss.aipe.sem.path(model=model, Sigma=Sigma, desired.width=desired.width, 
    which.path=which.path, conf.level=conf.level, assurance=NULL)$sample.size
    gamma<- assurance
    signal <- TRUE
    N<- N0
    #print(cat(N, N0))
    while(signal){
        c<- N0*qchisq(p=gamma, df=N-1)
        N.1 <- (1+sqrt(1+4*c))/2
        if (abs(N.1-N)>1e-5) N<- N.1
        else signal<-FALSE
        #print(cat(N, N.1, "\n"))
        }# end of while()     
    }# end of if(!is.null(assurance))

N <- ceiling(N)
result<- list()
result$parameters <- names(result.sem$coeff)
result$path.index <- j
result$sample.size <-N
result$obs.vars <- result.sem$var.names[1:p]
result$var.theta.j <- h.jj/N
return(result)
}# end of function()

