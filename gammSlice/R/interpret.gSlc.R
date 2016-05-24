### This code is copied partly and then modified from the function
### interpret.gam of the package mgcv (Simon Wood).

interpret.gSlc <- function(gf){ 
    
   p.env <- environment(gf) # enviroment of formula
   tf <- terms.formula(gf,specials= "s") # specials attribute indicates which terms are smooth and which term is random factor.
  
   terms <- attr(tf,"term.labels") # labels of the model terms
   nt <- length(terms) # how many terms?
  
   if (attr(tf,"response")>0){
      response <- as.character(attr(tf,"variables")[2])
      pf <- rf <- paste(response,"~",sep="")
   } else pf <- rf <- "~"

   sp <- attr(tf,"specials")$s # array of indices of smooth terms
  
  
   vtab <- attr(tf,"factors") # cross tabulation of vars to terms
  
   if (length(sp)>0)  for (i in 1: length(sp)) {
	ind <- (1:nt)[as.logical(vtab[sp[i],])]
      sp[i] <- ind # the term that smooth relates to
   } 
  
   k <- kp <- 1
   len.sp <- length(sp)
   smooth.spec <- list()
   smooth.var <- c()
   linear.var <- c()
   num.basis <- c()   
   
   av <- rep("",0)

   if (nt)
      for (i in 1:nt) # work through all terms
      { if (k<= len.sp && sp[k] == i)
        { str.i <- terms[i]
          if (length(grep("nBasis", str.i)) > 0) {
              str.i <- sub("nBasis", "k", str.i)
          }
          st <- eval(parse(text= str.i),envir = p.env)
	     smooth.spec[[k]] <- st
          smooth.var <- c(smooth.var,smooth.spec[[k]]$term)
          num.basis <- c(num.basis , smooth.spec[[k]]$bs.dim) 
          k <- k + 1 # counts smooth terms    
        } else # paramatric 
	  { if (kp>1) pf <- paste(pf,"+",terms[i],sep="") # add to parametric formula
          else pf <- paste(pf,terms[i],sep="")
          av[kp] <- terms[i]
          kp <- kp + 1 # counts parametric terms
        } 
      }
  
  if (length(av)) {
   
        pred.formula <- as.formula(paste("~", paste(av,collapse = "+")))
        pav <- all.vars(pred.formula)
        pred.formula <- reformulate(pav)

  } else {pred.formula <- ~1; pav <- rep("",0)}
  
  linear.var <- pav
  
  var.list <- c(linear.var, smooth.var)
  ret<-list(smooth.spec=smooth.spec, svar =smooth.var, lvar=linear.var, varlist = var.list, response = response,nbas = num.basis)
  class(ret)<-"gamm.formula"
  return(ret)
}
