# hazard funtion for gamlss.families
# Mikis Stasinopoulos Bob Rigby and Vlasios Voudouris
#----------------------------------------------------------------------------------------
hazardFun <-function(family = "NO",  ...)
  {
   fname <- family
   if (mode(family) != "character" && mode(family) != "name")
   fname <- as.character(substitute(family))
 distype <- eval(call(family))$type
    dfun <- paste("d",fname,sep="")
    pfun <- paste("p",fname,sep="")
     pdf <- eval(parse(text=dfun))
     cdf <- eval(parse(text=pfun))
if ( distype == "Discrete")
   {
 fun <- function(x, log = FALSE, ...)
        {
        hfun <- pdf(x,...)/cdf(x, lower.tail=FALSE, ...)
        hfun <- ifelse(hfun>1, 1, hfun)
        hfun
        }
   }
else
   {     
 fun <- function(x, log = FALSE, ...)
        {
        #dfun <- ifelse(cdf(x,...)>1e-10, pdf(x,log = TRUE,...)-log(1-cdf(x,...)), pdf(x,log = TRUE,...)+cdf(x,...))# this is not working
        #dfun <- pdf(x,log = TRUE,...)-log(1-cdf(x,...)) # the original
        # dfun <- pdf(x,...)/(1-cdf(x,...)) # version 2 still with problem if pdf is very very small 
        # d1<- ifelse(pdf(x,...)<1e-15, NA, pdf(x,...)) # it works for the film data  but has to checked in general 
        # sofx<- ifelse(cdf(x,...)==1, NA, cdf(x,lower.tail=FALSE,...))
# old version hfun <- pdf(x,...)/ ifelse(cdf(x,...)==1, NA, cdf(x,lower.tail=FALSE,...))
          hfun <- pdf(x, log=TRUE, ...)- cdf(x, lower.tail=FALSE,log.p=TRUE,...)
          hfun <- exp(hfun)
        #dfun <- if (log == TRUE) dfun else exp(dfun)# part of the oroginal
        hfun
        }
   }
  fun
  }
#  Hazard=d/(1-p)
#  log(Hazard)= log(d) - log(1-p)
#For p very small,
#   log(1-p) ~ -p
#and so
#  log(Hazard) ~ log(d) + p
#----------------------------------------------------------------------------------------
gen.hazard <-function (family = "NO", ...)
{
#------------------------------------------
     fam  <- as.gamlss.family(family) # family 
    fname <- fam$family[[1]] 
  # dorfun <- paste("d",fname,sep="")
  # porfun <- paste("p",fname,sep="")
     hfun <- paste("h",fname,sep="")
  #   sfun <- paste("s",fname,sep="")
 #   eval(dummy <- survivalFun(family = fname, ...))
 #   eval(call("<-",as.name(sfun),dummy), envir=sys.frame(0))#
    eval(dummy <- hazardFun(par, family = fname,  ...))
    eval(call("<-",as.name(hfun),dummy), envir=sys.frame(0))#  
     cat("A hazard function for",  fname, "is generated \n", 
  "and saved under the name: ", "\n",paste(hfun), "\n")#
}
#---------------------------------------------------------------------------------------- 
