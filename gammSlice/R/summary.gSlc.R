summary.gSlc <-
function(object,pages = 1,...) {

   #Obtain information from object

   formula.infor <- object$formulaInfor
   xList <- list(t(object$summ))
   random.factor <- object$randomFactor
   num.smooth.funcs <- length(formula.infor$svar)

   app.names <- function(ind, vec)
    lapply(ind,
           function(i, vec) substitute(vec[i],
                                       list(i=i, vec=vec)),
           vec=as.name(vec)) 

   linear.beta.names <- app.names(formula.infor$lvar,"beta")
   smooth.edf.names <- app.names(formula.infor$svar,"edf")
   if (num.smooth.funcs > 0){
       summ.names <- c(linear.beta.names, smooth.edf.names)
   } else {summ.names <- c(expression(beta[0]),linear.beta.names)}

   if(random.factor) {summ.names <- c(summ.names, expression(sigma^2)) }

   numrows <- length(summ.names)
   numrowspage <- (numrows-1) %/% pages + 1  
   pages <- (numrows-1)%/% numrowspage + 1 
   
   xList.i <- xList
   ihigh=0

   for (i in 1: pages) {
       ilow   <- ihigh + 1
       ihigh  <- min((ilow - 1 + numrowspage), numrows)  

       xList.i[[1]]   <- xList[[1]][ ,ilow:ihigh ]
       summ.names.i   <- summ.names[ilow : ihigh]
       summMCMC(xList = xList.i, parNames = summ.names.i)
       par(ask = TRUE)
   }
   par(ask = FALSE)
}
