`D1Amat` <-
function(J, filter.number = 10., family = "DaubLeAsymm", tol = 
9.9999999999999977e-101, verbose = FALSE)
{
if(verbose == TRUE)
   cat("Computing ipndacw\n")
now <- proc.time()[1.:2.]
if(J >= 0.)
   stop("J must be negative integer")
if(J - round(J) != 0.)
#
   stop("J must be an integer")
#
#See if matrix already exists. If so, return it
#
rmnorig <- rmname(J = J, filter.number = filter.number, family = family)
rm.there <- rmget(requestJ =  - J, filter.number = filter.number, 
family = family)
#
#
#See if partially computed matrix exists. If so, use it.
#
if(!is.null(rm.there)) {
   if(verbose == TRUE)
      cat("Returning precomputed version: using ", rm.there,"\n")
   speed <- proc.time()[1.:2.] - now
   if(verbose == TRUE)
      cat("Took ", sum(speed), " seconds\n")
   rmnexists <- rmname(J =  - rm.there, filter.number = 
   filter.number, family = family)
   tmp <- get(rmnexists,envir=DWEnv)[1.:( - J), 1.:( - J)]
   assign(rmnorig, tmp, envir=DWEnv)
   return(tmp)
}
#
#
#Otherwise have to compute whole matrix
#
if(J != -1.) {
   for(j in (1. + J):(-1.)) {
      rmn <- rmname(J = j, filter.number = filter.number,
      family = family)
         if(exists(rmn)) {
            if(verbose == TRUE) {
               cat("Partial matrix: ", rmn, " exists (")
               cat(paste(round(100. - (100. * (j *j))/(J * J), digits = 1.),
               "% left to do)\n", sep = ""))
            }
         fmat <- rep(0., J * J)
         H <- filter.select(filter.number = 
         filter.number, family = family)$H
         error <- 0.
         answer <- .C("rainmatPARTIAL",
         J = as.integer( - J),
         j = as.integer( - j),
         H = as.double(H),
         LengthH = as.integer(length(H)),
         fmat = as.double(fmat),
         tol = as.double(tol),
         error = as.integer(error), PACKAGE = "LS2W")
         if(answer$error != 0.)
           stop(paste("Error code was ", answer$error))
         m <- matrix(answer$fmat, nrow =  - J)
         m[1.:( - j), 1.:( - j)] <- get(rmn,envir=DWEnv)
         nm <- as.character(-1.:J)
         dimnames(m) <- list(nm, nm)
         speed <- proc.time()[1.:2.] - now
         if(verbose == TRUE)
            cat("Took ", sum(speed), " seconds\n")
         assign(rmnorig, m, envir=DWEnv)
         return(m)
         }
      }
   }
fmat <- rep(0., J * J)
H <- filter.select(filter.number = filter.number, family = family)$
H
error <- 0.
answer <- .C("rainmatPARENT",
J = as.integer( - J),
H = as.double(H),
LengthH = as.integer(length(H)),
fmat = as.double(fmat),
tol = as.double(tol),
error = as.integer(error),
PACKAGE = "LS2W")
if(answer$error != 0.)
stop(paste("Error code was ", answer$error))
speed <- proc.time()[1.:2.] - now
if(verbose == TRUE)
cat("Took ", sum(speed), " seconds\n")
m <- matrix(answer$fmat, nrow =  - J)
nm <- as.character(-1.:J)
dimnames(m) <- list(nm, nm)
assign(rmnorig, m, pos = DWEnv)
m
}

