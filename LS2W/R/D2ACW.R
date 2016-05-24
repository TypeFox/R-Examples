"D2ACW" <-
function(J, filter.number = 1., family = "DaubExPhase", switch = "direction",
tol = 1e-100, OPLENGTH = 2000., verbose = FALSE)
{
#
# Program to compute the 2-D discrete autocorrelation wavelets.
# This program basically calls PsiJ and PhiJ and goes from there.
#
#
# 
# Some quick checks first ...
#
if(verbose) {
   cat("Computing The two-dimensional (discrete) autocorrelation coefficients: \n")
   cat("\n")
}
if(switch != "direction" && switch != "level")
   stop("\nSorry, switch can only take the values direction and level!\n")
#  Two structuring options are available: we deal with my preferred one first:
#
if(J >= 0.) #
   stop("\n J must be a negaitve integer!!!\n")
#
#
# Now for the other method of structuring the list
#
if(switch == "direction") {
   now <- proc.time()[1.:2.]
#
# See if the list already exists.
#
   Wavorig2 <- # 
   Psi2Dname(J = J, filter.number = filter.number, family = family,
   switch = switch)
#
# If the list doesn't exist ... go make one!
#
   if(exists(Wavorig2)) {
      if(verbose)
         cat("Returning precomputed version \n")
      speed <- proc.time()[1.:2.] - now
      if(verbose)
         cat("Took ", sum(speed), " seconds \n")
      return(get(Wavorig2,envir=DWEnv))
   }
   temp1 <- PsiJ(J = J, filter.number = filter.number, family = 
   family, tol = tol, OPLENGTH = OPLENGTH, verbose = verbose)
   temp2 <- PhiJ(J = J, filter.number = filter.number, family = 
   family, tol = tol, OPLENGTH = OPLENGTH, verbose = 
   verbose)
   mat <- vector("list", -3. * J)
   for(j in 1.:( - J)) {
#
#
# Vertical autocorrelation wavelets:
#
#
# Horizontal autocorrelation wavelets:
#
      mat[[j]] <- #
      outer(temp1[[j]], temp2[[j]])
#
# Diagonal autocorrelation wavelets:
#
      mat[[ - J + j]] <- #
      outer(temp2[[j]], temp1[[j]])
      mat[[-2. * J + j]] <- outer(temp1[[j]], temp1[[j]])
   }
   assign(Wavorig2, mat, envir=DWEnv)
}
if(switch == "level") {
   now <- proc.time()[1.:2.]
#
# See if it already exists.
#
   Wavorig2 <- # 
   Psi2Dname(J = J, filter.number = filter.number, family = family,
   switch = switch)
   if(exists(Wavorig2)) {
      if(verbose)
         cat("Returning precomputed version \n") 
      speed <- proc.time()[1.:2.] - now
      if(verbose)
         cat("Took ", sum(speed), " seconds \n")
      return(get(Wavorig2,envir=DWEnv))
   }
   temp1 <- PsiJ(J = J, filter.number = filter.number, family = 
   family, tol = tol, OPLENGTH = OPLENGTH, verbose = 
   verbose)
   temp2 <- PhiJ(J = J, filter.number = filter.number, family = 
   family, tol = tol, OPLENGTH = OPLENGTH, verbose = 
   verbose)
   mat <- vector("list", -3. * J)
   for(j in 1.:( - J)) {
#
# Vertical autocorrelation wavelets:
#
#
# Horizontal autocorrelation wavelets:
#
      mat[[3. * j - 2.]] <- #
      outer(temp1[[j]], temp2[[j]])
#
# Diagonal autocorrelation wavelets:
#
      mat[[3. * j - 1.]] <- #
      outer(temp2[[j]], temp1[[j]])
      mat[[3. * j]] <- outer(temp1[[j]], temp1[[j]])
   }
   assign(Wavorig2, mat, envir=DWEnv)
}
else if(switch != "level" && switch != "direction") {
   stop("FOOL - switch must either be direction or level!!!")
}
mat
}

