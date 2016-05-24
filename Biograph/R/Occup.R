Occup <-
function (Bdata) {
# A. State occupancies at exact ages 
# Input  = Bdata.dat
#          iagelow, iagehigh (global variables, defined in MAIN programme)
#          nsample, numstates, namstates (global variables, produced in StateSpace)
#          ages  (global variable, produced in AgeTrans )
#          iagelow = lowest age
#          iagehigh = highest age
#          ageprint: selection of ages to be printed (1 = age 0)
#          subjects: list of subjects for which data should be printed
# Output:
#   a. state_occup : state occupancies at consecutive ages
#   b.  tsjt_p : sojourn time 
#  ageentry[i]
#  agecens[i]
# METHOD:
#    With each element of agelist(continuous stay in state) is associated a state (st[ix])
#
# B. Sojourn times by age and state are grouped over individuals => tsjt
# st_agelist[i,ix] = state occupied at exact age ix by subject i
# st_age_1[i,ixx] = state occupied at exact age ixx by subject i (ixx = one-year age groups)
#        After censoring: state occupied is numstates + 1 
# sjt[ixx,j]  = number of years lived in state j during age interval (ix,ix+1) by given individual
# sjt_age_1[i,ixx,j]  = months spent in state j by individual i between ages ixx and ixx+1
# tsjt  = total sojourn time (entire population) by age and state
#     Note: tsjt[1:numstates+1,] with tsjt[numstates+1,] = person years after censoring
#            Hence: apply(tsjt_p[,-(numstates+1)],1,sum) is person years in observation window 
#                   For og98: apply(tsjt_p[,-6],1,sum)
# NOTE: tstate differs a little from TAB.OUT in Fortran. If indiv makes transition at exact age ix,
#    he occupies the OLD state in the fortran programme and the NEW state in the R programme
 z<- check.par (Bdata)
 iagehigh <- attr(Bdata,"param")$iagehigh
 iagelow <- attr(Bdata,"param")$iagelow
numstates <- attr(Bdata,"param")$numstates
nsample <- nrow(Bdata)
if (is.null(attr(Bdata,"format.date"))) stop ('in Occup: attr(Bdata,"format.date") missing.')
format.in <- attr(Bdata,"format.date")
nage <- iagehigh - iagelow + 1 
if (nage > 150) warning("Occup: Please check ages. The number of age groups may be numeric(0) or exceed 150. Check Attributes of data object (include Param) and run Parameters to be certain.")
if (nage > 150) warning("Occup: Please check ages. ") else print (paste ("The number of age groups is ",iagehigh - iagelow + 1 ,sep=""),quote=FALSE)
param <- Parameters(Bdata)
print (". . . .  Running Occup . . . . ")
#agetrans <- AgeTrans(Bdata) # at exit: timeunit is year!!
# z<- state_age (Bdata,age=iagelow:iagehigh,ID=Bdata$ID)
y <- state_time (Bdata,ID=Bdata$ID)


zzz <- y$state.n[,-1]
colnames(zzz)[numstates+1] <- "Censored"
class(zzz) <- "occup.S"
return (list ( state_occup = zzz,
               st_age_1 = y$state,
               sjt_age_1 = y$sjt_age_1,
               tsjt = y$tsjt[,2:(numstates+2)]))
}
