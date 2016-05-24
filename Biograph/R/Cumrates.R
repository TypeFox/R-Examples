Cumrates <-
function (irate,Bdata)
{ # Calculate MSLT from original data (Bdata) and radix
  # irate = 1 estimate Nelson-Aalenniter5
  # irate = 2 estimate occurrence-exposure rates
  # irate = 3 both
  # NOTE: if there is an absorbing state and st.absorb=NULL, then
  # calling mvna produces and error: "there are undefined transitions". 
  # In that case, remove the transition from aborbing state to censoring.
  if (missing(irate)) irate=3
  z<- check.par (Bdata)
  namstates <- attr(Bdata,"param")$namstates
  numstates <- length (namstates)
  namage <- attr(Bdata,"param")$namage
iagehigh <- attr(Bdata,"param")$iagehigh
# require (mvna)
print (". . . . . .  Removing intrastate transitions . . . . . ")
locpat <- locpath(Bdata)
removed <- Remove.intrastate(Bdata)
Bdata <- removed
z<- StateSpace (Bdata)
st.absorb <- z$absorbstates
# ========  Step 1: ESTIMATE TRANSITION RATES ========
print (". . . . . . Estimating rates . . . . .")
Lambda <- 0
Lambda.oe <- 0 
M.oe <- 0
namstates2 <- vector (mode="numeric",length=length(namstates))
for (i in 1:length(namstates))
{ namstates2[i] <- grep(namstates[i],namstates)} 
namstates2 <- as.character(1:length(namstates))
niter5 <- ifelse (irate==2,1,3)
if (irate %in% c(1,3)) # ===  mvna estimates transition rates  ====
  { print ("  ")
  	print ("Cumrates: calls function Biograph.mvna  . . . ")
  	# convert from months to years
  	Bdata <- date_b (Bdata=Bdata,format.out="age",covs=NULL)
  	#Bdata <- CMC.ages(Bdata)  # NO COVARIATES INCLUDED
  	print (Bdata[1:5,])
  	Dmvna <- Biograph.mvna (Bdata)
   print ("Cumrates 88")
  	# see GLHS_mvna.r
  	D2 <- Dmvna$D
    tra <- Dmvna$par$trans_possible
    cens <- Dmvna$cens
    #D3 <- data.frame(cbind(D2$id,D2$from,as.character(D2$to),D2$entry,D2$exit,D2$time))
    D3 <- data.frame(id=D2$id,from=D2$from, 
                              to=D2$to,
                              entry=round(D2$entry,2),
                              exit=round(D2$exit,2))
    D3$id <- as.numeric(D3$id)
    D3$entry <- as.numeric(D3$entry)
    D3$exit <- as.numeric(D3$exit)
    D3$exit <- ifelse (D3$exit<=D3$entry,D3$entry+1,D3$exit) # CORRECT if CMC exit = CMC entry
# ------------------------------------------------
    print (". . . . . . . . . . . . ")
    print (". . . . Running mvna . . . . . . ")
    # Absorbing states determined in StateSpace (letter)
    if (is.null(st.absorb)) zz2 <- D3 else 
         {  st.absorb2 <- which (namstates==st.absorb)
         	zz2 <- subset (D3,D3$from!=st.absorb) } 
         	       # NOTE: st.absorb used because name of state is given (nog mumber)
    na <- mvna(data=zz2,state.names=namstates,tra=attr(Dmvna$D,"param")$trans_possible,cens.name=Dmvna$cens)
    cumh <- predict (na,times=seq(0,iagehigh,by=1))
   #	print ("cumrates test78")   
# see MSLT.mvna.r
  #  Lambda = cumulative hazard by age, destination, origin
   	 nage = nrow(cumh[[1]])
   Lambda <- array (NA,dim=c(nage,numstates,numstates,niter5))
   if (length(namage)!=nage) 
      {  print ("WARNING")
      	 print  ("Names of ages groups (namage) not consistent with number of age groups")
      	 print ("Names are inferred. Please check the result. ")
      	 namage <- 0:(nage-1)}
   dimnames(Lambda) <- list(age=namage[1:nage],destination=namstates,origin=namstates,variant=c("Expected","Upper","Lower"))
   
   	for (i in 1:attr(removed,"param")$ntrans)
   	{ des <- as.numeric(as.character(attr(removed,"param")$transitions$DES[i]))
   	  or <- as.numeric(as.character(attr(removed,"param")$transitions$OR[i]))
   	  Lambda[,des,or,1] <- cumh[[i]]$na
      Lambda[,des,or,2] <- cumh[[i]]$upper
      Lambda[,des,or,3] <- cumh[[i]]$lower
   	}

    for (iter in 1:3)
    { for (ix in 1:nage)
      {  diag(Lambda[ix,,,iter]) <- - apply(Lambda[ix,,,iter],2,function(x) sum(x,na.rm=TRUE)) # column sum
        for (i in 1:numstates) {for (j in 1:numstates) 
        	 Lambda[ix,j,i,iter] <- ifelse (is.na(Lambda[ix,j,i,iter]),0,Lambda[ix,j,i,iter]) }
      }
    } 
  # Compute age-specific rates
  astr <- array (NA,dim=c(nage,numstates,numstates,niter5))
  dimnames (astr) <- dimnames (Lambda)
  astr[1,,,] <- Lambda[1,,,]
  for (ix in 2:nage)
  { astr[ix-1,,,] <- Lambda[ix,,,]-Lambda[(ix-1),,,]
  }
  astr[nage,,,] <- 0
  }
  if (irate %in% c(2,3))
  { # === occurrence-exposure rates ====
  	print ("Computing occurrence-exposure rates",quote=FALSE)
  	print ("Running Sequences.ind",quote=FALSE)  
  	ist <- Sequences.ind (Bdata$path,namstates)
  	print ("Running Occup",quote=FALSE)    
  	occup <- Occup (Bdata)
  	print ("Running Trans",quote=FALSE)  
  	trans <- Trans (Bdata)
  	print ("Running RateTable",quote=FALSE)  
    ratetable <- RateTable (Bdata,occup,trans)
    pc_ac <- 2    # age-cohort rates  ( 1 - age-cohort rates)
    print ("Running Rates",quote=FALSE)  
    M <- Rates.ac(ratetable$Stable) 
# NOTE: cum oe rates: Lambda.oe[14,,] = sum of rates to and including age 12
   Lambda.oe<- M$Mcum
   M.oe <- M$M
   # compare Lambda[51,,] and Lambda.oe[50,,]
  }
  
if (!exists("astr")) astr <- NA
if (!exists("Lambda")) Lambda <- NA
if (!exists("cumh")) cumh <- NA
if (!exists("Lambda.oe")) Lambda.oe <- NA
if (!exists("M.oe")) M.oe <- NA

 cum <- list(D=removed,
              irate=irate,
              NeAa = Lambda,
              predicted=cumh,
              astr = astr,
              oeCum = Lambda.oe,
              oe=M.oe)
 class(cum) <- 'cumrates'

 return (cum)

#   attr(cum$D$D,"param")$namstates
}
