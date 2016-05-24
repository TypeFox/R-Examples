"computeBounds" <-
function(n,drift,alpha,phi,t,t2,OneOrTwoSidedBounds,whatSpendingFunctionIsUsed,boundsTruncation)
{
###########################################################################
######################### INITIALIZE VARIABLES ############################
###########################################################################
  probExit<-0 # probExit is a vector of exit probabilities
  probDifference<-0 # probDifference(i) = probExit(i)-probExit(i-1) where probExit is a vector of exit probabilities.
  standardDeviation<-0 #  the standard deviation of the process increment
  standardDeviationProcess<-0 #the standard deviation of the process.
  lowerBounds <-0 # lowerBounds is the vector of lower standardized boundaries
  upperBounds <-0 # upperBounds is the vector of upper standardized boundaries
  gridSize<-0.05 # the grid size for numerical integration by trapezoidal rule
  lastGrid<-0 # the grid of the joint density at the last analysis.
  numberOfIntegrationIntervalls<-0 #"integer" the number of intervals for numerical integration

  lowerIntegrationLimit<-0 # the vector of lower integration limits.
  upperIntegrationLimit<-0 # the vector of upper integration limits.
  noError<-TRUE

  # qnorm is the inverse standard normal cdf
  # pnorm is the standard normal cdf



###########################################################################
############################### START #####################################
###########################################################################

  ## Negative infinity is set to -8 by default
  negInf<-(-8)  #negInf is "negative infinity" for the program.

########## COMPUTE PROBABILITIES ACCORDING TO USE FUNCTION ################

  ## Therefore function 'alphaByUseFunction(...)' is called ###
  probExit <- alphaByUseFunction(whatSpendingFunctionIsUsed,n,alpha,phi,OneOrTwoSidedBounds,t)

  ## Compute 'probDifference' - the change in type I error to spent
  toleranceProbDiff <- 1.0e-13

  for (i in 1:n)
  {
    ##catch first loop
    if (i==1)
    {
      probDifference[1]<-probExit[1]
    }
    ##do the rest
    else
    {
      probDifference[i]<-probExit[i]-probExit[i-1]
    }

    ##Check type I error to spend
    if (probDifference[i]<0 || probDifference[i]>1)
    {
      probDifference[i]<-min(1,probDifference[i])
      probDifference[i]<-max(0,probDifference[i])
      cat("\n")
      cat(" Error in spending function at interim time:",i,"\n")
      cat(" Calculated probabilites are:",probExit,"\n")
      cat(" Calculated function is not increasing strictly or out of range!","\n")
      cat(" the differences intype I error spent between analyses","\n")
      cat(" at this point will be set to:",probDifference[i],"\n")

    }

    if (probDifference[i]<toleranceProbDiff)
    {
      cat("\n")
      cat(" Type I error spent too small at interim time:",i,"\n")
      cat(" Zero used as approximation for:","\n")
      print(probDifference[i],digits=22)
      cat("\n")
    }
   }#end <--*for*


   ###--- Calculate standard deviations of increments and process ---###
   ##function stdDeviations(...) returns (standardDeviation, standardDeviationProcess)
   stdDeviationVector<-stdDeviations( n, t2 )
   standardDeviation <- stdDeviationVector[[1]]
   standardDeviationProcess <- stdDeviationVector[[2]]

###########################################################################
################### BEGIN CALCULATING BOUNDARIES ##########################
###########################################################################

##--------------------------------------------------------##
##-- Direct calculations can be made for first analysis --##
##--------------------------------------------------------##

  ###Check type I error to spend
  if (probDifference[1]<0 || probDifference[1] >1)
  {
    print("",quote=FALSE)
    print(" Error in spending function - alpha is not in [0,1].",quote=FALSE)
    probDifference[1]<-min(1,probDifference[1])
    probDifference[1]<-max(0,probDifference[1])
  }

  ##Spending probability is zero (or less - so it was set to zero)
  if (probDifference[1]==0)
  {
    ##if boundaries should be truncated set upperBounds[1] accordingly
    ##instead of -negInf that is "negative infinity" for the programm
    upperBounds[1]<- -negInf
    if (upperBounds[1] > boundsTruncation)
    {
      upperBounds[1] <- boundsTruncation
      probDifference[1] <- OneOrTwoSidedBounds* ( 1-pnorm(upperBounds[1]) )
      probExit[1] <- probDifference[1]

      ##Difference between second and first exit probability
      if (n > 1)
      {
        probDifference[2] <- probExit[2] - probExit[1]
      }
    }
    upperIntegrationLimit[1] <- upperBounds[1]*standardDeviation[1]
  }#end <--*if*

  ##Spending probability is one (or more - so it was set to one)
  else if (probDifference[1]==1)
       {
         upperBounds[1] <- 0
         upperIntegrationLimit[1] <- upperBounds[1]*standardDeviation[1]
       }

       else
       {
         ##First bound based on normal distribution.
         upperBounds[1] <- qnorm( 1- (probDifference[1]/OneOrTwoSidedBounds) )

         ##check whether bound must be truncated
         if (upperBounds[1]>boundsTruncation)
         {
           upperBounds[1] <- boundsTruncation
           probDifference[1] <- OneOrTwoSidedBounds * ( 1 - pnorm(upperBounds[1]) )
           probExit[1] <- probDifference[1]

           ##Difference between second and first exit probability
           if (n > 1)
           {
             probDifference[2] <- probExit[2] - probExit[1]
             cat("probExit: ",probExit,"\n")
             cat("probDifference: ",probDifference[2],"\n")
           }
         }
         upperIntegrationLimit[1] <- upperBounds[1]*standardDeviation[1]
       }#end <--*else*
  ###end Checking type I error to spend


  ##Lower bound is either "negative infinity" (one-sided test)
  ##or -upperIntegrationLimit (two-sided test )
  if (OneOrTwoSidedBounds==1)
  {
    lowerBounds[1] <- negInf
    lowerIntegrationLimit[1] <- lowerBounds[1]*standardDeviation[1]
  }
  else
  {
    lowerBounds[1] <- -upperBounds[1]
    lowerIntegrationLimit[1] <- -upperIntegrationLimit[1]
  }



  ##Number of intervals for numerical integration
  numberOfIntegrationIntervalls[1] <- trunc( (upperIntegrationLimit[1]-lowerIntegrationLimit[1]) / (gridSize*standardDeviation[1]) )

##------------------------------------------------##
##-- Calculations for second and later analyses --##
##------------------------------------------------##


  ##catch n<2 that is user choosed just one interim analysis
  if(n<2)
  {
    #do nothing
  }

  else #go on with further analysis
  {
    for (i in 2:n)
    {
      ##Calculate joint density for use in next step
      if (i==2)
      {
        ##call function 'jointDensity' with parameter==1,
        ##for computing joint Density in 1st analysis - look also function jointDensity
        lastGrid <- jointDensity(1, lowerIntegrationLimit[1], upperIntegrationLimit[1], standardDeviation[1], numberOfIntegrationIntervalls, lastGrid)
        #                    |
        #               parameter==1
      }

      ###Check type I error to spend
      if (probDifference[i]<0 || probDifference[i] >1)
      {
        print("",quote=FALSE)
        print(" Error in spending function - alpha is not in [0,1].",quote=FALSE)
        probDifference[i]<-min(1,probDifference[i])
        probDifference[i]<-max(0,probDifference[i])
      }

      ##Spending probability is zero (or less than tolerance)##
      if (probDifference[i] < toleranceProbDiff)
      {
        upperBounds[i] <- -negInf

        ##check whether bound must be truncated
        if (upperBounds[i]>boundsTruncation)
        {
          upperBounds[i] <- boundsTruncation
          probDifference[i]<-OneOrTwoSidedBounds*tailProbability(upperBounds[i]*standardDeviationProcess[i],lastGrid,numberOfIntegrationIntervalls[i-1],lowerIntegrationLimit[i-1],upperIntegrationLimit[i-1],standardDeviation[i])
          probExit[i] <- probDifference[i] + probExit[i-1]

          if(n>i)
          {
            probDifference[i+1] <- probExit[i+1]-probExit[i]
          }
        }
        upperIntegrationLimit[i] <- upperBounds[i]*standardDeviationProcess[i]

      }#end <--*if*


      ##Spending probability is one (or more).##
      else if (probDifference[i]==1)
           {
             upperBounds[i] <- 0
             upperIntegrationLimit[i] <- upperBounds[i]*standardDeviation[i] ## that is <- 0
           }

      ##-------------------------------------------------##
      ##-- Using a search algorithm to find the bounds --##
      ##-------------------------------------------------##

           ##Bounds are found using a search starting at the bound from
           ##the previous analysis
           else
           {
             upperIntegrationLimit[i] <- searchForBound(lastGrid,numberOfIntegrationIntervalls,i,gridSize, probDifference[i]/OneOrTwoSidedBounds, standardDeviation[i], lowerIntegrationLimit, upperIntegrationLimit,n)

             ##check if function searchForBound(...) worked correctly
             if (!is.numeric(upperIntegrationLimit[i]))
             {
               ##in this case something went wrong - function returned noError==FALSE
               print(" Error in function 'searchForBound' - analysis aborted!",quote=FALSE)
               noError <- FALSE
               break
             }
             else ##everything went ok - function returned a numeric
             {
               #standarize upper boundary
               upperBounds[i] <- upperIntegrationLimit[i]/standardDeviationProcess[i]
             }

             ##If a truncation point is used, check to see if it
             ##applies and recompute probabilities if necessary.
             if(upperBounds[i]>boundsTruncation)
             {
               upperBounds[i]<-boundsTruncation
               probDifference[i]<-OneOrTwoSidedBounds*tailProbability(upperBounds[i]*standardDeviationProcess[i],lastGrid,numberOfIntegrationIntervalls[i-1],lowerIntegrationLimit[i-1],upperIntegrationLimit[i-1],standardDeviation[i])

               if(n>i)
               {
                 probDifference[i+1] <- probExit[i+1]-probExit[i]
               }
             }
             upperIntegrationLimit[i] <- upperBounds[i]*standardDeviationProcess[i]

           }#end <--*else*


      ##Lower bound is either "negative infinity" (one-sided test)
      ##or -upperIntegrationLimit (two-sided test )
      if (OneOrTwoSidedBounds==1)
      {
        lowerIntegrationLimit[i] <- negInf*standardDeviationProcess[i]
        lowerBounds[i] <- negInf
      }
      else
      {
        lowerIntegrationLimit[i] <- -upperIntegrationLimit[i]
        lowerBounds[i] <- -upperBounds[i]
      }

      ##Number of intervals for numerical integration.
      numberOfIntegrationIntervalls[i] <- trunc( (upperIntegrationLimit[i]-lowerIntegrationLimit[i]) / (gridSize*standardDeviation[i]) )

      ##Calculate joint density for use in next step.
      if (i!=n)
      {
        lastGrid <- jointDensity(i, lowerIntegrationLimit, upperIntegrationLimit, standardDeviation[i], numberOfIntegrationIntervalls, lastGrid)
      }

      #goto next step i.e. to i-th step in the for-loop
    }#end <--*for (i in 2:n)*
  }#end <--*else #go on with further analysis*

  ##Return a list containing the vectors lowerBounds, upperBounds, probExit, probDifference
  toBeReturned<-list(lowerBounds=lowerBounds, upperBounds=upperBounds, exitProbabilities=probExit,differencesExitProbabilities=probDifference)
  return(toBeReturned)
}#end <--*function(...)*
