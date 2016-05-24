"alphaByUseFunction" <-
function(whatSpendingFunctionIsUsed,n,alpha,phi,OneOrTwoSidedBounds,t)
{
  ###Initialize variables###
  probExit<-0


  #calculate probabilities according to use function

  ###O'Brien and Fleming-type ###
  if (whatSpendingFunctionIsUsed==1)
  {
    for (k in 1:n)
    {
    probExit[k] <- OneOrTwoSidedBounds * asOBF( alpha, t[k], OneOrTwoSidedBounds )
    }
  }

    ### Pocock-type ###
    else if (whatSpendingFunctionIsUsed==2 || whatSpendingFunctionIsUsed==5)
         {
           for (k in 1:n)
           {
             probExit[k] <- OneOrTwoSidedBounds *  asPocock( alpha, t[k], OneOrTwoSidedBounds )
           }
         }

         ### Power Family ###
         else if (whatSpendingFunctionIsUsed==3)
              {
                for (k in 1:n)
                {
                  probExit[k] <- OneOrTwoSidedBounds * asPowerFamily( alpha, t[k], OneOrTwoSidedBounds, phi )
                }
              }

              ### Hwang-Shih-de Cani family ###
              else if (whatSpendingFunctionIsUsed==4)
                   {
                     for (k in 1:n)
                     {
                     probExit[k] <- OneOrTwoSidedBounds * asHSdeCani( alpha, t[k], OneOrTwoSidedBounds, phi )
                     }
                   }

                   ###-----------------------------------------------###
                   ##-- more spending functions could be added here --##
                   ###-----------------------------------------------###
                   else
                   {
                     #if this else-tree is reached something went wrong before
                     print(" Chosen Use Function does not exist", quote=FALSE)
                     print(" or is not implemented correctly yet!", quote=FALSE)
                   }


  ## returns the vector 'probExit' of type I error spent at each analysis
  ## with length 'n'
  return(probExit)


}#end <--*function(...)*
