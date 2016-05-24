#### Author: Adriano Zanin Zambom
#### contact: adriano.zambom@gmail.com
#### last modified: 2/Oct/2015
####
#### papers of reference: 
#### Zambom, A. Z. and Akritas, M. G. (2012). a) Nonparametric Model Checking and Variable Selection. Statistica Sinica, v. 24, pp. 1837
#### Zambom, A. Z. and Akritas, M. G. (2012). b) Signicance Testing and Group Variable Selection. Journal of Multivariate Analysis, v. 133, pp. 51

 
##########################################################################################################################
# Function: npvarselec
#
##########################################################################################################################

npvarselec <- function(X, Y, method = "backward", p = 7, degree.pol = 0, kernel.type = "epanech", bandwidth = "CV", gridsize = 10, dim.red = c(1,10)) UseMethod("npvarselec")


print.npvarselec <- function(x,...)
{

    cat("\n\nNumber of Covariates Selected: ")
    if (x$selected[1] == -1)
       cat("No significant covariates") else
    {
       cat(length(x$selected))
       cat("\n\nCovariate(s) Selected: \n")
       cat("---------------------------\n")
       cat("Covariate Index  |  p-value\n")
       for (i in 1:length(x$selected))
          if (x$selected[i] < 10)
             cat("         ",x$selected[i],"     | ",x$p_values[i],"\n") else
          if (x$selected[i] < 100)
             cat("        ",x$selected[i],"     | ",x$p_values[i],"\n") else
          if (x$selected[i] < 1000)
             cat("       ",x$selected[i],"     | ",round(x$p_values[i],digits = 2),"\n") else
             cat("      ",x$selected[i],"     | ",round(x$p_values[i],digits = 2),"\n") 

       cat("---------------------------\n")
    }
    
}



##########################################################################################################################
# Function: npvarselec
# Parameters Received:
#
# X = Matrix with columns being the variables and rows the observations
# Y = vector of observations
# method = backward, forward, forward2
# p = number of covariate values included in the window W_i, has to be an odd number!
# degree.pol = the degree of the polynomial to fit the nonparametric regression
# kernel.type = type of kernel (box, truncated normal, etc) for estimation of m_1
# bandwidth = "CV" for cross validation, "GCV" for GCV, "Adp" for adaptive, or any other number will be fixed
# dim.red: vector where 1st entry indicates 1 for SIR and 2 for SPC. Second entry indicates number of slices (if SIR) or number of components (if SPC)
#
# Values Returned:
# test = list 
#   test$selected = indices of covariates selected
#   test$p_values = of each selected covariate 
##########################################################################################################################

npvarselec.default <- function(X, Y, method = "backward", p = 7, degree.pol = 0, kernel.type = "epanech", bandwidth = "CV", gridsize = 10,  dim.red = c(1,10), ...)
{
    if ((!identical(method,"backward")) && (!identical(method,"forward")) && (!identical(method,"forward2")))
       stop("\n\nInvalid method: ",method,"\n\n") else
    if ((!identical(bandwidth,"CV")) && (!identical(bandwidth,"GCV")) && (!identical(bandwidth,"CV2")) && (!identical(bandwidth,"GCV2")))
       stop("\n\nInvalid bandwidth: npvarselec function only takes bandwidth = 'CV', or 'GCV', or 'CV2', or 'GCV2'. \n\n") 


   print_iteration <- function(iter, not_eliminated) ## just to print nicely on the screen
   {
         cat(iter) 
         if (iter < 10)
            cat("     | ") else
         if (iter < 100)
            cat("    | ") else
         if (iter < 1000)
            cat("   | ") else
            cat("  | ")      
            cat(not_eliminated)
            cat("\n")
   }


   options(digits=2)

   if (is.null(dim(X))) # only 1 covariate
   {
      ANOVA = anova_test_univariate(X, Y, p)
      if (ANOVA > .05)
          not_eliminated = -1 else
          not_eliminated = 1
   } else
   {
      n = dim(X)[1]
      d = dim(X)[2]
      
      cat("-------------------------------------------------------------\n")
      cat("Iter. | Variables in the model\n")

      if (identical(method,"backward"))  ##______________________________________ Backward Elimination
      {
         X_aux = X
         not_eliminated = seq(1,d)  
      
         iter = 1
         print_iteration(iter, not_eliminated)

         STOP = 0 
         while (STOP == 0)  
         {

            ANOVA = 0
            for (ind_test in 1:d) ## test each variable with all others in the model
               ANOVA[ind_test] = npmodelcheck(X_aux, Y, ind_test, p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value

            ordem = order(ANOVA)
            benjamini = 0
            for (ben in 1:d)
     	        if (ANOVA[ordem][ben] < ben*0.05/d/(sum(1/seq(1,d))))
	               benjamini = benjamini + 1

            if (benjamini == d)  ## all variables are significant with FDR correction
               STOP = 1 else
            {
               X_aux = X_aux[,-ordem[d]]
               ANOVA = ANOVA[-ordem[d]]
               not_eliminated = not_eliminated[-ordem[d]]
               ordem = ordem[-d]
               d = d-1

               iter = iter + 1
               print_iteration(iter, not_eliminated)
            
               if (d == 1)
               {
                  STOP = 1  
                  ordem = 1 
                  ANOVA = anova_test_univariate(X_aux, Y, p)
                  if (ANOVA > .05)
                     not_eliminated = -1
               }
            }
         }
         
      } else
      if (identical(method,"forward"))  ##______________________________________ Forward Selection
      {
          not_in_model = seq(1,d) 
          not_eliminated = 0
         
          ANOVA = 0
          for (ind_test in 1:d)
             ANOVA[ind_test] = npmodelcheck(X[,ind_test], Y, 1, p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value 
          
          ordem = order(ANOVA)
          if (ANOVA[ordem][1] > 0.05) ## no covariate is significant to enter at first, then stop here
             not_eliminated = -1
          else
          {
              not_eliminated = ordem[1]
              not_in_model = not_in_model[-ordem[1]]
              print_iteration(1, not_eliminated)
              iter = 2
              
              p_val = 0;
              STOP = 0 
              while (STOP == 0)  
              {
                  ANOVA = 0
                  for (ind_test in 1:length(not_in_model))
                     ANOVA[ind_test] = npmodelcheck(cbind(X[,not_eliminated],X[,not_in_model[ind_test]]), Y, (length(not_eliminated)+1), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
                  
                  ordem = order(ANOVA)
                  not_eliminated = c(not_eliminated,not_in_model[ordem[1]])             ## include the smallest p-value
                  not_in_model = not_in_model[-ordem[1]]

                  ANOVA = 0                                        ## test all in the model with the new included
                  for (ind_test in 1:length(not_eliminated))
                     ANOVA[ind_test] = npmodelcheck(X[,not_eliminated], Y, ind_test, p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize,  dim.red = dim.red)$p_value

                  ordem = order(ANOVA)
                  
                  benjamini = 0
                  for (ben in 1:length(not_eliminated))
                     if (ANOVA[ordem][ben] < ben*0.05/length(not_eliminated)/(sum(1/seq(1,length(not_eliminated)))))
                        benjamini = benjamini + 1
                  
                  
                  if (benjamini == length(not_eliminated))  ## all variables are significant with FDR correction
                  {
                      print_iteration(iter, not_eliminated)
                      iter = iter + 1
                      if (benjamini == d)
                         STOP = 1
                  } else
                  {
                      STOP = 1
                      not_eliminated = not_eliminated[-length(not_eliminated)]      
                  }
              }
          }
          
          
      } else
      if (identical(method,"forward2"))  ##______________________________________ "Forward type 2" Selection
      {
        not_in_model = seq(1,d) 
        not_eliminated = 0
         
          ANOVA = 0
          for (ind_test in 1:d)
             ANOVA[ind_test] = npmodelcheck(X[,ind_test], Y, 1, p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value 
          
          ordem = order(ANOVA)
          if (ANOVA[ordem][1] > 0.05) ## no covariate is significant to enter at first, then stop here
             not_eliminated = -1
          else
          {
              not_eliminated = ordem[1]
              not_in_model = not_in_model[-ordem[1]]
              print_iteration(1, not_eliminated)
              iter = 2
              
              p_val = 0;
              STOP = 0 
              while (STOP == 0)  
              {
                  best = c(-1,2)
                  qtos = length(not_in_model)
                  for (ind_test in 1:qtos)  ## choose the next predictor that gives smallest p-value of last in model with FDR
                  {
                     not_eliminated = c(not_eliminated,not_in_model[1])   ## this is always 1, cause when I return 
                     not_in_model = not_in_model[-1]                      ## it, it goes to the end of the vector

                     ANOVA = 0                                        ## test all in the model with the new included
                     for (ind_test2 in 1:length(not_eliminated))
                        ANOVA[ind_test2] = npmodelcheck(X[,not_eliminated], Y, ind_test2, p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize,  dim.red = dim.red)$p_value

                     ordem = order(ANOVA)   ## check if the whole model is signif
                     benjamini = 0
                     for (ben in 1:length(not_eliminated))
                        if (ANOVA[ordem][ben] < ben*0.05/length(not_eliminated)/(sum(1/seq(1,length(not_eliminated)))))
                           benjamini = benjamini + 1

                     if ((max(ANOVA) < best[2]) && (benjamini == length(not_eliminated)))
                         best = c(ind_test, max(ANOVA))

                     not_in_model = c(not_in_model, not_eliminated[length(not_eliminated)])
                     not_eliminated = not_eliminated[-length(not_eliminated)]
                  }

                  if (best[2] < 0.05/(sum(1/seq(1,length(not_eliminated)+1))))
                  {
                      not_eliminated = c(not_eliminated,not_in_model[best[1]])   ## this is always 1, cause when I return 
                      not_in_model = not_in_model[-best[1]]                      ## it, it goes to the end of the vector
                      print_iteration(iter, not_eliminated)
                      iter = iter + 1
                      if (length(not_eliminated) == d)
                         STOP = 1
                  } else
                  {
                      STOP = 1
                  }
              }
          }
      } 
      

   }   
      
   ANOVA = 0                                        ## test all in the model to return the p-values
   for (ind_test in 1:length(not_eliminated))
       ANOVA[ind_test] = npmodelcheck(X[,not_eliminated], Y, ind_test, p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize,  dim.red = dim.red)$p_value



   cat("-------------------------------------------------------------\n")


   test = list()
   test$selected = not_eliminated
   test$p_values = ANOVA

   result = test
   result$call <- match.call()

   class(result) <- "npvarselec"
   result
}





