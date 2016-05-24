#### Author: Adriano Zanin Zambom
#### contact: adriano.zambom@gmail.com
#### last modified: 2/Oct/2015
####
#### papers of reference: 
#### Zambom, A. Z. and Akritas, M. G. (2012). a) Nonparametric Model Checking and Variable Selection. Statistica Sinica, v. 24, pp. 1837.
#### Zambom, A. Z. and Akritas, M. G. (2012). b) Signicance Testing and Group Variable Selection. Journal of Multivariate Analysis, v. 133, pp. 51.


##########################################################################################################################
# Function: group.npvarselec
#
##########################################################################################################################

group.npvarselec <- function(X, Y, groups, method = "backward", p = 7,  fitSPC = TRUE, degree.pol = 0, kernel.type = "epanech", bandwidth = "CV", gridsize = 10, dim.red = c(1,10)) UseMethod("group.npvarselec")


print.group.npvarselec <- function(x,...)
{
    cat("\n\nNumber of Groups Selected: ")
    if (x$selected[1] == -1)
       cat("No significant groups") else
    {
       cat(length(x$selected))
       cat("\n\nGroups Selected: \n")
       cat("---------------------------\n")
       cat("Group Index  |  p-value\n")
       for (i in 1:length(x$selected))
          if (x$selected[i] < 10)
             cat("        ",x$selected[i],"  | ",x$p_values[i],"\n") else
          if (x$selected[i] < 100)
             cat("       ",x$selected[i],"  | ",x$p_values[i],"\n") else
          if (x$selected[i] < 1000)
             cat("      ",x$selected[i],"  | ",round(x$p_values[i],digits = 2),"\n") else
             cat("     ",x$selected[i],"  | ",round(x$p_values[i],digits = 2),"\n") 

       cat("---------------------------\n")
    }
    
}




##########################################################################################################################
# Function: group.npvarselec
# Parameters Received:
#
# X = Matrix with columns being the variables and rows the observations
# Y = vector of observations
# groups = vector with indices of groups to be tested
# method = backward(default), forward or forward2
# fitSPC = if true, m_1 is estimated with the first SPC of each group
# p = size of the window Wi (number of neighbors to include in the ANOVA cell)
# degree.pol = the degree of the polynomial to fit the nonparametric regression
# kernel.type = type of kernel for estimation of m_1
# bandwidth = "CV" for cross validation, "GCV" for GCV, "Adp" for adaptive, or any other number/vector will be fixed
# dim.red: vector where 1st entry indicates 1 for SIR and 2 for SPC. Second entry indicates number of slices (if SIR) or number of components (if SPC)
#
# Values Returned:
# test = list 
#   test$selected = indices of selected groups
#   test$p_values = of each group
##########################################################################################################################

group.npvarselec.default <- function(X, Y, groups, method = "backward", p = 7, fitSPC = TRUE, degree.pol = 0, kernel.type = "epanech", bandwidth = "CV", gridsize = 10,   dim.red = c(1,10), ...)
{
    
   if ((!identical(method,"backward")) && (!identical(method,"forward")) && (!identical(method,"forward2")))
      stop("\n\nInvalid method: ",method,"\n\n") else
   if (!is.vector(groups,"list"))
      stop("\n\ngroups must be a vector of type list.\n\n") else
   if ((bandwidth != "CV") && (bandwidth != "GCV") && (bandwidth != "CV2") && (bandwidth != "GCV2"))
      stop("\n\nInvalid bandwidth: npvarselec function only takes bandwidth = 'CV', or 'GCV', or 'CV2', or 'GCV2'. \n\n")


   print_iteration <- function(iter, selected) ## just to print nicely on the screen
   {
         cat(iter) 
         if (iter < 10)
            cat("     | ") else
         if (iter < 100)
            cat("    | ") else
         if (iter < 1000)
            cat("   | ") else
            cat("  | ")      
            cat(selected)
            cat("\n")
   }

    get_quais <- function(groups, group_number) ## function to get indices of covariates of X are in "group"
    {
        result = groups[[group_number[1]]]
        if (length(group_number) > 1)
        for (i in 2:length(group_number))
        result = c(result,groups[[group_number[i]]])
        
        return(result)     
    }

    
   options(digits=2)

   if (is.null(dim(X))) # only 1 covariate
   {
      ANOVA = anova_test_univariate(X, Y, p)
      if (ANOVA > .05)
          selected = -1 else
          selected = 1
   } else
   {
      n = dim(X)[1]
      d = length(groups) # number of groups
      
      quais = groups[[1]]                         ## quais is a vector with list of all covariates of all groups
      length_subgroups = length(groups[[1]])      ## this is the number of covariates in each group
      for (i in 2:d)
      {
         quais = c(quais,groups[[i]])
         length_subgroups[i] = length(groups[[i]])
      }
         

       if (identical(fitSPC,TRUE))   ## need to fit local poly with only 1st SPC of each group
       {
           X_SPC = matrix(0,n,d)
           for (i in 1:d)                 ## computes the 1st SPC of each group
           {
               X_I_aux = X[,get_quais(groups,i)]
               if (is.null(dim(X_I_aux)))
               X_SPC[,i] = X_I_aux else
               {
                   s0 = 0
                   for (prin in 1:(dim(X_I_aux)[2]))
                      s0[prin] = anova_test_univariate(X_I_aux[,prin], Y, p = p)
                   quais2 = which(s0 < 0.3)       #### theta = 0.3 threshold parameter for SPC
                   if (length(quais2) == 0) quais2 = 1    ### if no p-value is < theta, get the smallest
                   quais2 = order(s0)[1:(length(quais2))]   
                   Princ_Comp = princomp(X_I_aux[,quais2])
                   if (length(Princ_Comp$loadings[,1]) == 1)
                      X_SPC[,i] =  as.numeric(X_I_aux[,quais2]*(Princ_Comp$loadings[,1])) else
                      X_SPC[,i] =  as.numeric(X_I_aux[,quais2]%*%(Princ_Comp$loadings[,1]))                     
               }
           }
       }        

      cat("-------------------------------------------------------------\n")
      cat("Iter. | Variables in the model\n")

      if (identical(method,"backward"))  ##______________________________________ Group Backward Elimination
      {
         selected = seq(1,d)  
      
         iter = 1
         STOP = 0 
         if (identical(fitSPC,TRUE))   ## need to fit local poly with only 1st SPC of each group
         {
            while (STOP == 0)  # with 1st SPC
            {
                ANOVA = 0
                for (ind_test in 1:length(selected))
                   ANOVA[ind_test] = npmodelcheck(cbind(X[,get_quais(groups,selected[ind_test])],X_SPC[,selected[-ind_test]]), Y, seq(1,length(get_quais(groups,selected[ind_test]))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value

                ordem = order(ANOVA)
                
                benjamini = 0
                for (ben in 1:d)
                   if (ANOVA[ordem][ben] < ben*0.05/d/(sum(1/seq(1,d))))
                      benjamini = benjamini + 1
                
                if (benjamini == d)  ## all variables are significant with FDR correction
                   STOP = 1 else
                {
                    ANOVA = ANOVA[-ordem[d]]
                    selected = selected[-ordem[d]]
                    ordem = ordem[-d]
                    d = d-1
                    
                    print_iteration(iter, selected)
                    iter = iter + 1
                    
                    if (d == 1)
                    {
                        STOP = 1  
                        ordem = 1 
                        ANOVA = npmodelcheck(X[,get_quais(groups,selected[1])], Y, seq(1,length(get_quais(groups,selected[1]))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
                        if (ANOVA > .05)
                           selected = -1
                    }
                }
            }          
         } else  ## dont use 1st SPC of groups
            while (STOP == 0)  # with all covariates fo fit
            {
               ANOVA = npmodelcheck(X[,quais], Y, seq(1,length_subgroups[1]), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
               for (ind_test in 2:d) ## test each group with all others in the model
                  ANOVA[ind_test] = npmodelcheck(X[,quais], Y, seq(sum(length_subgroups[1:(ind_test-1)])+1,sum(length_subgroups[1:ind_test])), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value

               ordem = order(ANOVA)

               benjamini = 0
               for (ben in 1:d)
     	           if (ANOVA[ordem][ben] < ben*0.05/d/(sum(1/seq(1,d))))
	                  benjamini = benjamini + 1


               if (benjamini == d)  ## all variables are significant with FDR correction
                  STOP = 1 else
               {
                  ANOVA = ANOVA[-ordem[d]]
                  selected = selected[-ordem[d]]
                  if (ordem[d] == 1)
                     quais = quais[-seq(1,length_subgroups[1])] else
                     quais = quais[-seq(sum(length_subgroups[1:(ordem[d]-1)])+1,sum(length_subgroups[1:ordem[d]]))]
                  length_subgroups = length_subgroups[-ordem[d]]
                  ordem = ordem[-d]
                  d = d-1

                  print_iteration(iter, selected)
                  iter = iter + 1
            
                  if (d == 1)
                  {
                     STOP = 1  
                     ordem = 1 
                     ANOVA = npmodelcheck(X[,quais], Y, seq(1,length_subgroups[1]), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
                     if (ANOVA > .05)
                        selected = -1
                  }
               }
            }
          
      } else
      if (identical(method,"forward"))  ##______________________________________ Group Forward Selection
      {
         not_in_model = seq(1,d) 
                    
         ANOVA = 0 
         for (ind_test in 1:d) ## test each group with all others in the model
            ANOVA[ind_test] = npmodelcheck(X[,get_quais(groups,ind_test)], Y, seq(1,length(get_quais(groups,ind_test))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
          
          
          ordem = order(ANOVA)
          if (ANOVA[ordem][1] > 0.05) ## no covariate is significant to enter at first, then stop here
             selected = -1
          else
          {
              selected = ordem[1]
              not_in_model = not_in_model[-ordem[1]]
              print_iteration(1, selected)
              iter = 2
              
              p_val = 0;
              STOP = 0 
              if (identical(fitSPC,TRUE))   ## need to fit local poly with only 1st SPC of each group
              {
                  while (STOP == 0)  # with 1st SPC
                  {
                      ANOVA = 0
                      for (ind_test in 1:length(not_in_model))
                         ANOVA[ind_test] = npmodelcheck(cbind(X[,get_quais(groups,not_in_model[ind_test])],X_SPC[,selected]), Y, seq(1,length(get_quais(groups,not_in_model[ind_test]))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value

                      ordem = order(ANOVA)
                      selected = c(selected,not_in_model[ordem[1]])             ## include the smallest p-value
                      not_in_model = not_in_model[-ordem[1]]

                      ANOVA = 0                                        ## test all in the model with the new included
                      for (ind_test in 1:length(selected))
                         ANOVA[ind_test] = npmodelcheck(cbind(X[,get_quais(groups,selected[ind_test])],X_SPC[,selected[-ind_test]]), Y, seq(1,length(get_quais(groups,selected[ind_test]))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
                      
                      ordem = order(ANOVA)
                      
                      benjamini = 0
                      for (ben in 1:length(selected))
                         if (ANOVA[ordem][ben] < ben*0.05/length(selected)/(sum(1/seq(1,length(selected)))))
                            benjamini = benjamini + 1


                      if (benjamini == length(selected))  ## all variables are significant with FDR correction
                      {
                          print_iteration(iter, selected)
                          iter = iter + 1
                          if (benjamini == d)
                             STOP = 1
                      } else
                      {
                          STOP = 1
                          not_in_model = c(not_in_model,selected[length(selected)])
                          ANOVA = ANOVA[ordem]
                          ANOVA = ANOVA[-length(ANOVA)]
                          selected = selected[-length(selected)]      
                      }
                  }
              } else
                 while (STOP == 0)  
                 {
                     ANOVA = 0
                     for (ind_test in 1:length(not_in_model))
                        ANOVA[ind_test] = npmodelcheck(cbind(X[,get_quais(groups,not_in_model[ind_test])],X[,get_quais(groups,selected)]), Y, seq(1,length(get_quais(groups,not_in_model[ind_test]))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
                  

                     ordem = order(ANOVA)
                     selected = c(selected,not_in_model[ordem[1]])             ## include the smallest p-value
                     not_in_model = not_in_model[-ordem[1]]

                     ANOVA = 0                                        ## test all in the model with the new included
                     for (ind_test in 1:length(selected))
                        ANOVA[ind_test] = npmodelcheck(cbind(X[,get_quais(groups,selected[ind_test])],X[,get_quais(groups,selected[-ind_test])]), Y, seq(1,length(get_quais(groups,selected[ind_test]))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
                  
                     ordem = order(ANOVA)
                  
                     benjamini = 0
                     for (ben in 1:length(selected))
                        if (ANOVA[ordem][ben] < ben*0.05/length(selected)/(sum(1/seq(1,length(selected)))))
                           benjamini = benjamini + 1
                     if (benjamini == length(selected))  ## all variables are significant with FDR correction
                     {
                         print_iteration(iter, selected)
                         iter = iter + 1
                         if (benjamini == d)
                            STOP = 1
                     } else
                     {
                         STOP = 1
                         not_in_model = c(not_in_model,selected[length(selected)])
                         ANOVA = ANOVA[ordem]
                         ANOVA = ANOVA[-length(ANOVA)]
                         selected = selected[-length(selected)]      
                     }

                 }
          }

      } else
      if (identical(method,"forward2"))  ##______________________________________ Group Stepwise Selection
      {
         not_in_model = seq(1,d) 
                    
         ANOVA = 0 
         for (ind_test in 1:d) ## test each group with all others in the model
            ANOVA[ind_test] = npmodelcheck(X[,get_quais(groups,ind_test)], Y, seq(1,length(get_quais(groups,ind_test))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
          
          
          ordem = order(ANOVA)
          if (ANOVA[ordem][1] > 0.05) ## no covariate is significant to enter at first, then stop here
             selected = -1
          else
          {
              selected = ordem[1]
              not_in_model = not_in_model[-ordem[1]]
              print_iteration(1, selected)
              iter = 2
              
              p_val = 0;
              STOP = 0 
              if (identical(fitSPC,TRUE))   ## need to fit local poly with only 1st SPC of each group
              {
                  while (STOP == 0)  # with 1st SPC
                  {
                     best = c(-1,2)
                     qtos = length(not_in_model)
                     for (ind_test in 1:qtos)  ## choose the next predictor that gives smallest p-value of last in model with FDR
                     {
                        selected = c(selected,not_in_model[1])   ## this is always 1, cause when I return 
                        not_in_model = not_in_model[-1]                      ## it, it goes to the end of the vector

                        ANOVA = 0                                        ## test all in the model with the new included
                        for (ind_test2 in 1:length(selected))
                           ANOVA[ind_test2] = npmodelcheck(cbind(X[,get_quais(groups,selected[ind_test2])],X_SPC[,selected[-ind_test2]]), Y, seq(1,length(get_quais(groups,selected[ind_test2]))) , p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize,  dim.red = dim.red)$p_value

                        ordem = order(ANOVA)   ## check if the whole model is signif
                        benjamini = 0
                        for (ben in 1:length(selected))
                           if (ANOVA[ordem][ben] < ben*0.05/length(selected)/(sum(1/seq(1,length(selected)))))
                              benjamini = benjamini + 1
                         
                        if ((max(ANOVA) < best[2]) && (benjamini == length(selected)))
                            best = c(ind_test, max(ANOVA))
  
                        not_in_model = c(not_in_model, selected[length(selected)])
                        selected = selected[-length(selected)]
                     }

                     if (best[2] < 0.05/(sum(1/seq(1,length(selected)+1))))
                     {
                         selected = c(selected,not_in_model[best[1]])   ## this is always 1, cause when I return 
                         not_in_model = not_in_model[-best[1]]                      ## it, it goes to the end of the vector
                         print_iteration(iter, selected)
                         iter = iter + 1
                         if (length(selected) == d)
                            STOP = 1
                     } else
                     {
                         STOP = 1
                     }
                  }
              } else
              while (STOP == 0)  # not fitSPC
              {
                  best = c(-1,2)
                  qtos = length(not_in_model)
                  for (ind_test in 1:qtos)  ## choose the next predictor that gives smallest p-value of last in model with FDR
                  {
                      selected = c(selected,not_in_model[1])   ## this is always 1, cause when I return 
                      not_in_model = not_in_model[-1]                      ## it, it goes to the end of the vector
                      
                      ANOVA = 0                                        ## test all in the model with the new included
                      for (ind_test2 in 1:length(selected))
                         ANOVA[ind_test2] = npmodelcheck(cbind(X[,get_quais(groups,selected[ind_test2])],X[,get_quais(groups,selected[-ind_test2])]), Y, seq(1,length(get_quais(groups,selected[ind_test2]))) , p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize,  dim.red = dim.red)$p_value
                      
                      ordem = order(ANOVA)   ## check if the whole model is signif
                      benjamini = 0
                      for (ben in 1:length(selected))
                         if (ANOVA[ordem][ben] < ben*0.05/length(selected)/(sum(1/seq(1,length(selected)))))
                            benjamini = benjamini + 1
                      
                      if ((max(ANOVA) < best[2]) && (benjamini == length(selected)))
                         best = c(ind_test, max(ANOVA))
                      
                      not_in_model = c(not_in_model, selected[length(selected)])
                      selected = selected[-length(selected)]
                  }
                  
                  if (best[2] < 0.05/(sum(1/seq(1,length(selected)+1))))
                  {
                      selected = c(selected,not_in_model[best[1]])   ## this is always 1, cause when I return 
                      not_in_model = not_in_model[-best[1]]                      ## it, it goes to the end of the vector
                      print_iteration(iter, selected)
                      iter = iter + 1
                      if (length(selected) == d)
                         STOP = 1
                  } else
                  {
                      STOP = 1
                  }
              }
          }      
      }

   }

    # here I do the test again to return the right ANOVA
       ANOVA = 0                                        
    if (identical(fitSPC,TRUE))
       for (ind_test in 1:length(selected))
           ANOVA[ind_test] = npmodelcheck(cbind(X[,get_quais(groups,selected[ind_test])],X_SPC[,selected[-ind_test]]), Y, seq(1,length(get_quais(groups,selected[ind_test]))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value
    else
       for (ind_test in 1:length(selected))
           ANOVA[ind_test] = npmodelcheck(cbind(X[,get_quais(groups,selected[ind_test])],X[,get_quais(groups,selected[-ind_test])]), Y, seq(1,length(get_quais(groups,selected[ind_test]))), p = p, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize, dim.red = dim.red)$p_value


   cat("-------------------------------------------------------------\n")


   test = list()
   test$selected = selected
   test$p_values = ANOVA

   result = test
   result$call <- match.call()

   class(result) <- "group.npvarselec"
   result
}
