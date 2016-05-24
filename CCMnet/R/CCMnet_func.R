

CCMnet_constr <- function(Network_stats, Prob_Distr, Prob_Distr_Params, 
                         samplesize = 5000, burnin=1000, interval=1000,
                         statsonly=TRUE, P=NULL,
                         population, covPattern = NULL, remove_var_last_entry = FALSE) {
  
  if (length(population) == 1) {
    return(uni_modal_constr(Network_stats, Prob_Distr, Prob_Distr_Params, 
                            samplesize, burnin, interval,
                            statsonly, P,
                            population, covPattern, remove_var_last_entry)
           )
  } else if (length(population) == 2) {
    return(bi_modal_constr(Network_stats, Prob_Distr, Prob_Distr_Params, 
                           samplesize, burnin, interval,
                           statsonly, P,
                           population, covPattern, remove_var_last_entry)
           )
  }
}


uni_modal_constr <- function(Network_stats, Prob_Distr, Prob_Distr_Params, 
                          samplesize, burnin, interval,
                          statsonly, P,
                          population, covPattern, remove_var_last_entry) {

  error = 0

  if(is.null(covPattern)) {
    covPattern = rep(1,population)
  }
  G_max_degree_bool = FALSE
  ER_prob = .05
  
  if ((length(Network_stats) == 1) && (Network_stats == "DegreeDist")){
    max_degree_f = max_degree = length(Prob_Distr_Params[[1]][[1]])-1
  } else if  ((length(Network_stats) == 1) && (Network_stats == "Density")) {
    max_degree_f = max_degree = population - 1
    if (Prob_Distr == "NP") {
      ER_prob = (max(which(Prob_Distr_Params[[1]][[1]] > 0))-1)/choose(population,2) * .8
    }
  } else if ((length(Network_stats) == 2) && (Network_stats[1] == "DegreeDist") && (Network_stats[2] == "Mixing")) {
    max_degree_1 = length(Prob_Distr_Params[[1]][[1]][[1]])-1
    max_degree_2 = length(Prob_Distr_Params[[1]][[1]][[2]])-1
    max_degree = min(max_degree_1, max_degree_2)
    max_degree_f = max(max_degree_1, max_degree_2)
  } else if ((length(Network_stats) == 2) && (Network_stats[1] == "Mixing") && (Network_stats[2] == "DegreeDist")) {
    max_degree_1 = length(Prob_Distr_Params[[2]][[1]][[1]])-1
    max_degree_2 = length(Prob_Distr_Params[[2]][[1]][[2]])-1
    max_degree = min(max_degree_1, max_degree_2)
    max_degree_f = max(max_degree_1, max_degree_2)
  } else if ((length(Network_stats) == 1) && (Network_stats == "DegMixing")) {
    max_degree_f = max_degree = floor(sqrt(2*length(upper.tri(Prob_Distr_Params[[1]][[1]], diag = TRUE))))
  } else if  ((length(Network_stats) == 2) && (Network_stats[1] == c("DegMixing")) && (Network_stats[2] == c("Triangles"))) {
    max_degree_f = max_degree = floor(sqrt(2*length(upper.tri(Prob_Distr_Params[[1]][[1]], diag = TRUE))))
  } else if  ((length(Network_stats) == 2) && (Network_stats[1] == c("Triangles")) && (Network_stats[2] == c("DegMixing"))) {
    max_degree_f = max_degree = floor(sqrt(2*length(upper.tri(Prob_Distr_Params[[2]][[1]], diag = TRUE))))
  } else {
    max_degree_f = max_degree = population - 1
  }
  
  Gen_Net_counter = 1
  if (is.null(P)) {
    print("Generating Random Initial Network...")
    while (!G_max_degree_bool) { 
      if (Gen_Net_counter > 1) {
        print(paste("Try ", Gen_Net_counter, sep=""))
      }
      g = as.network(rgraph(n=population, m=1, tprob=ER_prob, mode="graph", diag=FALSE, replace=FALSE,
             tielist=NULL, return.as.edgelist=FALSE), directed = FALSE)
      g %v% "CovAttribute" = covPattern
      ER_prob = ER_prob/2
      G_max_degree_bool = max(degree(g)/2) <= max_degree
      Gen_Net_counter =   Gen_Net_counter + 1
    }
    print("COMPLETED: Generated Random Initial Network")
    P = g
  } else {
    g = P
  }
  
  max_degree = max_degree_f
  
  P_edge_mat = unlist(P$mel)
  dim(P_edge_mat) = c(3,network.edgecount(P))
  Trans_nedges = c(network.edgecount(P),0,0)
  Trans_networktails = P_edge_mat[1,]
  Trans_networkheads = P_edge_mat[2,]
  
  edge_mat = unlist(g$mel)
  dim(edge_mat) = c(3,network.edgecount(g))
  
  nedges = c(network.edgecount(g),0,0)
  tails = edge_mat[1,]
  heads = edge_mat[2,]
  Clist_n = network.size(g)
  
  if ((length(Network_stats) == 1) && (Network_stats == "DegreeDist")) {
    if (length(Prob_Distr_Params[[1]][[1]]) < 2) {
      print("Error: length of mean vector is less than 2")
      error = 1
    }
    if (Prob_Distr == "Normal") {
      mean_vector = Prob_Distr_Params[[1]][[1]]
      var_vector = Prob_Distr_Params[[1]][[2]]
   
      mean_vector = mean_vector / population
      var_vector = var_vector / population^2
      prob_type = c(1,0,0,0,1)
      if (dim(var_vector)[1] != dim(var_vector)[2]) {
        print("Error: Covariance matrix is not square")
        error = 1
      }
      if (dim(var_vector)[1] != length(mean_vector)) {
        print("Error: Dimension mismatch between covariance matrix and mean vector")
        error = 1
      }
      
      if (remove_var_last_entry == TRUE) {
        inverse_var_x = solve(var_vector[-length(mean_vector),-length(mean_vector)])
        inverse_var_x = rbind(inverse_var_x,0)
        inverse_var_x = cbind(inverse_var_x,0)
        var_vector = c(inverse_var_x)
      } else {
        var_vector = solve(var_vector)
      }      
      
    } else if (Prob_Distr == "NegBin") {
      mean_vector = Prob_Distr_Params[[1]][[1]] 
      var_vector = c(0,0)
      prob_type = c(2,0,0,0,1)
    } else if (Prob_Distr == "DirMult") {
      mean_vector = Prob_Distr_Params[[1]][[1]] 
      var_vector = c(0,0)
      prob_type = c(3,0,0,0,1)
    } else {
      print("Error: No such distribution for degree distribution currently implemented.")
      print("Email ravi.goyal@mail.harvard.edu to add feature.")
      error = 1
    }
    if (error == 0) {
      Clist_nterms = 2 #Number of different terms
      Clist_fnamestring = "edges degree"
      Clist_snamestring = "CCMnet CCMnet"
      inputs = c(c(0,1,0,0), length(mean_vector), length(mean_vector), c(0:(length(mean_vector)-1)))
      eta0 = rep(-999.5,length(c(nedges[1], mean_vector,0)))
      stats = c(nedges[1], tabulate(degree(g)/2 + 1),rep(0, length(mean_vector) - length(tabulate(degree(g)/2 + 1))))      
      MHproposal_name = "TNT"
      MHproposal_package = "CCMnet"
    }
  } else if ((length(Network_stats) == 1) && (Network_stats == "Density")) {
    if (Prob_Distr == "Normal") {
      prob_type = c(0,0,0,0,1)
      mean_vector = c(Prob_Distr_Params[[1]][[1]],Prob_Distr_Params[[1]][[1]])
      var_vector = c(Prob_Distr_Params[[1]][[2]], Prob_Distr_Params[[1]][[2]])
      if (length(Prob_Distr_Params[[1]][[1]]) != 1) {
        print("Error: mean value for network density is one positive value")
        error = 1
      }
      if (length(Prob_Distr_Params[[1]][[2]]) != 1) {
        print("Error: variance for network density is one positive value")
        error = 1
      }
    } else if (Prob_Distr == "NP") {
      prob_type = c(0,0,0,0,99)
      mean_vector = Prob_Distr_Params[[1]][[1]]
      var_vector = c(0,0)
    } else {
      print("Error: No such distribution for degree distribution currently implemented.")
      print("Email ravi.goyal@mail.harvard.edu to add feature.")
      error = 1
    }
    if (error == 0) {
      Clist_nterms = 2 #Number of different terms
      Clist_fnamestring = "edges nfstab"
      Clist_snamestring = "CCMnet CCMnet"
      inputs = c(0,1,0,0,1,0)
      eta0 = c(-999.5, -999.5)
      stats = c(nedges[1],nedges[1])
      MHproposal_name = "TNT"
      MHproposal_package = "CCMnet"
    }
  } else if ((length(Network_stats) == 1) && (Network_stats == "Mixing")) {
      print("Error: No such distribution for mixing currently implemented.")
      print("Email ravi.goyal@mail.harvard.edu to add feature.")
      error = 1
  } else if (((length(Network_stats) == 2) && (Network_stats[1] == "DegreeDist") && (Network_stats[2] == "Mixing")) ||
             ((length(Network_stats) == 2) && (Network_stats[1] == "Mixing") && (Network_stats[2] == "DegreeDist")))  {
    if (Network_stats[1] == "Mixing") { #swap prob_distr_params
      Prob_Distr_Params_temp = Prob_Distr_Params[[1]]
      Prob_Distr_Params[[1]] = Prob_Distr_Params[[2]]
      Prob_Distr_Params[[2]] = Prob_Distr_Params_temp
    }
    if (length(Prob_Distr_Params[[1]][[1]][[1]]) != length(Prob_Distr_Params[[1]][[1]][[2]])) {
      print("Error: Current limitation requires mean degree distributions to be of equal length.")
       error = 1
    }
    if (dim(Prob_Distr_Params[[1]][[2]][[1]])[1] != dim(Prob_Distr_Params[[1]][[2]][[2]])[1]) {
      print("Error: Current limitation requires covariance matrices to be of equal dimensions.")
      error = 1
    }
    if (dim(Prob_Distr_Params[[1]][[2]][[1]])[1] != dim(Prob_Distr_Params[[1]][[2]][[1]])[2]) {
      print("Error: Covariance matrix is not square.")
      error = 1
    }
    if (dim(Prob_Distr_Params[[1]][[2]][[2]])[1] != dim(Prob_Distr_Params[[1]][[2]][[2]])[2]) {
      print("Error: Covariance matrix is not square.")
      error = 1
    } 
    if ((Prob_Distr[1] == "Normal") && ((Prob_Distr[2] == "Normal"))) {          
        Clist_nterms = 3 #Number of different terms
        Clist_fnamestring = "edges degree_by_attr nodemix"
        Clist_snamestring = "CCMnet CCMnet CCMnet"
        MHproposal_name = "TNT"
        MHproposal_package = "CCMnet"
        covariate_list = covPattern
        
        inputs1 = c(rbind(c(0:(length(Prob_Distr_Params[[1]][[1]][[1]])-1)), rep(1,length(Prob_Distr_Params[[1]][[1]][[1]]))))  
        inputs2 = c(rbind(c(0:(length(Prob_Distr_Params[[1]][[1]][[2]])-1)), rep(2,length(Prob_Distr_Params[[1]][[1]][[2]]))))
        
        inputs = c(c(0,1,0,0), length(Prob_Distr_Params[[1]][[1]][[1]]) + length(Prob_Distr_Params[[1]][[1]][[2]]),
                   2*(length(Prob_Distr_Params[[1]][[1]][[1]])+length(Prob_Distr_Params[[1]][[1]][[2]])) + g$gal$n, inputs1, inputs2, covariate_list, c(4,2,4 + g$gal$n), c(1,2,2,2), covariate_list)
        eta0 = rep(-999.5,length(c(nedges[1], Prob_Distr_Params[[1]][[1]][[1]],Prob_Distr_Params[[1]][[1]][[2]],1,1)))
        
        mixing = c(0,0,0)
        edge_list = unlist(g$mel)
        dim(edge_list) = c(3,nedges[1])
        for (num_edge in c(1:nedges[1])) {
          if ((covariate_list[edge_list[1,num_edge]] == 1) && (covariate_list[edge_list[2,num_edge]] == 1)) {
            mixing[1] = mixing[1] + 1
          }
          if ((covariate_list[edge_list[1,num_edge]] == 1) && (covariate_list[edge_list[2,num_edge]] == 2)) {
            mixing[2] = mixing[2] + 1
          }
          if ((covariate_list[edge_list[1,num_edge]] == 2) && (covariate_list[edge_list[2,num_edge]] == 1)) {
            mixing[2] = mixing[2] + 1
          }
          if ((covariate_list[edge_list[1,num_edge]] == 2) && (covariate_list[edge_list[2,num_edge]] == 2)) {
            mixing[3] = mixing[3] + 1
          }
        }
        deg_dist_1 = tabulate(degree(g, gmode="graph")[which(covariate_list == 1)]+1)
        deg_dist_2 = tabulate(degree(g, gmode="graph")[which(covariate_list == 2)]+1)
        
        deg_dist_1 = c(tabulate(degree(g, gmode="graph")[which(covariate_list == 1)]+1), rep(0,max(0,length(Prob_Distr_Params[[1]][[1]][[1]])-length(deg_dist_1))))
        deg_dist_2 = c(tabulate(degree(g, gmode="graph")[which(covariate_list == 2)]+1), rep(0,max(0,length(Prob_Distr_Params[[1]][[1]][[2]])-length(deg_dist_2))))
        
        #Assume max degree of both node types is the same
        deg_dist_1 = c(deg_dist_1, rep(0,max(0,length(deg_dist_2)-length(deg_dist_1))))
        deg_dist_2 = c(deg_dist_2, rep(0,max(0,length(deg_dist_1)-length(deg_dist_2))))
        
        stats = c(nedges[1], deg_dist_1, deg_dist_2, mixing[c(2,3)])
        
        mean_vector = c(Prob_Distr_Params[[1]][[1]][[1]], Prob_Distr_Params[[1]][[1]][[2]],  Prob_Distr_Params[[2]][[1]])      

        if (remove_var_last_entry == TRUE) {
          inverse_var_x1 = solve(Prob_Distr_Params[[1]][[2]][[1]][-length(Prob_Distr_Params[[1]][[1]][[1]]),-length(Prob_Distr_Params[[1]][[1]][[1]])])
          inverse_var_x1 = rbind(inverse_var_x1,0)
          inverse_var_x1 = cbind(inverse_var_x1,0)
          
          inverse_var_x2 = solve(Prob_Distr_Params[[1]][[2]][[2]][-length(Prob_Distr_Params[[1]][[1]][[2]]),-length(Prob_Distr_Params[[1]][[1]][[2]])])
          inverse_var_x2 = rbind(inverse_var_x2,0)
          inverse_var_x2 = cbind(inverse_var_x2,0)
        } else {
          inverse_var_x1 = solve(Prob_Distr_Params[[1]][[2]][[1]])
          inverse_var_x2 = solve(Prob_Distr_Params[[1]][[2]][[2]])
        }   
        
        var_vector = c(c(inverse_var_x1),c(inverse_var_x2), Prob_Distr_Params[[2]][[2]])
        
        prob_type = c(1,1,0,0,1)
        
      } else if ((Prob_Distr[1] == "Tdist") && ((Prob_Distr[2] == "Tdist"))) {
        if (length(Prob_Distr_Params[[1]]) != 3) {
          
        }
        if (dim(Prob_Distr_Params[[1]][[3]][1]) > 0) {
          print("Error: Degrees of freedom are not greater than 0.")
          error = 1
        }
        if (dim(Prob_Distr_Params[[1]][[3]][2]) > 0) {
          print("Error: Degrees of freedom are not greater than 0.")
          error = 1
        } 
        Clist_nterms = 3 #Number of different terms
        Clist_fnamestring = "edges degree_by_attr nodemix"
        Clist_snamestring = "CCMnet CCMnet CCMnet"
        MHproposal_name = "TNT"
        MHproposal_package = "CCMnet"
        covariate_list = covPattern
        
        inputs1 = c(rbind(c(0:(length(Prob_Distr_Params[[1]][[1]][[1]])-1)), rep(1,length(Prob_Distr_Params[[1]][[1]][[1]]))))  
        inputs2 = c(rbind(c(0:(length(Prob_Distr_Params[[1]][[1]][[2]])-1)), rep(2,length(Prob_Distr_Params[[1]][[1]][[2]]))))
        
        inputs = c(c(0,1,0,0), length(Prob_Distr_Params[[1]][[1]][[1]]) + length(Prob_Distr_Params[[1]][[1]][[2]]),
                   2*(length(Prob_Distr_Params[[1]][[1]][[1]])+length(Prob_Distr_Params[[1]][[1]][[2]])) + g$gal$n, inputs1, inputs2, covariate_list, c(4,2,4 + g$gal$n), c(1,2,2,2), covariate_list)
        eta0 = rep(-999.5,length(c(nedges[1], Prob_Distr_Params[[1]][[1]][[1]],Prob_Distr_Params[[1]][[1]][[2]],1,1)))
        
        mixing = c(0,0,0)
        edge_list = unlist(g$mel)
        dim(edge_list) = c(3,nedges[1])
        for (num_edge in c(1:nedges[1])) {
          if ((covariate_list[edge_list[1,num_edge]] == 1) && (covariate_list[edge_list[2,num_edge]] == 1)) {
            mixing[1] = mixing[1] + 1
          }
          if ((covariate_list[edge_list[1,num_edge]] == 1) && (covariate_list[edge_list[2,num_edge]] == 2)) {
            mixing[2] = mixing[2] + 1
          }
          if ((covariate_list[edge_list[1,num_edge]] == 2) && (covariate_list[edge_list[2,num_edge]] == 1)) {
            mixing[2] = mixing[2] + 1
          }
          if ((covariate_list[edge_list[1,num_edge]] == 2) && (covariate_list[edge_list[2,num_edge]] == 2)) {
            mixing[3] = mixing[3] + 1
          }
        }
        deg_dist_1 = tabulate(degree(g, gmode="graph")[which(covariate_list == 1)]+1)
        deg_dist_2 = tabulate(degree(g, gmode="graph")[which(covariate_list == 2)]+1)
        
        deg_dist_1 = c(tabulate(degree(g, gmode="graph")[which(covariate_list == 1)]+1), rep(0,max(0,length(Prob_Distr_Params[[1]][[1]][[1]])-length(deg_dist_1))))
        deg_dist_2 = c(tabulate(degree(g, gmode="graph")[which(covariate_list == 2)]+1), rep(0,max(0,length(Prob_Distr_Params[[1]][[1]][[2]])-length(deg_dist_2))))
        
        #Assume max degree of both node types is the same
        deg_dist_1 = c(deg_dist_1, rep(0,max(0,length(deg_dist_2)-length(deg_dist_1))))
        deg_dist_2 = c(deg_dist_2, rep(0,max(0,length(deg_dist_1)-length(deg_dist_2))))
        
        stats = c(nedges[1], deg_dist_1, deg_dist_2, mixing[c(2,3)])
        
        mean_vector = c(Prob_Distr_Params[[1]][[1]][[1]], Prob_Distr_Params[[1]][[1]][[2]],  Prob_Distr_Params[[2]][[1]], Prob_Distr_Params[[1]][[3]][1], Prob_Distr_Params[[1]][[3]][2], Prob_Distr_Params[[2]][[3]])

        if (remove_var_last_entry == TRUE) {
          inverse_var_x1 = solve(Prob_Distr_Params[[1]][[2]][[1]][-length(Prob_Distr_Params[[1]][[1]][[1]]),-length(Prob_Distr_Params[[1]][[1]][[1]])])
          inverse_var_x1 = rbind(inverse_var_x1,0)
          inverse_var_x1 = cbind(inverse_var_x1,0)
          
          inverse_var_x2 = solve(Prob_Distr_Params[[1]][[2]][[2]][-length(Prob_Distr_Params[[1]][[1]][[2]]),-length(Prob_Distr_Params[[1]][[1]][[2]])])
          inverse_var_x2 = rbind(inverse_var_x2,0)
          inverse_var_x2 = cbind(inverse_var_x2,0)
        } else {
          inverse_var_x1 = solve(Prob_Distr_Params[[1]][[2]][[1]])
          inverse_var_x2 = solve(Prob_Distr_Params[[1]][[2]][[2]])
        } 
        
        var_vector = c(c(inverse_var_x1),c(inverse_var_x2), Prob_Distr_Params[[2]][[2]])
        
        prob_type = c(2,2,0,0,1)
          
      } else {
          print("Error: No such distribution for degree distribution and mixing currently implemented.")
          print("Email ravi.goyal@mail.harvard.edu to add feature.")
          error = 1
      }
  } else if ((length(Network_stats) == 1) && (Network_stats == "DegMixing")) {
    if (Prob_Distr == "Normal") {
      
      if (class(Prob_Distr_Params[[1]][[1]]) != "numeric") {
        print("Error: Mean degree mixing should be a vector representing upper triangle of degree mixing matrix.")
        error = 1
      }
      if (class(Prob_Distr_Params[[1]][[2]]) != "matrix") {
        print("Error: Covariance of degree mixing matrix should be a matrix.")
        error = 1
      }
      if (dim(Prob_Distr_Params[[1]][[2]])[1] != dim(Prob_Distr_Params[[1]][[2]])[2]) {
        print("Error: Covariance matrix is not square.")
        error = 1
      }
      if (length(Prob_Distr_Params[[1]][[1]]) != dim(Prob_Distr_Params[[1]][[2]])[2]) {
        print("Error: mean vector and covariance matrix are not similar dimensions.")
        error = 1
      }
      
      
      m1 = matrix(c(1:max_degree), nrow = max_degree, ncol = max_degree)
      m1 = m1[upper.tri(m1, diag = TRUE)]
      
      m2 = t(matrix(c(1:max_degree), nrow = max_degree, ncol = max_degree))
      m2 = m2[upper.tri(m2, diag = TRUE)] 
      
      inputs = c(c(0,1,0), c(((max_degree+1)*max_degree), ((max_degree+1)*max_degree*.5), (((max_degree+1)*max_degree)+1)))
      inputs = c(inputs, m1, m2, max_degree)
      
      eta0 = rep(-999.5,length(c(nedges[1])) + .5*((max_degree+1)*max_degree))
      
      g_dmm = matrix(0,  nrow = max_degree, ncol = max_degree)
      edge_list = unlist(g$mel)
      dim(edge_list) = c(3,nedges[1])
      g_degree = degree(g, gmode = "graph")
      for (num_edge in c(1:nedges[1])) {
        deg1 = g_degree[edge_list[1,num_edge]]
        deg2 = g_degree[edge_list[2,num_edge]]
        if ((deg1 <= max_degree) && (deg2 <= max_degree)) {
          g_dmm[deg1, deg2] = g_dmm[deg1,deg2] + 1
          if (deg1 != deg2) {
            g_dmm[deg2, deg1] = g_dmm[deg2,deg1] + 1
          }
        }
      }
      
      stats = c(nedges[1],g_dmm[upper.tri(g_dmm, diag = TRUE)] )
      prob_type = c(0,0,1,0,1)
      
      mean_vector = Prob_Distr_Params[[1]][[1]] 

      if (remove_var_last_entry == TRUE) {
        inverse_var_x = solve(Prob_Distr_Params[[1]][[2]] [-length(mean_vector),-length(mean_vector)])
        inverse_var_x = rbind(inverse_var_x,0)
        inverse_var_x = cbind(inverse_var_x,0)
      } else {
        inverse_var_x = solve(Prob_Distr_Params[[1]][[2]])
      } 
            
      var_vector = c(inverse_var_x)
    } else {
      print("Error: No such distribution for degree mixing currently implemented.")
      print("Email ravi.goyal@mail.harvard.edu to add feature.")
      error = 1
    }
    Clist_nterms = 2 #Number of different terms
    Clist_fnamestring = "edges degmix"
    Clist_snamestring = "CCMnet CCMnet"
    MHproposal_name = "TNT"
    MHproposal_package = "CCMnet"
  } else if ((length(Network_stats) == 2) && (Network_stats[1] == c("DegMixing")) && (Network_stats[2] == c("Triangles"))  ||
             (length(Network_stats) == 2) && (Network_stats[1] == c("Triangles")) && (Network_stats[2] == c("DegMixing"))
             ) {
    if (Network_stats[1] == "Triangles") { #swap prob_distr_params
      Prob_Distr_Params_temp = Prob_Distr_Params[[1]]
      Prob_Distr_Params[[1]] = Prob_Distr_Params[[2]]
      Prob_Distr_Params[[2]] = Prob_Distr_Params_temp
    }
    if (class(Prob_Distr_Params[[1]][[1]]) != "numeric") {
      print("Error: Mean degree mixing should be a vector representing upper triangle of degree mixing matrix.")
      error = 1
    }
    if (class(Prob_Distr_Params[[1]][[2]]) != "matrix") {
      print("Error: Covariance of degree mixing matrix should be a matrix.")
      error = 1
    }
    if (dim(Prob_Distr_Params[[1]][[2]])[1] != dim(Prob_Distr_Params[[1]][[2]])[2]) {
      print("Error: Covariance matrix is not square.")
      error = 1
    }
    if (length(Prob_Distr_Params[[1]][[1]]) != dim(Prob_Distr_Params[[1]][[2]])[2]) {
      print("Error: mean vector and covariance matrix are not similar dimensions.")
      error = 1
    }
    if (length(Prob_Distr_Params[[2]][[1]]) != 1) {
      print("Error: Mean Triangles such be a single positive value.")
      error = 1
    }
    if (length(Prob_Distr_Params[[2]][[2]]) != 1) {
      print("Error: Variance of Triangles such be a single positive value.")
      error = 1
    }
    if ((Prob_Distr[1] == "Normal") && (Prob_Distr[2] == "Normal")) {
      m1 = matrix(c(1:max_degree), nrow = max_degree, ncol = max_degree)
      m1 = m1[upper.tri(m1, diag = TRUE)]
    
      m2 = t(matrix(c(1:max_degree), nrow = max_degree, ncol = max_degree))
      m2 = m2[upper.tri(m2, diag = TRUE)] 
    
      inputs = c(c(0,1,0), c(((max_degree+1)*max_degree), ((max_degree+1)*max_degree*.5), (((max_degree+1)*max_degree)+1)))
      inputs = c(inputs, m1, m2, max_degree, c(0,1,0))
    
      eta0 = rep(-999.5,length(c(nedges[1])) + .5*((max_degree+1)*max_degree) + 1)
    
      g_dmm = matrix(0,  nrow = max_degree, ncol = max_degree)
      edge_list = unlist(g$mel)
      dim(edge_list) = c(3,nedges[1])
      g_degree = degree(g, gmode = "graph")
      for (num_edge in c(1:nedges[1])) {
        deg1 = g_degree[edge_list[1,num_edge]]
        deg2 = g_degree[edge_list[2,num_edge]]
        if ((deg1 <= max_degree) && (deg2 <= max_degree)) {
          g_dmm[deg1, deg2] = g_dmm[deg1,deg2] + 1
          if (deg1 != deg2) {
            g_dmm[deg2, deg1] = g_dmm[deg2,deg1] + 1
          }
        }
      }
      
      stats = c(nedges[1],g_dmm[upper.tri(g_dmm, diag = TRUE)], triad.census(dat=g, g=NULL, mode = "graph")[4])
      prob_type = c(0,0,1,1,1)
    
      mean_vector = c(Prob_Distr_Params[[1]][[1]], Prob_Distr_Params[[2]][[1]] )
    
      if (remove_var_last_entry == TRUE) {
        inverse_var_x = solve(Prob_Distr_Params[[1]][[2]][-length(mean_vector[-1]),-length(mean_vector[-1])])
        inverse_var_x = rbind(inverse_var_x,0)
        inverse_var_x = cbind(inverse_var_x,0)
      } else {
        inverse_var_x = solve(Prob_Distr_Params[[1]][[2]])
      } 
        
      inverse_var_x = rbind(inverse_var_x,0)
      inverse_var_x = cbind(inverse_var_x,0)
      inverse_var_x[dim(inverse_var_x)[1], dim(inverse_var_x)[1]] = 1/Prob_Distr_Params[[2]][[2]]
    
      var_vector = c(inverse_var_x)
    } else {
      print("Error: No such distribution for degree mixing currently implemented.")
      print("Email ravi.goyal@mail.harvard.edu to add feature.")
      error = 1
    }    
    Clist_nterms = 3 #Number of different terms
    Clist_fnamestring = "edges degmix triangle"
    Clist_snamestring = "CCMnet CCMnet CCMnet"
    MHproposal_name = "TNT"
    MHproposal_package = "CCMnet"

  } else {
    print("Error: No such distribution for mixing currently implemented.")
    print("Email ravi.goyal@mail.harvard.edu to add feature.")
    error = 1   
  }
  
  
  if (error == 0) {
    
    numnetworks = 0 #MCMC_wrapper required
    Clist_dir = FALSE
    Clist_bipartite = FALSE
    maxedges = 200001
    verbose = FALSE
    
    BayesInference = 0 #Required for Bayesian Inference
    TranNet = NULL
    P = P
    Ia = NULL
    Il = NULL
    R_times = NULL
    beta_a = NULL
    beta_l = NULL
    
    NetworkForecast = 0  #Required for Network Forecasting
    evolution_rate_mean = 0
    evolution_rate_var = 0 
    
    if (statsonly == FALSE) {
      samplesize_v = rep(1, samplesize)
      burnin_v = c(burnin, rep(interval, samplesize))
      interval_v = rep(interval, samplesize)
    } else {
      samplesize_v = samplesize
      burnin_v = burnin
      interval_v = interval
    }
    
    statsmatrix <- c()
    newnetwork = list()

    for (sample_net in c(1:length(samplesize_v))) {
      samplesize = samplesize_v[sample_net]
      burnin = burnin_v[sample_net]
      interval = interval_v[sample_net]
        
      z <- .C("MCMC_wrapper", as.integer(numnetworks), 
            as.integer(nedges), as.integer(tails), as.integer(heads), 
            as.integer(Clist_n), as.integer(Clist_dir), as.integer(Clist_bipartite), 
            as.integer(Clist_nterms), as.character(Clist_fnamestring), 
            as.character(Clist_snamestring), as.character(MHproposal_name), 
            as.character(MHproposal_package), as.double(inputs), as.double(eta0), as.integer(samplesize), 
            s = as.double(rep(stats, samplesize)), as.integer(burnin), 
            as.integer(interval), newnwtails = integer(maxedges), 
            newnwheads = integer(maxedges), as.integer(verbose), 
            as.integer(NULL), 
            as.integer(NULL), 
            as.integer(NULL), 
            as.integer(NULL), 
            as.integer(NULL), 
            as.integer(FALSE), 
            as.integer(0), 
            as.integer(maxedges), status = integer(1), 
            as.integer(prob_type),   ###MOD ADDED RAVI
            as.integer(max_degree),
            as.double(mean_vector),
            as.double(var_vector),
            as.integer(BayesInference),
            as.integer(Trans_nedges),
            as.integer(Trans_networktails),
            as.integer(Trans_networkheads),
            as.double(Ia), 
            as.double(Il), 
            as.double(R_times), 
            as.double(beta_a), 
            as.double(beta_l),
            as.integer(NetworkForecast),
            as.double(evolution_rate_mean),
            as.double(evolution_rate_var),
            PACKAGE = "CCMnet")   
      nw = g
      statsmatrix <- rbind(statsmatrix, matrix(z$s, nrow = samplesize, ncol = length(stats), byrow = TRUE))
      newnetwork[[sample_net]] <- newnw.extract(nw, z, output = "network")
      
      edge_mat = unlist(newnetwork[[sample_net]]$mel)
      dim(edge_mat) = c(3,network.edgecount(newnetwork[[sample_net]]))
      
      nedges = c(network.edgecount(newnetwork[[sample_net]]),0,0)
      tails = edge_mat[1,]
      heads = edge_mat[2,]
      stats = statsmatrix[sample_net,]
      
      gc()
    }

    if (!(is.matrix(statsmatrix))) {
      len_statsmatrix = length(statsmatrix)
      statsmatrix = as.matrix(statsmatrix)
      dim(statsmatrix) = c(1,len_statsmatrix)
    }
    
    if ((length(Network_stats) == 1) && (Network_stats == "DegreeDist")){
      statsmatrix = statsmatrix[,-1]
      colnames(statsmatrix) = paste("Degree", c(0:(dim(statsmatrix)[2]-1)), sep = " ")
    } else if  ((length(Network_stats) == 1) && (Network_stats == "Density")) {
      statsmatrix = statsmatrix[,1]/choose(population,2)
    } else if ((length(Network_stats) == 2) && (Network_stats[1] == "DegreeDist") && (Network_stats[2] == "Mixing")) {
      statsmatrix = statsmatrix[,-1]
    } else if ((length(Network_stats) == 2) && (Network_stats[1] == "Mixing") && (Network_stats[2] == "DegreeDist")) {
      statsmatrix = statsmatrix[,-1]
    } else if ((length(Network_stats) == 1) && (Network_stats == "DegMixing")) {
      statsmatrix = statsmatrix[,-1]
      colnames(statsmatrix) = paste("Edges", paste(m1, m2, sep = "-"), sep = " ")      
    } else if  ((length(Network_stats) == 2) && (Network_stats[1] == c("DegMixing")) && (Network_stats[2] == c("Triangles"))) {
      statsmatrix = statsmatrix[,-1]
      colnames(statsmatrix) = c(paste("Edges", paste(m1, m2, sep = "-"), sep = " "), "Triangles")
    } else if  ((length(Network_stats) == 2) && (Network_stats[1] == c("Triangles")) && (Network_stats[2] == c("DegMixing"))) {
      statsmatrix = statsmatrix[,-1]
      statsmatrix = cbind(statsmatrix[,dim(statsmatrix)[2]], statsmatrix[,-dim(statsmatrix)[2]])
      colnames(statsmatrix) = c("Triangles", paste("Edges", paste(m1, m2, sep = "-"), sep = " "))
    } else {
      statsmatrix = statsmatrix[,-1]
    }
    
    return(list(statsmatrix,newnetwork))
  } else {
    return(list(NULL, NULL))
  }
  
}

############################
############################
############################
############################
############################
############################
############################
############################
############################
############################

bi_modal_constr <- function(Network_stats, Prob_Distr, Prob_Distr_Params, 
                             samplesize, burnin, interval,
                             statsonly, P,
                             population, covPattern, remove_var_last_entry) {
  
  error = 0
  
  if(is.null(covPattern)) {
    covPattern = rep(1,population[1] + population[2])
  }
  
  G_max_degree_bool = FALSE
  ER_prob = .05
  
  if ((length(Network_stats) == 1) && (Network_stats == "DegreeDist")){
    max_degree_1 = length(Prob_Distr_Params[[1]][[1]])-1
    max_degree_2 = length(Prob_Distr_Params[[2]][[1]])-1
    max_degree = min(max_degree_1, max_degree_2)
    max_degree_f = max(max_degree_1, max_degree_2)
    
    len_deg_dist = length(Prob_Distr_Params[[1]][[1]])
    Prob_Distr_Params[[1]][[1]] = as.matrix(Prob_Distr_Params[[1]][[1]])
    dim(Prob_Distr_Params[[1]][[1]]) = c(1, len_deg_dist)
    
    len_deg_dist = length(Prob_Distr_Params[[2]][[1]])
    Prob_Distr_Params[[2]][[1]] = as.matrix(Prob_Distr_Params[[2]][[1]])
    dim(Prob_Distr_Params[[2]][[1]]) = c(1, len_deg_dist)      
  } else if ((length(Network_stats) == 2) && (Network_stats[1] == "DegreeDist") && (Network_stats[2] == "Mixing")) {
    max_degree_1 = dim(Prob_Distr_Params[[1]][[1]])[2]-1
    max_degree_2 = dim(Prob_Distr_Params[[2]][[1]])[2]-1
    max_degree = min(max_degree_1, max_degree_2)
    max_degree_f = max(max_degree_1, max_degree_2)
  } else if ((length(Network_stats) == 2) && (Network_stats[1] == "Mixing") && (Network_stats[2] == "DegreeDist")) {
    max_degree_1 = length(Prob_Distr_Params[[1]][[1]])-1
    max_degree_2 = length(Prob_Distr_Params[[2]][[1]])-1
    max_degree = min(max_degree_1, max_degree_2)
    max_degree_f = max(max_degree_1, max_degree_2)
  } else {
    max_degree_f = max_degree = sum(population) - 1
  }
  
  if ((length(Network_stats) == 1) && (Network_stats == "Density")){
    if (is.null(P)) {
      print("Generating Random Initial Network...")
      
      g <-network.initialize(population[1] + population[2],directed = FALSE, bipartite=population[1])
      g %v% "CovAttribute" = covPattern
      
      add.edges(g, tail = 1, head = population[1] + 1)
      
      print("COMPLETED: Generated Random Initial Network")
      P = g
    } else {
      g = P
    }
  }

  subpopulation = tabulate(covPattern[c(1:population[1])])
  Gen_Net_counter = 1
  if (is.null(P)) {
    print("Generating Random Initial Network...")
    
    g <-network.initialize(population[1] + population[2],directed = FALSE, bipartite=population[1])
    g %v% "CovAttribute" = covPattern
    
    current_node = 1
    for (i in c(1:(dim(Prob_Distr_Params[[1]][[1]])[2]))) { #Go through mean degree distribution - by column
      for (j in c(1:(dim(Prob_Distr_Params[[1]][[1]])[1]))) { #by rows
        num_nodes_deg_i = floor(Prob_Distr_Params[[1]][[1]][j,i]/sum(Prob_Distr_Params[[1]][[1]][j,]) * subpopulation[j]) #i - 1
        if (num_nodes_deg_i > 0) {
          for (k in c(1:num_nodes_deg_i)) { # by number of nodes
            if (i > 1)  {
              if (network.edgecount(g) > 0) {
                valid_population_2 = which(degree(g, gmode = 'graph')[c((population[1]+1):(population[1]+population[2]))] < max_degree)
              } else {
                valid_population_2 = c(1:population[2])
              }
              nodes_type_2 = sample(valid_population_2, size = (i-1), replace = FALSE, prob = NULL)
              add.edges(g, tail = rep(current_node, i-1), head = population[1] + nodes_type_2)
            }
            current_node = current_node + 1
          }           
        }
      }
    }
    print("COMPLETED: Generated Random Initial Network")
    P = g
  } else {
    g = P
  }
  
  if (is.null(P)) {
    P = g
  }
  
  #max(degree(g, gmode='graph'))
  
  max_degree = max_degree_f
  
  P_edge_mat = unlist(P$mel)
  dim(P_edge_mat) = c(3,network.edgecount(P))
  Trans_nedges = c(network.edgecount(P),0,0)
  Trans_networktails = P_edge_mat[1,]
  Trans_networkheads = P_edge_mat[2,]
  
  edge_mat = unlist(g$mel)
  dim(edge_mat) = c(3,network.edgecount(g))
  
  nedges = c(network.edgecount(g),0,0)
  tails = edge_mat[1,]
  heads = edge_mat[2,]
  Clist_n = network.size(g)
  
  if ((length(Network_stats) == 1) && (Network_stats == "DegreeDist")) {
    if (length(Prob_Distr_Params[[1]][[1]]) < 2) {
      print("Error: length of mean vector of type 1 is less than 2")
      error = 1
    }
    if (length(Prob_Distr_Params[[2]][[1]]) < 2) {
      print("Error: length of mean vector of type 2 is less than 2")
      error = 1
    }
    if (Prob_Distr == "Normal") {
      mean_vector = c(Prob_Distr_Params[[1]][[1]],Prob_Distr_Params[[2]][[1]])
      
      var_vector = c()
      inverse_var_x1 = solve(Prob_Distr_Params[[1]][[2]][-length(Prob_Distr_Params[[1]][[1]]),-length(Prob_Distr_Params[[1]][[1]])])
      inverse_var_x1 = rbind(inverse_var_x1,0)
      inverse_var_x1 = cbind(inverse_var_x1,0)
      
      var_vector = c(var_vector, c(inverse_var_x1))
      
      inverse_var_x2 = solve(Prob_Distr_Params[[2]][[2]][-length(Prob_Distr_Params[[2]][[1]]),-length(Prob_Distr_Params[[2]][[1]])])
      inverse_var_x2 = rbind(inverse_var_x2,0)
      inverse_var_x2 = cbind(inverse_var_x2,0)
      
      var_vector = c(var_vector, c(inverse_var_x2))
      
      prob_type = c(1,0,0,0,1)
      
    } else if (Prob_Distr == "DirMult") {
      mean_vector = c(Prob_Distr_Params[[1]][[1]],Prob_Distr_Params[[2]][[1]])
      var_vector = c(0,0)
      
      prob_type = c(2,0,0,0,1)
    } else if (Prob_Distr == "NP") {
      mean_vector = Prob_Distr_Params[[3]][[1]]
      var_vector = c(0,0)
      
      prob_type = c(99,0,0,0,1)
    } else {
      print("Error: No such distribution for degree distribution currently implemented.")
      print("Email ravi.goyal@mail.harvard.edu to add feature.")
      error = 1
    }
    if (error == 0) {
      Clist_nterms = 3 #Number of different terms
      Clist_fnamestring = "edges b1degree b2degree"
      Clist_snamestring = "CCMnet CCMnet CCMnet CCMnet"
      MHproposal_name = "TNT"
      MHproposal_package = "CCMnet"
      
      inputs =  c(c(0,1,0), 
                  c(0,length(Prob_Distr_Params[[1]][[1]]),length(Prob_Distr_Params[[1]][[1]]),c(0:(length(Prob_Distr_Params[[1]][[1]])-1))),
                  c(0,length(Prob_Distr_Params[[2]][[1]]),length(Prob_Distr_Params[[2]][[1]]),c(0:(length(Prob_Distr_Params[[2]][[1]])-1)))
      )
      
      stats = c(summary(g ~ edges), 
                as.numeric(summary(g ~ b1degree(c(0:max_degree)))),
                as.numeric(summary(g ~ b2degree(c(0:max_degree))))
      )
      
      eta0 = rep(-999.5, length(stats))
    }
  } else if ((length(Network_stats) == 1) && (Network_stats == "Density")) {
    if (Prob_Distr == "Normal") {
      prob_type = c(0,0,0,0,1)
      mean_vector = c(Prob_Distr_Params[[1]][[1]],Prob_Distr_Params[[1]][[1]])
      var_vector = c(Prob_Distr_Params[[1]][[2]], Prob_Distr_Params[[1]][[2]])
      if (length(Prob_Distr_Params[[1]][[1]]) != 1) {
        print("Error: mean value for network density is one positive value")
        error = 1
      }
      if (length(Prob_Distr_Params[[1]][[2]]) != 1) {
        print("Error: variance for network density is one positive value")
        error = 1
      }
    } else if (Prob_Distr == "NP") {
      prob_type = c(0,0,0,0,99)
      mean_vector = Prob_Distr_Params[[1]][[1]]
      var_vector = c(0,0)
    } else {
      print("Error: No such distribution for degree distribution currently implemented.")
      print("Email ravi.goyal@mail.harvard.edu to add feature.")
      error = 1
    }
    if (error == 0) {
      Clist_nterms = 2 #Number of different terms
      Clist_fnamestring = "edges nfstab"
      Clist_snamestring = "CCMnet CCMnet"
      inputs = c(0,1,0,0,1,0)
      eta0 = c(-999.5, -999.5)
      stats = c(nedges[1],nedges[1])
      MHproposal_name = "TNT"
      MHproposal_package = "CCMnet"
    }    
  } else if ((length(Network_stats) == 1) && (Network_stats == "Mixing")) {
    print("Error: No such distribution for mixing currently implemented.")
    print("Email ravi.goyal@mail.harvard.edu to add feature.")
    error = 1
  } else if (((length(Network_stats) == 2) && (Network_stats[1] == "DegreeDist") && (Network_stats[2] == "Mixing")) ||
               ((length(Network_stats) == 2) && (Network_stats[1] == "Mixing") && (Network_stats[2] == "DegreeDist")))  {
    if (Network_stats[1] == "Mixing") { #swap prob_distr_params
      Prob_Distr_Params_temp = Prob_Distr_Params[[1]]
      Prob_Distr_Params[[1]] = Prob_Distr_Params[[2]]
      Prob_Distr_Params[[2]] = Prob_Distr_Params_temp
      
      Prob_Distr_Params_temp = Prob_Distr_Params[[2]]
      Prob_Distr_Params[[2]] = Prob_Distr_Params[[3]]
      Prob_Distr_Params[[3]] = Prob_Distr_Params_temp
    }
    if ((Prob_Distr[1] == "Normal") && ((Prob_Distr[2] == "Normal"))) {

      mean_vector = c(c(t(Prob_Distr_Params[[1]][[1]])),c(t(Prob_Distr_Params[[2]][[1]])),c(t(Prob_Distr_Params[[3]][[1]])))
      
      var_vector = c()
      for (i in c(1:dim(Prob_Distr_Params[[1]][[1]])[1])) {
        inverse_var_x1 = solve(Prob_Distr_Params[[1]][[2]][[i]][-length(Prob_Distr_Params[[1]][[1]][i,]),-length(Prob_Distr_Params[[1]][[1]][i,])])
        inverse_var_x1 = rbind(inverse_var_x1,0)
        inverse_var_x1 = cbind(inverse_var_x1,0)
          
        var_vector = c(var_vector, c(inverse_var_x1))
      }
        
      for (i in c(1:dim(Prob_Distr_Params[[2]][[1]])[1])) {
        inverse_var_x2 = solve(Prob_Distr_Params[[2]][[2]][[i]][-length(Prob_Distr_Params[[2]][[1]][i,]),-length(Prob_Distr_Params[[2]][[1]][i,])])
        inverse_var_x2 = rbind(inverse_var_x2,0)
        inverse_var_x2 = cbind(inverse_var_x2,0)
        
        var_vector = c(var_vector, c(inverse_var_x2))
        
      }
        
      inverse_var_x3 = matrix(0, ncol = length(c(Prob_Distr_Params[[3]][[1]])), nrow =  length(c(Prob_Distr_Params[[3]][[1]])))
      for (i in c(1:dim(Prob_Distr_Params[[3]][[1]])[1])){
        r_index = c((1+(i-1)*dim(Prob_Distr_Params[[3]][[1]])[1]):(i*dim(Prob_Distr_Params[[3]][[1]])[1]))
        inverse_eta_i = Prob_Distr_Params[[3]][[2]][r_index,r_index] 
        inverse_eta_i = solve(inverse_eta_i[-dim(Prob_Distr_Params[[3]][[1]])[1],-dim(Prob_Distr_Params[[3]][[1]])[1]])
        inverse_eta_i = rbind(inverse_eta_i,0)
        inverse_eta_i = cbind(inverse_eta_i,0)
        inverse_var_x3[r_index,r_index] = inverse_eta_i
      }
      var_vector = c(var_vector, c(inverse_var_x3))    
        
      prob_type = c(1,1,0,0,1)
      
    } else if ((Prob_Distr[1] == "DirMult") && ((Prob_Distr[2] == "DirMult"))) {
      mean_vector = c(c(t(Prob_Distr_Params[[1]][[1]])),c(t(Prob_Distr_Params[[2]][[1]])),c(t(Prob_Distr_Params[[3]][[1]])))
      var_vector = rep(1, dim(Prob_Distr_Params[[1]][[1]])[2]^2 *  dim(Prob_Distr_Params[[1]][[1]])[1] * 2 +  dim(Prob_Distr_Params[[1]][[1]])[1]^4)
      prob_type = c(2,2,0,0,1)      
    } else {
      print("Error: No such distribution for degree distribution and mixing currently implemented.")
      print("Email ravi.goyal@mail.harvard.edu to add feature.")
      error = 1
    }

    if (error == 0) {
      Clist_nterms = 4 #Number of different terms
      Clist_fnamestring = "edges b1degree_by_attr b2degree_by_attr mix"
      Clist_snamestring = "CCMnet CCMnet CCMnet CCMnet"
      MHproposal_name = "TNT"
      MHproposal_package = "CCMnet"
      
      b1deg_by_attr_input = c()
      b2deg_by_attr_input = c()
      for (i in c(1:dim(Prob_Distr_Params[[1]][[1]])[1])) {
        b1deg_by_attr_input = c(b1deg_by_attr_input, c(rbind(c(0:(dim(Prob_Distr_Params[[1]][[1]])[2]-1)),i)))
        b2deg_by_attr_input = c(b2deg_by_attr_input, c(rbind(c(0:(dim(Prob_Distr_Params[[1]][[1]])[2]-1)),i)))
      }
      
      inputs =  c(c(0,1,0), 
                  c(0,dim(Prob_Distr_Params[[1]][[1]])[2]*dim(Prob_Distr_Params[[1]][[1]])[1],dim(Prob_Distr_Params[[1]][[1]])[2]*dim(Prob_Distr_Params[[1]][[1]])[1]*2+population[1]),
                  b1deg_by_attr_input, covPattern[c(1:population[1])],
                  c(0,dim(Prob_Distr_Params[[1]][[1]])[2]*dim(Prob_Distr_Params[[1]][[1]])[1],dim(Prob_Distr_Params[[1]][[1]])[2]*dim(Prob_Distr_Params[[1]][[1]])[1]*2+population[2]),
                  b2deg_by_attr_input, covPattern[c((population[1]+1):(population[1] + population[2]))], 
                  c(dim(Prob_Distr_Params[[1]][[1]])[1]^2, dim(Prob_Distr_Params[[1]][[1]])[1]^2, dim(Prob_Distr_Params[[1]][[1]])[1]^2*2 + sum(population)),
                  rep(c(1:dim(Prob_Distr_Params[[1]][[1]])[1]), dim(Prob_Distr_Params[[1]][[1]])[1]),
                  rep(c((dim(Prob_Distr_Params[[1]][[1]])[1]+1):(dim(Prob_Distr_Params[[1]][[1]])[1]*2)), each=dim(Prob_Distr_Params[[1]][[1]])[1]),
                  c(covPattern + c(rep(0,population[1]),rep(dim(Prob_Distr_Params[[1]][[1]])[1], population[2]) ))
      )
      
      stats = c(summary(g ~ edges), 
                as.numeric(summary(g ~ b1degree(c(0:max_degree), by = "CovAttribute"))),
                as.numeric(summary(g ~ b2degree(c(0:max_degree), by = "CovAttribute"))),
                as.numeric(summary(g ~ nodemix(attrname = "CovAttribute")))
      )
      
      eta0 = rep(-999.5, length(stats))
    }
    
  } else {
    print("Error: No such Network Statistics currently implemented.")
    print("Email ravi.goyal@mail.harvard.edu to add feature.")
    error = 1   
  }
  
  
  if (error == 0) {
    
    numnetworks = 0 #MCMC_wrapper required
    Clist_dir = FALSE
    Clist_bipartite = population[1]
    maxedges = 200001
    verbose = FALSE
    
    BayesInference = 0 #Required for Bayesian Inference
    TranNet = NULL
    P = P
    Ia = NULL
    Il = NULL
    R_times = NULL
    beta_a = NULL
    beta_l = NULL
    
    NetworkForecast = 0  #Required for Network Forecasting
    evolution_rate_mean = 0
    evolution_rate_var = 0 
    
    if (statsonly == FALSE) {
      samplesize_v = rep(1, samplesize)
      burnin_v = c(burnin, rep(interval, samplesize))
      interval_v = rep(interval, samplesize)
    } else {
      samplesize_v = samplesize
      burnin_v = burnin
      interval_v = interval
    }
    
    statsmatrix <- c()
    newnetwork = list()
    
    for (sample_net in c(1:length(samplesize_v))) {
      samplesize = samplesize_v[sample_net]
      burnin = burnin_v[sample_net]
      interval = interval_v[sample_net]
      
      z <- .C("MCMC_wrapper", as.integer(numnetworks), 
              as.integer(nedges), as.integer(tails), as.integer(heads), 
              as.integer(Clist_n), as.integer(Clist_dir), as.integer(Clist_bipartite), 
              as.integer(Clist_nterms), as.character(Clist_fnamestring), 
              as.character(Clist_snamestring), as.character(MHproposal_name), 
              as.character(MHproposal_package), as.double(inputs), as.double(eta0), as.integer(samplesize), 
              s = as.double(rep(stats, samplesize)), as.integer(burnin), 
              as.integer(interval), newnwtails = integer(maxedges), 
              newnwheads = integer(maxedges), as.integer(verbose), 
              as.integer(NULL), 
              as.integer(NULL), 
              as.integer(NULL), 
              as.integer(NULL), 
              as.integer(NULL), 
              as.integer(FALSE), 
              as.integer(0), 
              as.integer(maxedges), status = integer(1), 
              as.integer(prob_type),   ###MOD ADDED RAVI
              as.integer(max_degree),
              as.double(mean_vector),
              as.double(var_vector),
              as.integer(BayesInference),
              as.integer(Trans_nedges),
              as.integer(Trans_networktails),
              as.integer(Trans_networkheads),
              as.double(Ia), 
              as.double(Il), 
              as.double(R_times), 
              as.double(beta_a), 
              as.double(beta_l),
              as.integer(NetworkForecast),
              as.double(evolution_rate_mean),
              as.double(evolution_rate_var),
              PACKAGE = "CCMnet")   
      nw = g
      statsmatrix <- rbind(statsmatrix, matrix(z$s, nrow = samplesize, ncol = length(stats), byrow = TRUE))
      newnetwork[[sample_net]] <- newnw.extract(nw, z, output = "network")
      
      edge_mat = unlist(newnetwork[[sample_net]]$mel)
      dim(edge_mat) = c(3,network.edgecount(newnetwork[[sample_net]]))
      
      nedges = c(network.edgecount(newnetwork[[sample_net]]),0,0)
      tails = edge_mat[1,]
      heads = edge_mat[2,]
      stats = statsmatrix[sample_net,]
      
      gc()
    }
    
    if (!(is.matrix(statsmatrix))) {
      len_statsmatrix = length(statsmatrix)
      statsmatrix = as.matrix(statsmatrix)
      dim(statsmatrix) = c(1,len_statsmatrix)
    }
    
     if ((length(Network_stats) == 1) && (Network_stats == "DegreeDist")){
       statsmatrix = statsmatrix[,-1]
       colnames(statsmatrix) = c(paste("Type1_Degree", c(0:(dim(statsmatrix)[2]/2-1)), sep = " "),
                                 paste("Type1_Degree", c(0:(dim(statsmatrix)[2]/2-1)), sep = " ")
       )                                 
     } else if  ((length(Network_stats) == 1) && (Network_stats == "Density")) {
       #statsmatrix[,1] = statsmatrix[,1]/choose(population[1],2)
       #statsmatrix[,2] = statsmatrix[,2]/choose(population[2],2)
       statsmatrix = statsmatrix[,c(1:2)]
     } else if ((length(Network_stats) == 2) && (Network_stats[1] == "DegreeDist") && (Network_stats[2] == "Mixing")) {
       statsmatrix = statsmatrix[,-1]
     } else if ((length(Network_stats) == 2) && (Network_stats[1] == "Mixing") && (Network_stats[2] == "DegreeDist")) {
       statsmatrix = statsmatrix[,-1]
     } else {
       statsmatrix = statsmatrix[,-1]
     }
    
    return(list(statsmatrix,newnetwork))
  } else {
    return(list(NULL, NULL))
  }
  
}

