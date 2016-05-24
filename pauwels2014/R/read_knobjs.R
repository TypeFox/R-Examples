if(getRversion() >= "2.15.1")  utils::globalVariables(c("exps"))
read_knobjs <- function(objs){
  ## Arrange the results in a table
	lapply(objs, 
    FUN = function(obj){
      knobj <- eval(parse(text = obj))
      if( is.null(knobj$datas[[ length(knobj$datas) ]]$thetas) ){
        knobj$datas[[length(knobj$datas)]] <- NULL
      }
      risks_mean <- sapply(1:length(knobj$datas), FUN = function(k){
          thetas <- knobj$datas[[k]]$thetas
          theta <- apply(thetas, 2, mean)
          theta_trans <- knobj$transform_params(theta)
          risk_theta_fun(theta_trans, 
            knobj$global_parameters$true_params, 
            knobj$global_parameters$n_params
          )
        }
      )
      
      cost <- 0
      costs <- 0:(length(knobj$datas)-1)
      temp <- paste(exps[,3],exps[,1])
    
      for(i in 2:length(knobj$datas)){    
        cost <- cost + exps$Cost[temp == knobj$experiments[i]]
        costs[i] <- cost
      }
      rbind(risks_mean, costs)
      #risks_mean
    }
  )
}
