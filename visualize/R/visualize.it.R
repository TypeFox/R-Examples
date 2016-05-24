visualize.it <-
function(dist='norm', stat = c(0,1), params = list(mu = 0, sd = 1), section = "lower", strict = c(0,1)) {
  dist = visualize.distributions[[casefold(dist)]]
  if(is.null(dist)) stop("Distribution not found.\n")
  if(length(params) != dist$params) stop("Invalid amount of parameters provided.\n")
    
  if(length(stat)>1 & (section != "bounded" & section != "tails")){ 
    stop(paste('Supplied stat length > 1 and section="',section,'" requires one statistic. Please resubmit with stat=your_test_statistic.\n'))
  }
  else if(length(stat)<2 & (section == "bounded" | section == "tails")){ 
    stop(paste('Supplied stat length < 2 and section="',section,'" requires two statistics. Please resubmit with stat=c(lower,upper).\n'))
  }

  #distribution specific graphing call.
  if(dist$type == "continuous"){ visualize.continuous(dist, stat, params, section)}
  else{     
    
    #Ensures array is inbounds according to conditions
    inequality = if(strict[[1]] == 0) {"equal to"} else {"strict"}
    if(length(strict)<2 & (section == "bounded" | section == "tails")){ 
      strict = c(strict[[1]],strict[[1]]) 
      cat(paste("Supplied strict length < 2, setting inequalities to ", inequality, " inequality.\n"))
    }
    
    visualize.discrete(dist, stat, params, section, strict)
  } 
}
