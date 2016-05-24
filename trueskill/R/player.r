#' Player: class to hold the rank, skill and names of players
#' @param rank rank of player in the match outcome
#' @param skill skill of player represented by Gaussian object e.g. Gaussian(mu = 25, sigma = 25/3)
#' @param name name the player for display purposes
#' @param ... Reference Class inheritance
Player <- setRefClass("Player",
  fields = list(rank = "numeric", skill = "Gaussian", name = "character"),
  methods = list(
    initialize = function(rank = 1, skill = Gaussian$new(), name = "") {
      .self$skill <- skill
      .self$rank <- rank
      .self$name <- name
    }, 
    UpdateSkill = function(mu, sigma) {
      .self$skill <- Gaussian$new(mu = mu, sigma = sigma)
    },
    show = function() {
      print(sprintf("[rank, skill, player]: [%s, [(%s, %s), (%s, %s)], %s]", 
        rank, round(skill$MuSigma()[1], 3), round(skill$MuSigma()[2], 3), round(skill$pi, 3), round(skill$tau, 3), name))
    }    
  )
)

#' Adjust the skills of a list of players.
#' 
#' @param players is a list of player objects, for all the players who
#' participated in a single game.  A 'player object' is any object with
#' a 'skill' attribute (a Gaussian) and a 'rank' field.
#' Lower ranks are better; the lowest rank is the overall winner of the
#' game.  Equal ranks mean that the two players drew.
#' 
#' This function updates all the 'skill' attributes of the player
#' objects to reflect the outcome of the game.  The input list is not
#' altered.
#
#' Create all the variable nodes in the graph.  "Teams" are each a
#' single player; there's a one-to-one correspondence between players
#' and teams.  (It would be straightforward to make multiplayer
#' teams, but it's not needed for my current purposes.) 
AdjustPlayers = function(players, parameters = Parameters()) {
  
  # prior factor, gamma
  # likelihood factor, beta
  # truncate factor, epsilon
  
  SortPlayers = function(players) {
    GetRank = function(x) return(x$rank)
    sorted_players <- players[order(unlist(Map(GetRank, players)))] 
    return(sorted_players)
  }
  
  players <- SortPlayers(players)

  # skill, performance and team variable, as well as the difference variable

  # create var from player object, then map over list of players
  GenSkill <- function(player, varname) return(Variable$new(value = player$skill, name = paste(varname, player$name, sep = "")))
  GenVar <- function(player, varname) return(Variable$new(name = paste(varname, player$name, sep = "")))
                               
  # Create each layer of factor nodes.  At the top we have priors
  # initialized to the player's current skill estimate.
  
  GenPriorFactor <- function(skill_var, player, gamma) {
    new_sigma <- sqrt(player$skill$sigma()^ 2 + gamma ^ 2)
    param <- Gaussian$new(mu = player$skill$mu(), sigma = new_sigma)
    return(PriorFactor$new(variable = skill_var, param = param, name = paste("PF", player$name, sep = "")))
  }
  
  GenLikelihoodFactor <- function(skill_var, perf_var, player, beta) {
    return(LikelihoodFactor$new(skill_var, perf_var, beta ^ 2, name = paste("LF", player$name, sep = "")))
  }                                                               
  
  GenSumFactor <- function(team_var, perf_var, player) {
    return(SumFactor$new(sum_variable = team_var, term_variables = list(perf_var), coeff = list(1), name = paste("SF", player$name, sep = "")))
  }                                                                                    
  
  # Team Diff SumFactor
  GenTeamDiff <- function(diff_var, match_list) {
    match_name <- paste("SF", match_list[[1]]$name, "vs.", match_list[[2]]$name, sep = " ")
    return(SumFactor$new(sum_variable = diff_var, term_variables = match_list, coeff = list(1, -1), name = match_name))
  }
  
  # zip teams less last team, with teams less first team (t1, t2, t3) with (t2, t3, t4)
  GenTeamDiffList <- function(diff_vars, team_vars) {
    match_list <- mapply(list, team_vars[-length(team_vars)], team_vars[-1], SIMPLIFY = F)
    return(mapply(GenTeamDiff, diff_vars, match_list, SIMPLIFY = F))
  }

  GenTruncateFactor <- function(diff_var, player1, player2, epsilon) {
    
    if(player1$rank == player2$rank) { V <- Vdraw; W <- Wdraw }
    else { V <- Vwin ; W <- Wwin }
    
    return(TruncateFactor$new(diff_var, V, W, epsilon, name = paste("TF", player1$name, player2$name, sep = "_")))	  
  }
  
  skill_vars <- mapply(GenSkill, players, "skill")
  perf_vars <- mapply(GenVar, players, "perf") 
  team_vars <- mapply(GenVar, players, "team")         
  diff_vars <- mapply(GenVar, players[-length(players)], "diff")
                                               
  skill <- mapply(GenPriorFactor, skill_vars, players, parameters$gamma)
  skill_to_perf <- mapply(GenLikelihoodFactor, skill_vars, perf_vars, players, parameters$beta)
  perf_to_team <- mapply(GenSumFactor, team_vars, perf_vars, players)
  team_diff <- GenTeamDiffList(diff_vars, team_vars)

  # At the bottom we connect adjacent teams with a 'win' or 'draw'
  # factor, as determined by the rank values.
  
  trunc <- mapply(GenTruncateFactor, diff_vars, players[-length(players)], players[-1], parameters$epsilon)
  
  # Start evaluating the graph by pushing messages 'down' from the
  # priors.
  Map(function(x) x$Start(), skill)
  Map(function(x) x$UpdateValue(), skill_to_perf)
  Map(function(x) x$UpdateSum(), perf_to_team)
  
  # Because the truncation factors are approximate, we iterate,
  # adjusting the team performance (t) and team difference (d)
  # variables until they converge.  In practice this seems to happen
  # very quickly, so I just do a fixed number of iterations.
  #
  # This order of evaluation is given by the numbered arrows in Figure
  # 1 of the Herbrich paper.

  z <- 1
  while (z <= 10) {
    Map(function(x) x$UpdateSum(), team_diff)                        # arrows (1) and (4)
    Map(function(x) x$Update(), trunc)                               # arrows (2) and (5)
    Map(function(x) { x$UpdateTerm(1); x$UpdateTerm(2) }, team_diff) # arrows (3) and (6)
    z <- z + 1
  }
    		
  # Now we push messages back up the graph, from the teams back to the
  # player skills.

  Map(function(x) x$UpdateTerm(1), perf_to_team)
  Map(function(x) x$UpdateMean(), skill_to_perf)
  
  # Finally, the players' new skills are the new values of the s
  # variables.
                                 
  UpdateSkill <- function(player, skill_var) {
     player$UpdateSkill(skill_var$value$mu(), skill_var$value$sigma())     
  }
  
  mapply(UpdateSkill, players, skill_vars, SIMPLIFY = F)
  return(players)      
}

#' display list of players nicer
#' @param list a list of player objects
PrintList = function(list) {
  for(i in 1:length(list)) {
    print(list[[i]])
  }
}

