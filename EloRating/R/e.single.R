e.single <-
function(ELO1old, ELO2old, outcome, k=100) {
  
  # outcome must be one of the following:
  # "1" = first individual wins and second looses
  # "2" = second individual wins and first looses
  # "0" = interaction ends in a draw/tie (no winner and no looser) 
  
  # make sure ELO ratings are given as numerics
  ELO_1 <- as.numeric(as.character(ELO1old))
  ELO_2 <- as.numeric(as.character(ELO2old))
  
  # calculates the difference between the two ratings
  ELO_diff <- (ELO_1 - ELO_2)
  
  # z score based on fixed SD=200 (see Elo 1978)   
  z_score <- ELO_diff/(200*sqrt(2))  
  
  # calculates the winning probabilty
  p_win <- pnorm(z_score)
  
  # product of winning probabilty and k-factor
  kp <- k * p_win
  
  # the actual updating calculations
  if(outcome==1) { ELO1new <- ELO_1 - kp + k
                   ELO2new <- ELO_2 + kp - k
  }
  if(outcome==2) { ELO1new <- ELO_1 - kp
                   ELO2new <- ELO_2 + kp
  }
  if(outcome==0) { ELO1new <- ELO_1 - kp + 0.5 * k
                   ELO2new <- ELO_2 + kp - 0.5 * k
  }
  # returning the updated ratings
  # ratings are rounded to integers
  return(round(c(ELO1new, ELO2new), 0))
}
