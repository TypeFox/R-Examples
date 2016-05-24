SplitSteal <- rbind( 
  do(187) * data.frame(agegroup = "Under40", decision = "Split"),
  do(195) * data.frame(agegroup = "Under40", decision = "Steal"),
  do(116) * data.frame(agegroup = "Over40",  decision = "Split"),
  do(76)  * data.frame(agegroup = "Over40",  decision = "Steal")
  )

