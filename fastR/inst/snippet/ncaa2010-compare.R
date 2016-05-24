compareTeams <-  
  function(team1,team2,model,abilities=BTabilities(model)) { 
    a <- abilities[team1,1]
    b <- abilities[team2,1]
    return(ilogit(a-b))
} 
compareTeams('Kansas','Kentucky',ab=ratings)
compareTeams('Michigan St.','Butler',ab=ratings)
compareTeams('Butler','Duke',ab=ratings)
