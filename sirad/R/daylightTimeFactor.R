daylightTimeFactor <-
  function(lat,i) {
    ws <- acos(-tan(lat)*tan(solarDecl(i)))
    ws
  }

