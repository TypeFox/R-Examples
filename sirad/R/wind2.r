wind2 <-
  function(uz,meah) {
    u2 <- round(uz*4.87/log(67.8*meah-5.42),digits=3)
    u2
  }