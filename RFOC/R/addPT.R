`addPT` <-
function(MEC, pch=5)
  {

    focpoint(MEC$P$az, MEC$P$dip, pch=5, lab="P", UP=MEC$UP)
    focpoint(MEC$T$az, MEC$T$dip, pch=5, lab="T", UP=MEC$UP)

  }

