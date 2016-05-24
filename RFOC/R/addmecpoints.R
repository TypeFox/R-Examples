`addmecpoints` <-
function(MEC, pch=5)
  {

    focpoint(MEC$F$az, MEC$F$dip, pch=5, lab="F", UP=MEC$UP)
    focpoint(MEC$G$az, MEC$G$dip, pch=5, lab="G", UP=MEC$UP)

     focpoint(MEC$U$az, MEC$U$dip, pch=5, lab="U", UP=MEC$UP)
    focpoint(MEC$V$az, MEC$V$dip, pch=5, lab="V", UP=MEC$UP)

  
   focpoint(MEC$P$az, MEC$P$dip, pch=5, lab="P", UP=MEC$UP)
    focpoint(MEC$T$az, MEC$T$dip, pch=5, lab="T", UP=MEC$UP)

  }

