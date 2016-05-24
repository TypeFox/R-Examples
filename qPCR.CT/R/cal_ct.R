cal_ct <-
function(con.con, tr.con, con.tr, tr.tr)
{
  #get length of data
  lcc <- length(con.con)
  ltc <- length(tr.con)
  lct <- length(con.tr)
  ltt <- length(tr.tr)
  #check length
  if(lcc == 0 || ltc == 0 || lct == 0 || ltt ==0)
  {
    cat("Invalid dataset! Length of dataset is 0. Please check dataset!")
  } else if ((lcc != ltc) || (lcc != lct) || (lcc!=ltt) || (ltc != lct) || (lct != ltt)) {
    cat("Invalid dataset! The length of dataset is not equal.Please check dataset!")
  } else {
    mean.cc <- mean(con.con)
    mean.tc <- mean(tr.con)
    mean.ct <- mean(con.tr)
    mean.tt <- mean(tr.tr)
    dcon.tr <- con.tr - mean.cc
    dtr.tr  <- tr.tr  - mean.tc
    mean.dcon.tr <- mean(dcon.tr)
    mean.dtr.tr <- mean(dtr.tr)
    ddcon.tr <-  mean.dcon.tr - dcon.tr
    ddtr.tr  <-  mean.dcon.tr - dtr.tr
    two.ddcon.tr <- 2^ddcon.tr
    two.ddtr.tr <- 2^ddtr.tr
    mean.two.ddcon.tr <- mean(two.ddcon.tr)
    mean.two.ddtr.tr <- mean(two.ddtr.tr)
    final.con <- two.ddcon.tr / mean.two.ddcon.tr
    final.tr  <- two.ddtr.tr /  mean.two.ddcon.tr
    final <- data.frame(final.con, final.tr)
    return(final)
  }
}

