frmtMntse <-
function (inNum, mntse) {
  # function to control the significands of a number
return(as.numeric(format(round(inNum, mntse), nsmall=mntse)))
}
