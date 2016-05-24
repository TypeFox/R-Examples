## A QAbstractItemDelegate extension that formats numbers to text using R

qrTextFormattingDelegate <- function(parent = NULL)
{
  .Call("qt_qrTextFormattingDelegate", parent, PACKAGE="qtbase")
}
