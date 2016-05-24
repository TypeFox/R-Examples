#define function row.any and col.any
row.any<-function (x)
{
  apply(x, 1, any)
}
row.all<-function (x)
{
  apply(x, 1, all)
}

col.any<-function (x)
{
  apply(x, 2, any)
}
col.all<-function (x)
{
  apply(x, 2, all)
}

remove.na.rows <- function(x)
{
  x[!row.any(is.na(x)),]
}
