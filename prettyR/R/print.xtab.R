print.xtab<-function (x,col.width=7,or=TRUE,chisq=FALSE,phi=FALSE,
 rowname.width=NA,html=FALSE,bgcol="lightgray",...) {

 tdim <- dim(x$counts)
 ncols<-tdim[2]+2
 row.labels<-format(names(x$row.margin))
 if(is.na(rowname.width)) rowname.width<-max(nzchar(row.labels))
 if(any(nzchar(row.labels) > rowname.width))
  truncString(row.labels,rowname.width)
 rowname.space<-paste(rep(" ",rowname.width),sep="",collapse="")
 if(html)
  cat("<table border=0 relief=\"flat\" style=\"background-color:",bgcol,
   "\">\n<tr>\n<td colspan=",ncols,">\n",sep="")
 cat("Crosstabulation of",x$varnames[1],"by",x$varnames[2],"\n")
 if(html) cat("<tr><td><td colspan=",ncols-2," align=center>",sep="")
 else cat("\t")
 cat(x$varnames[2])
 col.labels<-names(x$col.margin)
 if(max(nzchar(col.labels)) > col.width)
  col.labels<-truncString(col.labels,maxlen=col.width,justify="left")
 if(html) cat("<tr><td>",x$varnames[1],sep="")
 else cat("\n",format(x$varnames[1],width=rowname.width),sep="")
 for(i in 1:length(col.labels)) {
  if(html) cat("<td>")
  cat(format(col.labels[i],width=col.width,justify="right"))
 }
 cat("\n")
 gt <- sum(x$counts)
 for(i in 1:tdim[1]) {
  if(html) cat("<tr><td>")
  cat(row.labels[i])
  if(html) {
   for(j in 1:tdim[2]) {
    if(html) cat("<td>")
    cat(x$counts[i,j])
    if(html) cat("<br>")
    cat(round(100*x$counts[i,j]/x$row.margin[i],2))
    if(html) cat("<br>")
    cat(round(100*x$counts[i,j]/x$col.margin[j],2))
   }
   if(html) cat("<td>")
   cat(x$row.margin[i])
   if(html) cat("<br>-<br>")
   cat(round(100*x$row.margin[i]/gt,2),"\n")
  }
  else {
   cat(format(c(x$counts[i,],x$row.margin[i]),width=col.width),
    "\n",sep="")
   cat(rowname.space,format(c(round(100*x$counts[i,]/x$row.margin[i],2),"-"),
    width=col.width,justify="right"),"\n",sep="")
   cat(rowname.space,format(round(c(100*x$counts[i,]/x$col.margin,
    100*x$row.margin[i]/gt),2),width=col.width),"\n\n",sep="")
  }
 }
 if(html) {
  cat("<tr><td>\n")
  for(i in 1:tdim[2]) {
   if(html) cat("<td>")
   cat(x$col.margin[i])
   if(html) cat("<br>")
   cat(round(100*x$col.margin[i]/gt,2))
  }
  cat("<td>",gt,"<br>100\n")
 }
 else {
  cat(rowname.space,format(c(x$col.margin,gt),width=col.width),"\n",sep="")
  cat(rowname.space,
   format(c(round(100*x$col.margin/gt,2),100),width=col.width),"\n",sep="")
 }
 if(chisq) {
  x2<-chisq.test(x$counts, ...)
  if(html) cat("<tr><td colspan=",ncols,">",sep="")
  cat("X2[",x2$parameter,"]=",round(x2$statistic,3),", p=",
   x2$p.value,"\n\n",sep = "")
 }
 if(tdim[1]==2 && tdim[2]==2) {
  logical.names<-function(x,names1=c("FALSE","0","N","NO"),
   names2=c("TRUE","1","Y","YES")) {
   for (nameindex in 1:length(names1))
    if(toupper(x[1]) == names1[nameindex] &&
     toupper(x[2]) == names2[nameindex]) 
     return(TRUE)
    return(FALSE)
  }
  if(or) {
   if(html) cat("<tr><td colspan=",ncols,">",sep="")
   cat("odds ratio =",
    round(x$counts[1,1]*x$counts[2,2]/(x$counts[1, 2]*x$counts[2,1]),2),"\n")
   if(logical.names(names(x$col.margin))) {
    if(html) cat("<tr><td colspan=",ncols,">",sep="")
    cat("relative risk (",x$varnames[1],"-",names(x$row.margin)[2], 
     ") = ",round((x$counts[2,2]/x$row.margin[2])/
     (x$counts[1,2]/x$row.margin[1]),2),"\n",sep="")
  
   }
   if(phi) {
    num<-x$counts[1,1]*x$counts[2,2]-x$counts[1,2]*x$counts[2,1]
    denom<-sqrt(as.numeric(x$row.margin[1])*
     as.numeric(x$row.margin[2])*as.numeric(x$col.margin[1])*
     as.numeric(x$col.margin[2]))
    if(html) cat("<tr><td colspan=",ncols,">",sep="")
    cat("phi =",round(num/denom,3),"\n")
   }
  }
  else if(phi) cat("phi coefficient only valid for 2x2 table\n")
 }
 if(html) cat("</table>\n")
}
