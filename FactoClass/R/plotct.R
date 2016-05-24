#---------------------------------------------------------
# Función para graficar perfiles fila y columna de tablas
# de contingencia
# versión inicial de Camilo Torres
# Modificaciones julio 15 de 2010. CE Pardo
# modificación julio 30 de 2015: agregar salida de tablas 
# de perfiles y de contingencia con marginales
#---------------------------------------------------------
plotct <- function(x,profiles="both",legend.text=TRUE,tables=FALSE,nd=1,... )
{
 x <- as.matrix( x )
 total <- sum( x )
 f.marginal <- colSums( x ) / total
 c.marginal <- rowSums( x ) / total
 f.perfil <- rbind( prop.table( x, 1 ), marg=f.marginal )
 c.perfil <- cbind( prop.table( x, 2 ), marg=c.marginal )
 # graficas con leyenda
 if (legend.text==TRUE)
   {
      if (profiles=="both" | profiles=="row")
       barplot( t(f.perfil), legend.text=legend.text, beside=FALSE, horiz=TRUE,
               las=1, xlim=c(0,1.5),
               xaxp=c(0,1,5),
               args.legend = list( x = "right"), ... )
     if (profiles=="both") dev.new()
     if (profiles=="both" | profiles=="col")
       barplot( c.perfil, legend.text=legend.text, beside=FALSE, las=2,
               xlim=c(0,ncol(x)+6.5),
               args.legend = list( x = "right" ), ... )
  }
  # graficas sin leyenda
 if (legend.text==FALSE)
   {
      if (profiles=="both" | profiles=="row")
       barplot( t(f.perfil),beside=FALSE, horiz=TRUE,las=1, ... )
     if (profiles=="both") dev.new()
     if (profiles=="both" | profiles=="col")
       barplot( c.perfil,beside=FALSE, las=2,... )
 }
 # adicionado por CEPT jul 30 2015
 if (tables) {
    tcm <- cbind(x,marR=rowSums(x))
    tcm <- rbind(tcm,marC=colSums(tcm))
    tab<-NULL
    tab$ctm=tcm; tab$perR<-round(f.perfil*100,nd); tab$perC<-round(c.perfil*100,nd)
    return(tab)
 }    
}

