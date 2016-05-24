victory_int<-function(x=x,y=y,col=c("black","red","green"),alpha=0.2,what=which.min){
   quimin=apply(y,1,what)
   lim <- par("usr")
   minloc=quimin[1]
   xdeb=x[1]
   for(i in 2:length(x)){
      if(quimin[i]!=minloc & i!=length(x)){
         couleur_loc=c(t(col2rgb(col[quimin[i-1]]))/255)
         couleur_loc=rgb(red=couleur_loc[1],green=couleur_loc[2],blue=couleur_loc[3],alpha)
         rect(xdeb, lim[3]-1, (x[i-1]+x[i])/2, lim[4]+1, border = couleur_loc, col = couleur_loc)
         xdeb=(x[i-1]+x[i])/2
         minloc=quimin[i]
      }else if(quimin[i]!=minloc & i==length(x)){#changement et on est au dernier, on ferme le precedent et on ajoute le rectange final
         couleur_loc=c(t(col2rgb(col[quimin[i-1]]))/255)
         couleur_loc=rgb(red=couleur_loc[1],green=couleur_loc[2],blue=couleur_loc[3],alpha)
         rect(xdeb, lim[3]-1, (x[i-1]+x[i])/2, lim[4]+1, border = couleur_loc, col = couleur_loc)
         xdeb=(x[i-1]+x[i])/2
         minloc=quimin[i]
         couleur_loc=c(t(col2rgb(col[quimin[i]]))/255)
         couleur_loc=rgb(red=couleur_loc[1],green=couleur_loc[2],blue=couleur_loc[3],alpha)
         rect(xdeb, lim[3]-1, x[i], lim[4]+1, border = couleur_loc, col = couleur_loc)
      }else if(i==length(x)){#fini donc on trace
         couleur_loc=c(t(col2rgb(col[quimin[i]]))/255)
         couleur_loc=rgb(red=couleur_loc[1],green=couleur_loc[2],blue=couleur_loc[3],alpha)
         rect(xdeb, lim[3]-1, x[i], lim[4]+1, border = couleur_loc, col = couleur_loc)
      }
   }
}