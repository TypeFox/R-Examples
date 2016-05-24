a.getcol.ade <-
function(n, alpha=1, type='normal'){

#########################################
# make single values
if(length(n)>1){
y<-as.list(n)
for(i in 1:length(n)){
y[[i]] <- a.getcol.ade(n[i], alpha=alpha, type=type)
}
return(y)
}
#########################################

if(length(n)==1){
if(type=='normal' | type=='n'){
if(n==1) color <- 1
if(n>1 & n<10) color <- c('blue','red','green','cyan','magenta','yellow','thistle3', 'orange', 'brown')[1:n]
if(n>=10)       color <- rainbow(n+round(n/10))[1:n]
}

if(type=='pastel' | type=='p'){
if(n==1) color <- 'honeydew4'
if(n>1 & n<10) color <- c('steelblue2','brown1','palegreen3','turquoise2','orchid2','khaki2','thistle3', 'sienna1', 'salmon4')[1:n]
if(n>=10)       color <- terrain.colors(n+1)[1:n]
}

color<-a.alpha.ade(color, alpha)

return(color)
}
}
