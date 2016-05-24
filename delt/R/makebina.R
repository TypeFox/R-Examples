makebina<-function(et){

#source("~/delt/R/makebina.R")
#list(frame)
#frame:
#var  splitting variable  or "<leaf>"
#n    number of cases
#dev  deviance of the node
#yval fitted value
#split/splits two-column matrix of the label for the left and right splits
#splits.cutleft  splits.cutright

len<-length(et$split)

var<-matrix("",len,1)
split<-matrix("",len,2)
n<-matrix(0,len,1)
dev<-matrix(0,len,1)
yval<-matrix(0,len,1)

leve<-matrix(0,len,1)    #level for each node
runnum<-matrix(0,len,1)  #ordinal number in the (imaginary) full bin tree 
runredu<-matrix(0,len,1)

leve[1]<-1
runnum[1]<-1

pino<-matrix(0,len,1)
pinoind<-1
pino[1]<-1
laskuri<-0
while (pinoind>0){
   cur<-pino[pinoind]
   pinoind<-pinoind-1
      
   laskuri<-laskuri+1
   n[laskuri]<-et$nelem[cur]                #+1
   dev[laskuri]<-et$ssr[cur]                #abs(et$ssr)
   yval[laskuri]<-et$mean[cur]
   if (et$left[cur]==0){ 
          var[laskuri]<-"<leaf>"
   } 
   else{ 
          var[laskuri]<-paste("x",as.character(et$direc[cur]))
   }
   #var[i]<-as.character(et$direc[i])
   if (et$left[cur]!=0){        #(!is.na(et$split[cur])){
     #split[laskuri,1]<-as.character(et$split[cur])
     #split[laskuri,2]<-as.character(et$split[cur])
     split[laskuri,1]<-format(et$split[cur],digits=2,nsmall=1)
     split[laskuri,2]<-format(et$split[cur],digits=2,nsmall=1)
   }
   runredu[laskuri]<-runnum[cur]
   
   while (et$left[cur]>0){
       #
       # laita oikea pinoon, laske "leve" ja "runnum"
       #
       oikea<-et$right[cur]
       #
       pinoind<-pinoind+1
       pino[pinoind]<-oikea
       #
       leve[oikea]<-leve[cur]+1
       runnum[oikea]<-childcode(leve[cur],runnum[cur])$right
       #
       # count "leve" ja "runnum", go to et$left, update
       #
       vasen<-et$left[cur]
       #
       leve[vasen]<-leve[cur]+1
       runnum[vasen]<-childcode(leve[cur],runnum[cur])$left
       #
       cur<-vasen
       #
       laskuri<-laskuri+1
       n[laskuri]<-et$nelem[cur]                #+1
       dev[laskuri]<-et$ssr[cur]                #abs(et$ssr)
       yval[laskuri]<-et$mean[cur]
       if (et$left[cur]==0){ 
           var[laskuri]<-"<leaf>" 
       } 
       else{ 
           var[laskuri]<-paste("x",as.character(et$direc[cur])) 
       }
       #var[i]<-as.character(et$direc[i])
       if (et$left[cur]!=0){    #if not child
         #split[laskuri,1]<-as.character(et$split[cur])
         #split[laskuri,2]<-as.character(et$split[cur])
         split[laskuri,1]<-format(et$split[cur],digits=2,nsmall=1)
         split[laskuri,2]<-format(et$split[cur],digits=2,nsmall=1)
       }
       runredu[laskuri]<-runnum[cur]
       #
   }
}
var<-var[1:laskuri]
split<-split[1:laskuri,]
n<-n[1:laskuri]
dev<-dev[1:laskuri]
yval<-yval[1:laskuri]
#
runredu<-runredu[1:laskuri]
#
#
frame<-data.frame(var=var,n=n,dev=dev,yval=yval) 
row.names(frame)<-runredu
frame$splits<-array(split,c(laskuri,2),
                    list(character(0), c("cutleft", "cutright")))
#
where<-matrix(1,len,1)
terms<-""
call<-""
y<-FALSE
weigths<-FALSE
#
bintree<-list(frame=frame,where=where,terms=terms,call=call,y=y,weigths=weigths)
attr(bintree,"class")<-"tree" 
return(bintree)
}





