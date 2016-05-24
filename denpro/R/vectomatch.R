vectomatch<-function(vec1,vec2)
{

d<-dim(vec1)[2]
prenum<-dim(vec1)[1]
curnum<-dim(vec2)[1]
parento<-matrix(0,curnum,1)

smallernum<-min(prenum,curnum)
greaternum<-max(prenum,curnum)

dista<-matrix(NA,smallernum,greaternum)
for (ap in 1:smallernum){
    for (be in 1:greaternum){
           if (d==1){
               if (prenum<=curnum){
               precenter<-vec1[ap]
               curcenter<-vec2[be]
               }
               else{
               precenter<-vec2[ap]
               curcenter<-vec1[be]
               }
           }
           else{
               if (prenum<=curnum){
               precenter<-vec1[ap,]
               curcenter<-vec2[be,]
               }
               else{
               precenter<-vec2[ap,]
               curcenter<-vec1[be,]
               }
           }
           dista[ap,be]<-etais(curcenter,precenter)
    }
}

match<-matrix(0,smallernum,1)  #for each mode the best match
findtie<-TRUE

# find the best match for all and check whether there are ties
match<-matrix(0,smallernum,1)
for (bm in 1:smallernum){
      minimi<-min(dista[bm,],na.rm=TRUE)
      match[bm]<-which(minimi==dista[bm,])[1]
}
findtie<-FALSE
bm<-1
while ((bm<=smallernum) && (findtie==FALSE)){
         koe<-match[bm]
         bm2<-bm+1
         while (bm2<=smallernum){
            if (koe==match[bm2]){
                  findtie<-TRUE
            }
            bm2<-bm2+1
         }
         bm<-bm+1
}
    
onkayty<-FALSE
tiematch<-matrix(0,smallernum,1)

while (findtie){

      onkayty<-TRUE
      
      # find the best match for all
      bestmatch<-matrix(0,smallernum,1)
      for (bm in 1:smallernum){
          allna<-TRUE
          am<-1
          while ((am<=greaternum) && (allna)){
             if (!is.na(dista[bm,am])) allna<-FALSE
             am<-am+1
          }
          if (!(allna)){
             minimi<-min(dista[bm,],na.rm=TRUE)
             bestmatch[bm]<-which(minimi==dista[bm,])[1]
          }
          else bestmatch[bm]<-tiematch[bm]
      }

      # find the first tie
      findtie<-FALSE

      tieset<-matrix(0,smallernum,1)
      bm<-1
      while ((bm<=smallernum) && (findtie==FALSE)){
         koe<-bestmatch[bm]
         bm2<-bm+1
         while (bm2<=smallernum){
            if (koe==bestmatch[bm2]){
                  findtie<-TRUE
                  tieset[bm]<-1
                  tieset[bm2]<-1
            }
            bm2<-bm2+1
         }
         bm<-bm+1
      }

      # solve the first tie
      if (findtie==TRUE){
         numofties<-sum(tieset)
         kavelija<-0
         tiepointer<-matrix(0,numofties,1) 
         # find the second best
         secondbest<-matrix(0,smallernum,1)
         for (bm in 1:smallernum){
            if (tieset[bm]==1){
               redudista<-dista[bm,]
               redudista[bestmatch[bm]]<-NA
               minimi<-min(redudista,na.rm=TRUE)
               secondbest[bm]<-which(minimi==redudista)[1]

               kavelija<-kavelija+1
               tiepointer[kavelija]<-bm
            }
         }
         # try different combinations       
         # try all subsets of size 2 from the set of ties
         numofsubsets<-choose(numofties,2)
            #gamma(numofties+1)/gamma(numofties-2+1)
         valuelist<-matrix(0,numofsubsets,1)
         vinnerlist<-matrix(0,numofsubsets,1)
         matchlist<-matrix(0,numofsubsets,1)
         runneri<-1
         eka<-1
         while (eka<=numofties){
            ekapo<-tiepointer[eka]
            toka<-eka+1
            while (toka<=numofties){
               tokapo<-tiepointer[toka]
               # try combinations for this subset (there are 2)
               # 1st combination
               fvinner<-ekapo
               fvinnermatch<-bestmatch[fvinner]
               floser<-tokapo
               flosermatch<-secondbest[floser]
               fvalue<-dista[fvinner,fvinnermatch]+dista[floser,flosermatch]
                # 2nd combination
               svinner<-tokapo
               svinnermatch<-bestmatch[svinner]
               sloser<-ekapo
               slosermatch<-secondbest[sloser]
               svalue<-dista[svinner,svinnermatch]+dista[sloser,slosermatch]
               # tournament
               if (fvalue<svalue){
                   valuelist[runneri]<-fvalue
                   vinnerlist[runneri]<-fvinner
                   matchlist[runneri]<-fvinnermatch
               }
               else{ 
                   valuelist[runneri]<-svalue
                   vinnerlist[runneri]<-svinner
                   matchlist[runneri]<-svinnermatch
               }
               runneri<-runneri+1 
               # 
               toka<-toka+1
            }
            eka<-eka+1
         }
         minimi<-min(valuelist,na.rm=TRUE)
         bestsub<-which(minimi==valuelist)[1]
         vinnerson<-vinnerlist[bestsub]
         matcherson<-matchlist[bestsub]

         tiematch[vinnerson]<-matcherson
         dista[vinnerson,]<-NA
         dista[,matcherson]<-NA

      }

}  #while (findtie)

if (onkayty){  #there was one tie
          
          for (sepo in 1:smallernum){
               if (tiematch[sepo]!=0) match[sepo]<-tiematch[sepo]
               else match[sepo]<-bestmatch[sepo]
          }
}

newnode<-matrix(0,curnum,1)
if (prenum>curnum) parento<-match
else{
    for (i in 1:prenum){   #kaannetaan linkit
        linko<-match[i]
        parento[linko]<-i
    }    
    for (j in 1:curnum){
        if (parento[j]==0){   #jos ei linkkia, haetaan lahin vanhemmaksi
             newnode[j]<-1    #we label the rest-labels
             distvec<-dista[,j]  #sarake antaa etaisyyden
             minimi<-min(distvec,na.rm=TRUE)
             parento[j]<-which(minimi==distvec)[1]
        }
    }
}

return(list(parent=parento,newnode=newnode))
}



















