coloallo<-function(mt,paletti=NULL)
{
# fast allocation of colors (matching of modes)
# mt is mode tree
# paletti gives a list of colors

if (is.null(paletti))
paletti<-c("red","blue","green","turquoise","orange","navy",
"darkgreen","orchid",colors()[50:100])

d<-dim(mt$xcoor)[2]

snum<-0
for (i in 1:length(mt$mlabel)){
  if (mt$mlabel[i]==1) snum<-snum+1
}

xcoor<-mt$xcoor
ycoor<-mt$ycoor
mlabel<-mt$mlabel
lenni<-length(ycoor)

colot<-matrix("",lenni,1)

# find the locations for the information for each h

low<-matrix(0,snum,1)
upp<-matrix(0,snum,1)
low[1]<-1
glob<-2
while ((glob<=lenni) && (mlabel[glob]!=1)){
       glob<-glob+1
}
upp[1]<-glob-1
# now glob is at the start of new block
i<-2
while (i<=snum){
   low[i]<-glob
   glob<-glob+1
   while ((glob<=lenni) && (mlabel[glob]!=1)){
       glob<-glob+1
   }
   upp[i]<-glob-1
   i<-i+1
}

# first we allocate colors for the largest h

run<-1  #low[1]
while (run<=upp[1]){
   colot[run]<-paletti[run]
   run<-run+1
}

firstnewcolo<-run

i<-2

while (i<=snum){
   prenum<-upp[i-1]-low[i-1]+1
   curnum<-upp[i]-low[i]+1

   smallernum<-min(prenum,curnum)
   greaternum<-max(prenum,curnum)

   if (prenum==smallernum){   
        bases<-i
        compa<-i-1
   }
   else{
        bases<-i-1
        compa<-i
   }

      dista<-matrix(NA,smallernum,greaternum)
      for (ap in low[bases]:upp[bases]){
         for (be in low[compa]:upp[compa]){
           if (d==1){
               curcenter<-xcoor[ap]
               precenter<-xcoor[be]
           }
           else{
               curcenter<-xcoor[ap,]
               precenter<-xcoor[be,]
           }
           dista[be-low[compa]+1,ap-low[bases]+1]<-etais(curcenter,precenter)
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

      while (findtie){

      onkayty<-TRUE
      tiematch<-matrix(0,smallernum,1)
      
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
          else bestmatch[bm]<-match[bm]
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

      # finally allocate colors
      run<-1
      while (run<=smallernum){
          
          if (prenum==smallernum){
             xind<-run
             yind<-match[xind]
          }
          else{
             yind<-run
             xind<-match[yind]
          }

          colot[low[i]+yind-1]<-colot[low[i-1]+xind-1]    
          run<-run+1
      }
                    
      if (prenum<greaternum){

        run<-low[bases]
        while (run<=upp[bases]){
            if (colot[run]==""){
               colot[run]<-paletti[firstnewcolo]
               firstnewcolo<-firstnewcolo+1   
            }
            run<-run+1
        }

     }

     i<-i+1
}

return(colot)
}


















