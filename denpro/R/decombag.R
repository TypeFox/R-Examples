decombag<-function(numofall,levseq,
left,right,val,vec,infopointer,parent,
nodenumOfTree,terminalnum,
value,low,upp,nodefinder,
d)
{
#Makes density tree
#returns list(parent,level,component)
#  -component is pointer to AtomlistAtom, AtomlistNext, ei tarvita

#apu2<-0

levnum<-length(levseq)

AtomlistNext<-matrix(0,numofall,1)
AtomlistAtom<-matrix(0,numofall,1) #point to value,..: values in 1,...,atomnum

componum<-numofall
parentLST<-matrix(0,componum,1)
level<-matrix(0,componum,1)
component<-matrix(0,componum,1)    

pinoComponent<-matrix(0,componum,1)   #pointer to component, level,...
pinoTaso<-matrix(0,componum,1)        #ordinal of level (pointer to levseq)

# Initilize the lists

AtomlistAtom[1:terminalnum]<-seq(1,terminalnum)
AtomlistNext[1:(terminalnum-1)]<-seq(2,terminalnum)
AtomlistNext[terminalnum]<-0
listEnd<-terminalnum
beg<-1

# Let us divide the lowest level set to disconnected parts

begi<-1
tobehandled<-matrix(0,nodenumOfTree,1)
atto<-AtomlistAtom[begi]
while (begi>0){
   if (value[atto]>0){
      node<-nodefinder[atto]
      tobehandled[node]<-1
   }    
   begi<-AtomlistNext[begi]
   atto<-AtomlistAtom[begi]
}
curterminalnum<-terminalnum

dld<-declevgen(tobehandled,curterminalnum,
left,right,val,vec,infopointer,parent,
low,upp,
d)
 
totbegSepary<-dld$totbegSepary
begsSepaNext<-dld$begsSepaNext
begsSepaBegs<-dld$begsSepaBegs
atomsSepaNext<-dld$atomsSepaNext
atomsSepaAtom<-dld$atomsSepaAtom

lc<-listchange(AtomlistAtom,AtomlistNext,totbegSepary,
begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom,
curterminalnum,beg)
begs<-lc$begs
AtomlistNext<-lc$AtomlistNext
AtomlistAtom<-lc$AtomlistAtom

koko<-length(begs)
# Talletetaan osat 
component[1:koko]<-begs
level[1:koko]<-levseq[1]       #arvo toistuu 
efek<-koko                     #kirjataan uusien osien lkm  ????? jos vain yksi
# Laitetaan kaikki osat pinoon
pinoComponent[1:koko]<-seq(1,koko)      #1,2,...,koko
pinoTaso[1:koko]<-1     #kaikki osat kuuluvat alimpaan tasojoukkoon
pinind<-koko            #indeksi pinoon

if (levnum>1){  while (pinind>=1){
  # Take from stack
  ind<-pinoComponent[pinind]      #indeksi tasoon
  levind<-pinoTaso[pinind]        #ko tason korkeus
  pinind<-pinind-1                #otettu pinosta
  partlevsetbeg<-component[ind]  
  higlev<-levseq[levind+1]
  
  # Make intersection with the curr. component and higher lev.set
  PrevlistEnd<-listEnd
  addnum<-0      #num of atoms in the intersection
  removenum<-0   #num of atoms which have to be removed to get intersec.
  runner<-partlevsetbeg
  origiListEnd<-listEnd
  while ((runner>0) && (runner<=origiListEnd)){
      atom<-AtomlistAtom[runner]
      arvo<-value[atom]
      if (arvo>=higlev){
          listEnd<-listEnd+1    
          AtomlistAtom[listEnd]<-atom
          AtomlistNext[listEnd]<-listEnd+1 
          addnum<-addnum+1     
      }
      else{           
          removenum<-removenum+1
      }                                
      runner<-AtomlistNext[runner] 
  }
  AtomlistNext[listEnd]<-0      # we have to correct the end to terminate
  if (addnum>0){
      AtomlistNext[PrevlistEnd]<-0
      beghigher<-PrevlistEnd+1 
  }
  if (removenum==0){  #jos leikkaus ei muuta, niin tasoj sailyy samana  
      level[ind]<-levseq[levind+1]  #remove lower part
          #component and parentLST stay same, it is enough to change level
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pinoComponent[pinind+1]<-ind     
        pinoTaso[pinind+1]<-levind+1 #tasojouk taso on levind+1 
        pinind<-pinind+1       
      }
  }
  else if (addnum>0){     #leikkaus ei tyhja
      beg<-beghigher

      begi<-beghigher
      tobehandled<-matrix(0,nodenumOfTree,1)
      while (begi>0){
          atto<-AtomlistAtom[begi]
          node<-nodefinder[atto]
          tobehandled[node]<-1
          begi<-AtomlistNext[begi]
      }
      curterminalnum<-addnum

      dld<-declevgen(tobehandled,curterminalnum,
      left,right,val,vec,infopointer,parent,
      low,upp,
      d)
 
      totbegSepary<-dld$totbegSepary
      begsSepaNext<-dld$begsSepaNext
      begsSepaBegs<-dld$begsSepaBegs
      atomsSepaNext<-dld$atomsSepaNext
      atomsSepaAtom<-dld$atomsSepaAtom

#      apu2<-apu2+1
#if (apu2==1) an<-atomsSepaNext     

      lc<-listchange(AtomlistAtom,AtomlistNext,totbegSepary,
      begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom,
      curterminalnum,beg)
      begs<-lc$begs
      AtomlistNext<-lc$AtomlistNext
      AtomlistAtom<-lc$AtomlistAtom

      koko<-length(begs)    #jos vain yksi ?????????
      #
      # paivitetaan kumu tulokseen
      level[(efek+1):(efek+koko)]<-levseq[levind+1]   #arvo toistuu  
      parentLST[(efek+1):(efek+koko)]<-ind
      component[(efek+1):(efek+koko)]<-begs
      efek<-efek+koko
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pinoComponent[(pinind+1):(pinind+koko)]<-seq(efek-koko+1,efek)
        pinoTaso[(pinind+1):(pinind+koko)]<-levind+1  #tasjouk tas on levind+1 
        pinind<-pinind+koko       
      }
  }
}} 
level<-t(level[1:efek])
parentLST<-t(parentLST[1:efek])
component<-t(component[1:efek])
#
return(list(level=level,parent=parentLST,
component=component,AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext))
}









