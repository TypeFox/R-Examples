listchange<-function(AtomlistAtom,AtomlistNext,totbegSepary,
begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom,
terminalnum,beg){
#
#create begs: beginnings of lists of atoms
#beg is index to AtomlistAtom/Next
#totbegsepary is index to begsSepaBegs/Next
#
begs<-matrix(0,terminalnum,1)
#
runnerBegs<-totbegSepary  #total beginning of list is at the root of the tree
runnerOrigi<-beg
runnerOrigiprev<-beg
sepalkm<-0
while (runnerBegs>0){
   sepalkm<-sepalkm+1
   runnerAtoms<-begsSepaBegs[runnerBegs]
   begs[sepalkm]<-runnerOrigi
   #  first step (in order to get also runnerOrigiprev to play)
   AtomlistAtom[runnerOrigi]<-atomsSepaAtom[runnerAtoms]
   runnerOrigiprev<-runnerOrigi
   runnerOrigi<-AtomlistNext[runnerOrigi]
   runnerAtoms<-atomsSepaNext[runnerAtoms]
   while (runnerAtoms>0){
       AtomlistAtom[runnerOrigi]<-atomsSepaAtom[runnerAtoms]
       runnerOrigiprev<-runnerOrigi
       runnerOrigi<-AtomlistNext[runnerOrigi]
       runnerAtoms<-atomsSepaNext[runnerAtoms]
   }
   AtomlistNext[runnerOrigiprev]<-0      #mark the end of the list
   runnerBegs<-begsSepaNext[runnerBegs]
}
#
begs<-begs[1:sepalkm]
#
return(list(begs=begs,AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext))
}
