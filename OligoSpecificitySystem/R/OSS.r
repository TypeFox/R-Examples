##########################################
### Oligo Specificity System
##########################################
#This function is based on ARB software output. It allows to
#1/ Calculate the number on sequences matched by a degenerated primer (with on ly one base pair difference).
#2/ Calculate the real efficiency of PCR system based on the number of sequences amplified by both primers.
#####################################################
### How it woks 
#1/ perform the first primer match (first primer of the PCR system or first sequence of the degenerated primer) and print it in ASCII .TXT  using ARB
#2/ perform the second primer match (second primer of the PCR system or second sequence of the degenerated primer) and print it in ASCII .TXT  using ARB
#3/ Copy and paste the function in the prompt of R
#4/ start the function using OSS() in the promt
##########################################################################

"OSS"<-function()
{ 
  tclRequire("BWidget") 
  tt <- tktoplevel()
  tkwm.title(tt,"Oligo Specificity System (OSS)")
  tkgrid(tklabel(tt,text="                                                               ")) 
  
  closeAllConnections()
  ONDB1name<<-"Not loaded"
  ONDB2name<<-"Not loaded"
  ONDB3name<<-"Not loaded. Loading of this data base is optionnal"
  ONDB4name<<-"Not loaded. Optionnal"
  .write<-function(){
          ####OligoNucleotide data base 1
          tkconfigure(ONDB1txt, state="normal")
          tkdelete(ONDB1txt,"0.0","100000.0")
          tkinsert(ONDB1txt, "end", ONDB1name)
          tkconfigure(ONDB1txt, state="disabled")
          ####OligoNucleotide data base 2
          tkconfigure(ONDB2txt, state="normal")
          tkdelete(ONDB2txt,"0.0","100000.0")
          tkinsert(ONDB2txt, "end", ONDB2name)
          tkconfigure(ONDB2txt, state="disabled")
          ####OligoNucleotide data base 3
          tkconfigure(ONDB3txt, state="normal")
          tkdelete(ONDB3txt,"0.0","100000.0")
          tkinsert(ONDB3txt, "end", ONDB3name)
          tkconfigure(ONDB3txt, state="disabled")
          ####OligoNucleotide data base 4
          tkconfigure(ONDB4txt, state="normal")
          tkdelete(ONDB4txt,"0.0","100000.0")
          tkinsert(ONDB4txt, "end", ONDB4name)
          tkconfigure(ONDB3txt, state="disabled")}
  run <- function(){
          code <- tclvalue(tkget(txt,"0.0","end"))
          e <- try(parse(text=code))
          if (inherits(e, "try-error")){
          tkmessageBox(message="Syntax error",icon="error")
          return()}
          .write()
          print(eval(e))
          .write()}
  
  example1<-function(){optxt(n=1)
                      OnOk()
                      print("Please wait while loading")
                      ONDB1name<<-fil
                      txt1<-scan(file=fil,what="character",fill=TRUE,quiet=TRUE)
                      a<-which(txt1==separator)
                      b<-vector(length=length(a))
                      for (i in 1:c(length(a)-1)){
                      b[i]<-paste(txt1[c(a[i]+1):c(a[i+1]-1)],collapse="")}
                      b[length(a)]<-paste(txt1[c(a[length(a)]+1):length(txt1)],collapse="")
                      txt1<<-b
                      print(paste("Your oligonucleotide database has",length(b),"sequences"))
                      print("DataBase 1 successfully imported")}
  example2<-function(){optxt(n=2)
                      OnOk()
                      print("Please wait while loading")
                      ONDB2name<<-fil
                      txt2<-scan(file=fil,what="character",fill=TRUE,quiet=TRUE)
                      a<-which(txt2==separator)
                      b<-vector(length=length(a))
                      for (i in 1:c(length(a)-1)){
                      b[i]<-paste(txt2[c(a[i]+1):c(a[i+1]-1)],collapse="")}
                      b[length(a)]<-paste(txt2[c(a[length(a)]+1):length(txt2)],collapse="")
                      txt2<<-b
                      print(paste("Your oligonucleotide database has",length(b),"sequences"))
                      print("DataBase 2 successfully imported")}
  example3<-function(){optxt(n=3)
                      OnOk()
                      print("Please wait while loading")
                      ONDB3name<<-fil
                      txt3<-scan(file=fil,what="character",fill=TRUE,quiet=TRUE)
                      a<-which(txt3==separator)
                      b<-vector(length=length(a))
                      for (i in 1:c(length(a)-1)){
                      b[i]<-paste(txt3[c(a[i]+1):c(a[i+1]-1)],collapse="")}
                      b[length(a)]<-paste(txt3[c(a[length(a)]+1):length(txt3)],collapse="")
                      txt3<<-b
                      print(paste("Your oligonucleotide database has",length(b),"sequences"))
                      print("DataBase 3 successfully imported")}
  export1<-function(){exportglobal(n=1)}
  export2<-function(){exportglobal(n=2)}
  export3<-function(){exportglobal(n=3)}
  
  reset<-function(){txt1<-0;txt1<<-txt1;txt2<-0;txt2<<-txt2;txt3<-0;txt3<<-txt3;txt4<-0;txt4<<-txt4}
  reset1<-function(){txt1<-0;txt1<<-txt1;ONDB1name<<-"Not loaded"}
  reset2<-function(){txt2<-0;txt2<<-txt2;ONDB2name<<-"Not loaded"}
  reset3<-function(){txt3<-0;txt3<<-txt3;ONDB3name<<-"Not loaded. Loading of this data base is optionnal"}
  reset4<-function(){txt4<-0;txt4<<-txt4;ONDB4name<<-"Not loaded. Optionnal"}
          reset()
  OnOk<- function(){
	        separator <- tclvalue(Name)
	        separator <<-separator
          nbsequence <- tclvalue(Namenbseq)
	        nbsequence <<-nbsequence}
  calc<-function(){
          security()  
          OnOk()
          if(length(txt4)==1 & length(txt3)==1) pcrcal2()
          if(length(txt4)==1 & length(txt3)!=1) pcrcal3()
          if(length(txt4)!=1) pcrcal4()}
  dprim<-function(){
          security()
          OnOk()
          if(length(txt4)==1 & length(txt3)==1) degencal2()
          if(length(txt4)==1 & length(txt3)!=1) degencal3()
          if(length(txt4)!=1) degencal4()}
  dprimconc<-function(){
          security()
          OnOk()
          concatenate()}
  effic<-function(){
          security()
          OnOk()
          if(length(txt4)==1 & length(txt3)==1) glob2d()
          if(length(txt4)==1 & length(txt3)!=1) glob3d()
          if(length(txt4)!=1) tkmessageBox(message="Sorry this function is not available for more than 3 oligonucleotide databases")}
  impr<-function(){
          security()
          OnOk()
          if(length(txt4)==1 & length(txt3)==1) impr2d()
          if(length(txt4)==1 & length(txt3)!=1) impr3d()
          if(length(txt4)!=1) tkmessageBox(message="Sorry this function is not available for more than 3 oligonucleotide databases")}


#### bouton aide à améliorer
fontHeading <- tkfont.create(family="times",size=9,weight="bold")
fontHeading2 <- tkfont.create(family="times",size=8,weight="bold")
tkgrid(tkbutton(tt,text="H E L P",font=fontHeading,command=oppdf),sticky="e")
####OligoNucleotide data base 1
tt1<- tkframe(tt,relief="groove",borderwidth=3)
labtt1<-tklabel(tt1,text="DATA  MANAGEMENT",font=fontHeading)
ONDB1<- tkframe(tt1)
ONDB1lab<-tklabel(ONDB1,text="OligoNucleotide database 1: ")
ONDB1b1<-tkbutton(ONDB1,text="Load",command=dbb1)
ONDB1b2<-tkbutton(ONDB1,text="Reset",command=reset1)
ONDB1b3<-tkbutton(ONDB1,text="Export",command=export1)
ONDB1b4<-tkbutton(ONDB1,text="Example database 1",command=example1)
ONDB1txt <- tktext(ONDB1,bg="#d8d8d8", width=81,height=1,fg="dark green")
tkpack(ONDB1lab,ONDB1b1,ONDB1b2,ONDB1b3,ONDB1txt,ONDB1b4,side="left")
####OligoNucleotide data base 2
ONDB2<- tkframe(tt1)
ONDB2lab<-tklabel(ONDB2,text="OligoNucleotide database 2: ")
ONDB2b1<-tkbutton(ONDB2,text="Load",command=dbb2)
ONDB2b2<-tkbutton(ONDB2,text="Reset",command=reset2)
ONDB2b3<-tkbutton(ONDB2,text="Export",command=export2)
ONDB2b4<-tkbutton(ONDB2,text="Example database 2",command=example2)
ONDB2txt <- tktext(ONDB2,bg="#d8d8d8", width=81,height=1,fg="dark green")
tkpack(ONDB2lab,ONDB2b1,ONDB2b2,ONDB2b3,ONDB2txt,ONDB2b4,side="left")
####OligoNucleotide data base 3
ONDB3<- tkframe(tt1)
ONDB3lab<-tklabel(ONDB3,text="OligoNucleotide database 3: ")
ONDB3b1<-tkbutton(ONDB3,text="Load",command=dbb3)
ONDB3b2<-tkbutton(ONDB3,text="Reset",command=reset3)
ONDB3b3<-tkbutton(ONDB3,text="Export",command=export3)
ONDB3b4<-tkbutton(ONDB3,text="Example database 3",command=example3)
ONDB3txt <- tktext(ONDB3,bg="#d8d8d8", width=81,height=1,fg="dark green")
tkpack(ONDB3lab,ONDB3b1,ONDB3b2,ONDB3b3,ONDB3txt,ONDB3b4,side="left")
####OligoNucleotide more data bases
ONDB4<- tkframe(tt1)
ONDB4lab<-tklabel(ONDB4,text="         ")
ONDB4b1<-tkbutton(ONDB4,text="Add more oligonucleotide databases",command=dbb4)
ONDB4b2<-tkbutton(ONDB4,text="Reset",command=reset4)
ONDB4lab2<-tklabel(ONDB4,text="",width=16)
ONDB4txt <- tktext(ONDB4,bg="#d8d8d8", width=81,height=1,fg="dark green")
tkpack(ONDB4lab,ONDB4b1,ONDB4b2,ONDB4txt,ONDB4lab2,side="left")
                                     
#####field separator
fieldsep<- tkframe(tt1)
Name <<- tclVar("*")
fieldseptxt<-tklabel(fieldsep,text="Sequence separator in oligonucleotide databases ")
fieldseptxt1<-tklabel(fieldsep,text="  ' * ' is the default of ARB ASCII output")
entry.Name <-tkentry(fieldsep,width="12",textvariable=Name)
tkpack(fieldseptxt,entry.Name,fieldseptxt1,side="left")
tkpack(labtt1,ONDB1,ONDB2,ONDB3,ONDB4,fieldsep,side="top")
tkgrid(tt1)
#####number of sequence
nbseq<- tkframe(tt1)
Namenbseq <<- tclVar("optional")
nbseqtxt<-tklabel(nbseq,text="Number of sequences in the targeted database ")
nbseqtxt1<-tklabel(nbseq,text="                                                             ")
nbseq.entry.Name <-tkentry(nbseq,width="12",textvariable=Namenbseq)
tkpack(nbseqtxt,nbseq.entry.Name,nbseqtxt1,side="left")

tkpack(labtt1,ONDB1,ONDB2,ONDB3,ONDB4,fieldsep,nbseq,side="top")
tkgrid(tt1)
######Choice pcr
tkgrid(tklabel(tt,text="  "))
tt2<- tkframe(tt)
choicePCR<-tkframe(tt2,relief="groove",borderwidth=3)
choicePCRlab1<-tklabel(choicePCR,text="OLIGONUCLEOTIDES SYSTEM",font=fontHeading)
choicePCRlab5<-tklabel(choicePCR,text="based on PCR or hybridization technologies",font=fontHeading2)

choicePCRlab2<-tklabel(choicePCR,text="")
choicePCRlab3<-tklabel(choicePCR,text="")
PCRb1<-tkbutton(choicePCR,text="Global efficiency",command=calc)
PCRb2<-tkbutton(choicePCR,text="All information on oligonucleotide databases",command=effic)
tkpack(choicePCRlab1,choicePCRlab5,choicePCRlab2,PCRb1,choicePCRlab3,PCRb2,side="top")
######Separation
separation<-tkframe(tt2,relief="groove")
separationlab1<-tklabel(separation,text="                              ")
tkpack(separationlab1)
######Choice degenrate primers
choiceDEG<-tkframe(tt2,relief="groove",borderwidth=3)
choiceDEGlab1<-tklabel(choiceDEG,text="DEGENERATE OLIGONUCLEOTIDE WITH MISMATCHES",font=fontHeading)
choiceDEGlab2<-tklabel(choiceDEG,text="")
DEGb1<-tkbutton(choiceDEG,text="Global efficiency",command=dprim)
DEGb3<-tkbutton(choiceDEG,text="All information on oligonucleotide databases",command=effic)
DEGb2<-tkbutton(choiceDEG,text="Concatenate databases of degenerate oligonucleotide",command=dprimconc)
DEGb4<-tkbutton(choiceDEG,text="Improvement of degeneracies",command=impr)
tkpack(choiceDEGlab1,DEGb1,DEGb4,DEGb3,choiceDEGlab2,DEGb2,side="top")
tkpack(choiceDEG,separation,choicePCR,side="left")
tkgrid(tt2)
#######logo et contact
tkgrid(tklabel(tt,text="  "))
im1<-tkframe(tt)  
  zz<-file.path(paste(.libPaths(), "/OligoSpecificitySystem/Rlogo.GIF",sep=""))
  icn<-tkimage.create("photo", file = zz)
  Rlabel <- tklabel(im1, image = icn)
  zzz<-file.path(paste(.libPaths(), "/OligoSpecificitySystem/tcltk.GIF",sep=""))
  icnn<-tkimage.create("photo", file = zzz)
  tcltklab <- tklabel(im1, image = icnn)
  cont<-tklabel(im1,text="Contact : laurent.cauquil@toulouse.inra.fr")
  kk<-tklabel(im1,text="")
  tkpack(Rlabel,tcltklab,cont,kk,side="left")
tkgrid(im1)
tkgrid(tklabel(tt,text="                                                               ")) 
  




tkfocus(tt)

tkbind(ONDB1txt, "<Motion>",.write)
tkbind(ONDB1txt, "<Control-Return>",run)
tkbind(ONDB2txt, "<Motion>",.write)
tkbind(ONDB2txt, "<Control-Return>",run)
tkbind(ONDB3txt, "<Motion>",.write)
tkbind(ONDB3txt, "<Control-Return>",run)
tkbind(ONDB4txt, "<Motion>",.write)
tkbind(ONDB4txt, "<Control-Return>",run)
tkbind(entry.Name, "<Return>",OnOk)
tkbind(nbseq.entry.Name, "<Return>",OnOk)
}


"Oss"   <-function(){OSS()}                                
"OligoSpecificitySystem"   <-function(){OSS()}
"oligospecificitysystem"   <-function(){OSS()}
"Oligospecificitysystem"   <-function(){OSS()}
"oss"   <-function(){OSS()}



 
 
 