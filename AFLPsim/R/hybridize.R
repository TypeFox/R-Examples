hybridize <-
function(pa,pb,Nf1,Nbxa=Nf1, Nbxb=Nf1, Nf2=Nf1, type="selection",hybrid="all",Nsel=Nmarker*0.1,S=0){

   Nmarker=ncol(pa) 
   if (ncol(pa)!= ncol(pb))
     stop("You need the same number of markers in parentals")
   if (Nsel==0){ type="neutral"}
   if (Nf1==0)
    stop("at least 2 F1 hybrids")  
   if (ncol(pa)!=ncol(pb))
    stop("different number of fragments")
   if (type=="neutral"){
    cat("########Neutral hybridization########")
    Nsel=Nmarker
    S=0
  }
   markers<-paste("M",seq(1,Nmarker,1),sep="")
   PA<-paste("PA",seq(1,nrow(pa),1),sep="")
   PB<-paste("PB",seq(1,nrow(pb),1),sep="")
   
   colnames(pa)<-markers
   rownames(pa)<-PA
   colnames(pb)<-markers
   rownames(pb)<-PB
   
   simulation<-list(PA= pa,PB= pb)

  #Se calculas las frecuencias de lso parentales
  pa<-as.matrix(pa)
  pb<-as.matrix(pb)
  fa1000<-apply(pa,2,mean)
  fb1000<-apply(pb,2,mean)
  pa1000<-1-sqrt(1-fa1000)
  pb1000<-1-sqrt(1-fb1000)
  
  sel<-sample(1:Nmarker,Nsel)
  sel<-sort(sel)
  
 if (any(hybrid == "all" | hybrid=="F1"))
    {
    #El modelo de introgresion neutra
    neutralmodel<-function(x,y){x+y-(x*y)}
    
    ff1<-neutralmodel(pa1000,pb1000)
    
    progeny<-matrix(1,1000,Nmarker)
    
    for (i in 1:Nmarker)
      progeny[,i]<-rbinom(1000,1,ff1[i])
    
    # Aqui es donde metemos seleccion. Para ello hay que definir sel 
    # que es un vector de longitud Nsel que equivale 
    # al numero de fragmentos a seleccionarse positivamente
    ## Nsel<=Nmarker
  
    #Ya tenemos los fragmentos a seleccionar. Ahora vamos a seleccionar los fragmentos sel con una coeficiente de selecci??n S [1-10]
    y2 <-(1-pa1000)*(1-pb1000)
    y2[y2==0]<-0.00000000000000000000000000000000001
    y2<-y2[sel]
    selection<-function(z){
      ((S+1)*z)/(((S+1)*z)+ y2)}
    
    fS<-selection(ff1[sel])
    
    selected<-matrix(1,1000,Nsel)
    for (i in 1: Nsel)
      selected[,i]<-rbinom(1000,1,fS[i])
    progeny[,sel]<-selected
    progeny[progeny=="NaN"]<-0
    
    #Submuestreamos al azar un Numero Nf1 de hibridos con seleccion sobre esos marker
    f1sel<-progeny[sample(nrow(progeny), Nf1), ]
    
    #La ponemos bonita nombres de alelos e individuos

    F1<-paste("F1",seq(1,nrow(f1sel),1),sep="_")
    colnames(f1sel)<-markers
    rownames(f1sel)<-F1
    simulation <- c(simulation, list(F1= f1sel))
  }
  else {simulation <- c(simulation, list(F1= NA))}
  
  if (any(hybrid == "all" | hybrid=="BxA")){
    
    #El modelo de introgresion neutra
    neutralmodelbxa<-function(x,y){
      ((3*x)+y-(x^2)-(x*y))/2}
    
    fbxa<-neutralmodelbxa(pa1000,pb1000)
    
    progeny<-matrix(1,1000,Nmarker)
    
    for (i in 1:Nmarker)
      progeny[,i]<-rbinom(1000,1,fbxa[i])
    
    # Aqui es donde metemos seleccion. Para ello hay que definir sel 
    # que es un vector de longitud Nsel que equivale 
    # al numero de fragmentos a seleccionarse positivamente
    ## Nsel<=Nmarker

    #Ya tenemos los fragmentos a seleccionar. Ahora vamos a seleccionar los fragmentos sel con una coeficiente de selecci??n S [1-10]
    y2 <-((pa1000^2)-(3*(pa1000))- pb1000 + (pa1000*pb1000) + 2 )/2
    y2[y2==0]<-0.00000000000000000000000000000000001
    y2<-y2[sel]
    selection<-function(z){
      ((S+1)*z)/(((S+1)*z)+ y2)}
    
    fS<-selection(fbxa[sel])
    
    selected<-matrix(1,1000,Nsel)
    for (i in 1: Nsel)
      selected[,i]<-rbinom(1000,1,fS[i])
    progeny[,sel]<-selected
    progeny[progeny=="NaN"]<-0
    
    #Submuestreamos al azar un Numero Nbxa de hibridos con seleccion sobre esos marker
    bxasel<-progeny[sample(nrow(progeny), Nbxa), ]
    
    #La ponemos bonita nombres de alelos e individuos

    BxA<-paste("BxA",seq(1,nrow(bxasel),1),sep="_")

    colnames(bxasel)<-markers
    rownames(bxasel)<-BxA
    simulation <- c(simulation, list(BxA= bxasel))
    }
  else {simulation <- c(simulation, list(BxA= NA))}
  
  if (any(hybrid == "all" | hybrid=="BxB")){
    
    #El modelo de introgresion neutra
    neutralmodelbxb<-function(x,y){
      ((3*x)+y-(x^2)-(x*y))/2}
    
    fbxb<-neutralmodelbxb(pb1000,pa1000)
    
    progeny<-matrix(1,1000,Nmarker)
    
    for (i in 1:Nmarker)
      progeny[,i]<-rbinom(1000,1,fbxb[i])
    
    # Aqui es donde metemos seleccion. Para ello hay que definir sel 
    # que es un vector de longitud Nsel que equivale 
    # al numero de fragmentos a seleccionarse positivamente
    ## Nsel<=Nmarker
    sel<-sample(1:Nmarker,Nsel)
    sel<-sort(sel)
    #Ya tenemos los fragmentos a seleccionar. Ahora vamos a seleccionar los fragmentos sel con una coeficiente de selecci??n S [1-10]
    y2 <-((pb1000^2)-(3*(pb1000))- pa1000 + (pa1000*pb1000) + 2 )/2
    y2[y2==0]<-0.00000000000000000000000000000000001
    y2<-y2[sel]
    selection<-function(z){
      ((S+1)*z)/(((S+1)*z)+ y2)}
    
    fS<-selection(fbxb[sel])
    
    selected<-matrix(1,1000,Nsel)
    for (i in 1: Nsel)
      selected[,i]<-rbinom(1000,1,fS[i])
    progeny[,sel]<-selected
    progeny[progeny=="NaN"]<-0
    
    #Submuestreamos al azar un Numero Nbxa de hibridos con seleccion sobre esos marker
    bxbsel<-progeny[sample(nrow(progeny), Nbxb), ]
    
    #La ponemos bonita nombres de alelos e individuos

    BxB<-paste("BxB",seq(1,nrow(bxbsel),1),sep="_")
    
    colnames(bxbsel)<-markers
    rownames(bxbsel)<-BxB
    
    simulation<-c(simulation,list(BxB=bxbsel))
  }
  else {simulation <- c(simulation, list(BxB= NA))}
  
  if (any(hybrid == "all" | hybrid=="F2")){
    
    #El modelo de introgresion neutra
    neutralmodel2<-function(x,y){
      x + y - (0.5*(x*y)) - (0.25*(x^2+y^2))}
      
          
    f2<-neutralmodel2(pa1000,pb1000)
    
    progeny<-matrix(1,1000,Nmarker)
    
    for (i in 1:Nmarker)
      progeny[,i]<-rbinom(1000,1,f2[i])
    
    # Aqui es donde metemos seleccion. Para ello hay que definir sel 
    # que es un vector de longitud Nsel que equivale 
    # al numero de fragmentos a seleccionarse positivamente
    ## Nsel<=Nmarker
    sel<-sample(1:Nmarker,Nsel)
    sel<-sort(sel)
    #Ya tenemos los fragmentos a seleccionar. Ahora vamos a seleccionar los fragmentos sel con una coeficiente de selecci??n S [1-10]
    y2 <-((pa1000^2)-(3*(pa1000))- pb1000 + (pa1000*pb1000) + 2 )/2
    y2[y2==0]<-0.00000000000000000000000000000000001
    y2<-y2[sel]
    selection<-function(z){
      ((S+1)*z)/(((S+1)*z)+ y2)}
    
    fS<-selection(f2[sel])
    
    selected<-matrix(1,1000,Nsel)
    for (i in 1: Nsel)
      selected[,i]<-rbinom(1000,1,fS[i])
    progeny[,sel]<-selected
    progeny[progeny=="NaN"]<-0
    
    #Submuestreamos al azar un Numero Nbxa de hibridos con seleccion sobre esos marker
    f2sel<-progeny[sample(nrow(progeny), Nf2), ]
    
    #La ponemos bonita nombres de alelos e individuos

    F2<-paste("F2",seq(1,nrow(f2sel),1),sep="_")
    
    colnames(f2sel)<-markers
    rownames(f2sel)<-F2
    
    simulation<-c(simulation,list(F2=f2sel))
  }
  else {simulation <- c(simulation, list(F2= NA))}
  
    
  if (S==0) {sel= NA}
  simulation<-c(simulation,list(SelMarkers=sel, S = S))
  class(simulation) <- "hybridsim"
  simulation
}
