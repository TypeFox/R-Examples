hybridsim <-function(Nmarker,Na,Nb,Nf1,Nbxa=Nf1, Nbxb=Nf1, Nf2=Nf1, type="selection",hybrid="all",Nsel=Nmarker*0.1,S=0,apa=0.5,apb=0.5){
  if (Na==0 |Nb==0)
    stop("at least 1 individuals in each parental population")
  if (Na==0 |Nb==0)
    stop("at least 1 individuals in each parental population")
    #if (Nf1==0 && any(hybrid=="F1" | hybrid=="all"))
   # stop("at least 2 F1 hybrids")
  if (type=="neutral"){
    cat("########Neutral hybridization########")
    Nsel=Nmarker
    S=0}
   #Crea matriz de individuos con todo 1
  pa<-matrix( 1,1000,Nmarker)
  #Crea la matriz con alelos siguiendo la probabilida de la beta
  p1<-rbeta(Nmarker,apa,apa)
  fa<-p1*(2-p1)
  for (i in 1:Nmarker)
  pa[,i]<-rbinom(1000,1,fa[i])
  pasel<-pa[sample(nrow(pa), Na), ]
  
  pb<-matrix(1,1000,Nmarker)
  p2<-rbeta(Nmarker,apb,apb)
  fb<-p2*(2-p2)
  for (i in 1:Nmarker)
    pb[,i]<-rbinom(1000,1,fb[i])
    pbsel<-pb[sample(nrow(pb), Nb), ]
  
  #Creamos los individuos F1 siguiendo un modelo neutral pero queremos seleccionar algunos alelos. Para ello fabricamos un numero grande de progenie y luego recogemos los individuos que tengan los fragmentos outlier donde outlier es un vector de alelos (M1,M2,...)
  #Se calculas las frecuencias de los parentales
  fa1000<-apply(pa,2,mean)
  fb1000<-apply(pb,2,mean)
  pa1000<-1-sqrt(1-fa1000)
  pb1000<-1-sqrt(1-fb1000)
  #Markers under selection
  sel<-sample(1:Nmarker,Nsel)
  sel<-sort(sel)
    
  markers<-paste("M",seq(1,Nmarker,1),sep="")
  PA<-paste("PA",seq(1,nrow(pasel),1),sep="")
  PB<-paste("PB",seq(1,nrow(pbsel),1),sep="")
  
  colnames(pasel)<-markers
  rownames(pasel)<-PA
  colnames(pbsel)<-markers
  rownames(pbsel)<-PB
  simulation<-list(PA= pasel,PB= pbsel)
  
  if (any(hybrid == "all" | hybrid=="F1")){
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
    

    #Ya tenemos los fragmentos a seleccionar. Ahora vamos a seleccionar los fragmentos sel con una coeficiente de selecci??n S [1-10]
    ######Errrorrr
    y2 <-(1-(((pa1000)+(pb1000))/2))^2
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
