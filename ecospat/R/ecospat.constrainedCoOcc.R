
###########################################################
###########################################################
#########        	CO-OCCURRENCE ANALYSIS			     ########	
#########  					         &         					   ########		
#########  ENVIRONMENTALLY CONSTRAINED NULL MODELS ########
#########                  CtRA1                   ######## 
###########################################################
###########################################################
# 
# Adapted by Anne Dubuis from P.R Peres-Neto codes used for: 
#
# Peres-Neto P.R., Olden J.D. and Jackson D.A. (2001) Environmentally constrained null models: sites suitablity as occupancy criterion. 
# Oikos 93:110-120
#
# and modified by Manuela D'Amen
#
# Constrperm() produces a null matrix of 0/1 in which the 1 are weighted by probability values
# from the db "proba". The sum column is costant i.e. null species have the same occurrences as the observed ones
# while row sum is free: richness is variable in the null community sites.
#
# SpeciesCooccurrenceStats() produces a classical C-score index for any confusion matrix: this is used to calculate
# both the C-score for the observed dataset and for the null matrices from  the former functions
#
# NullModels() uses the former two functions and returns the C-score index for the observed community (ObsCscoreTot), the mean of C-score for the simulated communities (SimCscoreTot)
# p.value (PValTot) and standardized effect size (SES.Tot). It saves also a table in the path specified where the same 
# metrics are calculated for each species pair (only the table with species pairs with significant p.values is saved in this version)
# # NOTE: a SES that is greater than 2 or less than -2 is statistically significant with a tail probability of less than 0.05 (Gotelli & McCabe 2002 - Ecology)
#
# presence: presence absence table
# pred: table with values of probability of presence for each species in each site
# nbpermut: number of permutation in the null model
# outpath: path to specify where to save the tables
# NB: Format required for imput databases: a plots (rows) x species (columns) matrix
# Input matrices should have column names (species names) and row names (sampling plots)


##################################################################################################################################################

# NullModels(presence, pred, 1000, outpath) # launch command
#library(ade4)

#############################
## Constrained permutation ##
#############################
Constrperm<-function(presence, pred) 
{
  nbsps<-ncol(presence) # species number
  nbsites<-nrow(presence) # sites number
  
  nbocc<- as.vector(apply(presence,MARGIN=2,sum)) # occurrences number for each species
  sumprob<-as.vector(apply(pred,MARGIN=2,sum)) # occurrence probability sum for each species
  
  matsum<-matrix(0,ncol=nbsps,nrow=nbsites)
  
  for(i in 1:nbsps)			#loop per ciascuna specie
  {
    matsum[,i]<-sumprob[i] #riempie la matrice "specie X siti" della somma delle prob di ciascuna specie (quindi le colonne contengono valori tutti uguali)
  }
  
  # relative probability matrix
  pred<-pred/matsum  # divisione tra matrici per avere probabilita relative ogni valore e diviso per la somma totale delle probablita per ogni specie 
  
  
  transpo<-t(as.matrix(presence))  # trasposizione della matrice 0/1 
  sps.names<-row.names(transpo) # salva i nomi delle specie 
  sites<-row.names(presence) # salva i nomi dei siti 
  noms<-list(sites, sps.names) # lista dei 2 vettori con i nomi
  
  nullmat<-matrix (0,nrow=nbsites,ncol=nbsps,dimnames=noms)	# matrice vuota specie X sito 
  randr<-matrix(runif(nbsites*nbsps),nbsites)	# random number matrix - stesse dimensioni di proba # QUESTO E' IL PASSAGGIO DI RANDOMIZZAZIONE!
  randr<-randr*pred #matrice con valori random ma pesati per le probabilita di ogni specie in ogni sito	
  
  for (i in 1: nbsps) # matrix randomisation
  {		
    a<-as.vector(randr[,i]) # per ogni specie un vettore 
    names(a)<-(1:nbsites) # nomi del vettore: un altro vettore di numeri continui
    a<-sort(a)  # il vettore viene ordinato secondo probabilita crescente e con esso i numeri che sono i nomi di riga         
    b<-as.numeric(names(a))	# un secondo vettore fatto solo dei nomi resi numerici e ordinati secondo le probabilita crescente di a 
    x<-nbsites-nbocc[i]+1  # numero di siti non occupati dalla specie 
    
    for(j in x:nbsites){   # mette tanti 1 quanti erano nella distribuzione originale per i valori piu alti di probabilita
      r<-b[j]
      nullmat[r,i]<-1
    }    
  }
  return (nullmat)
}

### END function ###


##############################
## Co-occurrence statistics ##                      
## 		calculations	    ##	

##############################
SpeciesCooccurrenceStats<-function (presence)  
{
  nbsps<-ncol(presence) # species number
  nbsites<-nrow(presence) # sites number
  nbocc<- as.vector(apply(presence,MARGIN=2,sum))	# occurrences number for each species
  
  presence<-as.matrix(presence)       # sps trasformato in matrice 
  coocc<-t(presence)%*%presence	      # produit matricielle de la matrice presence abscence par sa transposee donne le nb de checkboard unit pour chaque paire d'sps
  nbspec=dim(coocc)[1]			# le nombre d'espece= premiere dimension de la matrice obtenue
  mat1<-array(apply(presence,MARGIN=2,sum),dim=c(nbsps,nbsps)) # creation d'une matrice carree de dimension n.spec avec les nb occ pr chaque sps (le mm chiffre sur tt la ligne)
  mat2<-t(array(apply(presence,MARGIN=2,sum),dim=c(nbsps,nbsps))) # transpo de mat1
  ### le calcul matriciel permet d'ameliorer enormement la rapidite de calcul par l'ordi, par contre ca ne rend pas la comprehension tres aisee:)
  Cscoreperspecies <- (mat1 - coocc)*(mat2 - coocc) # c-score per species
  
  df.synthesis1 <- data.frame(Col = rep(1:ncol(Cscoreperspecies),each=ncol(Cscoreperspecies)),  # dataframe creation to summarize results
                              Row = rep(1:nrow(Cscoreperspecies),nrow(Cscoreperspecies)),Sps1 = rep(colnames(Cscoreperspecies),
                                                                                                    each=ncol(Cscoreperspecies)), Sps2 = rep(rownames(Cscoreperspecies),nrow(Cscoreperspecies)),CScore= c(Cscoreperspecies)) 
  v.diago.inf <- c(rownames(df.synthesis1)[df.synthesis1[,1]>df.synthesis1[,2]],rownames(df.synthesis1)[df.synthesis1[,1]==df.synthesis1[,2]])
  df.synthesis <- df.synthesis1[-as.numeric(v.diago.inf),]
  
  ### ci-dessus, creation d'un dataframe en collant les numeros des colonnes et des lignes, les noms des especes ainsi que les valeurs de scores
  ### puis on enleve les paire d'especes a double ainsi que la diagonale (2 mm sps)
  
  Cscore<-mean(df.synthesis[,5]) # c-score on all matrix   c-score sur toute la matrice = moyenne des c-score pour chaque paire d'sps
  
  return(list(coocc=Cscoreperspecies, synth=df.synthesis1, synth_lt=df.synthesis, Cscore=Cscore))
  
}### END function ###


#################
## Null models ##
#################
#nbpermut=100

ecospat.cons_Cscore<- function(presence,pred,nbpermut,outpath)
{ 	
  cat("Computing observed co-occurence matrix", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  
  nbsps<-ncol(presence) 
  nbsites<-nrow(presence) 
  
  Obs<-SpeciesCooccurrenceStats(presence) # cooccurrence statistics on observed matrix
  
  CooccObs<-Obs$coocc # double matrix of species with co-occurrence values for each couple
  synthObs<-Obs$synth_lt # summary table with co-occurrence values
  CscoreTot<-Obs$Cscore
  
  CooccProb<-matrix(0,nrow=nrow(synthObs),ncol=nbpermut, dimnames = list(c(paste(synthObs[,3],synthObs[,4])),c(1:nbpermut)) )
  
  cat("Computing permutations", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  
  if (nbpermut>0){
    for (z in 1:nbpermut){	# here in the Anne's version (z in 2:nbpermut)..why 2??
      #cat(z, "\n",append =T)
      degenerated<-TRUE
      while (degenerated==TRUE){ # looking for a matrix without degenerated sites
        random.distrib.matrix<-Constrperm(presence,pred) # null matrix (0/1) column sum is the true species n obs
        sumSites<- as.vector(apply(random.distrib.matrix,MARGIN=2,sum)) # row sum: sites richness 
        for (i in 1:nbsites){
          ifelse(sumSites[i]==0,degenerated<-TRUE,degenerated<-FALSE)		
        }
      }
      
      Rnd<-SpeciesCooccurrenceStats(random.distrib.matrix) # cooccurrence statistics on randomized matrix
      
      synthRnd<-Rnd$synth_lt	
      
      CooccProb[,z] <- synthRnd[,5] #qui si crea la matrice con i valori di c-score delle coppie di specie delle comunita nulle 
    }
    
    ## for the whole community
    vec.CScore.tot<-as.vector(apply(CooccProb,MARGIN=2,mean)) # C-score for all null communities (mean on the columns)
    SimulatedCscore<-mean(vec.CScore.tot) # mean of Simulation C-score: Simulated C-score
    sd.SimulatedCscore<-sd(vec.CScore.tot) # standard deviation of null communities
    ses<-(CscoreTot-SimulatedCscore)/sd.SimulatedCscore # standardized effect size:  It scales the results in units of standard deviations, 
    # which allows for meaningful comparisons among different tests
    
    randtest.less<-as.randtest(vec.CScore.tot, CscoreTot, alter="less")
    pval.less<-randtest.less$pvalue
    randtest.greater<-as.randtest(vec.CScore.tot, CscoreTot, alter="greater")
    pval.greater<-randtest.greater$pvalue
    plot(randtest.greater, xlab= "Simulated C-scores",main=paste("", sep=""))
    
    ## for species pairs
    vec.Cscore.pairs<-as.vector(apply(CooccProb,MARGIN=1,mean)) # Mean of null communities C-score for any pair of species
    mat_pval <- matrix(0,nrow(synthRnd),5) # matrice per raccogliere i valori di probabilita per ogni coppia di C-score 
    mat_pval[,1] <- synthObs[,5] # Observed C-score in the first column
    mat_pval[,2] <- vec.Cscore.pairs # Mean of simulated C-scores in the 2nd column
    
    for (i in 1:nrow(CooccProb))
    {		  
      randtest.less<-as.randtest(CooccProb[i,], mat_pval[i,1], alter="less")
      mat_pval[i,3]<-randtest.less$pvalue
      randtest.greater<-as.randtest(CooccProb[i,], mat_pval[i,1], alter="greater")
      mat_pval[i,4]<-randtest.greater$pvalue
      mat_pval[i,5]<-(mat_pval[i,1]-mean(CooccProb[i,]))/(sd(CooccProb[i,])) #ses for each couple
    }    
  }	
  
  cat("Permutations finished",date(), "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  cat("Exporting dataset", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  cat(".............", "\n",append = F)
  
  df.synthesis <- data.frame(Col = rep(1:ncol(CooccObs),each=ncol(CooccObs)),Row = rep(1:nrow(CooccObs),
                                                                                       nrow(CooccObs)),Sps1 = rep(colnames(CooccObs),each=ncol(CooccObs)), 
                             Sps2 = rep(rownames(CooccObs),nrow(CooccObs)),Cooc_score = c(CooccObs)) 	
  v.diago.inf <- c(rownames(df.synthesis)[df.synthesis[,1]>df.synthesis[,2]],rownames(df.synthesis)[df.synthesis[,1]==df.synthesis[,2]])
  df.synthesis <- df.synthesis[-as.numeric(v.diago.inf),]
  df.synthesis<- cbind(df.synthesis,mat_pval[,2:5])
  names(df.synthesis)[5:9]<-c("C-scoreObs","C-scoreExp","p.less","p.greater","ses")
  
  hist(as.vector(df.synthesis[,9]), xlab="ses", main = paste("Histogram of standardized effect size"))
  abline(v=c(2,-2),col = "red")
  #write.table(df.synthesis,file=paste(outpath,"Results_const_Cscores.txt", sep=""),sep="\t",append=F,row.names=F,col.names=T,quote=F)
  
  #estrazione solo delle righe con p <=0.05
  tab<-df.synthesis
  v<-c(0)
  for (i in 1:nrow(tab)){
    if (tab[i,7]<=0.05||tab[i,8]<=0.05){
      v<-c(v,i)
    }
  }
  m<-data.frame()
  for(j in 1:length(v)){
    m<-rbind(m,tab[v[j],])
  }
  
  m1<-na.omit(m)
  
  write.table(m1,file=paste(outpath,"/Signific_const_Cscores.txt",sep=""),sep="\t",append=F,row.names=F,col.names=T,quote=F)
  
  l<-list(ObsCscoreTot=CscoreTot, SimCscoreTot=SimulatedCscore, PVal.less=pval.less,PVal.greater=pval.greater,SES.Tot=ses)
  return(l)
  
}
### END function ###








