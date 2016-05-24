class.percent <-
function(data, mode="both", empty.col=FALSE, lang="en-US")
{
 clz <- as.numeric(data[1,])
 if (sum(clz) > 550){
   clz<-(-log2(clz/1000))
 }

 tab <- data[2:nrow(data),]
 tab.res <- as.data.frame(tab[,1],row.names=rownames(tab)) 
 m<-0 

 if (mode=="both")
 {
  tab.res <- as.data.frame(cbind(tab.res,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")
   colnames(tab.res) <- c("Cascalho","Areia","Silte","Argila","Matac\u00E3o","Calhau","Seixo","Gr\u00E2nulo","Areia.Muito.Grossa",
"Areia.Grossa","Areia.M\u00E9dia","Areia.Fina","Areia.Muito.Fina","Silte.Grosso","Silte.M\u00E9dio","Silte.Fino","Silte.Muito.Fino","Argilas")
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")
   colnames(tab.res) <- c("Gravel","Sand","Silt","Clay","Boulder","Cobble","Pebble","Granules","Very.Coarse.Sand",
"Coarse.Sand","Medium.Sand","Fine.Sand","Very.Fine.Sand","Coarse.Silt","Medium.Silt","Fine.Silt","Very.Fine.Silt","Clays")
 }

 if (mode=="total")
 {
  tab.res <- as.data.frame(cbind(tab.res,NA,NA,NA))
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p") 
   colnames(tab.res) <- c("Cascalho","Areia","Silte","Argila")
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e") 
   colnames(tab.res) <- c("Gravel","Sand","Silt","Clay")
 }

 if (mode=="classes")
 {
  tab.res <- as.data.frame(cbind(tab.res,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p") 
   colnames(tab.res) <- c("Matac\u00E3o","Calhau","Seixo","Gr\u00E2nulo","Areia.Muito.Grossa",
"Areia.Grossa","Areia.M\u00E9dia","Areia.Fina","Areia.Muito.Fina","Silte.Grosso","Silte.M\u00E9dio","Silte.Fino","Silte.Muito.Fino","Argilas")
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e") 
   colnames(tab.res) <- c("Boulder","Cobble","Pebble","Granules","Very.Coarse.Sand",
"Coarse.Sand","Medium.Sand","Fine.Sand","Very.Fine.Sand","Coarse.Silt","Medium.Silt","Fine.Silt","Very.Fine.Silt","Clays")
 }

 for (j in 1:nrow(tab)) 
 {  
  vbz <- as.numeric(tab[j,])
  soma <- sum(vbz) 
  vfz <- vbz*100/soma 
  r<-length(clz)

  p.sand <- 0
  p.gravel <- 0 
  p.silt <- 0 
  p.clay <- 0 
  p.boulder <- 0
  p.cobble <- 0 
  p.pebble <- 0 
  p.granules <- 0 
  p.very.coarse.sand <- 0 
  p.coarse.sand <- 0 
  p.medium.sand <- 0
  p.fine.sand <- 0 
  p.very.fine.sand <- 0 
  p.coarse.silt <- 0 
  p.medium.silt <- 0 
  p.fine.silt <- 0 
  p.very.fine.silt <- 0 
  p.clays <- 0 

  for (i in 2:r) 
  {
   if (clz [i] <= -8) p.boulder <- p.boulder + vfz[i]   
   if ((clz [i] > -8) & (clz [i] <= -6)) p.cobble <- p.cobble + vfz[i]
   if ((clz [i] > -6) & (clz [i] <= -2)) p.pebble <- p.pebble +vfz[i]
   if ((clz [i] > -2) & (clz [i] <= -1)) p.granules <- p.granules + vfz[i]
   if ((clz [i] > -1) & (clz [i] <= 0)) p.very.coarse.sand <- p.very.coarse.sand + vfz[i]
   if ((clz [i] > 0) & (clz [i] <= 1)) p.coarse.sand <- p.coarse.sand + vfz[i]
   if ((clz [i] > 1) & (clz [i] <= 2)) p.medium.sand <- p.medium.sand + vfz[i]
   if ((clz [i] > 2) & (clz [i] <= 3)) p.fine.sand <- p.fine.sand +vfz[i]
   if ((clz [i] > 3) & (clz [i] <= 4)) p.very.fine.sand <- p.very.fine.sand + vfz[i]
   if ((clz [i] > 4) & (clz [i] <= 5)) p.coarse.silt <- p.coarse.silt + vfz[i]
   if ((clz [i] > 5) & (clz [i] <= 6)) p.medium.silt <- p.medium.silt + vfz[i]
   if ((clz [i] > 6) & (clz [i] <= 7)) p.fine.silt <- p.fine.silt + vfz[i]
   if ((clz [i] > 7) & (clz [i] <= 8)) p.very.fine.silt <- p.very.fine.silt + vfz[i]
   if (clz [i] > 8) p.clays <- p.clays + vfz[i]
  

   p.sand <- p.very.coarse.sand + p.coarse.sand + p.medium.sand + p.fine.sand + p.very.fine.sand   
   p.silt <- p.coarse.silt + p.medium.silt + p.fine.silt + p.very.fine.silt  
   p.clay <- p.clays   
   p.gravel <- 100 -(p.sand + p.silt + p.clay) 
   if (p.gravel < 0.0001) p.gravel <-0   
  }

  if (mode=="both") 
  {
   tab.res[j,1] <- p.gravel
   tab.res[j,2] <- p.sand
   tab.res[j,3] <- p.silt
   tab.res[j,4] <- p.clay
   tab.res[j,5] <- p.boulder
   tab.res[j,6] <- p.cobble
   tab.res[j,7] <- p.pebble
   tab.res[j,8] <- p.granules
   tab.res[j,9] <- p.very.coarse.sand
   tab.res[j,10] <- p.coarse.sand
   tab.res[j,11] <- p.medium.sand
   tab.res[j,12] <- p.fine.sand
   tab.res[j,13] <- p.very.fine.sand
   tab.res[j,14] <- p.coarse.silt
   tab.res[j,15] <- p.medium.silt
   tab.res[j,16] <- p.fine.silt
   tab.res[j,17] <- p.very.fine.silt
   tab.res[j,18] <- p.clays
  } 

  if (mode=="classes") 
  {
   tab.res[j,1] <- p.boulder
   tab.res[j,2] <- p.cobble
   tab.res[j,3] <- p.pebble
   tab.res[j,4] <- p.granules
   tab.res[j,5] <- p.very.coarse.sand
   tab.res[j,6] <- p.coarse.sand
   tab.res[j,7] <- p.medium.sand
   tab.res[j,8] <- p.fine.sand
   tab.res[j,9] <- p.very.fine.sand
   tab.res[j,10] <- p.coarse.silt
   tab.res[j,11] <- p.medium.silt
   tab.res[j,12] <- p.fine.silt
   tab.res[j,13] <- p.very.fine.silt
   tab.res[j,14] <- p.clays
  }

  if (mode=="total") 
  {
   tab.res[j,1] <- p.gravel
   tab.res[j,2] <- p.sand
   tab.res[j,3] <- p.silt
   tab.res[j,4] <- p.clay
  }
 }

 if (empty.col==FALSE) 
 {
   for (k in 1:length(tab.res)) 
   {
     if (colnames(tab.res[k])=="Sand" | colnames(tab.res[k])=="Silt" | colnames(tab.res[k])=="Clay" | colnames(tab.res[k])=="Areia" | colnames(tab.res[k])=="Silte" | colnames(tab.res[k])=="Argila"){
       m<-m+1   
       output[m]<-tab.res[k]
       next
     }
     if (sum(tab.res[k]) > 0.0001)
     {
       m<-m+1
       if (m==1) {output<-tab.res[k]} else {output[m]<-tab.res[k]}
     }
   }
 }
 else {output<-tab.res} 

 return(output) 

}

