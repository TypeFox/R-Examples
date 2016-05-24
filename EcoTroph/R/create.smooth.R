create.smooth <-
function(tab_input,smooth_type=NULL,sigmaLN_cst=NULL,pas=NULL,shift=NULL,smooth_param=NULL)
{
shift_init <- shift
if (is.null(shift)) shift <- 1.8               ##default values
if (is.null(pas)) pas <- 0.1
if (is.null(smooth_param)) smooth_param <- 0.07  
if (is.null(sigmaLN_cst)) sigmaLN_cst <- 0.12 
if (is.null(smooth_type)) smooth_type <- 1 

Troph_round <- tab_input$TL
for (i in 1:length(tab_input$group_name)){
tmp_tl=tab_input[i,]$TL
Troph_round[i] <- seq(0,7,pas)[as.numeric(cut(tmp_tl+0.0000000000001,seq(0,7,pas),right=FALSE))]  ## assignment of a trophic class to the trophic groups
}

TL_out <- c(1,seq(from=2, to=7, by=pas))
tab_smooth <- array(dim=c(length(TL_out),length(tab_input$group_name)),dimnames=list(TL_out,Troph_round))

toto <- c(1,2,3)
if (is.na(match(smooth_type,toto))){
cat(paste("You didn't choose a right value for smooth_type. Type 1, 2 or 3.\n"))
}

else{
if (smooth_type==1){                                  ## constant sigma
sigmaLN <- rep(sigmaLN_cst,length(Troph_round))
}
if (smooth_type==2){
sigmaLN <- smooth_param*log(Troph_round-0.05)        ##lognormal sigma 
if (is.null(shift_init)) shift <- 0.95               ##default values
}
if (smooth_type==3){
sigmaLN <- tab_input$OI/Troph_round
if (is.null(shift_init)) shift <- 0  
}

#Handling of the zero-values in sigmaLN (notably for sigmaLN=OI)
for (i in 1:length(tab_input$group_name)){
if(sigmaLN[i]==0){
sigmaLN[i]<- 0.01
cat(paste("the value of the sigmaLN was 0 for the group", tab_input$group_name[i],". We change it for a value of 0.01. If you want to correct it, change your dataset (fix() or else).\n"))
}}

for (i in 1:length(tab_input$group_name))
for (j in 2:length(TL_out))
if (colnames(tab_smooth)[i]>=2)
tab_smooth[j,i]<-exp(-1/2*((log(TL_out[j]-shift)-log(Troph_round[i]-shift))/sigmaLN[i])^2)/((TL_out[j]-shift)*sigmaLN[i]*((2*pi)^(1/2)))


for(i in 1:length(Troph_round))
if (colnames(tab_smooth)[i]==1)
tab_smooth[1,i]<-1

for (i in 1:length(Troph_round))
if (colnames (tab_smooth)[i]>1 & colnames (tab_smooth)[i]<2){
tab_smooth[,i] <- 0
tab_smooth[1,i] <- 2- as.numeric(colnames(tab_smooth)[i])
tab_smooth[2,i] <- as.numeric(colnames(tab_smooth)[i])-1
}

tab_smooth[is.na(tab_smooth)]<-0

sumcol <- colSums(tab_smooth)
for (i in 1:length(tab_input$group_name))
tab_smooth[,i]<- tab_smooth[,i]/sumcol[i] 
class(tab_smooth)<-"smooth"
return(tab_smooth)
}
}

