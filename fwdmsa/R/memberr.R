memberR <- function(X,minsize){
#membership of observations to bins based on restscores
#this function is needed for estimating MSA during the forward search
#needs library(mokken)
X.tot <- apply(X,1,sum)
N <- dim(X)[1]
J <- dim(X)[2]
R <- as.matrix(X) %*% (matrix(1, J, J) - diag(J))    #restscores voor alle persone voor alle items
mono <- check.monotonicity(X)
lo.mono <- matrix(,round(N/minsize)+1,J)
hi.mono <- matrix(,round(N/minsize)+1,J)
n.mono <- matrix(,round(N/minsize)+1,J)
for(j in 1:J){
lo.mono[,j] <- c(length(mono$results[[j]][[2]][,2]),mono$results[[j]][[2]][,2],rep(NA,round(N/minsize)-length(mono$results[[j]][[2]][,2])))
hi.mono[,j] <- c(length(mono$results[[j]][[2]][,3]),mono$results[[j]][[2]][,3],rep(NA,round(N/minsize)-length(mono$results[[j]][[2]][,2])))
n.mono[,j] <- c(length(mono$results[[j]][[2]][,4]),mono$results[[j]][[2]][,4],rep(NA,round(N/minsize)-length(mono$results[[j]][[2]][,4])))}
#n.mono: eerste rij geeft aantal restscoregroepen voor dat item weer

#in welke restscore group zitten de personen
group.member <- list()
for(j in 1:J){
        #binning group membership:
        group.member[[j]] <- list()
        for(l in 1:dim(mono$results[[j]][[2]])[1]){
            group.member[[j]][[l]] <- which(R[,j]>=mono$results[[j]][[2]][l,2] & R[,j]<=mono$results[[j]][[2]][l,3])
            }
    }
#group.member[[1]][[3]] #dit welke personen zitten in van item 1 in restscoregroep 3

gr.member <- matrix(,N,J)
for(j in 1:J){
for(i in 1:N){
for(bin in 1:dim(mono$results[[j]][[2]])[1]){
gr.member[group.member[[j]][[bin]],j]=bin}}}

member.list <- list(lo.mono=lo.mono,hi.mono=hi.mono,n.mono=n.mono,group.member=group.member,gr.member=gr.member)
return(member.list)
}
