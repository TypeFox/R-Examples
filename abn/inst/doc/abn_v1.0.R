### R code from vignette source 'abn_v1.0.Rnw'

###################################################
### code chunk number 1: abn_v1.0.Rnw:171-176 (eval = FALSE)
###################################################
## library( abn)    # Load the library
## 
## mydat <- pigs.vienna[,-11]  # Get data, drop batch variable
## 
## mydat[ complete.cases(mydat),]


###################################################
### code chunk number 2: abn_v1.0.Rnw:180-181 (eval = FALSE)
###################################################
## mydat[,1] <- as.factor(mydat[,1])


###################################################
### code chunk number 3: abn_v1.0.Rnw:184-185 (eval = FALSE)
###################################################
## levels( mydat[,1]) 


###################################################
### code chunk number 4: abn_v1.0.Rnw:196-197 (eval = FALSE)
###################################################
## ban <- matrix( rep(0,dim(mydat)[2]^2),ncol=dim(mydat)[2])


###################################################
### code chunk number 5: abn_v1.0.Rnw:200-204 (eval = FALSE)
###################################################
## colnames( ban) <- rownames(ban) <- names(mydat)
## 
## retain <- matrix( rep(0,dim(mydat)[2]^2),ncol=dim(mydat)[2])
## colnames( retain) <- rownames(retain) <- names(mydat)


###################################################
### code chunk number 6: abn_v1.0.Rnw:207-212 (eval = FALSE)
###################################################
## mydists <- list( PC="binomial", PT="binomial", MS="binomial", 
##                  HS="binomial", TAIL="binomial", 
##                  Abscess="binomial", Pyaemia="binomial", 
##                  EPcat="binomial", PDcat="binomial", 
##                  plbinary="binomial")


###################################################
### code chunk number 7: abn_v1.0.Rnw:215-218 (eval = FALSE)
###################################################
## mycache <- buildscorecache( data.df=mydat, 
##            data.dists=mydists, dag.banned=ban, 
##            dag.retained=retain, max.parents=max.par) 


###################################################
### code chunk number 8: abn_v1.0.Rnw:221-227 (eval = FALSE)
###################################################
## mp.dag <- mostprobable( score.cache=mycache)
## 
## fabn <- fitabn( dag.m=mp.dag,data.df=mydat,
##                 data.dists=mydists) 
## 
## datadir <- tempdir()


###################################################
### code chunk number 9: abn_v1.0.Rnw:230-232 (eval = FALSE)
###################################################
## save( mycache, mp.dag, fabn, file=
##       paste(datadir,"mp",max.par,".RData",sep=""))


###################################################
### code chunk number 10: abn_v1.0.Rnw:246-253 (eval = FALSE)
###################################################
## load( "RData/Rout1.RData") 
## ml1 <- fabn$mlik;  
## 
## tographviz( dag.m=mp.dag,data.df=mydat,data.dists=mydists,
##             outfile="DAG_cycle.dot")
## system( "dot -Tpdf -o DAG_cycle.pdf DAG_cycle.dot")
## system( "evince DAG_cycle.pdf&")


###################################################
### code chunk number 11: abn_v1.0.Rnw:256-257 (eval = FALSE)
###################################################
## mp.mlik <- c( -44711.62,-44685.53,-44684.64,-44684.64,-44684.64)


###################################################
### code chunk number 12: abn_v1.0.Rnw:283-285 (eval = FALSE)
###################################################
## marg.f <- fitabn( dag.m=mydag,data.df=mydat,data.dists=mydists,
##                   compute.fixed=TRUE,n.grid=1000)


###################################################
### code chunk number 13: abn_v1.0.Rnw:293-307 (eval = FALSE)
###################################################
## library( Cairo) # Available from CRAN
## CairoPDF( "MargPlots_PigsData.pdf")
## for( i in 1:length(marg.f$marginals)){
##   cat( "processing marginals for node:",
##   nom1 <- names( marg.f$marginals)[i],"\n")
##   cur.node <- marg.f$marginals[i]
##   # Get the marginal for current node, this is a matrix [x,f(x)]
##   cur.node <- cur.node[[1]]
##   # This is always [[1]] for models without random effects
##   for( j in 1:length(cur.node)){
##   cat( "processing parameter:",nom2<- names(cur.node)[j],"\n")
##   cur.param <- cur.node[[j]]
##   plot( cur.param,type="l",main=paste(nom1,":",nom2))}}
## dev.off()


###################################################
### code chunk number 14: abn_v1.0.Rnw:321-325 (eval = FALSE)
###################################################
## marnew <- marg.f$marginals[[1]]
## 
## for( i in 2: length(marg.f$marginals)){
##   marnew <- c(marnew, marg.f$marginals[[i]])}


###################################################
### code chunk number 15: abn_v1.0.Rnw:329-337 (eval = FALSE)
###################################################
## myarea <- rep( NA,length( marnew))
## names( myarea) <- names( marnew)
## for( i in 1:length(marnew)){
##   tmp <- spline(marnew[[i]]) 
##   # Spline just helps make the estimation smoother
##   myarea[i] <- sum(diff(tmp$x)*tmp$y[-1])}
##   # Just width x height of rectangles
## cbind( myarea)


###################################################
### code chunk number 16: abn_v1.0.Rnw:340-356 (eval = FALSE)
###################################################
## library( Cairo)
## mycols <- rep("green",length(marnew))
## mycols[1:2] <- "red"; # PC 
## mycols[3] <- "orange"; # PT 
## mycols[4:5] <- "yellow"; # MS
## mycols[6:7] <- "blue"; # HS
## mycols[8:9] <- "lightblue"; # TAIL 
## mycols[10] <- "mistyrose"; # Abscess 
## mycols[11:12] <- "green"; # Pyaemia, lightcyan 
## mycols[13:14] <- "lavender"; # EPcat
## mycols[15:18] <- "cornsilk"; # PDcat 
## mycols[19:22] <- "brown" ; # plbinary
## CairoPNG( "Area_PigsData.png",pointsize=10,width=720,height=640)
## par( las=2, mar=c(8.1,4.1,4.1,2.1))
## barplot( myarea,ylab="Area under Density",ylim=c(0,2), col=mycols)
## dev.off()


###################################################
### code chunk number 17: abn_v1.0.Rnw:359-388 (eval = FALSE)
###################################################
## print( names(marnew))
## # want to bind all the marginals the same nodes into a matrix
## m <- marnew; 
## # PT -> PC, 1 parent, 2 params;
## PC.p <- cbind( m[["PC|(Intercept)"]],  m[["PC|PT"]])
## # PC --> NO PARENTS, 1 params;
## PT.p <- cbind( m[["PT|(Intercept)"]]) 
## # HS -> MS, 1 parent, 2 params; 
## MS.p <- cbind( m[["MS|(Intercept)"]], m[["MS|HS"]]) 
## # EPcat -> HS, 1 parent, 2 params;
## HS.p <- cbind( m[["HS|(Intercept)"]],  m[["HS|EPcat"]]) 
## # PDcat -> TAIL, 1 parent, 2 params;
## TAIL.p <- cbind( m[["TAIL|(Intercept)"]],  m[["TAIL|PDcat"]]) 
## # Abscess --> NO PARENTS, 1 param;
## Abscess.p <- cbind( m[["Abscess|(Intercept)"]]) 
## # TAIL -> Pyaemia, 1 parent, 2 params;
## Pyaemia.p <- cbind(m[["Pyaemia|(Intercept)"]],m[["Pyaemia|TAIL"]]) 
## # plbinary -> EPcat, 1 parent, 2 params;
## EPcat.p <- cbind( m[["EPcat|(Intercept)"]], m[["EPcat|plbinary"]]) 
## # MS, EPcat, plbinary --> PDcat, 3 parents, 4 params;
## PDcat.p <- cbind( m[["PDcat|(Intercept)"]], m[["PDcat|MS"]], 
##                   m[["PDcat|EPcat"]], m[["PDcat|plbinary"]]) 
## # PC, PT, Abscess --> plbinary, 3 parents, 4 params;
## plbinary.p <- cbind(m[["plbinary|(Intercept)"]],m[["plbinary|PC"]], 
##                     m[["plbinary|PT"]], m[["plbinary|Abscess"]]) 
## 
## dump( c( "PC.p","PT.p","MS.p","HS.p","TAIL.p","Abscess.p", 
##          "Pyaemia.p","EPcat.p","PDcat.p","plbinary.p"), 
##           file=paste('Jags/',"pigs_post_params.R",sep=""))


###################################################
### code chunk number 18: abn_v1.0.Rnw:400-401 (eval = FALSE)
###################################################
## orig.data <- pigs.vienna[,-11] 


###################################################
### code chunk number 19: abn_v1.0.Rnw:405-406 (eval = FALSE)
###################################################
## system("jags jags_pigs_script.R")


###################################################
### code chunk number 20: abn_v1.0.Rnw:410-417 (eval = FALSE)
###################################################
## library( coda)
## boot.data <- read.coda("out1chain1.txt","out1index.txt")
## 
## boot.data <- as.data.frame(boot.data)
## for( j in 1:dim(orig.data)[2]){if(is.factor(orig.data[,j]))
##         { boot.data[,j]<- as.factor(boot.data[,j])
##           levels(boot.data[,j])<- levels(orig.data[,j])}}


###################################################
### code chunk number 21: abn_v1.0.Rnw:420-433 (eval = FALSE)
###################################################
## ban <- matrix( rep(0,dim(orig.data)[2]^2),
##                ncol=dim(orig.data)[2])
## colnames( ban) <- rownames(ban) <- names(orig.data)
## 
## retain <- matrix( rep(0,dim(orig.data)[2]^2),
##                   ncol=dim(orig.data)[2])
## colnames( retain) <- rownames(retain) <- names(orig.data)
## 
## mydists <- list( PC="binomial",  PT="binomial", MS="binomial",   
##                  HS="binomial",  TAIL="binomial", 
##                  Abscess="binomial", Pyaemia="binomial", 
##                  EPcat="binomial", PDcat="binomial", 
##                  plbinary="binomial")


###################################################
### code chunk number 22: abn_v1.0.Rnw:436-437 (eval = FALSE)
###################################################
## max.par <- 3;


###################################################
### code chunk number 23: abn_v1.0.Rnw:440-443 (eval = FALSE)
###################################################
## boot1.cache <- buildscorecache( data.df=boot.data, 
##                data.dists=mydists, max.parents=max.par, 
##                dag.banned=ban, dag.retained=retain) 


###################################################
### code chunk number 24: abn_v1.0.Rnw:446-450 (eval = FALSE)
###################################################
## boot1.mp <- mostprobable( score.cache=boot1.cache)
## 
## datadir <- tempdir()
## save( boot1.mp,file=paste( datadir,"boot1run.RData",sep=''))


###################################################
### code chunk number 25: abn_v1.0.Rnw:513-542 (eval = FALSE)
###################################################
## mydata <- pigs.vienna[,-11] 
## 
## N <- 10240;
## # Write out manually, clearer than using rep() 
## mydag <- matrix(c( 
## # PC  PT  MS  HS TAIL Abscess Pyaemia EPcat PDcat plbinary
##   0,  1,  0,  0,  0,  0,       0,     0,    0,     0, 
##   0,  0,  0,  0,  0,  0,       0,     0,    0,     0, 
##   0,  0,  0,  1,  0,  0,       0,     0,    0,     0, 
##   0,  0,  0,  0,  0,  0,       0,     1,    0,     0, 
##   0,  0,  0,  0,  0,  0,       0,     0,    1,     0, 
##   0,  0,  0,  0,  0,  0,       0,     0,    0,     0, 
##   0,  0,  0,  0,  1,  0,       0,     0,    0,     0, 
##   0,  0,  0,  0,  0,  0,       0,     0,    0,     1, 
##   0,  0,  1,  0,  0,  0,       0,     1,    0,     1, 
##   1,  1,  0,  0,  0,  1,       0,     0,    0,     0),
##   byrow=TRUE,ncol=10)
## colnames(mydag) <- rownames(mydag) <- names(mydata)
## sum( mydag) # 12 arcs, as in original model, Figure 2
## 
## mydists <- list( PC="binomial", PT="binomial", MS="binomial", 
##                  HS="binomial", TAIL="binomial", 
##                  Abscess="binomial", Pyaemia="binomial", 
##                  EPcat="binomial", PDcat="binomial", 
##                  plbinary="binomial")
## 
## # Use fitabn to check mydag is correct (no typos mlik = -44684.64)
## print( fitabn(dag.m=mydag,data.df=mydata,data.dists=mydists)$mlik)
## bestdag <- mydag


###################################################
### code chunk number 26: abn_v1.0.Rnw:552-564 (eval = FALSE)
###################################################
## boot.dags <- list()
## these <- grep("mp10Kboot\\d+.RData", dir()) 
## num <- 1
## for( i in dir()[these]){# Load each file
##   load(i) # Provides dags, a list
##   tmp <- dags[which(unlist(lapply(dags,sum))>0)]
##   # Get valid entries in dags but as a list
##   for( j in 1:length(tmp)){
##   # For each entry copy into boot.dags, and increment counter
##     boot.dags[[num]]<- tmp[[j]]; num <- num+1 }
##   rm( dags)
## }


###################################################
### code chunk number 27: abn_v1.0.Rnw:567-579 (eval = FALSE)
###################################################
## if( FALSE){
##   scores <- rep(0,length(boot.dags))
##   for(i in 1:length(boot.dags)){
##     scores[i] <- fitabn(dag.m=boot.dags[[i]],data.df=mydata, 
##   	  									data.dists=mydists)$mlik
##   }
##   scores.b <- scores[-which(scores< -N)] 
##   orig.score <- fitabn(dag.m=bestdag,data.df=mydata, 
##   										 data.dists=mydists)$mlik
##   plot(density(scores.b,from=min(scores.b),to=max(scores.b)))
##   abline(v=orig.score,lwd=2,col="blue")
## }


###################################################
### code chunk number 28: abn_v1.0.Rnw:583-607 (eval = FALSE)
###################################################
## boot.dags.trim <- boot.dags
## for( i in 1:length(boot.dags)){
##    boot.dags.trim[[i]] <- boot.dags.trim[[i]]*bestdag  }
## 
## arc.freq <- lapply(boot.dags.trim,sum)
## arc.freq <- table(unlist(arc.freq))
## library( Cairo)
## CairoPNG("PigsFreqBootRes.png",pointsize=10,width=720,height=700)
## par(las=1, mar=c(6.1,6.1,4.1,2.1))
## barplot( arc.freq,ylab="",xlab="",col="skyblue",
##          names.arg=names(arc.freq), ylim=c(0,2500))
## par( las=1)
## mtext( "Number of arcs in bootstrap DAG",1,line=3,cex=1.5)
## par( las=3)
## mtext( "Frequency out of 10240",2,line=4,cex=1.5)
## dev.off()
## 
## total.dag <- matrix(rep(0,dim(bestdag)[2]^2),ncol=dim(bestdag)[2])
## colnames(total.dag) <- rownames(total.dag)<- colnames(bestdag)
## # Get the support for each arc, total.dag:
## for( i in 1:length(boot.dags)){
##   if(sum(boot.dags[[i]])>0){total.dag <- total.dag+boot.dags[[i]]}}  
## total.dag <- total.dag*bestdag # We only want arcs in the best DAG
## total.dag 


###################################################
### code chunk number 29: abn_v1.0.Rnw:613-616 (eval = FALSE)
###################################################
## f <- function(val,limit){ if(val<limit){return(0)} 
##                           else {return(1)}}
## bestdag.trim <- apply( total.dag,c(1,2),FUN=f,limit=N/2)


###################################################
### code chunk number 30: abn_v1.0.Rnw:620-640 (eval = FALSE)
###################################################
## bestdag.trim.nodir <- bestdag
## bestdag.trim.nodir[,] <- 0 # Set zero
## child <- NULL; parent <- NULL
## for( i in 1:dim(total.dag)[1]){
##   for( j in 1:dim(total.dag)[2]){
##      if(i>j){ # Get most supported direction
##               if(total.dag[i,j]>total.dag[j,i]){
##                 m.i <- i; m.j <- j;} 
##               else {m.i <- j; m.j <- i;}
##               # Does arc quality - exceed threshold of support
##               if(total.dag[i,j]+total.dag[j,i]>N/2){
##               # We want this as more than 5000
##               bestdag.trim.nodir[m.i,m.j] <- 1}}}}
## 
## tographviz( dag.m=bestdag.trim,data.df=mydata,
##             data.dists=mydists, outfile="postbootpigs.dot")
## system( "dot -Tpdf -o postbootpigs.pdf postbootpigs.dot")
## system( "evince postbootpigs.pdf&")
## 
## save( bestdag.trim,file=paste("bestdagpigs_trim.RData",sep=''))


###################################################
### code chunk number 31: abn_v1.0.Rnw:656-694 (eval = FALSE)
###################################################
## mydata <- pigs.vienna[,-11] 
## 
## mydists <- list( PC="binomial", PT="binomial", MS="binomial",  
##                  HS="binomial", TAIL="binomial",  
##                  Abscess="binomial", Pyaemia="binomial",   
##                  EPcat="binomial", PDcat="binomial",
##                  plbinary="binomial")
## 
## marg.f <- fitabn(dag.m=bestdag.trim,data.df=mydata,
##                data.dists=mydists,compute.fixed=TRUE,
##                n.grid=1000)
## 
## library( Cairo)
## CairoPDF("Pigs_PostBootPlots.pdf")
## for( i in 1:length(marg.f$marginals)){
##   cat( "processing marginals for node:",
##   nom1 <- names(marg.f$marginals)[i],"\n")
##   cur.node <- marg.f$marginals[i]
##   # Get marginal for current node, this is a matrix [x,f(x)]
##   cur.node <- cur.node[[1]]
##   # This is always [[1]] for models without random effects
##    for( j in 1:length(cur.node)){
##    cat("processing parameter:",
##    nom2 <- names(cur.node)[j],"\n")
##    cur.param <- cur.node[[j]]
##    plot( cur.param,type="l",main=paste(nom1,":",nom2))}}
## dev.off()
## 
## marnew <- marg.f$marginals[[1]]
## for(i in 2: length(marg.f$marginals)){
##   marnew <- c(marnew, marg.f$marginals[[i]])}
## 
## margs <- marnew;
## mymat <- matrix(rep(NA,length(margs)*3),ncol=3)
## rownames(mymat) <- names(margs)
## colnames(mymat) <- c("2.5%","50%","97.5%")
## ignore.me <- union(grep("\\(Int",names(margs)),
##                    grep("prec",names(margs)))


###################################################
### code chunk number 32: abn_v1.0.Rnw:697-710 (eval = FALSE)
###################################################
## comment <- rep("",length(margs))
## for(i in 1:length(margs)){
##   tmp <- margs[[i]]
##   tmp2 <- cumsum(tmp[,2])/sum(tmp[,2])
##   mymat[i,] <- c(tmp[which(tmp2>0.025)[1]-1,1], 
##                  tmp[which(tmp2>0.5)[1],1],
##                  tmp[which(tmp2>0.975)[1],1])
##   myvec <- mymat[i,]
##   if( !(i%in%ignore.me) &&  (myvec[1]<0 && myvec[3]>0)){
##       comment[i] <- "not sig. at 5%"} 
##   # Truncate for printing
##   mymat[i,] <- as.numeric(formatC(mymat[i,],digits=3,format="f"))}
## cbind( mymat)


###################################################
### code chunk number 33: abn_v1.0.Rnw:792-798
###################################################
library( abn)
bin.nodes <- c( 1,3,4,6,9,10,11,12,15,18,19,20,21,26,27,28,32) 
var33.cat <- var33[,bin.nodes] #Categorical nodes only

dag33 <- matrix( 0, 17, 17) 
colnames( dag33) <- rownames( dag33) <- names( var33.cat)#Set names


###################################################
### code chunk number 34: abn_v1.0.Rnw:801-802
###################################################
dag33["v11","v12"] <- 0; dag33["v11","v10"]<- 0; dag33["v4","v3"]<- 0;


###################################################
### code chunk number 35: abn_v1.0.Rnw:805-813
###################################################
mydists.cat <- list( v1 ="binomial", v3 = "binomial",
       v4 = "binomial",  v6 = "binomial",  v9 = "binomial",
      v10 = "binomial", v11 = "binomial", v12 = "binomial",
      v15 = "binomial", v18 = "binomial", v19 = "binomial",
      v20 = "binomial", v21 = "binomial", v26 = "binomial",
      v27 = "binomial", v28 = "binomial", v32 = "binomial")
ind.mod.cat <- fitabn( data.df=var33.cat, dag.m=dag33,
                      data.dists=mydists.cat, verbose=FALSE) 


###################################################
### code chunk number 36: abn_v1.0.Rnw:819-820
###################################################
ind.mod.cat$mlik


###################################################
### code chunk number 37: abn_v1.0.Rnw:825-830
###################################################
dag33["v11","v12"] <- 1;
dag33["v11","v10"] <- 1;
dag33["v4","v3"] <- 1;
dep.mod.cat <- fitabn( data.df=var33.cat, dag.m=dag33,
                       data.dists=mydists.cat, verbose=FALSE)


###################################################
### code chunk number 38: abn_v1.0.Rnw:833-834
###################################################
dep.mod.cat$mlik


###################################################
### code chunk number 39: abn_v1.0.Rnw:837-839
###################################################
tographviz( dag=dag33, data.df=var33.cat, data.dists=mydists.cat,
           outfile="mydagcat.dot", directed=TRUE) #Create file


###################################################
### code chunk number 40: abn_v1.0.Rnw:853-856
###################################################
var33.cts <- var33[,-bin.nodes] 
dag33 <- matrix( 0, 16, 16)
colnames( dag33) <- rownames( dag33) <- names( var33.cts)


###################################################
### code chunk number 41: abn_v1.0.Rnw:859-865
###################################################
mydists.cts <- list(       v2 = "gaussian",  v5 = "gaussian",
         v7 = "gaussian",  v8 = "gaussian", v13 = "gaussian", 
        v14 = "gaussian", v16 = "gaussian", v17 = "gaussian",
        v22 = "gaussian", v23 = "gaussian", v24 = "gaussian",
        v25 = "gaussian", v29 = "gaussian", v30 = "gaussian",
        v31 = "gaussian", v33 = "gaussian")


###################################################
### code chunk number 42: abn_v1.0.Rnw:868-870
###################################################
ind.mod.cts <- fitabn( data.df=var33.cts, dag.m=dag33,
               data.dists=mydists.cts, verbose=FALSE)


###################################################
### code chunk number 43: abn_v1.0.Rnw:876-877
###################################################
ind.mod.cts$mlik


###################################################
### code chunk number 44: abn_v1.0.Rnw:882-887
###################################################
dag33["v33","v31"] <- 1;
dag33["v24","v23"] <- 1;
dag33["v14","v13"] <- 1;
dep.mod.cts <- fitabn( data.df=var33.cts, dag.m=dag33,
               data.dists=mydists.cts, verbose=FALSE)


###################################################
### code chunk number 45: abn_v1.0.Rnw:890-894
###################################################
dep.mod.cts$mlik

tographviz( dag=dag33, data.df=var33.cts, data.dists=mydists.cts,
            outfile="mydagcts.dot", directed=TRUE) #Create file


###################################################
### code chunk number 46: abn_v1.0.Rnw:907-909
###################################################
dag33 <- matrix( 0, 33, 33) 
colnames( dag33) <- rownames( dag33)<- names( var33)


###################################################
### code chunk number 47: abn_v1.0.Rnw:912-924
###################################################
mydists.mix <- list(      v1 = "binomial",  v2 = "gaussian",
        v3 = "binomial",  v4 = "binomial",  v5 = "gaussian",
        v6 = "binomial",  v7 = "gaussian",  v8 = "gaussian",
        v9 = "binomial", v10 = "binomial", v11 = "binomial",
       v12 = "binomial", v13 = "gaussian", v14 = "gaussian",
       v15 = "binomial", v16 = "gaussian", v17 = "gaussian",
       v18 = "binomial", v19 = "binomial", v20 = "binomial",
       v21 = "binomial", v22 = "gaussian", v23 = "gaussian",
       v24 = "gaussian", v25 = "gaussian", v26 = "binomial",
       v27 = "binomial", v28 = "binomial", v29 = "gaussian",
       v30 = "gaussian", v31 = "gaussian", v32 = "binomial",
       v33 = "gaussian")


###################################################
### code chunk number 48: abn_v1.0.Rnw:927-929
###################################################
ind.mod <- fitabn( data.df=var33, dag.m=dag33,
          data.dists=mydists.mix, verbose=FALSE)


###################################################
### code chunk number 49: abn_v1.0.Rnw:932-933
###################################################
ind.mod$mlik


###################################################
### code chunk number 50: abn_v1.0.Rnw:936-962
###################################################
dag33[2,1] <- 1;
dag33[4,3] <- 1;
dag33[6,4] <- 1; dag33[6,7] <- 1;
dag33[5,6] <- 1;
dag33[7,8] <- 1;  
dag33[8,9] <- 1;
dag33[9,10] <- 1;
dag33[11,10] <- 1; dag33[11,12] <- 1; dag33[11,19] <- 1;
dag33[14,13] <- 1;
dag33[17,16] <- 1;dag33[17,20] <- 1;
dag33[15,14] <- 1; dag33[15,21] <- 1;
dag33[18,20] <- 1;
dag33[19,20] <- 1;
dag33[21,20] <- 1;
dag33[22,21] <- 1;
dag33[23,21] <- 1;
dag33[24,23] <- 1;
dag33[25,23] <- 1; dag33[25,26] <- 1;
dag33[26,20] <- 1;
dag33[33,31] <- 1;
dag33[33,31] <- 1;
dag33[32,21] <- 1; dag33[32,31] <- 1; dag33[32,29] <- 1;    
dag33[30,29] <- 1;
dag33[28,27] <- 1; dag33[28,29] <- 1; dag33[28,31] <- 1;       
dep.mod <- fitabn( data.df=var33, dag.m=dag33, 
           data.dists=mydists.mix, verbose=FALSE)


###################################################
### code chunk number 51: abn_v1.0.Rnw:966-969
###################################################
dep.mod$mlik
tographviz( dag=dag33, data.df=var33, data.dists=mydists.mix,
            outfile="mydag_all.dot", directed=TRUE)#Create file


###################################################
### code chunk number 52: abn_v1.0.Rnw:991-995
###################################################
bin.nodes <- c( 1,3,4,6,9,10,11,12,15,18,19,20,21,26,27,28,32) 
var33.cat <- var33[,bin.nodes] #Categorical nodes only
dag33 <- matrix( 0, 17, 17) 
colnames(dag33) <- rownames(dag33) <- names(var33.cat)#Set names


###################################################
### code chunk number 53: abn_v1.0.Rnw:998-1002
###################################################
banned.cat <- matrix( 0, 17, 17)
colnames(banned.cat) <- rownames(banned.cat) <- names(var33.cat)
retain.cat <- matrix( 0, 17, 17)
colnames(retain.cat) <- rownames(retain.cat) <- names(var33.cat)


###################################################
### code chunk number 54: abn_v1.0.Rnw:1005-1011
###################################################
mydists.cat <- list(     v1 = "binomial",  v3 = "binomial",
       v4 = "binomial",  v6 = "binomial",  v9 = "binomial",
      v10 = "binomial", v11 = "binomial", v12 = "binomial",
      v15 = "binomial", v18 = "binomial", v19 = "binomial",
      v20 = "binomial", v21 = "binomial", v26 = "binomial",
      v27 = "binomial", v28 = "binomial", v32 = "binomial")


###################################################
### code chunk number 55: abn_v1.0.Rnw:1014-1017
###################################################
mycache.cat <- buildscorecache( data.df=var33.cat,
              data.dists=mydists.cat, dag.banned=banned.cat, 
              dag.retained=retain.cat, max.parents=1)


###################################################
### code chunk number 56: abn_v1.0.Rnw:1028-1031 (eval = FALSE)
###################################################
## heur.res.cat <- search.hillclimber( score.cache=mycache.cat,
##                 num.searches=1, seed=0, verbose=FALSE,
##                 trace=FALSE, timing.on=FALSE)


###################################################
### code chunk number 57: abn_v1.0.Rnw:1038-1045
###################################################
var33.cts <- var33[,-bin.nodes] #Drop categorical nodes
dag33 <- matrix( 0, 16, 16) 
colnames(dag33) <- rownames(dag33) <- names(var33.cts) #Set names
banned.cts <- matrix( 0, 16, 16)
colnames(banned.cts) <- rownames(banned.cts) <- names(var33.cts)
retain.cts <- matrix( 0, 16, 16)
colnames(retain.cts) <- rownames(retain.cts) <- names(var33.cts)


###################################################
### code chunk number 58: abn_v1.0.Rnw:1048-1054
###################################################
mydists.cts <- list(       v2 = "gaussian",  v5 = "gaussian",
         v7 = "gaussian",  v8 = "gaussian", v13 = "gaussian", 
        v14 = "gaussian", v16 = "gaussian", v17 = "gaussian",
        v22 = "gaussian", v23 = "gaussian", v24 = "gaussian",
        v25 = "gaussian", v29 = "gaussian", v30 = "gaussian",
        v31 = "gaussian", v33 = "gaussian")


###################################################
### code chunk number 59: abn_v1.0.Rnw:1058-1061
###################################################
mycache.cts<- buildscorecache( data.df=var33.cts,
             data.dists=mydists.cts, dag.banned=banned.cts,
             dag.retained=retain.cts, max.parents=1)


###################################################
### code chunk number 60: abn_v1.0.Rnw:1064-1067 (eval = FALSE)
###################################################
## heur.res.cts<- search.hillclimber( score.cache=mycache.cts,
##               num.searches=1, seed=0, verbose=FALSE,
##               trace=FALSE, timing.on=FALSE)


###################################################
### code chunk number 61: abn_v1.0.Rnw:1075-1077
###################################################
dag33 <- matrix( 0, 33, 33) 
colnames(dag33) <- rownames(dag33) <- names(var33)#Set names


###################################################
### code chunk number 62: abn_v1.0.Rnw:1080-1084
###################################################
banned.mix <- matrix( 0, 33, 33)
colnames(banned.mix) <- rownames(banned.mix) <- names(var33)
retain.mix<- matrix( 0, 33, 33)
colnames(retain.mix) <- rownames(retain.mix) <- names(var33)


###################################################
### code chunk number 63: abn_v1.0.Rnw:1087-1099
###################################################
mydists.mix <- list(       v1 = "binomial",  v2 = "gaussian",
         v3 = "binomial",  v4 = "binomial",  v5 = "gaussian",
         v6 = "binomial",  v7 = "gaussian",  v8 = "gaussian",
         v9 = "binomial", v10 = "binomial", v11 = "binomial",
        v12 = "binomial", v13 = "gaussian", v14 = "gaussian",
        v15 = "binomial", v16 = "gaussian", v17 = "gaussian",
        v18 = "binomial", v19 = "binomial", v20 = "binomial",
        v21 = "binomial", v22 = "gaussian", v23 = "gaussian",
        v24 = "gaussian", v25 = "gaussian", v26 = "binomial",
        v27 = "binomial", v28 = "binomial", v29 = "gaussian",
        v30 = "gaussian", v31 = "gaussian", v32 = "binomial",
        v33 = "gaussian")


###################################################
### code chunk number 64: abn_v1.0.Rnw:1103-1106
###################################################
mycache.mix <- buildscorecache( data.df=var33,
             data.dists=mydists.mix, dag.banned=banned.mix,
             dag.retained=retain.mix, max.parents=1)


###################################################
### code chunk number 65: abn_v1.0.Rnw:1109-1112
###################################################
heur.res.mix <- search.hillclimber( score.cache=mycache.mix,
              num.searches=1, seed=0, verbose=FALSE,
              trace=FALSE, timing.on=FALSE)


###################################################
### code chunk number 66: abn_v1.0.Rnw:1125-1127
###################################################
dag33 <- matrix( 0, 33, 33) 
colnames(dag33) <- rownames(dag33) <- names(var33) #Set names


###################################################
### code chunk number 67: abn_v1.0.Rnw:1130-1134
###################################################
banned.mix <- matrix( 0, 33, 33)
colnames(banned.mix)<- rownames(banned.mix)<- names(var33)
retain.mix <- matrix( 0, 33, 33)
colnames(retain.mix) <- rownames(retain.mix)<- names(var33)


###################################################
### code chunk number 68: abn_v1.0.Rnw:1137-1150
###################################################
mydists.mix <- list(      v1 = "binomial",  v2 = "gaussian",
        v3 = "binomial",  v4 = "binomial",  v5 = "gaussian",
        v6 = "binomial" , v7 = "gaussian",  v8 = "gaussian",
        v9 = "binomial", v10 = "binomial", v11 = "binomial",
       v12 = "binomial", v13 = "gaussian", v14 = "gaussian",
       v15 = "binomial", v16 = "gaussian", v17 = "gaussian",
       v18 = "binomial", v19 = "binomial", v20 = "binomial",
       v21 = "binomial", v22 = "gaussian", v23 = "gaussian",
       v24 = "gaussian", v25 = "gaussian", v26 = "binomial",
       v27 = "binomial", v28 = "binomial", v29 = "gaussian",
       v30 = "gaussian", v31 = "gaussian", v32 = "binomial",
       v33 = "gaussian")
n.searches <-  10 


###################################################
### code chunk number 69: abn_v1.0.Rnw:1154-1155
###################################################
max.par <- 1


###################################################
### code chunk number 70: abn_v1.0.Rnw:1159-1161
###################################################
mycache.mix <- buildscorecache( data.df=var33, data.dists=mydists.mix,
dag.banned=banned.mix, dag.retained=retain.mix, max.parents=max.par)


###################################################
### code chunk number 71: abn_v1.0.Rnw:1165-1168
###################################################
myres.mlp <- search.hillclimber(score.cache=mycache.mix,
             num.searches=n.searches, seed=0, verbose=FALSE,
           trace=FALSE, timing.on=FALSE)


###################################################
### code chunk number 72: abn_v1.0.Rnw:1180-1182
###################################################
tographviz( dag= myres.mlp$consensus, data.df=var33,
data.dists=mydists.mix, outfile="dagcon.dot") #Create file


