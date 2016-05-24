trans.ASCII <-
function (path, mz) 
{
flagmz=FALSE   
  
	if(missing(mz)||mz==""){
		cat("WARNING, Stop! mz values needed \n")
		flagmz=TRUE
	}
	if(flagmz==FALSE){
	#Rprof()
	 if(missing(path)||path==""){
		require("tcltk")
		path=tclvalue(tkchooseDirectory(title="Please, select your directory containing raw ASCII files"))
	}
    st <- strsplit(date(), " ")[[1]]
    stBis <- strsplit(st[4], ":")[[1]]
    Hour <- paste(stBis[1], stBis[2], stBis[3], sep = "-")
    Date <- paste(st[1], st[2], st[3], Hour, sep = "_")
    PathDate <- paste("output", "_",Date, sep = "")
    dir.create(PathDate)
    L_analyses <- dir(path)
    print(L_analyses)
    for (N_analyses in 1:length(L_analyses)) {
	#change file.path instead of paste(path..
       # an_tempo <- scan(paste(paste(path, "/", sep = ""), L_analyses[N_analyses], 
       #     sep = ""), blank.lines.skip = FALSE)
			an_tempo <- scan(file.path(path, L_analyses[N_analyses]), blank.lines.skip = FALSE)
        chrom <- c(an_tempo[1], an_tempo[2])
        v <- vector()
        v[1] <- 4
        k <- 1
        for (i in 1:length(an_tempo)) {
            if (is.na(an_tempo[i]) == TRUE) 
                if (is.na(an_tempo[i + 1]) == FALSE) 
                  if (is.na(an_tempo[i + 3]) == TRUE) 
                    if (is.na(an_tempo[i + 6]) == FALSE) {
                      chrom <- rbind(chrom, c(an_tempo[i + 1], 
                        an_tempo[i + 2]))
                      k <- k + 1
                      v[k] <- i + 1
                    }
        }
        
        On<-chrom[,2]!=0
        chrom<-chrom[On,]
        v<-v[On]
        MS <- list()
        l <- 0
        
        for (nu in 1:(length(v)-1)) 
        {
            MS[[nu]] <- c(an_tempo[v[nu] + 3], an_tempo[v[nu] + 
                4])
            lo <- (v[nu + 1] - 2) - (v[nu] + 3) + 1
            Se <- seq(from = 2, to = (lo - 2), by = 2)
            if (lo > 2) 
                for (k in Se) MS[[nu]] <- rbind(MS[[nu]], c(an_tempo[v[nu] + 
                  3 + k], an_tempo[v[nu] + k + 4]))
        }
        
        MS <- MS[4:length(MS)]
        for (long in 1:length(MS)) {
            if (is.matrix(MS[[long]]) == TRUE) {
                MS[[long]][, 1] <- round(MS[[long]][, 1])
            }
            else {
                MS[[long]][1] <- round(MS[[long]][1])
            }
        }
        
        an <- matrix(nrow = length(MS), ncol = length(mz))
        colnames(an) <- mz
        for (i in 1:length(MS)) {
            if (is.matrix(MS[[i]]) == TRUE) {
                Ma <- match(mz, MS[[i]][, 1])
                for (p in 1:length(Ma)) if (is.na(Ma[p]) == FALSE) 
                  an[i, p] <- MS[[i]][Ma[p], 2]
                else (an[i, p] <- 0)
            }
            else {
                for (p in 1:length(mz)) an[i, p] <- 0
            }
        }
        
        chrom <- chrom[4:(dim(chrom)[1]-1), ]
        an <- cbind(chrom, as.data.frame(an))
		#change here file.path instead of paste(...
        write.table(an, file = file.path(PathDate, L_analyses[N_analyses]), row.names = FALSE)
    }
    #Rprof(NULL)
    #summaryRprof(filename = "Rprof.out")$sampling.time
	print(paste("Your working directory is :", eval(getwd())))
	print(paste("Transformed data files are in the folder:", PathDate, cat("\n")))
}
}
