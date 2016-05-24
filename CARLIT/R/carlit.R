carlit <-
function(X, List, EQR_r){
	
	EQ <- rep(NA, length(unique(paste(X[,1], X[,2], X[,3], sep="_"))))
L <- rep(NA, length(unique(paste(X[,1], X[,2], X[,3], sep="_"))))

EQssi <- data.frame(Sit=rep(NA, length(EQ)) , Situation=as.character(unique(paste(X[,1], X[,2], X[,3], sep="_"))), eq=rep(NA, length(unique(paste(X[,1], X[,2], X[,3], sep="_")))), l=L)

for (i in 1:length(unique(paste(X[,1], X[,2], X[,3], sep="_")))){
	
	Site <- X[paste(X[,1], X[,2], X[,3], sep="_")==as.character(unique(paste(X[,1], X[,2], X[,3], sep="_")))[i],]

	Num_par <- rep(NA, length(unique(Site[,5])))
	L_par <- rep(NA, 1)
	
		for (j in 1:length(unique(List[,1]))){
	
				Num_par[j] <- sum(as.integer(paste(Site[as.character(Site[,5])==as.character(List[j,1]),4])))*as.integer(paste(List[j,2]))
				
				Lpar <- sum(as.integer(paste(Site[,4])))
			   
	
											}
				L[i] <- Lpar
				Num <- rep(NA, length(unique(paste(X[,1], X[,2], X[,3], sep="_"))))
				Num[i] <- sum(Num_par)
				Den <- rep(NA, length(unique(paste(X[,1], X[,2], X[,3], sep="_"))))
				Den[i] <- sum(as.integer(paste(Site[,4])))
				
				EQ[i] <- as.numeric(Num[i])/as.numeric(Den[i])
				
				EQssi[i,1] <- as.character(unique(paste(Site[,1])))
				
				EQssi[i,3] <- as.numeric(EQ[i])
				
				EQssi[i,4] <- L[i]
	
									}
rm(EQ)
colnames(EQssi) <- c("Site", "GRS", "EQ.ssi", "l")

EQR <- rep(NA, length(EQssi[,1]))

for (l in 1:length(EQssi[,1])){
	
	Sub <- EQssi[as.character(EQssi[,1])==as.character(unique(EQssi[,1])[l]), ]
	
	EQR_par_A <- numeric()
	EQR_par_B <- numeric()
	EQR_par <- numeric()
	
		for (m in 1:length(unique(Sub[,2]))){
		
		EQR_par_A[   m   ] <- ((Sub[  m  ,3]) / (EQR_r[as.character(paste(rep(unique(Sub[,1]), each=length(EQR_r[,1])),paste(EQR_r[,2], EQR_r[,3], sep='_'), sep='_')) == unique(as.character(paste(Sub[,2]))[   m   ])    , 4]))
		
		
			if (paste(EQR_par_A[  m  ])>1) {
			EQR_par_A[  m  ] <- 1
			} else {
			EQR_par_A[  m  ] <- as.numeric(paste(EQR_par_A[  m  ]))
			}
			
		EQR_par_B[m] <-  Sub[  m  ,4]
		EQR_par[m] <- EQR_par_A[  m  ]* EQR_par_B[  m  ]
		
		
		}
	
	EQR <- c(EQR, rep(unique(sum(EQR_par)/sum(Sub[,4])), length(Sub[,1])))  
		EQR <- EQR[!is.na(EQR)]
}

EQssi[, 5] <- EQR
colnames(EQssi)[5] <- "EQR"
EQssi
}
