`DemoBoot` <-
function (stagedat = NULL, fruitdat = NULL, seeddat = NULL, n = 1) 
{
    nL <- nrow(stagedat)
    nF <- nrow(fruitdat)
    nS <- nrow(seeddat)
    out <- lapply(1:n, function(i) {
        stageR <- stagedat[sample(1:nL, nL, replace = TRUE), 
            ]
        fruitR <- fruitdat[sample(1:nF, nF, replace = TRUE), 
            ]
        seedR <- as.data.frame(seeddat[sample(1:nS, nS, replace = TRUE), 
            ])
        matR <- ProjMat(stagedat = stageR, fruitdat = fruitR, 
            seeddat = seedR)
        DemoInfo(matR)
    })
    return(out)
}
