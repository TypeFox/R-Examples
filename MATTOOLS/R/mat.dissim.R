 
mat.dissim <- function(inFossil, inModern, llMod=c(), modTaxa=c(), llFoss=c(), fosTaxa=c(), numAnalogs = 1, counts=T,
        sitenames = 1:length(inFossil[, 1]), dist.method = "euclidean")
{

#Get arguments for use in final reconstruction
        modchar=deparse(substitute(inModern))
        fosschar=deparse(substitute(inFossil))

        if (length(modTaxa) != length(fosTaxa)) stop("Number of taxa in modern sample does not equal number of taxa in fossil sample")
        if (counts) {

                inModern[,modTaxa]=inModern[,modTaxa]/rowSums(inModern[,modTaxa])
                inFossil[,fosTaxa]=inFossil[,fosTaxa]/rowSums(inFossil[,fosTaxa])

        }

        nColsMatrix = length(inFossil[, 1])
        LocMinRow = matrix(NA, nrow = numAnalogs, ncol = nColsMatrix)
        PosMinRow = matrix(NA, nrow = numAnalogs, ncol = nColsMatrix)
        DistMinRow <- matrix(NA, nrow = numAnalogs, ncol = nColsMatrix)
        DirMinRow <- matrix(NA, nrow = numAnalogs, ncol = nColsMatrix)
        CompxMinRow <- matrix(NA, nrow = numAnalogs, ncol = nColsMatrix)
        CompyMinRow <- matrix(NA, nrow = numAnalogs, ncol = nColsMatrix)
        colNames = vector("character")
        modMatrix = as.matrix(inModern[, modTaxa])
        fossilMatrix = sqrt(as.matrix(inFossil[, fosTaxa]))
        fossilLongLatm=as.matrix(inFossil[,llFoss])
        modernLongLatm=as.matrix(inModern[, llMod])
        dimnames(fossilLongLatm)=NULL
        dimnames(modernLongLatm)=NULL
        dimnames(modMatrix) = NULL
        dimnames(fossilMatrix) = NULL
        modMatrix = sqrt(t(modMatrix))
        for(i in 1:length(inFossil[, 1])) {
                currSpectrum = fossilMatrix[i,  ]
                sqdistVec = (currSpectrum - modMatrix)
                sqdistVec= sqdistVec*sqdistVec
                x = colSums(sqdistVec)
                y = rank(x,ties.method="first")
                ysubset = y <= numAnalogs
                xysubset = x[ysubset]
                zorder = order(xysubset)
                x = sort(xysubset)
                tseq=1:length(y)
                topn = tseq[ysubset][zorder]

                if(dist.method == "spherical") {
                        DistMinRow[, i]  <- apply(rbind(modernLongLatm[topn, ]), 1, great.circle.distance.f, fossilLongLatm[i, ])
                        DirMinRow[, i] <- rep(NA, numAnalogs)#<- apply(rbind(modernLongLatm[topn, ]), 1, spherical.direction.f, fossilLongLatm[i, ])
                        #CompxMinRow[, i] <- rep(NA, numAnalogs)
                        #CompyMinRow[, i] <- rep(NA, numAnalogs)
                }
                else if(dist.method == "euclidean") {
                        #ddiff=fossilLongLatm[i, ]-t(modernLongLatm[topn, ])#sqrt(colSums(ddiff*ddiff)) #
                        DistMinRow[, i] <- apply(rbind(modernLongLatm[topn, ]), 1, euclidean.distance.f, fossilLongLatm[i, ])
                        DirMinRow[, i] <- apply(rbind(modernLongLatm[topn, ]), 1, euclidean.direction.f, fossilLongLatm[i, ])
                        CompxMinRow[, i] <- apply(rbind(modernLongLatm[topn, ]), 1, euclidean.compx.f, fossilLongLatm[i, ])
                        CompyMinRow[, i] <- apply(rbind(modernLongLatm[topn, ]), 1, euclidean.compy.f, fossilLongLatm[i, ])
                }

                LocMinRow[, i] = x
                PosMinRow[, i] = topn
        }

        colNames = sitenames
        dimnames(LocMinRow) = list(NULL, colNames)
        dimnames(PosMinRow) = list(NULL, colNames)
        dimnames(DistMinRow) = list(NULL, colNames)
        dimnames(DirMinRow) = list(NULL, colNames)

        

        return(list(
        x = inFossil[, llFoss[1]], 
        y = inFossil[, llFoss[2]], 
        sqdist =  LocMinRow, 
        position = PosMinRow, 
        distance = DistMinRow, 
        direction = DirMinRow, 
        xcomponent = CompxMinRow, 
        ycomponent = CompyMinRow,
        inModern=modchar,
        inFossil=fosschar,
        llmod=llMod,
        modTaxa=modTaxa,
        counts=counts
        
        ))

  }



