##**** need some more examples/test cases
##**** add standard, grid versions

##**** Redo with 21 cases, breaking face ambiguity by always cutting
##**** off high vertices.
##**** Allow processing of one slice at a time (maybe even for multiple
##**** contours?)
##**** Need more complete documentation/commenting

PreProcessing <- local({
    explode <- function(x)
        floor(((x - 1) %% 2^(1:8))/2^(0:7))

    BasicRotation <-
        matrix(c(1,2,3,4,5,6,7,8,5,6,2,1,8,7,3,4,8,7,6,5,4,3,2,1,
                 4,3,7,8,1,2,6,5,2,6,7,3,1,5,8,4,6,5,8,7,2,1,4,3,
                 5,1,4,8,6,2,3,7,4,1,2,3,8,5,6,7,3,4,1,2,7,8,5,6,
                 2,3,4,1,6,7,8,5,6,7,3,2,5,8,4,1,7,8,4,3,6,5,1,2,
                 8,5,1,4,7,6,2,3,7,3,2,6,8,4,1,5,4,8,5,1,3,7,6,2,
                 3,2,6,7,4,1,5,8,2,1,5,6,3,4,8,7,1,4,8,5,2,3,7,6,
                 1,5,6,2,4,8,7,3,5,8,7,6,1,4,3,2,8,4,3,7,5,1,2,6,
                 3,7,8,4,2,6,5,1,7,6,5,8,3,2,1,4,6,2,1,5,7,3,4,8),
               ncol=8, byrow=TRUE)

    CaseRotation <-
        matrix(c(1,24,2,19,2,17,3,17,2,24,4,24,3,24,6,10,2,15,3,19,
                 4,17,6,9,3,9,6,8,6,1,9,23,2,20,3,18,4,7,6,16,5,24,
                 7,5,7,24,12,9,4,20,6,22,8,24,10,24,7,9,15,24,13,20,
                 6,20,2,21,4,6,3,16,6,4,4,16,8,23,6,14,10,23,5,21,
                 7,10,7,16,15,9,7,2,13,8,12,23,6,6,3,6,6,17,6,18,
                 9,18,7,4,13,17,15,18,6,13,7,6,12,16,13,18,6,2,11,24,
                 7,3,7,12,3,12,2,23,5,23,4,23,7,1,3,14,7,14,6,21,
                 15,23,4,15,7,19,8,19,13,23,6,11,12,17,10,19,6,23,4,12,
                 7,18,8,22,13,16,7,13,11,23,13,21,7,15,8,21,13,22,14,24,
                 8,15,13,11,7,7,8,12,4,22,3,23,7,23,6,24,12,18,6,7,
                 13,19,9,24,6,19,7,21,11,18,13,24,7,20,15,16,7,22,6,15,
                 3,22,6,3,15,17,10,22,6,12,12,24,7,11,6,5,3,15,13,10,
                 7,8,8,20,4,9,7,17,5,22,4,18,2,22,2,22,4,18,5,22,
                 7,17,4,9,8,20,7,8,13,10,3,15,6,5,7,11,12,24,6,12,
                 10,22,15,17,6,3,3,22,6,15,7,22,15,16,7,20,13,24,11,18,
                 7,21,6,19,9,24,13,19,6,7,12,18,6,24,7,23,3,23,4,22,
                 8,12,7,7,13,11,8,15,14,24,13,22,8,21,7,15,13,21,11,23,
                 7,13,13,16,8,22,7,18,4,12,6,23,10,19,12,17,6,11,13,23,
                 8,19,7,19,4,15,15,23,6,21,7,14,3,14,7,1,4,23,5,23,
                 2,23,3,12,7,12,7,3,11,24,6,2,13,18,12,16,7,6,6,13,
                 15,18,13,17,7,4,9,18,6,18,6,17,3,6,6,6,12,23,13,8,
                 7,2,15,9,7,16,7,10,5,21,10,23,6,14,8,23,4,16,6,4,
                 3,16,4,6,2,21,6,20,13,20,15,24,7,9,10,24,8,24,6,22,
                 4,20,12,9,7,24,7,5,5,24,6,16,4,7,3,18,2,20,9,23,
                 6,1,6,8,3,9,6,9,4,17,3,19,2,15,6,10,3,24,4,24,
                 2,24,3,17,2,17,2,19,1,24),
               ncol=2,byrow=TRUE)

    CaseRotationFlip <-
        cbind(CaseRotation,
              c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,
                1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,-1,-1,1,1,
                1,1,1,1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,1,
                -1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,1,1,1,1,
                1,-1,1,-1,1,1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,-1,-1,-1,-1,-1,
                -1,-1,
                -1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,1,1,1,1,-1,
                -1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,-1,
                -1,-1,
                -1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,
                -1,-1,
                1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,
                -1,-1,1,
                1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,
                -1,-1,-1,
                -1,-1,-1,-1,-1,-1))
    EdgePoints <-
        matrix(c(1,1,2,2,2,3,3,3,4,4,4,1,5,5,6,6,6,7,
                 7,7,8,8,8,5,9,1,5,10,2,6,11,3,7,12,4,8,13,9,9),
               ncol=3, byrow=TRUE)
    BasicEdges <- list(c(1,4,9),
                       c(2,4,9,10),
                       c(1,4,5,6,9,10),
                       c(1,4,6,7,9,11),
                       c(1,4,10,11,12),
                       c(2,4,6,7,9,10,11),
                       c(1,2,5,6,7,8,9,10,11,13),
                       c(9,10,11,12),
                       c(1,2,7,8,9,11),
                       c(1,2,3,4,5,6,7,8,13),
                       c(1,2,6,7,9,12),
                       c(1,4,5,8,9,10,11,12,13),
                       c(1,2,3,4,5,6,7,8,9,10,11,12,13),
                       c(1,4,7,8,10,11))
    EdgeSequence1 <-
        list(c(1,2,3),
             c(1,2,4,2,3,4),
             list(c(1,2,5,3,4,6),
                  c(1,2,6,2,4,6,2,3,4,2,5,3)),
             list(c(1,2,5,3,4,6),
                  c(1,3,5,3,4,5,2,5,4,1,2,6,1,6,3,2,4,6)),
             c(1,3,2,2,3,5,3,4,5),
             list(c(1,2,6,2,5,6,3,4,7),
                  c(1,2,7,2,4,7,2,5,4,3,4,5,3,5,6,1,6,3,1,3,7,7,3,4),
                  c(1,2,7,2,4,7,2,5,4,3,4,5,3,5,6)),
             list(c(1,8,2,3,7,6,4,5,9),
                  c(1,8,2,3,9,4,3,7,9,5,9,7,5,7,6),
                  c(3,7,6,1,8,4,1,4,5,1,9,2,1,5,9),
                  c(4,5,9,1,7,6,1,6,2,2,6,3,2,3,8),
                  c(1,10,2,2,10,9,5,9,10,5,10,6,6,10,7,3,7,10,3,10,4,4,10,8,1,
                    8,10),
                  c(1,10,2,1,10,7,6,10,7,5,10,6,5,9,10,4,9,10,3,10,4,3,8,10,2,
                    10,8),
                  c(5,9,10,2,10,9,1,10,2,1,7,10,6,10,7,3,10,6,3,8,10,4,10,8,4,
                    5,10),
                  c(1,7,6,1,6,9,1,9,2,5,9,6,3,8,4),
                  c(1,7,8,3,8,7,3,7,6,3,6,5,3,5,4,4,5,9,4,9,8,2,8,9,1,8,2)),
             c(1,2,3,1,3,4),
             c(1,2,6,1,3,5,1,6,3,3,4,5),
             list(c(1,4,5,2,6,3,3,6,7,4,8,5),
                  c(1,4,3,1,3,2,3,4,8,3,8,7,5,8,6,6,8,7,1,2,6,1,6,5),
                  c(1,2,9,1,9,5,5,9,8,4,8,9,3,4,9,3,9,7,6,7,9,2,6,9),
                  c(5,9,6,1,9,5,1,4,9,4,8,9,7,9,8,3,9,7,2,9,3,2,6,9),
                  c(1,2,5,2,6,5,3,4,8,3,8,7)),
             c(1,2,3,1,3,6,1,6,5,3,4,6),
             list(c(1,6,2,2,6,8,3,5,4,6,7,8),
                  c(1,5,3,1,3,6,3,6,7,3,7,4,2,4,8,2,5,4,1,5,2,4,7,8),
                  c(1,9,2,2,9,5,3,5,9,4,9,8,3,9,4,7,8,9,6,7,9,1,9,6),
                  c(4,9,5,1,5,9,1,9,2,2,9,8,7,8,9,6,9,7,3,6,9,3,9,4),
                  c(1,5,2,3,8,4,3,6,8,6,7,8)),
             list(
                  ##13.1
                  c(7,8,12,2,3,11,1,4,9,5,6,10),
                  ##13.2
                  c(2,  3, 11,  7,  8, 12,  9,  5,  4,  5,  4,  6,  4,  6,  1,
                    6,  1, 10),
                  c(1,  4,  9,  7,  8, 12, 10,  2,  5,  2,  5,  3,  5,  3,  6,
                    3,  6, 11),
                  c(5,  6, 10,  1,  4,  9, 11,  7,  2,  7,  2,  8,  2,  8,  3,
                    8,  3, 12),
                  c(5,  6, 10,  2,  3, 11, 12,  4,  7,  4,  7,  1,  7,  1,  8,
                    1,  8,  9),
                  c(5,  6, 10,  7,  8, 12,  2, 11,  1, 11,  1,  9, 11,  9,  3,
                    9,  3,  4),
                  c(2,  3, 11,  4,  1,  9,  5, 10,  8, 10,  8, 12, 10, 12,  6,
                    12,  6,  7),
                  ##13.3
                  c(7, 8, 12, 13, 3, 11, 13, 11, 6, 13,  6, 5, 13, 5,  9, 13,
                    9,  4, 13,  4,  1, 13, 1,10,13,10, 2,13,2,3),
                  c(2, 3, 11, 13, 6, 10, 13, 10, 1, 13,  1, 4, 13, 4, 12, 13,
                    12,  7, 13,  7,  8, 13, 8, 9,13,9, 5, 13,5,6),
                  c(7, 8, 12, 13, 6,  5, 13,  5, 9, 13,  9, 4, 13, 4,  3, 13,
                    3, 11, 13, 11,  2, 13, 2, 1,13,1,10,13,10,6),
                  c(2, 3, 11, 13,  4,  1, 13,  1, 10, 13, 10, 6, 13, 6, 7,13,
                    7, 12, 13, 12,  8, 13, 8, 5,13,5, 9,13, 9,4),
                  c(1, 4, 9,  13,  8, 12, 13, 12,  3, 13,  3, 2, 13, 2, 10, 13,
                    10, 5, 13,  5, 6, 13, 6,11,13,11,7,13, 7,8),
                  c(7, 8, 12, 13,  5,  6, 13,  6, 11, 13, 11, 3, 13, 3,  4, 13,
                    4, 9, 13,  9, 1, 13, 1, 2,13,2,10,13,10,5),
                  c(1, 4,  9, 13,  3,  2, 13,  2, 10, 13, 10, 5, 13, 5,  8, 13,
                    8, 12,13, 12, 7, 13, 7, 6,13,6,11,13,11,3),
                  c(5, 6, 10, 13,  1,  9, 13,  9,  8, 13,  8, 7, 13, 7, 11, 13,
                    11,  2,13,  2,  3, 13, 3, 12, 13, 12, 4,13,  4, 1),
                  c(5, 6, 10, 13,  8,  7, 13,  7, 11, 13, 11, 2, 13, 2,  1, 13,
                    1,  9,13,  9,  4, 13, 4,  3, 13, 3, 12,13, 12, 8),
                  c(1, 4,  9, 13,  2,  3, 13,  3, 12, 13, 12, 8, 13, 8,  5, 13,
                    5, 10, 13, 10, 6, 13, 6,  7, 13, 7, 11,13, 11, 2),
                  c(5, 6, 10, 13,  7,  8, 13,  8,  9, 13,  9, 1, 13, 1,  2, 13,
                    2, 11, 13, 11, 3, 13, 3,  4, 13, 4, 12,13, 12, 7),
                  c(2, 3, 11, 13,  1,  4, 13,  4, 12, 13, 12, 7, 13, 7,  6, 13,
                    6, 10, 13, 10, 5, 13, 5,  8, 13, 8, 9, 13,  9, 1),
                  ##13.4
                  c(13, 3, 11, 13, 11, 6, 13, 6, 7, 13, 7, 12, 13, 12, 8, 13,
                    8, 5, 13, 5, 9, 13, 9, 4,
                    13, 4, 1,  13, 1, 10, 13,10, 2, 13, 2,  3),
                  c(13, 4, 12, 13, 12, 7, 13, 7, 8, 13, 8,  9, 13,  9, 5, 13,
                    5, 6, 13, 6, 10, 13, 10, 1, 13,
                    1,  2, 13,  2, 11, 13, 11,  3, 13,  3,  4),
                  c(13, 2, 10, 13, 10, 5, 13, 5, 6, 13, 6, 11, 13, 11, 7, 13,
                    7, 8, 13, 8, 12, 13,12, 3, 13,
                    3,  4, 13,  4,  9, 13,  9,  1, 13,  1,  2),
                  c(13, 1,  9, 13,  9, 8, 13, 8, 5, 13, 5, 10, 13, 10, 6, 13,
                    6, 7, 13, 7, 11, 13,11, 2, 13,
                    2,  3, 13,  3, 12, 13, 12,  4, 13,  4,  1),
                  ##13.5.1
                  c(7,  8, 12,  2,  1, 10,  3,  4, 11,  4, 11,  6,  4,  6,  9,
                    6,  9,  5),
                  c(3,  2, 11,  8,  5,  9,  4,  1, 12,  1, 12,  7,  1,  7, 10,
                    7, 10,  6),
                  c(1,  4,  9,  6,  7, 11,  2,  3, 10,  3, 10,  5,  3,  5, 12,
                    5, 12,  8),
                  c(5,  6, 10,  4,  3, 12,  1,  2,  9,  2,  9,  8,  2,  8, 11,
                    8, 11,  7),
                  ##13.5.2
                  c(1, 2, 10, 8, 5, 9,  8, 9,  4, 8, 4,12, 4,12, 3,12, 3,11,
                    12,11, 7,11, 7, 6, 7, 6, 8, 6, 8, 5),
                  c(8, 5,  9, 3, 4, 12, 3, 12, 7, 3, 7,11, 7,11, 6,11, 6,10,
                    11,10, 2,10, 2, 1, 2, 1, 3, 1, 3, 4),
                  c(6, 7, 11, 1, 2, 10, 1, 10, 5, 1, 5, 9, 5, 9, 8, 9, 8,12,
                    9,12, 4,12, 4, 3, 4, 3, 1, 3, 1, 2),
                  c(3, 4, 12, 6, 7, 11, 6, 11, 2, 6, 2, 10,2, 10, 1,10,1, 9,
                    10,9, 5, 5, 9, 8, 5, 8, 6, 8, 6,  7)),
             c(1,4,2,1,6,4,1,5,6,3,6,4))

    switch23 <- function(x){
        num <- length(x) / 3
        temp <- x[0: (num-1)*3 + 2]
        x[0: (num-1)*3 + 2] <- x[0: (num-1)*3+ 3]
        x[0: (num-1)*3 + 3] <- temp
        x
    }

    SwitchSeq <- function(ed){
        if (is.list(ed)){
            lapply(1:length(ed), function(x) SwitchSeq(ed[[x]]))
        }
        else
            switch23(ed)
    }

    EdgeSequence2 <- SwitchSeq(EdgeSequence1)

    getedge <- function(x){
        case <- x[1]
        rotation <- x[2]
        map <- rep(0,8)
        for(i in 1:8){
            temp <- as.integer(BasicRotation[rotation,][i])
            map[temp] <- i
        }
        sapply(BasicEdges[[case-1]], function(x){
            if (x!=13){
                EndP1 <- EdgePoints[x,2]
                EndP2 <- EdgePoints[x,3]
                newEdge <- EdgePoints[(EdgePoints[,2]==map[EndP1]
                                       &EdgePoints[,3]==map[EndP2])|
                                      (EdgePoints[,3]==map[EndP1]
                                       &EdgePoints[,2]==map[EndP2]),][1]
             }
            else  newEdge <- 13
            newEdge})
    }

    GetEdges <- local({
        Edges <- apply(CaseRotationFlip[-c(1,256),], 1, getedge)
        Case <- cbind(seq(1:256), CaseRotationFlip[,c(1,3)])
        Edges <- apply(Case[-c(1,256),], 1, function(x){
            case <- x[2]-1
            EdgeNum <- x[1]-1
            if (x[3]==1)
                sapply(EdgeSequence1[[case]], function(x) Edges[[EdgeNum]][x])
            else sapply(EdgeSequence2[[case]], function(x) Edges[[EdgeNum]][x])
        })
        Edges
    })

    BasicFace <- list(c(0),c(0),c(1),c(7),c(0),c(2,7),c(1,2,6,7),c(0),c(0),
                      c(5,6,7),c(0), c(1,4,7),c(1,2,3,4,5,6,7),c(0))
    FacePoints <- matrix(c(seq(1,6),1,2,4,1,1,5,6,7,7,8,3,7,2,3,3,4,2,
                           6,5,6,8,5,4,8),
                         ncol=5)
    FacePoints <- cbind(FacePoints, apply(FacePoints[,2:5],1,prod))

    getface <- function(x) {
        case <- x[1]
        rotation <- x[2]
        map <- rep(0,8)
        for(i in 1:8){
            temp <- as.integer(BasicRotation[rotation,][i])
            map[temp] <- i
        }
        sapply(BasicFace[[case-1]], function(x){
            EndP <- rep(0,4)
            if (x==0) newFace <- 0
            else if (x==7) newFace <- 7
            else {
                for (i in 1:4){
                     point <- FacePoints[x,i+1]
                     EndP[i] <- map[point]
                }
                newFace<- FacePoints[FacePoints[,6]==prod(EndP[1:4]),][1]
            }
            newFace})
    }

    flipface <- function(case, face){
        if (face!=0){
            index <- explode(case+1)
            if (sum(index) > 4)
                index <- ifelse(index==0,1,0)
            if (face!=7 && index[FacePoints[face,2]]==0)
                face <- -face
            else if (face==7){
                tcase <- CaseRotationFlip[case+1,1]-1
                if ((tcase == 4 || tcase==6 || tcase==10 ||tcase==12)
                    && !(index[1]+index[7]==2) && !(index[3]+index[5]==2))
                    face <- -face
                else if (tcase==7
                         && !(index[1]+index[7]==0) && !(index[3]+index[5]==0))
                    face <- -face
            }
        }
        face
    }

    GetFaces <- local({
        Faces <- apply(CaseRotationFlip[-c(1,256),], 1, getface)
        for (i in 1:254)
            for(j in 1:length(Faces[[i]]))
                Faces[[i]][j] <- flipface(i, Faces[[i]][j])
        Faces
    })

    ## special
    ## name : name of the case
    ## nface: how many cases need to be checked
    ## sev: whether face 7 need to be checked
    ## nedge: total number of edges in the lookuptable
    ## ind: the index needed to check the lookuptable
    ## position: the corresponding positions in the lookuptable.
    special <- list(name =  c(3, 4, 6, 7,  10,12,13),
                    nface=  c(1, 1, 2, 4,  3, 3, 7),
                    sev =   c(0, 1, 1, 1,  1, 1, 1),
                    nedge = c(18,24,48,177,96,96,816),
                    ind = list(c(0,1), c(0,1), c(0,2,1,3),
                               c(0,8,4,12,2,10,1,9,6,14,5,13,3,11,15,7),
                               c(0,4,1,5,2,6,3,7),c(0,4,2,6,1,5,3,7),
                               c(0,1,2,4,8,16,32,3,9,17,33,6,18,34,12,20,36,24,
                                 40,35,25,22,44,19,41,38,28,83,105,102,92)),
                    position=list(list(c(1:6),c(7:18)),
                                  list(c(1:6),c(7:24)),
                                  list(c(1:9),c(10:33),c(34:48),c(34:48)),
                                  list(c(1:9),c(1:9),c(10:24),c(10:24),
                                       c(25:39),c(25:39),c(40:54), c(40:54),
                                       c(55:81), c(55:81),c(82:108),c(82:108),
                                       c(109:135),c(109:135),c(136:150),
                                       c(151:177)),
                                  list(c(1:12),c(13:36),c(37:60),c(37:60),
                                       c(61:84), c(61:84),c(85:96),c(85:96)),
                                  list(c(1:12), c(13:36), c(37:60), c(37:60),
                                       c(61:84), c(61:84),c(85:96),c(85:96)),
                                  list(c(1:12),
                                       c(13:30),c(31:48),c(49:66),c(67:84),
                                       c(85:102),c(103:120),
                                       c(121:150),c(151:180),c(181:210),
                                       c(211:240),c(241:270),c(271:300),
                                       c(301:330),
                                       c(331:360),c(361:390),c(391:420),
                                       c(421:450),c(451:480),
                                       c(481:516),c(517:552),c(553:588),
                                       c(589:624),
                                       c(625:642),c(643:660),c(661:678),
                                       c(679:696),
                                       c(697:726),c(727:756),c(757:786),
                                       c(787:816))
                                 )
                    )

    list(Edges = GetEdges, Faces = GetFaces, EdgePoints = EdgePoints,
         FacePoints = FacePoints, CaseRotationFlip = CaseRotationFlip,
         special = special)
})

Faces <- PreProcessing$Faces
Edges <- PreProcessing$Edges
EdgePoints <- PreProcessing$EdgePoints
FacePoints <- PreProcessing$FacePoints
CaseRotationFlip <- PreProcessing$CaseRotationFlip
special <- PreProcessing$special

fgrid <- function(fun, x, y, z) {
    g <- expand.grid(x = x, y = y, z = z)
    array(fun(g$x, g$y, g$z), c(length(x), length(y), length(z)))
}

faceType <- function(v, nx, ny, level, maxvol) {
    ## the following line replaces: v <- ifelse(v > level, 1, 0)
    if(level==maxvol)
      p <- v >= level
    else p <- v > level
    v[p] <- 1; v[! p] <- 0
    v[-nx, -ny] + 2 * v[-1, -ny] + 4 * v[-1, -1] + 8 * v[-nx, -1]
}

levCells <- function(v, level, maxvol) {
    nx <- dim(v)[1]
    ny <- dim(v)[2]
    nz <- dim(v)[3]
    cells <- vector("list", nz - 1)
    types <- vector("list", nz - 1)

    bottomTypes <- faceType(v[,,1], nx, ny, level, maxvol)
    for (k in 1 : (nz - 1)) {
        topTypes <- faceType(v[,, k + 1], nx, ny, level, maxvol)
        cellTypes <- bottomTypes + 16 * topTypes
        contourCells <- which(cellTypes > 0 & cellTypes < 255)
        cells[[k]] <- contourCells + (nx - 1) * (ny - 1) * (k - 1)
        types[[k]] <- as.integer(cellTypes[contourCells])
        bottomTypes <- topTypes
    }
    cells <- unlist(cells)
    i <- as.integer((cells - 1) %% (nx - 1) + 1)
    j <- as.integer(((cells - 1) %/% (nx - 1)) %% (ny - 1) + 1)
    k <- as.integer((cells - 1) %/% ((nx - 1) * (ny - 1)) + 1)
    t <- unlist(types)
    list(i = i, j = j, k = k, t = t)
}

CalPoint <- function(x1,x2,y1,y2,z1,z2,v1,v2){
    s <- v1 / (v1-v2)
    x <- x1+s*(x2-x1)
    y <- y1+s*(y2-y1)
    z <- z1+s*(z2-z1)
    c(x,y,z)
}

GetPoints<-function(edge, p1, info){
    ##**** need better name than info
    ## info is the output from GetBasic()
    x1 <- EdgePoints[edge,2]
    x2 <- EdgePoints[edge,3]
    c((1-floor(x1/9))*info[p1+x1-1,1]+floor(x1/9)*info[p1,1],
      (1-floor(x1/9))*info[p1+x2-1,1]+floor(x1/9)*info[p1+1,1],
      (1-floor(x1/9))*info[p1+x1-1,2]+floor(x1/9)*info[p1+1,2],
      (1-floor(x1/9))*info[p1+x2-1,2]+floor(x1/9)*info[p1+2,2],
      (1-floor(x1/9))*info[p1+x1-1,3]+floor(x1/9)*info[p1+1,3],
      (1-floor(x1/9))*info[p1+x2-1,3]+floor(x1/9)*info[p1+5,3],
      (1-floor(x1/9))*info[p1+x1-1,4]+floor(x1/9)*(0*info[p1+1,3]+1),
      (1-floor(x1/9))*info[p1+x2-1,4]+floor(x1/9)*(0*info[p1+1,3]-1))
}

FaceNo7 <- function(faces, p1, info){
    ##**** need better name than info
    ## info is the output from GetBasic()
    index <- ifelse(faces > 0, 1, -1)
    faces <- abs(faces)
    e1 <- FacePoints[faces,2]
    e2 <- FacePoints[faces,3]
    e3 <- FacePoints[faces,4]
    e4 <- FacePoints[faces,5]
    A <- info[p1+e1-1,4]
    B <- info[p1+e2-1,4]
    C <- info[p1+e3-1,4]
    D <- info[p1+e4-1,4]
    index <- index*ifelse (A*B-C*D > 0, 1, -1)
    ifelse(index==1, 1, 0)
}

Face7 <- function(faces, p1, info){
    ## info is the output from GetBasic()
    index <- ifelse(faces > 0, 1, -1)
    A0 <- info[p1,4];   B0 <- info[p1+3,4]
    C0 <- info[p1+2,4]; D0 <- info[p1+1,4]
    A1 <- info[p1+4,4]; B1 <- info[p1+7,4]
    C1 <- info[p1+6,4]; D1 <- info[p1+5,4]
    a <- (A1 - A0)*(C1 - C0) - (B1 - B0)*(D1 - D0)
    b <- C0*(A1 - A0) + A0*(C1 - C0) - D0*(B1 - B0) - B0*(D1 - D0)
    c <- A0*C0 - B0*D0
    tmax <- -b/(2*a)
    maximum <- a*tmax^2 + b*tmax + c
    maximum <- ifelse(maximum=="NaN",-1,maximum)
    cond1 <- ifelse (a < 0, 1 ,0)
    cond2 <- ifelse (tmax > 0, 1 ,0)
    cond3 <- ifelse (tmax < 1, 1, 0)
    cond4 <- ifelse (maximum >0, 1, 0)
    totalcond <- cond1 * cond2 * cond3 * cond4
    index <- index*ifelse(totalcond==1, 1, -1)
    ifelse(index==1, 1, 0)
}


## GetBasic()--- The output matrix "information" consists of 4 columns
## and #-of-cubes*8 rows.  The first 3 columns tell the
## coordinate(x,y,z) of each vertex of the cube, the 4th gives the
## intensity minus the threshold, which actually makes the threshold
## eaqual to 0. This is convenient for further judgment of subcases

GetBasic <- function(R, vol, level, v) {
    cube.1 <- cbind(v$i[R], v$j[R], v$k[R])
    index <- matrix(c(0,1,1,0,0,1,1,0,
                      0,0,1,1,0,0,1,1,
                      0,0,0,0,1,1,1,1),
                    nrow=8)
    ax.inc <- c(1,1,1)

    ver.inc <- t(apply(index,1, function(x) x*ax.inc))
    cube.co <-
        kronecker(rep(1,nrow(cube.1)),ver.inc) + kronecker(cube.1,rep(1,8))

    value <- vol[cube.co] - level
    information <- cbind(cube.co, value)
    information <- rbind(information, rep(0, 4))
    p1 <- (1:length(R) - 1) * 8 + 1
    cases <- v$t[R]
    list(information=information, p1 = p1, cases=cases)
}

PreRender <- function(edges, p1, type, info) {
    if(type==1){
        if (typeof(edges)=="list"){
            count <- sapply(edges, function(x) length(x))
            edges <- cbind(unlist(edges), rep(p1,count))
        }
        else{
            count <- nrow(edges)
            edges <- cbind(as.vector(t(edges)), rep(p1,each=count))
        }
    }
    else{
        if (is.vector(edges))
            edges <- matrix(edges, ncol = length(edges))
        p1 <- edges[, 1]
        count <- ncol(edges) - 1
        edges <- cbind(as.vector(t(edges[, -1])), rep(p1, each = count))
    }
    ##The output of GetPoints() are coordinates of cubes.
    info <- GetPoints(edges[,1],edges[,2], info)
    info <- matrix(info,ncol=8)
    ##The output of CalPoint() are coordinates of triangles.
    info <- CalPoint(info[,1],info[,2],info[,3],info[,4],
                            info[,5],info[,6],info[,7],info[,8])
    matrix(info,ncol=3)
}

rescale <- function(i, x) {
    nx <- length(x)
    low <- pmin(pmax(1, floor(i)), nx - 1)
    x[low] + (i - low) * (x[low + 1] - x[low])
}

computeContour3d <- function (vol, maxvol = max(vol), level,
                              x = 1:dim(vol)[1],
                              y = 1:dim(vol)[2],
                              z = 1:dim(vol)[3], mask) {

    nx <- length(x)
    ny <- length(y)
    nz <- length(z)

    if (missing(mask)) mask <- NULL
    if (is.function(mask)) mask <- fgrid(mask, x, y, z)
    if (! all(mask)) vol[! mask] <- NA

    v <- levCells(vol, level, maxvol)
    tcase <- CaseRotationFlip[v$t+1,1]-1

    R <- which(tcase %in% c(1,2,5,8,9,11,14))
    if (length(R) > 0){
        Basics <- GetBasic(R, vol, level, v)
        information <- Basics$information
        p1 <- Basics$p1
        cases <- Basics$cases
        edges <- Edges[cases]
        triangles <- PreRender(edges, p1,type=1, information)
    }
    else triangles <- matrix(0, nrow=0,ncol=3) # emty contour, e.g.

    for (i in 1:length(special$name)){
        R <- which(tcase == special$name[i])
        if (length(R) > 0) {
            Basics <- GetBasic(R, vol, level, v)
            information <- Basics$information
            p1 <- Basics$p1
            cases <- Basics$cases

            nface <- special$nface[i]
            nedge <- special$nedge[i]
            faces <- matrix(unlist(Faces[cases]), ncol = nface, byrow = TRUE)

            if (i==1)
                index <- FaceNo7(faces[, 1], p1, information)
            else if (i==2)
                index <- Face7(faces[, 1], p1, information)
            else{
                index <-  Face7(faces[, nface], p1, information)*2^(nface-1)
                for(j in 1:(nface-1)){
                    temp <-  FaceNo7(faces[, j], p1, information)
                    index <- index + temp * 2^(j-1)
                }
            }
            edges <- matrix(unlist(Edges[cases]), ncol = nedge, byrow = TRUE)
            edges <- cbind(edges, p1, index)
            ind <- special$ind[[i]]
            position <- special$position[[i]]

            for (j in 1:length(ind)){
                ed <- edges[which(index == ind[j]), c(nedge+1, position[[j]])]
                if (length(ed) > 0) {
                    prtri <- PreRender(ed,nedge+1,type=2, information)
                    triangles <- rbind(triangles, prtri)
                }
            }
        }
    }

    if (! identical(x, 1 : nx)) triangles[,1] <- rescale(triangles[,1], x)
    if (! identical(y, 1 : ny)) triangles[,2] <- rescale(triangles[,2], y)
    if (! identical(z, 1 : nz)) triangles[,3] <- rescale(triangles[,3], z)

    triangles
}

contourTriangles <- function(vol, maxvol, level,
                             x = 1:dim(vol)[1],
                             y = 1:dim(vol)[2],
                             z = 1:dim(vol)[3],
                             mask = NULL, color = "white", color2 = NA,
                             alpha = 1, fill = TRUE,
                             col.mesh = if (fill) NA else color,
                             material = "default", smooth = 0) {
    if (length(level) > 1) {
        val <- vector("list", length(level))
        for (i in seq(along = level)) {
            m <- if (is.list(mask)) mask[[i]] else mask
            col <- if (length(color) > 1) color[[i]] else color
            col2 <- if (length(color2) > 1) color2[[i]] else color2
            a <- if (length(alpha) > 1) alpha[[i]] else alpha
            fl <- if (length(fill) > 1) fill[[i]] else fill
            cm <- if (length(col.mesh) > 1) col.mesh[[i]] else col.mesh
            mat <- if (length(material) > 1) material[[1]] else material
            sm <- if (length(smooth) > 1) smooth[[1]] else smooth
            val[[i]] <- contourTriangles(vol, maxvol, level[i], x, y, z, m,
                                         col, col2, a, fl, cm, mat, sm)
        }
        val
    }
    else makeTriangles(computeContour3d(vol, maxvol, level, x, y, z, mask),
                       color = color, color2 = color2, alpha = alpha,
                       fill = fill, col.mesh = col.mesh,
                       material = material, smooth = smooth)
}

contour3d <- function(f, level,
                      x = 1:dim(f)[1], y = 1:dim(f)[2], z = 1:dim(f)[3],
                      mask = NULL, color = "white", color2 = NA, alpha = 1,
                      fill = TRUE, col.mesh = if (fill) NA else color,
                      material = "default", smooth = 0,
                      add = FALSE, draw = TRUE, engine = "rgl",
                      separate=FALSE,...){

    if (! all(is.finite(x), is.finite(y), is.finite(z)))
        stop("'x', 'y', and 'z' values must be finite and non-missing")
    if (is.function(f) || is.array(f) && length(dim(f)) == 3){
        if (is.function(f)){
            if (length(formals(f)) < 3)
                stop("The function must have at least 3 arguments.")
            vol <- fgrid(f, x, y, z)
          }
        else{
          if (dim(f)[1] != length(x) || dim(f)[2] != length(y) ||  dim(f)[3] != length(z))
            stop("dimensions of f do not match x, y, or z")
          vol <- f
        }

        maxvol <- max(vol)
        minvol <- min(vol)
        #cat("The range of 'f' is between ", round(minvol,2), " and ", round(maxvol,2), ".\n", sep="")
        con <- which(! level <= maxvol & level >= minvol)
        if (length(con) == length(level))
            stop(paste("The 'level' has to be within the range of 'f' (between ", round(minvol, 2), " and ", round(maxvol, 2),").\n", sep=""))
        else if (length(con) > 0){
            warning(paste("The 'level' outside the range of 'f' (between ", round(minvol, 2), " and ", round(maxvol, 2), ") has been removed. \n", sep=""))
            level <- level[-con]
            if (is.list(mask)) mask <- mask[-con]
            if (length(color) > 1) color <- color[-con]
            if (length(color2) > 1) color2 <- color2[-con]
            if (length(alpha) > 1) alpha <- alpha[-con]
            if (length(fill) > 1) fill <- fill[-con]
            if (length(col.mesh) > 1) col.mesh <- col.mesh[-con]
            if (length(material) > 1) material <- material[-con]
            if (length(smooth) > 1) smooth <- smooth[-con]
        }

      }

    else stop("vol has to be a function or a 3-dimensional array")



    scene <- contourTriangles(vol, maxvol, level, x, y, z, mask, color, color2,
                              alpha, fill, col.mesh, material, smooth)
    if (! draw || engine == "none"){
        if (! any(separate))
            scene
        else{
            if (length(level)==1){
                newScene <- separateTriangles(scene)
                cat("Triangles are separated into ", length(newScene),
                          " chunks.", "\n", sep="")
            }
            else{
                if (length(separate) < length(level))
                    separate <- c(separate, rep(FALSE, length(level)-length(separate)))

                newScene <- NULL
                for (i in 1:length(level)){
                    if (separate[i]){
                        new <- separateTriangles(scene[[i]])
                        newScene <- c(newScene, new)
                        cat("Triangles from level ", level[i],
                                  " are separated into ", length(new), " chunks.",
                                  "\n", sep="")
                }
                    else
                        newScene <- c(newScene, list(scene[[i]]))
                }
            }
            newScene
        }
    }
    else {
        scene <- colorScene(scene)
        if (engine == "rgl")
            drawScene.rgl(scene, add = add, ...)
        else if (engine %in% c("standard", "grid"))
            drawScene(scene, add = add, engine = engine, ...)
        else stop(paste("unknown rendering engine:", engine))
    }
}
