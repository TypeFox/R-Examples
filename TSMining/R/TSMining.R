#' A function to perform symbolic approximation aggregate (SAX) for time series data
#' 
#' The function create SAX symbols for a univariate time series. The details of this method can be referred to J. Lin, E. Keogh, L. Wei, S. Lonardi. Experiencing SAX: a novel symbolic representation of time series
#'
#' @import foreach
#' @import plyr
#' 
#' @param x is a numeric vector representing the univariate time series
#' @param w is the word size and should be an integer
#' @param a is the alphabet size and should be an integer
#' @param eps is the minimum threshold for variance in x and should be a numeric value. If x has a smaller variance than eps, it will represented as a word using the middle alphabet.
#' @param norm is a logical value deciding whether standardization should be applied to x. If True, x is standardized using mean and standard deviation
#'
#' @return The function returns a SAX representation of x
#'
#' @examples
#' x <- runif(n = 20, min = 0, max = 20)
#' Func.SAX(x = x, w = 5, a = 5, eps = .01, norm = TRUE)
#'
#' @export

Func.SAX  <- function(x, w, a, eps, norm) {  
    i=NULL
    if(sd(x) <= eps) {sym <- rep(letters[round((1+a)/2, digits = 0)], w)} else {
        #Normalize the data to have 0 mean and 1 standard deviation before piecewise aggregation
        if(norm==TRUE) {data.nor <- (x-mean(x))/sd(x)} else {
            data.nor <- x
        }  
        
        #Perform the piecewise aggregation
        ind <- round(seq(from = 1, to = length(data.nor), length.out = w+1), digits = 0)
        pieces <- foreach::foreach(i=1:(length(ind)-1), .combine = c) %do% {
            if(i!=(length(ind)-1)) {piece <- data.nor[ind[i]:(ind[i+1]-1)]} else {piece <- data.nor[ind[i]:ind[i+1]]}
            return(mean(piece, na.rm = T))   
        }
        
        #Perform alphabet assignment
        let <- letters[1:a]
        #Create breaks points based on Gaussian normal distribution
        bks <- round(qnorm(p = seq(from = 0, to = 1, length.out = a+1)), digits = 2)
        sym <- foreach::foreach(i=1:length(pieces), .combine = c) %do% {
            obs <- pieces[i]
            let[max(which(bks<obs))]
        }    
    }
    
    return(sym)
}

#' A function to create the distance matrix for alphabets
#'
#' This function create a distance matrix for alphabets used for SAX transformation
#'
#' @param a is an integer specifying the alphabet size.
#'
#' @return The function returns a matrix showing the distance between alphabets
#'
#' @examples
#'
#' Func.matrix(a=5)
#'
#' @export

Func.matrix <- function(a) {
    i=j=NULL
    let <- letters[1:a]
    bks <- round(qnorm(p = seq(from = 0, to = 1, length.out = (a+1))), digits = 2)
    bks.upd <- bks[-c(1,length(bks))]
    
    #Create the matrix for distance calculation
    dist.m <- matrix(data = NA, nrow = a, ncol = a, dimnames = list(let, let))
    for(i in 1:dim(dist.m)[1]) {
        for(j in 1:dim(dist.m)[2]) {
            dist.m[i,j] <- ifelse(abs(i-j)<=1, 0, bks.upd[max(c(i,j))-1]-bks.upd[min(c(i,j))])
        }
    }
    return(dist.m)
}

#' A function to calculate the distance between two SAX representations
#'
#' This function calculates the distance between two SAX representations
#'
#' @import foreach
#' @import plyr
#'
#' @param x is a SAX representations. 
#' @param y is a SAX representations. It should have the same length as x.
#' @param mat is the distance matrix created by \code{Func.matrix}
#' @param n is the length of the original time series before the SAX transformation
#'
#' @return The function returns a numeric value, which is the distance between two SAX representations
#'
#' @examples
#' #Assuming the original time series has a length of 20, n=20
#' #Assuming the time series is transformed into SAX representations using w=4 and a=4
#' #Assuming one is a,b,c,d and the other is d,b,c,d
#' Func.dist(x=c("a","b","c","d"), y=c("d","b","c","d"), mat=Func.matrix(a=4), n=20) 
#' 
#' @export

Func.dist <- function(x, y, mat, n) {
    i=NULL
    w <- length(x)
    d <- foreach::foreach(i=1:length(x), .combine = c) %do% {
        mat[which(rownames(mat)==x[i]), which(colnames(mat)==y[i])]
    }
    return(sqrt(sum(d^2))*sqrt(n/w))
}

#' A function implementing the univaraite motif discovery algorithm using random projection
#' 
#' The function implements the univariate motif discovery algorithm proposed in B. Chiu, E. Keogh, S. Lonardi. Probabilistic discovery of time series motifs. ACM SIGKDD, Washington, DC, USA, 2003, pp. 493-498.
#' 
#' @import foreach
#' @import plyr
#' 
#' @param ts is a numeric vector representing the univarate time series
#' @param global.norm is a logical value specifying whether global standardization should be used for the whole time series
#' @param local.norm is a logical value specifying whether local standardization should be used for each subsequences
#' @param window.size is a integer which defines the length of the sliding window used to create subsequences
#' @param overlap is a numeric value ranging from 0 to 1. It defines the percentage of overlapping when using sliding window to create subsquences. 0 means subsequences are created without overlaps. 1 means subsequences are created with the maximum overlap possible.
#' @param w is an integer which defines the word size used for SAX transformation 
#' @param a is an integer which defines the alphabet size used for SAX transformation
#' @param eps is the minimum threshold for variance in subsequence and should be a numeric value. If the subsequence considered has a smaller variance than eps, it will be represented as a word using the middle alphabet. The default value is 0.1
#' @param mask.size is the mask size used for random projection. It should be an integer ranging from 1 to the word size w
#' @param iter is an integer which specifies the iteration number in random projection, default value is 25
#' @param max.dist.ratio is a numeric value used to add other possible members to a motif candidate. Default value is 1.2. Each motif candidate has two subsequences. The distance between these two candidates are calculated as a baseline, denoted as BASE. Any subsequence, whose distance to the motif candidate is smaller than max.dist.ratio*BASE, is considered as a member of that motif candidate.  
#' @param count.ratio.1 defines the ratio between the iteration number and the minimum value in the collision matrix to be considered as motif candidate. Default value is 1.5. For instance, if the iter is 100, any pair of subsequence, which results in a value larger than 67 in the collision matrix, is considered as a motif candidate.
#' @param count.ratio.2 defines the ratio between the maximum counts in  the collision matrix and any other count values that will be considered as potential members to a motif candidate
#' 
#' @return The function returns a list of 6 elements. The first element is Subs, which is a data frame containing all the subsequences in original data formatThe second element is Subs.SAX, which is a data frame containing all the subsequences in SAX representations. The third element is Motif.raw, which is a list showing the motifs discovered in original data format.The fourth element is Motif.SAX, which is a list showing the motifs discovered in SAX representations. The fifth element is Collision.matrix, which is matrix containing the results of random projection. The sixth element is Indices, which is a list showing the starting positions of subsequences for each motif discovered.
#' 
#' @examples
#' #Perform the motif discovery for the first time series in the example data
#' data(test)
#' res.1 <- Func.motif(ts = test$TS1, global.norm = TRUE, local.norm = FALSE, 
#' window.size = 10, overlap = 0, w = 5, a = 3, mask.size = 3, eps = .01)
#' #Check the number of motifs discovered
#' length(res.1$Indices)
#' #Check the starting positions of subsequences of each motif discovered
#' res.1$Indices
#' 
#' @export

Func.motif <- function(ts, global.norm, local.norm, window.size, overlap, w, a, mask.size,
                       eps=0.1, iter=25, max.dist.ratio=1.2, count.ratio.1=1.5, count.ratio.2=1.2) {
    i=j=q=m=k=l=g=h=u=NULL
    #Perform uniform normalization 
    if(global.norm==TRUE) {ts.nor <- (ts - mean(ts))/sd(ts)} else {ts.nor <- ts}
    
    #Create the subsequence data frame
    b <- ifelse(overlap==1, yes = 1, no = round((1-overlap)*window.size, digits = 0))
    ts.subs <- foreach::foreach(i=seq(from = 1, to = length(ts), by = b), .combine = rbind) %do% {
        c(i,subs.temp <- ts.nor[i:(i+window.size-1)])
    }
    ts.subs <- na.omit(ts.subs)
    
    #Local normalization if needed
    ts.sax <- foreach::foreach(i=1:dim(ts.subs)[1], .combine = rbind) %do% {
        if (sd(ts.subs[i,-1])<=eps) {sax.temp <- rep(letters[round(a+1/2)], times = w)} else {
            sax.temp <- Func.SAX(x = ts.subs[i,-1], w = w, a = a, eps = eps, norm = local.norm)
        } 
        c(ts.subs[i,1], sax.temp)
    }
    
    ts.sax <- as.data.frame(ts.sax, stringsAsFactors = FALSE)
    colnames(ts.sax) <- c("StartP", 1:w)
    ts.sax$StartP <- as.numeric(ts.sax$StartP)
    
    #Perform the random projection
    col.mat <- matrix(data = 0, nrow = dim(ts.sax)[1], ncol = dim(ts.sax)[1])
    for(i in 1:iter) {
        col.pos <- sort(sample(x = 2:dim(ts.sax)[2], size = mask.size, replace = F), decreasing = F)
        sax.mask <- ts.sax[,col.pos]
        
        unique.lab <- unique(sax.mask)
        
        mat <- foreach::foreach(j = 1:dim(unique.lab)[1], .combine = rbind) %do% {
            indices <- foreach::foreach(q = 1:dim(sax.mask)[1], .combine = c) %do% {
                identical(as.character(sax.mask[q,]), as.character(unique.lab[j,]))
            }
        }
        
        for(m in 1:dim(mat)[1]) {
            if(length(which(mat[m,]==TRUE))>1) {
                com <- t(combn(x = which(mat[m,]==TRUE), m = 2))
                col.mat[com] <- col.mat[com] + 1    
            }
        }   
    }
    
    #Extact the tentative motif pair
    counts <- sort(col.mat, decreasing=TRUE)
    counts.sel <- counts[which(counts>=(iter/count.ratio.1))]
    
    motif.pair <- foreach::foreach(k = 1:length(unique(counts.sel)), .combine = rbind) %do% {
        arrayInd(which(col.mat==unique(counts.sel)[k]), .dim = dim(col.mat))    
    }
    
    indices <- foreach::foreach(l = 1:dim(motif.pair)[1]) %do% {
        pair <- c(ts.sax[motif.pair[l,1],1], ts.sax[motif.pair[l,2],1])
        
        cand.1 <- ts.subs[motif.pair[l,1],-1]
        cand.2 <- ts.subs[motif.pair[l,2],-1]
        dist.raw <- sqrt(sum((cand.1 - cand.2)^2))
        
        col.no <- col.mat[motif.pair[l,1],]
        ind.cand <- which(col.no > (max(col.no)/count.ratio.2))
        if(length(ind.cand)>1) {
            ind.temp <- ind.cand[-which(ind.cand == motif.pair[l,2])]
            if(length(ind.temp)==1) {
                df.cand.sel <- as.matrix(ts.subs[ind.temp,-1])
                dist.res <- plyr::aaply(.data = df.cand.sel, .margins = 2, .fun = function(x) sqrt(sum((cand.1-x)^2)))
                ind.final <- ts.sax[ind.temp[which(dist.res<=max.dist.ratio*dist.raw)],1]   
            } else {
                df.cand.sel <- ts.subs[ind.temp,-1]
                dist.res <- plyr::aaply(.data = df.cand.sel, .margins = 1, .fun = function(x) sqrt(sum((cand.1-x)^2)))
                ind.final <- ts.sax[ind.temp[which(dist.res<=max.dist.ratio*dist.raw)],1]  
            }} else {
                ind.final <- NULL
            }
        
        pair.final <- c(pair, ind.final)
    }
    
    #Combine the indices if there is any overlap
    vec.subset <- rep(0, length(indices))
    foreach::foreach(g = 1:(length(indices)-1), .combine = rbind) %do% {
        for (h in (g+1):length(indices)) {
            if(length(which(indices[[g]] %in% indices[[h]]))>0) {
                indices[[h]] <- unique(c(indices[[g]], indices[[h]]))
                vec.subset[g] <- 1
            }
        }
    }
    indices <- indices[vec.subset==0]
    
    motif.raw <- foreach::foreach(u = 1:length(indices)) %do% {
        ts.subs[which(ts.subs[,1] %in% indices[[u]]),]
    }
    
    motif.sax <- foreach::foreach(u = 1:length(indices)) %do% {
        ts.sax[which(ts.sax[,1] %in% indices[[u]]),]
    }
    
    return(list(Subs=ts.subs, Subs.SAX=ts.sax, Motif.raw=motif.raw, Motif.SAX=motif.sax, Collision.matrix=col.mat, Indices=indices))
}

#' A function to implment the multivariate motif discovery
#' 
#' This function implements the multivariate motif discovery method proposed in A. Vahdatpour, N. Amini, M. Sarrafzadeh. Towards unsupervised activity discovery using multi-dimensional motif detection in time series. IJCAI 2009 21st International Joint Conference on Artificial Intelligence.
#' 
#' @import foreach
#' @import plyr
#' 
#' @param motif.list is a list of lists, each contains the univariate motifs discovered in a univariate time series. The component of motif.list is the results of Func.motif()$Indices, which store the starting position of subsequences of each univariate motif 
#' @param window.sizes is a vector containing the length of motifs in each univariate time series. It should have the same order as components in motif.list.
#' @param alpha is a numeric ranging from 0 to 1. It specifies the minimum correlation between two univariate motifs before considered as a multivariate motifs
#' 
#' @return The function returns a list containing two elements. The first element is Motif, which is a list containing the univarate motif IDs for different multivariate motifs. e.g., if there are two univariate time series and each has 3 motifs, then univariate ID is from 1 to 6. The second element is Info, which is a list storing the information of univariate motifs for different multivariate motifs
#' 
#' @examples
#' data(test)
#' #Perform univariate motif discovery for each dimension in the example data
#' res.1 <- Func.motif(ts = test$TS1, global.norm = TRUE, local.norm = FALSE, 
#' window.size = 10, overlap = 0, w = 5, a = 3, mask.size = 3, eps = .01)
#' res.2 <- Func.motif(ts = test$TS2, global.norm = TRUE, local.norm = FALSE, 
#' window.size = 20, overlap = 0, w = 5, a = 3, mask.size = 3, eps = .01)
#' #Perform multivaraite motif discovery
#' res.multi <- Func.motif.multivariate(motif.list = list(res.1$Indices, res.2$Indices), 
#' window.sizes = c(10,20), alpha = .8)
#' 
#' @export

Func.motif.multivariate <- function(motif.list, window.sizes, alpha) {
    i=j=q=p=t=o=f=x=z=NULL
    #Get the total motif.no
    tot.no <- sum(plyr::laply(.data = motif.list, .fun = length))
    
    #Create the weight matrix
    w.mat <- matrix(data = 0, nrow = tot.no, ncol = tot.no)
    
    #Get the characteristics of motifs
    info <- foreach::foreach(i=1:length(motif.list), .combine = rbind) %do% {
        info.sub <- foreach::foreach(j=1:length(motif.list[[i]]), .combine = rbind) %do% {
            info.sub.sub <- foreach::foreach(q=1:length(motif.list[[i]][[j]]), .combine = rbind) %do% {
                c(i,j,q, motif.list[[i]][[j]][q], motif.list[[i]][[j]][q]+window.sizes[i]-1)   
            }
        }
    }
    rownames(info) <- 1:dim(info)[1]
    info <- as.data.frame(info)
    colnames(info) <- c("Variable","Motif","Member","StartP","EndP")
    info$Lab <- as.numeric(as.character(factor(info$Variable*100 + info$Motif, levels = unique(info$Variable*100 + info$Motif), labels = 1:length(unique(info$Variable*100 + info$Motif)))))
    info.ori <- info
    
    #Generate the weights
    pb <- txtProgressBar(min = 0, max = length(unique(info$Lab)))
    temp <- foreach::foreach (i = 1:length(unique(info$Lab))) %do% {
        setTxtProgressBar(pb = pb, value = i)
        mot.con <- info[which(info$Lab==i),]
        lab.con <- unique(mot.con$Lab)
        n <- dim(mot.con)[1]
        variable.con <- unique(mot.con$Variable)
        mot.com.all <- info[-which(info$Variable==variable.con),]
        
        temp.ind <- foreach::foreach (j = 1:length(unique(mot.com.all$Lab)), .combine = rbind) %do% {
            mot.com <- mot.com.all[which(mot.com.all$Lab==unique(mot.com.all$Lab)[j]),]
            lab.com <- unique(mot.com$Lab)
            count <- 0
            temp.ind.ind <- foreach::foreach (p = 1:dim(mot.con)[1], .combine = c) %do% {
                res.temp <- foreach::foreach (x = 1:dim(mot.com)[1], .combine = c) %do% {
                    out <- ifelse(mot.con[p,"StartP"]>mot.com[x,"EndP"] | mot.con[p,"EndP"]<mot.com[x,"StartP"], yes = 0, no = 1)
                    return(out)
                }
                if(sum(res.temp)!=0) {count <- count+1}
                return(ifelse(sum(res.temp)!=0, yes = 1, no = 0))
            }
            w.mat[lab.con ,lab.com] <- round(count/n, 2)
            return(c(lab.con, lab.com, temp.ind.ind))
        } 
        return(temp.ind)
    }
    w.ori <- w.mat
    
    #Perform grouping 
    #Get the occurrence of each motifs
    occurence <- plyr::daply(.data = info, .variables = "Lab", .fun = function(x) dim(x)[1])
    ord <- order(occurence, decreasing = T)
    
    re.temp <- foreach::foreach(t = 1:length(ord)) %do% {
        w <- w.mat[ord[t],]    
        ind.add <- which(w>alpha)
        s <- c(ord[t], ind.add)
        
        #Remove some of the occurences of motif j when included by motif i
        if(length(ind.add)>0) {
            out <- foreach::foreach(o = 1:length(ind.add), .combine = c) %do% {
                obs.con <- info[which(info$Lab==ind.add[o]),]
                mat <- temp[[ind.add[o]]]
                w.mat[ind.add[o], ord[t]] <- 0
                oo <- which(info$Lab==ind.add[o] & info$Member %in% c(which(mat[which(mat[,2]==ord[t]),-c(1,2)]==1)))
                return(oo)
            }
            info <- info[-unique(out),]
            
            #Update the weights associated with motif j
            #Generate new the weights for the jth row
            for (f in 1:length(ind.add)) {
                if(length(which(info$Lab==ind.add[f]))!=0) {
                    mot.con <- info[which(info$Lab==ind.add[f]),]
                    lab.con <- unique(mot.con$Lab)
                    n <- dim(mot.con)[1]
                    variable.con <- unique(mot.con$Variable)
                    mot.com.all <- info[-which(info$Variable==variable.con),]
                    
                    for (j in 1:length(unique(mot.com.all$Lab))) {
                        mot.com <- mot.com.all[which(mot.com.all$Lab==unique(mot.com.all$Lab)[j]),]
                        lab.com <- unique(mot.com$Lab)
                        count <- 0
                        for (p in 1:dim(mot.con)[1]) {
                            res.temp <- foreach::foreach (x = 1:dim(mot.com)[1], .combine = c) %do% {
                                out <- ifelse(mot.con[p,"StartP"]>mot.com[x,"EndP"] | mot.con[p,"EndP"]<mot.com[x,"StartP"], yes = 0, no = 1)
                                return(out)
                            }
                            if(sum(res.temp)!=0) {count <- count+1}
                        }
                        w.mat[lab.con ,lab.com] <- round(count/n, 2)
                    } 
                }
            }
            return(sort(s, decreasing = F))
        }  
    }
    ind.null <- which(plyr::laply(.data = re.temp, .fun = length)>0)
    mot.final <- re.temp[ind.null] 
    
    indices <- foreach::foreach(z = 1:length(mot.final)) %do% {
        info.ori[which(info.ori$Lab %in% mot.final[[z]]),]
    }
    
    return(list(Motif=mot.final, Info=indices))
}

#' A function to prepare the dataset for visualizing the univaraite motifs discovered
#' 
#' This function create a data set for the use of visualizing the univariate motifs discovered
#' 
#' @import foreach
#' @import plyr
#' @import reshape2
#' @import ggplot2
#' 
#' @param single.ts is a numeric vector used to represent the univariate time series
#' @param window.size is the window size used to create subsequences. It is also the length of univariate motifs
#' @param motif.indices is the results of Func.motif()$Indices, which store the starting position of subsequences for each univaraite motifs
#' 
#' @return The function returns a list of three elements. The first element is data.1, which can be used to show the whole time series with motifs identifed highlighted. The second element is data.2, which can be used to visualize the members of each motif. It is a list containing data frames. Each data frame is designed to visualize the members in each motif.
#' 
#' @examples
#' data(test)
#' #Perform univariate motif discovery for the first dimension data in the example data
#' res.1 <- Func.motif(ts = test$TS1, global.norm = TRUE, local.norm = FALSE, 
#' window.size = 10, overlap = 0, w = 5, a = 3, mask.size = 3, eps = .01)
#' data.vis <- Func.visual.SingleMotif(single.ts=test$TS1, window.size=10, motif.indices=res.1$Indices)
#' #To visualize general information of motifs discovered on the whole time series
#' library(ggplot2)
#' ggplot(data = data.vis$data.1) + 
#'  geom_line(aes(x = 1:dim(data.vis$data.1)[1], y = X)) + 
#'  geom_point(aes(x = 1:dim(data.vis$data.1)[1], y = X, color=Y))
#' #To visualize the detailed information of the 1st motif
#' ggplot(data = data.vis$data.2[[1]]) + geom_line(aes(x = Time, y = Value, linetype=Instance))
#' 
#' @export 

Func.visual.SingleMotif <- function(single.ts, window.size, motif.indices) {
    i=j=NULL
    pos <- foreach::foreach(i=1:length(motif.indices)) %do% {
        pos.ind <- foreach::foreach(j=1:length(motif.indices[[i]]), .combine = c) %do% {
            motif.indices[[i]][j]:(motif.indices[[i]][j]+window.size-1)
        }
        return(pos.ind)
    }
    temp.1 <- data.frame(X=single.ts, Y="Ref", stringsAsFactors = F)
    for(m in 1:length(pos)) {
        temp.1[pos[[m]],"Y"] <- paste0("Motif.",m)
    }
    temp.1$Y <- factor(temp.1$Y)
    
    temp.2 <- foreach::foreach(i=1:length(motif.indices)) %do% {
        pos.ind <- foreach::foreach(j=1:length(motif.indices[[i]]), .combine = rbind) %do% {
            single.ts[motif.indices[[i]][j]:(motif.indices[[i]][j]+window.size-1)]
        }
        pos.ind <- cbind(pos.ind, ID=1:dim(pos.ind)[1])
        
        pos.trans <- reshape2::melt(data = pos.ind, id.vars = "ID")
        pos.trans <- pos.trans[-which(pos.trans$Var2=="ID"),]
        pos.trans$Var2 <- as.numeric(rep(1:window.size, each = dim(pos.ind)[1]))
        pos.trans$Var1 <- as.numeric(pos.trans$Var1)
        pos.trans <- pos.trans[with(pos.trans, order(Var1, decreasing=F)),]
        colnames(pos.trans) <- c("Instance","Time","Value")
        pos.trans[,1] <- factor(pos.trans[,1])
        
        return(pos.trans)
    }
    
    return(list(data.1=temp.1, data.2=temp.2))
}

#' A function to prepare the data for the visualization of multivariate motifs discovered
#' 
#' This function prepares the data used for visualizing multivariate motifs. 
#' 
#' @import foreach
#' @import plyr
#' @import ggplot2
#' 
#' @param data is a data frame containing the multivariate time series data. Each column represents a time series.
#' @param multi.motifs is the result of Func.motif.multivariate 
#' @param index is an integer which specifies the No. of multivariate motif to be plotted
#' 
#' @return The function returns a data frame for the ease of visualizing multivariate motif discovered 
#' 
#' @examples
#' data(test)
#' #Perform univariate motif discovery
#' res.1 <- Func.motif(ts = test$TS1, global.norm = TRUE, local.norm = FALSE, 
#' window.size = 10, overlap = 0, w = 5, a = 3, mask.size = 3, eps = .01)
#' res.2 <- Func.motif(ts = test$TS2, global.norm = TRUE, local.norm = FALSE, 
#' window.size = 20, overlap = 0, w = 5, a = 3, mask.size = 3, eps = .01)
#' res.multi <- Func.motif.multivariate(motif.list = list(res.1$Indices, res.2$Indices), 
#' window.sizes = c(10,20), alpha = .8)
#' #Use the function to prepare the data frame for visualizing the first multivariate motifs identified
#' data.multi <- Func.visual.MultiMotif(data = test, multi.motifs = res.multi, index = 1)
#' #Make the plot using ggplot2
#' library(ggplot2)
#' ggplot(data = data.multi) + 
#'  geom_line(aes(x = T, y = X)) +
#'  geom_point(aes(x = T, y = X, col=Lab, shape=Lab)) + facet_grid(Facet~.)
#' 
#' @export

Func.visual.MultiMotif <- function(data, multi.motifs, index) {
    i=q=j=NULL
    data.trans <- foreach::foreach(q = 1:dim(data)[2], .combine = c) %do% {
        as.numeric(data[,q])
    }
    data.trans <- cbind.data.frame(X=data.trans, T=rep(1:dim(data)[1], times = dim(data)[2]), 
                                   Facet=rep(1:dim(data)[2], each = dim(data)[1]))
    data.trans$Facet <- factor(data.trans$Facet)
    
    con <- multi.motifs$Info[[index]]
    subs <- foreach::foreach(i = 1:length(unique(con$Variable)), .combine = rbind.data.frame) %do% {
        con.sub <- con[which(con$Variable==unique(con$Variable)[i]),]
        xs <- foreach::foreach(j = 1:dim(con.sub)[1], .combine = c) %do% {
            c(con.sub[j,"StartP"]:con.sub[j,"EndP"])   
        }
        
        xs.upd <- data.frame(X=xs, Lab=rep(con.sub$Lab, each=con.sub[1,"EndP"]-con.sub[1,"StartP"]+1), Var.no=unique(con.sub$Variable))
        return(xs.upd)
    }
    
    subs$Index.adj <- (subs$Var.no-1)*dim(data)[1]+subs$X
    data.trans$Index.adj <- 1:dim(data.trans)[1]
    
    data.upd <- plyr::join(x = data.trans, y = subs[,-1])
    data.upd[is.na(data.upd$Lab),"Lab"] <- 0
    data.upd[is.na(data.upd$Var.no),"Var.no"] <- as.numeric(as.character(data.upd[is.na(data.upd$Var.no),"Facet"]))
    data.upd$Lab <- factor(data.upd$Lab)
    data.upd$Var.no <- NULL
    
    return(data.upd)
}

#' An example data set for univariate motif discovery
#' 
#' The data is a data frame containing 100 observations and 2 variables. 
#' The first time series is denoted as TS1. It is created in such a way that two motifs are embedded, 
#' each with two appearances and a length of 10. The rest are randomly generated.
#' The second time series is denoted as TS2. It is created in such a way that one motif with three appearances are embedded. 
#' It has a length of 20. The rest are randomly generated.
#' This synthetic data set is used as examples for motif discovery
#' @docType data
#' 
#' @usage data(test)
#' 
#' @format A data frame with 100 rows and 2 variables
#' 
#' @examples
#' library(ggplot2)
#' data(test)
#' ggplot(data = test, aes(x = 1:dim(test)[1], y = TS1)) + geom_line() + geom_point()
#' ggplot(data = test, aes(x = 1:dim(test)[1], y = TS2)) + geom_line() + geom_point()

"test"

#' An example 1-week data set containing the power consumption data of two building services sub-systems
#' 
#' The data contain 7-day power consumption data of two building services sub-systems, i.e., water-cooled chillers (WCC) and air handling unit (AHU)
#' The collection interval is 15-minute. The data have 672 observations and 6 variables, i.e., Month, Day, Hour, Minute, WCC and AHU
#' @docType data
#' 
#' @usage data(BuildOperation)
#' 
#' @format A data frame with 672 rows and 6 variables
#' 
#' @examples
#' library(ggplot2)
#' data(BuildOperation)
#' ggplot(data = BuildOperation, aes(x = 1:dim(BuildOperation)[1], y = WCC)) + 
#'  geom_line() + geom_point()
#' ggplot(data = BuildOperation, aes(x = 1:dim(BuildOperation)[1], y = AHU)) + 
#'  geom_line() + geom_point()

"BuildOperation"
