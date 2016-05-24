                                        # Algorithm for j0 estimation based on the method of Hall and Schimek (2012)
                                        # The algorithm consists of an ordered sequence of "test stages" s1, s2, ....
                                        # Stage sk is associated with an integer Jsk, which when k is odd, is a potential lower bound to j0
                                        # Idata - input data is a vector of 0's and 1's

calculate.maxK <- function(lists, L, d, v, threshold=50) {
  compared.lists <- list() #contains all pairwise compared lists (structure for aggmap)
  info <- matrix(ncol = 0, nrow = 3) #contains information about list names
  rownames(info) <- c("listname", "original listname", "ref-list or trunc-list")
  grayL <- list() #contains information which object in the list has to be gray-shaded
  grayO <- c() #contains all gray-shaded objects
  temp.sumtrunclists <- list() #contains the summarized truncated lists (number of lists = L)
  summarytable <- matrix(nrow = 0, ncol = (3 + L)) #contains the summary-table
  venn.values <- list() #contains the Venn-lists for the Venn-diagram and the Venn-table (a Venn-diagram in table form)
  cg_temp = c()
  
                                        #first step: estimate the j_0 (and so k)
  res.j0.temp <- j0.multi(lists, d, v)
  res.temp <- as.matrix(res.j0.temp$L)
                                        #print(res.temp)
  maxK <- res.j0.temp$maxK
  
  if(sum(is.na(res.temp[,5]))<nrow(res.temp)){
    tl = as.list(lists[1:maxK, ])
    tli = as.character(unique(unlist(tl)))
	sp = rep(list(tli),L)
    resS = CEMC(input = tl, space = sp, k = maxK)
    temp = tapply(as.numeric(res.temp[, 5]), res.temp[, 1], function(x) max(x, na.rm = TRUE))
    
    if (sum(temp!="-Inf")>1){
      if (sum(is.na(res.temp[,5]))>0){
        res.temp2 = res.temp[-which(is.na(res.temp[,5])),]
      } else
        {res.temp2 = res.temp}
    } else {
      if (sum(is.na(res.temp[,5]))>0)
        {res.temp2 = t(as.matrix(res.temp[-which(is.na(res.temp[,5])),]))
       } else {
         res.temp2 = t(as.matrix(res.temp))
       }
    }
    
    temp2 = tapply(as.numeric(res.temp2[,5]), res.temp2[,1], function(x) max(x, na.rm = TRUE))
    
    list_t = list()
    for (i in unique(res.temp2[,1])){
      res.temp.temp = res.temp2[res.temp2[,1]==i,]
      if(!is.matrix(res.temp.temp)){
        res.temp.temp = t(as.matrix(res.temp.temp))
      }
      list_t[[i]] = cbind(res.temp.temp[,2][order(res.temp.temp[,5], decreasing=T)], res.temp.temp[,5][order(res.temp.temp[,5],decreasing=T)])
    }
    
                                        #Calculating block order(bo)
    
    bo = names(sort(unlist(lapply(list_t, FUN=function(x) max(as.numeric(x[,2]), na.rm=T))), decreasing=T))
    
    ilor = lapply(list_t, FUN=function(x) x[,1][order(as.numeric(x[,2]), decreasing=T)]) #inblock list order

    if (length(bo)>1)
      {
        ilor_final = list() # inblock list order final 
        ilor_final[[bo[1]]] = ilor[[bo[1]]]
        bo_temp = c()
        
        for (i in c(2:length(bo)))
          {
            bo_temp = c(bo_temp, bo[[i-1]]);
            ilor_final[[bo[i]]] = setdiff(ilor[[bo[i]]], bo_temp)
          }
      }else{ilor_final = ilor}
    
    
    bo_final = names(ilor_final)[which(lapply(ilor_final, length)!=0)]
    
##### condition if maximal estimated j0 is NA, then warning is returned #####

    if (maxK!="-Inf")
      {	
        
###### building plotflow#####
        
        crl <- 0 #current reference list
                                        #iterate over all blocks (a block is a reference list with the corresponding truncation lists)
        
                                        # fln - first list name
        for (fln in bo_final) {  
                                        # selecting block
          temp2 = list_t[[fln]]
          rownames(temp2) = temp2[,1]

                                        #get the objects of the current reference list
          gnp <- as.vector(lists[, fln][1:maxK]) # gene names to plot
                                        #gnp <- as.vector(gnp)
          
          crl <- crl + 1
          compared.lists[[paste("R", crl, sep = "")]] <- gnp

          temp.sumtrunclists[[fln]] <- gnp
          info <- cbind(info, c(paste("R", crl, sep = ""), fln, "R"))

          ctr <- 0 #current truncated lists
          grayL[[paste("R", crl, sep = "")]] <- rep(FALSE, length(gnp))		
          temp.countgray <- matrix(ncol = 2, nrow = length(gnp), data = 0)
          
                                        #iterate over the truncated lists of the current block
          for (l in c(1:length(ilor_final[[fln]]))) {
            n.genes.to.plot <- as.numeric(temp2[ilor_final[[fln]][l],2])
            
            temp.distances = c(1:length(gnp)) - match(gnp, as.character(lists[, ilor_final[[fln]][l]]))
            temp.countgray[,2] = temp.countgray[,2]+c(abs(temp.distances)<=d)

            ##temp.sumtrunclists[[ilor_final[[fln]][l]]] <- lists[, ilor_final[[fln]][l]][1:(as.numeric(temp2[ilor_final[[fln]][l],2]))]
            
            temp.sumtrunclists[[ilor_final[[fln]][l]]] <- lists[, ilor_final[[fln]][l]][1:length(resS$TopK)]
            ##check for gray-shade of an object in the truncated list
            temp.grayshade = abs(temp.distances)<=d
            ##add the truncated list
            ctr <- ctr + 1
            compared.lists[[paste("R", crl, "_T", ctr, sep = "")]] <- temp.distances
            grayL[[paste("R", crl, "_T", ctr, sep = "")]] <- temp.grayshade
            info <- cbind(info, c(paste("R", crl, "_T", ctr, sep = ""), ilor_final[[fln]][l], "T"))
            
          }# end for l
          
                                        #calculate if respective object of the reference list has to be gray-shaded
          temp.percentage = apply(as.matrix(temp.countgray[,-1]),1,sum)/c(length(ilor_final[[fln]]))*100
          grayL[[paste("R", crl, sep = "")]] = temp.percentage >= threshold
                                        #add object to a new list which contains all gray-shaded objects (add only if it is not already in the list)
          grayO = union(grayO, lists[, fln][which(temp.percentage >= threshold)])
        }# end for fln
        


                                        #having all the necessary information, calculate the summary-table
        colnames(summarytable) <- c(names(lists), "Rank sum", "Freq in input lists", "Freq in truncated lists")

        for (j in 1:length(tli)) {
          cg <- tli[j] #current genesymbol
          
                                        #get the positions of the current object in the input lists
          temp.positions <- rep(NA,L)
          for (q in 1:L) {
            temp.positions[q] <- match(cg, lists[,names(lists)[q]])
          }#end for q
          
#### positions = apply(lists,2,FUN=function(x) match(tli,x))
          
          
                                        #calculate the rank sum
          temp.ranksum <- 0
          temp.missingvalues <- 0
          for (q in 1:L) {
            if (is.na(temp.positions[q])) {
              temp.missingvalues <- temp.missingvalues + 1
            }
          }#end for q
                                        #if one or more rank positions are unavaliable (e.g. an object does not exist in a list), interpolate the rank sum
          if (temp.missingvalues > 0) {
            temp.meanrank <- 0
            temp.partialranksum <- 0
                                        #get the rank sum of the valid rank-positions
            for (q in 1:L) {
              if (!is.na(temp.positions[q])) {
                temp.partialranksum <- temp.partialranksum + temp.positions[q]
              }
            }#end for q
                                        #calculate the mean rank-position for the unavailable rank-positions
            temp.meanrank <- round(temp.partialranksum / (L - temp.missingvalues))
            temp.ranksum <- temp.partialranksum + (temp.meanrank * temp.missingvalues)
          } else {
            temp.ranksum <- sum(temp.positions)
          }# end if temp.missingvalues
          
                                        #calculate the frequency in the input lists
          temp.freqinput <- L - temp.missingvalues
          
                                        #calculate the frequency in the summarized truncated lists
          temp.freqtrunc <- 0	
          for (curr.listname in names(temp.sumtrunclists)) {
            if (cg %in% temp.sumtrunclists[[curr.listname]]) {
              temp.freqtrunc <- temp.freqtrunc + 1
            }
          }

                                        #add the calculated row (of the current object) to the summary-table
          cg_temp = c(cg_temp, cg)
          summarytable <- rbind(summarytable, c(temp.positions, temp.ranksum, temp.freqinput, temp.freqtrunc))
        }# end for j
        
        ## calculate Venn table

                                        # changing into vectors of characters
        temp.sumtrunclists.vect <- sapply(temp.sumtrunclists,as.vector, simplify=FALSE)
        names(temp.sumtrunclists.vect) <- names(temp.sumtrunclists)
        
                                        #check for each object if entry is present in lists
        venn.table <- data.frame(do.call(cbind, lapply(na.exclude(names(temp.sumtrunclists.vect)), function(nn){
          ifelse(tli %in% as.vector(temp.sumtrunclists.vect[[nn]]),nn,NA)
        })), stringsAsFactors=FALSE)
        rownames(venn.table) <- tli
        venn.table$listname <-apply(venn.table, 1, function(x){paste(sort(x[!is.na(x)]), sep="", collapse="_")})
        venn.temp <- split(rownames(venn.table),venn.table$listname)
                                        # adding stars to those that were in the CEMC final list
        venn.list <- lapply(venn.temp, function(x)  {a = rep("*",length(x))[match(x,resS$TopK)>0]; a[which(is.na(a))]=""; paste(a,x,sep="")})
        venn.list <- venn.list[order(-sapply(names(venn.list), nchar), names(venn.list))]
        venntable <- data.frame(t(sapply(names(venn.list), function(nn){
          data.frame(intersection=nn,objects=paste(sort(venn.list[[nn]]), sep="",collapse=", "), stringsAsFactors=FALSE)
        })))
        rownames(venntable) <- NULL





                                        #conversion of the summary table into a data frame so that the rankings are given as numbers, not as characters, otherwise the ordering is wrong
        
    summarytable.temp = data.frame(Object=cg_temp,summarytable, stringsAsFactors=FALSE)
	TopKranks = data.frame(TopK=resS$TopK, ranks=1:length(resS$TopK))
	CEMCres =  rep("",nrow(summarytable.temp))
	CEMCres[match(TopKranks$TopK,summarytable.temp$Object)] = "YES"
	summarytable.temp2 = data.frame(Final.selection.CEMC = CEMCres, summarytable.temp)
	topk.table = summarytable.temp2[which(summarytable.temp2$Final.selection.CEMC=='YES'),]
	ranks = numeric()
	ranks[match(TopKranks$TopK,topk.table$Object)] = TopKranks$ranks
	topk.table.rank = topk.table[order(ranks),]
	rest.table = summarytable.temp2[which(summarytable.temp2$Final.selection.CEMC==''),]
	rest.table = rest.table[order(rest.table[,L+3]),]
	summarytable.final = rbind(topk.table.rank, rest.table)


        
                                        #calculate the Venn-lists (to view the Venn-diagram) and the Venn-table
                                        #the calculation takes place only for L between 2 and 4 (a Venn-diagram for L > 4 cannot be properly arranged)
        ##combine all necessary objects into one single list
        truncated.lists <- list()
        truncated.lists$comparedLists <- compared.lists
        truncated.lists$info <- info
        truncated.lists$grayshadedLists <- grayL
        truncated.lists$summarytable <- summarytable.final
        truncated.lists$vennlists <- temp.sumtrunclists
        truncated.lists$venntable <- venntable
        truncated.lists$v <- v
        truncated.lists$Ntoplot<-sum(unlist(lapply(ilor_final,length)))+sum(unlist(lapply(ilor_final,length))>0)
        truncated.lists$Idata <- res.j0.temp$Idata
        truncated.lists$d <- d
        truncated.lists$threshold <- threshold
        truncated.lists$L <- L
        truncated.lists$N <- nrow(lists)
        truncated.lists$lists <- lists
        truncated.lists$maxK <-maxK
        truncated.lists$topkspace <-resS$TopK
        return(truncated.lists)
      }
  } else {
    message(paste("!!!...For selected delta, the top L list cannot be estimated (little or no overlap)!!!", "\n"))
    return(truncated.lists=NULL)
  } # end if if (maxK)
}# end of function calculate.maxK

compute.stream<-function(Idata, const=0.251, v, r=1.2)
{
  if(sum(Idata, na.rm=T)==length(Idata)) 
    {
      return(list(j0_est=length( ), reason.break="Idata is identity", Js=NA, v=v))

    }else
      {
        zv 		= .moderate.deviation(const, v)
        Js 		= c()
        k 		= 1
        pj.plus 	= c()
        pj.minus 	= c()
        h		= 0
        reason.break="NA"

###########
        repeat
          {
            if (floor(k/2)<(k/2))
              {
                if (k==1)
                  {	h = 1;
                        j=h-1;
                        v.last=v[k]
                        if (j< length(Idata))
                          {
                            repeat 
                              {
                                j = j+1
                                
                                if ((j+v[k]-1) >= length(Idata)) {#print("(j+v[k]-1) >= length(Idata)")
                                  break}
                                
                                        # computing pjplus
                                pj.plus = .pjplus(Idata, v[k], j)
                                        #print(pj.plus)
                                        #print(Idata[j:(j+v[k]-1)])
                                
                                        # testing
                                if ((pj.plus-0.5)<=zv) {#print("pj.plus-0.5<=zv")
                                  break}
                              }#end repeat

                            Js[k] = j	
                            Js	
                            if(reason.break!="NA"){break;}				

                          }# end if j<length(Idata)

                      } else{
                        h = Js[k-1]-trunc(r*v.last)+1;
                        j = h-1;
                                        #print(paste("j=",j, sep=""))
                        if (j< length(Idata))
                          {
                            repeat 
                              {
                                j = j+1
                                
                                if ((j+v.last-1) >= length(Idata)) {break}
                                
                                        # computing pjplus
                                pj.plus = .pjplus(Idata, v.last, j)
                                
                                        # testing
                                if ((pj.plus-0.5)<=zv) {break}
                              }
                            if(reason.break!="NA"){break;}
                            Js[k] = j	
                            Js					
                                        # breaking the repeat loop condition 2
                            if ((k-3)>=1) {if (Js[k-2]==Js[k] & Js[k-1]==Js[k-3]) {reason.break="Js[i-2]==Js[i] & Js[i-1]==Js[i-3]";
                                                                                   break}
			 		 }
                          }# end if

                      }#end else
              } # end is.odd

            if ((floor(k/2)==(k/2)))
              {
                h = Js[k-1]+trunc(r*v.last)
                j = h

                if (j>1)
                  {
                    repeat
                      {
                        j = j-1

                                        # computing pjminus
                        if ((j-v.last+1) <= 0){break}
                        
                        pj.minus 	= .pjminus(Idata, v.last, j)
			
                                        # testing
                        if ((pj.minus-0.5)>zv) {break}	
                      }# end repeat

                    Js[k] = j
                    
                                        # breaking the repeat loop condition 3
                    if (Js[k]-trunc(r*v.last)+1<=1) {  v.temp=v.last; 
                                                       while((Js[k]-r*v.temp<=1) & v.temp>0){v.temp=v.temp-1}
                                                       v.last=v.temp}
                    if (v.last==1){reason.break = "v converged to 1"; v[k+1]=v.last ;break; }
                    if (v.last==0){reason.break = "v converged to 0"; v[k+1]=v.last ;break; }


                                        # breaking the repeat loop condition 1
                    if (is.infinite(Js[k-1]) & is.infinite(Js[k])) 
                      {reason.break="is.infinite(Js[i-1]) & is.infinite(Js[i])"; break}

                  } # end if
                if(reason.break!="NA"){break;}
              }

            k=k+1
            v[k]=v.last

          }#end repeat
        if(v.last<=1 | reason.break!="Js[i-2]==Js[i] & Js[i-1]==Js[i-3]"){j0_est=NA} else {if ((floor(k/2)<(k/2)) & k>2) {j0_est = ceiling(Js[k-2]+0.5*v[k-2])}else if ((floor(k/2)==(k/2)) & k>1){j0_est = ceiling(Js[k-1]+0.5*v[k-1])} else{j0_est=NA}}
        return(list(j0_est=j0_est+1, k=j0_est, reason.break=reason.break, Js=Js+1, v=v))
      }# end if sum(Idata)==length(Idata)
}

                                        #Method of Hall and Schimek (2012) adapted to the situation of multiple ranked lists 
                                        #Inputs:
                                        #Data - input matrix, each column represents one list
                                        #delta - the maximal distance between the rank positions of an object in a pair of lists
                                        #v - value of tuning parameter nu used for j0 estimation
                                        #Outputs: maximal estimated j0 based on all combinations of 2 lists

                                        #Data = lists
                                        #delta between 0 and 20 suggested for short lists with N<200
                                        #v=10 suggested for short lists with N<200

j0.multi<-function(lists,d,v) {
  if(is.null(colnames(lists))){
    colnames(lists)<-paste("L",1:ncol(lists),sep="")
    warning("No colnames given in lists. Replaced with default values: L1, L2...")
  }
  maxK=0
  L = c()
  Idata_ID = c()
  names_idata = c()
  for (i in 1:ncol(lists)){
    for (j in 1:ncol(lists)){
      if (i!=j) {
	ID = prepare.idata(lists[,c(i,j)],d=d)
   	Idata_ID=cbind(Idata_ID, ID$Idata)
	names_idata = c(names_idata, paste(colnames(lists)[i],"_",colnames(lists)[j],sep=""))
	J = compute.stream(ID$Idata,v=v)$j0_est
	k = compute.stream(ID$Idata,v=v)$k
	L = rbind(L, cbind(colnames(lists)[i], colnames(lists)[j], v,J,k, d))
      }# end for if
    }# end for j
  }# end for i
  colnames(Idata_ID) = names_idata
  L = data.frame(L, stringsAsFactors=F)
  names(L) = c("list1", "list2", "v", "j0_est","k","delta")
  if(sum(is.na(L$k))<nrow(L)){maxK = max(as.numeric(L$k), na.rm=T)}else{maxK=NA}
  return(list(maxK=maxK,L=L, Idata=Idata_ID))
}

                                        # Zv moderate deviation

.moderate.deviation <- function(C, v)
{
  zv = sqrt((C*log(v, base=10))/v)
  return(zv)
}

.pjminus <- function(Idata, v, j)
{
  pj.minus = 1/v*(sum(Idata[(j-v+1):j], na.rm = TRUE))
  return(pj.minus)
}

.pjplus <- function(Idata, v, j)
{
  pj.plus = (1/v)*(sum(Idata[j:(j+v-1)], na.rm = TRUE))
  return(pj.plus)
}

                                        # Function to prepare Idata for multiple rankings from several assessors allowing for missing values
                                        # x - data matrix, where columns represent the rank order of objects from two different assessors and the rows represent object names
                                        # delta - the maximal distance between the rank positions of an object in a pair of lists
                                        # num.omit - the maximal number of ommited objects from the analysis 
                                        #
                                        # The result is an object of type "Idata", which is a list containing Idata and the information on the distance delta.

prepare.idata <- function(x, d)
{
  if(ncol(x)<2){
    warning("You need a minimum of two lists to compare. Execution halted.","\n")
  }else{
    if(ncol(x)>2){warning("The data matrix you have submitted contains more than two lists (columns). Only first two will be used.","\n")}
    rank.diff = c(1:nrow(x))-match(x[,1],x[,2])
    Idata = as.numeric(abs(rank.diff)<=d)
    Idata[is.na(Idata)] = 0

    return(list(Idata = Idata, d = d))
  }
}

