permutation <-
function (nperm = 1000, write.output = TRUE) 
{
Class<-NULL
load("permubiome.RData")
df_norm<-df_norm
classes<-levels(df_norm$Class)
group1<-subset(df_norm, Class == classes[1])
group2<-subset(df_norm, Class == classes[2])
categories<-colnames(df_norm)
size1 <- nrow(group1)
size2 <- nrow(group2)
size <- size1 + size2
pvalue_matrix<-matrix(, nrow = ncol(df_norm)-2, ncol = 5, byrow=T)
colnames(pvalue_matrix)<-c("Category", paste(classes[1], "(median)") , paste(classes[2], "(median)"), "p.value","p.adjust (fdr)")
for (i in 3:(ncol(df_norm))){
category<-categories[i]
diff <- median(group1[,i]) - median(group2[,i])
x <- c(group1[,i], group2[,i])
y <- array(, nperm)
for (j in 1:nperm) {
set <- sample(size, size2, replace = FALSE)
diff_iter <- median(x[set]) - median(x[-set])
y[j] <- diff_iter
ref_score <- (diff-mean(y))/sd(y)
}
if (ref_score > 0){pvalue.i <- 1-pnorm(ref_score)
} else {pvalue.i <- pnorm(ref_score)
}
padjust.i <- p.adjust(pvalue.i, method = "fdr", n = nrow(pvalue_matrix))
    if (pvalue.i <= 1) {print(paste(category, signif(pvalue.i, 4), signif(padjust.i, 4)))
    } else {print(paste(category, "1.000", signif(padjust.i, 4 ) ))
    }
    if (pvalue.i != 0) {pvalue_matrix[(i-2),1] <- category
    }
if (pvalue.i != 0) {pvalue_matrix[(i-2),2] <- round(median(group1[,i]), digits = 0)
    }
if (pvalue.i != 0) {pvalue_matrix[(i-2),3] <- round(median(group2[,i]), digits = 0)
    }
if (pvalue.i <= 1) {pvalue_matrix[(i-2),4] <- format(pvalue.i, digits = 7, scientific=F)
} else {pvalue_matrix[(i-2),2] <- 1.000
}
    if (pvalue.i != 0) {pvalue_matrix[(i-2),5] <- format(padjust.i, digits = 7, scientific=F)
    }
    invisible()
    
}
if (write.output == TRUE){
all <- readline("Do you want to include all fetures in the output? (yes/no) : ")
if (substr(all, 1, 1) == "n"){
select <- as.numeric(readline("Features under what level of significance do you want to retrieve (i.e. 0.2) : "))
significant<-subset(pvalue_matrix, pvalue_matrix[,5] <= select)
ordered<-significant[order(significant[,5]),]
} else {
significant<-subset(pvalue_matrix, pvalue_matrix[,5] <= 1)
ordered<-significant[order(significant[,5]),]
}
write.table(ordered, file="permutation.output", quote=F, row.names=F, sep="\t")
print(paste("Permutation test done and output table printed!"))
} else {
significant<-subset(pvalue_matrix, pvalue_matrix[,4] <= 1)
ordered<-significant[order(significant[,4]),]
ordered
print(paste("Permutation test done!"))
}
}
