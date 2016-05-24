require(kinship2)


data(sample.ped)


df2 <- sample.ped[sample.ped$ped==2,]
relate2 <- matrix(c(210,211,1,
                   212,213,3), nrow=2, byrow=TRUE)

ped2 <- pedigree(df2$id, df2$father, df2$mother, 
       df2$sex,
       affected=cbind(df2$affected, df2$avail), 
       relation=relate2)




unavail2 <- findUnavailable(ped2, avail=ped2$affected[,2])
unavail2
## 205, 210, 213

ped2.trim.unavail <- pedigree.trim(unavail2, ped2)
ped2.sub <- ped2[-unavail2]
ped2df1 <- as.data.frame.pedigree(ped2.trim.unavail)
ped2df2 <- as.data.frame.pedigree(ped2.sub)

ped2.trim210 <- pedigree.trim(210, ped2)
ped2.sub210 <- ped2[-10]

ped2.trim210$id
ped2.trim210$relation

ped2.sub210$id
ped2.sub210$relation
##  indx1 indx2    code
##2    11    12 UZ twin

as.data.frame(ped2.trim210)
#    id dadid momid    sex affected.1 affected.2
#1  201     0     0   male          1          1
#2  202     0     0 female         NA          0
#3  203     0     0   male          1          1
#4  204   201   202 female          0          1
#5  205   201   202   male         NA          0
#6  206   201   202 female          1          1
#7  207   201   202 female          1          1
#8  208   201   202 female          0          0
#9  209     0     0   male          0          0
#10 211   203   204   male          0          1
#11 212   209   208 female          0          1
#12 213   209   208   male          0          0
#13 214   209   208   male          1          1



df1<- sample.ped[sample.ped$ped==1,]

ped1 <- pedigree(df1$id, df1$father, df1$mother, 
       df1$sex,
       affected=cbind(df1$affected, df1$avail))


unavail1 <- findUnavailable(ped1, avail=df1$avail)
unavail1 ## 101 102 107 108 111 121 122 123 131 132 134 139
ped1.avail <- pedigree.trim(unavail1, ped1)

as.data.frame(ped1.avail)
#    id dadid momid    sex affected.1 affected.2
#1  103   135   136   male          1          0
#2  104     0     0 female          0          0
#3  105     0     0   male         NA          0
#4  106     0     0 female         NA          0
#5  109     0     0 female          0          1
#6  110   103   104   male          1          1
#7  112   103   104   male          1          0
#8  113   103   104 female          0          1
#9  114   103   104   male          1          0
#10 115   105   106 female          0          0
#11 116   105   106 female          1          1
#12 117     0     0   male          1          0
#13 118   105   106 female          1          1
#14 119   105   106   male          1          1
#15 120     0     0 female          0          0
#16 124   110   109   male          1          1
#17 125   112   118 female          0          1
#18 126   112   118 female          0          1
#19 127   114   115   male          1          1
#20 128   114   115   male          1          1
#21 129   117   116   male          0          1
#22 130   119   120   male          0          1
#23 133   119   120 female          0          1
#24 135     0     0   male         NA          0
#25 136     0     0 female         NA          0
#26 137     0     0   male         NA          0
#27 138   135   136 female         NA          0
#28 140   137   138 female          0          1
#29 141   137   138 female          0          1

