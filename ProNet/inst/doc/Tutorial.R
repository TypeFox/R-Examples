### R code from vignette source 'Tutorial.Rnw'

###################################################
### code chunk number 1: Tutorial.Rnw:32-34
###################################################
options(width=100)

###################################################
### code chunk number 2: Tutorial.Rnw:51-53
###################################################
install.packages("ProNet")

###################################################
### code chunk number 3: Tutorial.Rnw:57-59
###################################################
library("ProNet")

###################################################
### code chunk number 4: Tutorial.Rnw:63-66
###################################################
nodes<-data.frame(c("1855","1856","1857"))
network<-construction(input=nodes,db="Biogrid",species="human",ID.type="Entrez Gene",hierarchy=1)

###################################################
### code chunk number 5: Tutorial.Rnw:70-74
###################################################
net1<-extraction(network, mode="sample", sample.number=20)
net2<-extraction(network, mode="exact", nodes=1:20)
net3<-assemble(net1, net2, mode="union")

###################################################
### code chunk number 6: Tutorial.Rnw:78-80
###################################################
visualization(network, layout="fruchterman.reingold",node.size=8,node.fill.color="red",node.border.color="red",node.label.color="yellow")

###################################################
### code chunk number 7: Tutorial.Rnw:84-87
###################################################
topology(network, simple.parameters=TRUE, degree.distribution=TRUE,clustering.coefficient=TRUE)
net.comparing(net1,net2,topology.parameters=TRUE)

###################################################
### code chunk number 8: Tutorial.Rnw:90-92
###################################################
cluster(network, method="MCODE", plot=TRUE, layout="fruchterman.reingold")

###################################################
### code chunk number 9: Tutorial.Rnw:95-98
###################################################
enrichment.annotation(network, onto="MF", pvalue=0.05)
go.profiles(V(net1)$name, V(net2)$name, onto="MF", mode="frequency", plot=TRUE)

###################################################
### code chunk number 10: Tutorial.Rnw:110-115
###################################################
library("ProNet")
iavPath <-file.path(system.file("example",package="ProNet"),"iav.txt")
iav <- read.table(iavPath, header=TRUE, sep="\t")
head(iav)

###################################################
### code chunk number 11: Tutorial.Rnw:118-124
###################################################
g1 <- construction(iav[,c("Gene_name_1","Gene_name_2")],local.net=TRUE)
sp <- unique(cbind(c(as.vector(iav[,"Gene_name_1"]),as.vector(iav[,"Gene_name_2"])),
                   c(as.vector(iav[,"Adscription_1"]),as.vector(iav[,"Adscription_2"]))))
V(g1)$species <- sp[,2]
summary(g1)

###################################################
### code chunk number 12: Tutorial.Rnw:128-134
###################################################  
hostPath <-file.path(system.file("example",package="ProNet"),"host.txt")
host <- read.table(hostPath, header=TRUE, sep="\t")
g2 <- construction(input=as.data.frame(unique(host[,"Protein.name"])),
                   hierarchy=1,db="HPRD",species="human",ID.type="Gene symbol")
summary(g2)

###################################################
### code chunk number 13: Tutorial.Rnw:138-144
###################################################  
hprd <- construction(db="HPRD",ID.type= c("Gene symbol"))
id <- match(unique(c(V(g1)$name,V(g2)$name)),V(hprd)$name)
gtemp <- induced.subgraph(hprd, id[!is.na(id)])
g3 <- assemble(g1,gtemp,mode="union")
summary(g3)

###################################################
### code chunk number 14: Tutorial.Rnw:149-157
###################################################  
color <- rep(1,vcount(g3))
color[V(g3)$species=="DHP of IAV"] <- "red"
color[V(g3)$species=="IAV protein"] <- "black"
color[is.na(V(g3)$species)] <- "green"
visualization(g3,node.size=3,node.fill.color=color,node.label="",edge.color="gray")
legend("topleft",col=c("black","red","green"),
       legend=c("virus","human_direct","human_indirect"),pch=19)
  
###################################################
### code chunk number 15: Tutorial.Rnw:168-172
###################################################
V(g3)$expression<-rexp(vcount(g3),1)
location(g3,species=c("human"),vertex.size=3,vertex.label.cex=0.5,
         vertex.color="expression",xlim=c(-1,1),ylim=c(-1,1))

###################################################
### code chunk number 16: Tutorial.Rnw:185-187
###################################################
topology(g3,simple.parameters=TRUE)

###################################################
### code chunk number 17: Tutorial.Rnw:191-194
###################################################
tp <- topology(g2,degree.distribution=TRUE)
head(as.data.frame(tp))

###################################################
### code chunk number 18: Tutorial.Rnw:198-201
###################################################
tp <- topology(g2,shortest.paths=TRUE)
head(as.data.frame(tp))

###################################################
### code chunk number 19: Tutorial.Rnw:216-218
###################################################
net.comparing(g3,hprd,topology.parameters=TRUE)

###################################################
### code chunk number 20: Tutorial.Rnw:222-224
###################################################
net.comparing(g3,hprd,topology.parameters=FALSE,degree=TRUE)

###################################################
### code chunk number 21: Tutorial.Rnw:235-237
###################################################
net.comparing(g3,hprd,topology.parameters=FALSE,degree=TRUE)

###################################################
### code chunk number 22: Tutorial.Rnw:241-243
###################################################
comp.rand.subnet(g3,hprd,nsim=10000,ave.path.len=TRUE)

###################################################
### code chunk number 23: Tutorial.Rnw:260-266
###################################################
result <- cluster(g3, method="FN")
clusters <- rep(1, vcount(g3))
for(i in 1:vcount(g3)){clusters[i] <- result[[i]]}
clusters <- as.factor(clusters)
table(clusters)

###################################################
### code chunk number 24: Tutorial.Rnw:269-273
###################################################
result <- mcode(g3,vwp=0.05,haircut=TRUE,fluff=FALSE,fdt=0.8,loops=FALSE)
summary(result$COMPLEX)
result$score

###################################################
### code chunk number 25: Tutorial.Rnw:277-281
###################################################
cluster1<-induced.subgraph(g3,result$COMPLEX[[1]])
summary(cluster1)
visualization(cluster1,node.size=4,node.label=V(cluster1)$name,node.label.color="blue")

###################################################
### code chunk number 26: Tutorial.Rnw:294-299
###################################################
idPath <-file.path(system.file("example",package="ProNet"),"hprd.id.txt")
id <- read.table(idPath, header=FALSE, sep="\t")
colnames(id) <- c("hprd_id","geneSymbol","nucleotide_accession","protein_accession","entrezgene_id","omim_id","swissprot_id","main_name")
head(id)

###################################################
### code chunk number 27: Tutorial.Rnw:302-307
###################################################
index1 <- match(V(cluster1)$name, as.vector(id$geneSymbol), nomatch=0)
entrez1 <- as.vector(id$entrezgene_id[index1])
go.mf <- enrichment.annotation(entrez1, onto="MF", pvalue=0.05)
head(go.mf[,c("GO_ID","GO_term","Evidence","p.value")])

###################################################
### code chunk number 28: Tutorial.Rnw:311-313
###################################################
go.profiles(entrez1, onto="MF",main="cluster1")

###################################################
### code chunk number 29: Tutorial.Rnw:324-329
###################################################
cluster2<-induced.subgraph(g3,result$COMPLEX[[2]])
index2 <- match(V(cluster2)$name, as.vector(id$geneSymbol), nomatch=0)
entrez2 <- as.vector(id$entrezgene_id[index2])
go.profiles(entrez1,entrez2,onto="MF",main=c("cluster1 vs 2"))
