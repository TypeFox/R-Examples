# <span style="font-weight:bold; color:#F87217; text-decoration:underline">To answer this question, without loss of generality we formulate it as an often-encountered problem: from the whole human gene network, how to derive a subnetwork that revolves around a gene 'NANOG' and its direct interacting neighbors.</span>

# 1) load the whole human gene network. The network is extracted from the STRING database (version 9.1). Only those associations with medium confidence (score>=400) are retained
network <- dRDataLoader(RData = "org.Hs.string")
# If required, can only focus on those edges with confidence score >= 900
network <- igraph::subgraph.edges(network, eids=E(network)[combined_score>=900])
network
# 2) find the gene 'NANOG' and its direct neighbors
ind <- match(V(network)$symbol, 'NANOG')
nodes_query <- V(network)$name[!is.na(ind)]
g <- dNetInduce(g=network, nodes_query=nodes_query, knn=1, remove.loops=F, largest.comp=T)
g
# 3) visualise the extracted subnetwork
visNet(g, vertex.label=V(g)$symbol, vertex.shape="sphere")
# 4) highlight the query gene in red, and its neighours in black
vcolors <- rep('black', vcount(g))
vcolors[match('NANOG', V(g)$symbol)] <- 'red'
visNet(g, vertex.label=V(g)$symbol, vertex.shape="sphere", vertex.color=vcolors)
# 5) output the detailed information for genes in the subnetwork onto the file called 'NANOG.subnetwork.txt'
out <- data.frame(GeneID=V(g)$geneid, Symbol=V(g)$symbol, Description=V(g)$description)
write.table(out, file="NANOG.subnetwork.txt", quote=F, row.names=F,col.names=T,sep="\t")
out

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">A more complicated question could be: the interaction neighbors are also required to be functionally related to 'stem'.</span> For this, we can further use terms from Gene Ontology (and its subontology: Biological Process).
keyword <- 'stem'
org.Hs.egGOBP <- dRDataLoader(RData = "org.Hs.egGOBP")
ind <- grep(keyword, org.Hs.egGOBP$set_info$name, ignore.case=T, perl=T)
genes <- unique(unlist(org.Hs.egGOBP$gs[ind]))

# Next, identify genes in the subnetwork that are functionally related to 'stem'
ind <- match(V(g)$geneid, genes)
intersect_genes <- V(g)$symbol[!is.na(ind)]
intersect_genes

# Now, highlight the query gene in red, its 'stem'-related neighours in blue, and the rest in black
vcolors <- rep('black', vcount(g))
vcolors[match(intersect_genes, V(g)$symbol)] <- 'blue'
vcolors[match('NANOG', V(g)$symbol)] <- 'red'
visNet(g, vertex.label=V(g)$symbol, vertex.shape="sphere", vertex.color=vcolors)

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">Alternatively, remain only genes related to 'stem' and the query gene 'NANOG'.</span>
both_genes <- c(intersect_genes, 'NANOG')
both_genes

# extract the subnetwork containing the gene 'CDKN2D' and its differentiation-related neighbors
ind <- match(V(g)$symbol, both_genes)
nodes_query <- V(g)$name[!is.na(ind)]
g1 <- dNetInduce(g, nodes_query=nodes_query, knn=0, remove.loops=F, largest.comp=T)
vcolors <- rep('blue', vcount(g1))
vcolors[match('NANOG', V(g1)$symbol)] <- 'red'
visNet(g1, vertex.label=V(g1)$symbol, vertex.shape="sphere", vertex.color=vcolors)

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">Also, only edges linked to the query gene 'NANOG' are remained.</span>
ind <- match('NANOG', V(g1)$symbol)
eids <- E(g1)[from(ind)]
g2 <- subgraph.edges(g1, eids=eids)
vcolors <- rep('blue', vcount(g2))
vcolors[match('NANOG', V(g2)$symbol)] <- 'red'
visNet(g2, vertex.label=V(g2)$symbol, vertex.shape="sphere", vertex.color=vcolors)

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">In addition to the direct interacting neighbors, the indirect neighbors can be specifed via the parameter 'knn' in the function dNetInduce.</span>
