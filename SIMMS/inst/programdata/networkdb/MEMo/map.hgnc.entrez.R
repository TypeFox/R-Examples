tmp <- read.table("pathway_based_sub_networks__HGNC.txt", header = F, comment.char = "", sep = "\t", stringsAsFactors = FALSE)

map <- read.table("HGNC_Entrez_map.txt", row.names = 1, header = T, sep = "\t", stringsAsFactors = FALSE)

for (i in 1:nrow(tmp)) {
	if (length(grep("^NAME=", tmp[i, 2])) > 0) {
		cat("\nignoring you: ", i, tmp[i, 2]);
		}
	else {
		tmp[i, 2] <- map[tmp[i, 2], 1];
		tmp[i, 3] <- map[tmp[i, 3], 1];
		}
	}

write.table(
	tmp,
	"pathway_based_sub_networks.txt",
	row.names = F,
	col.names = F,
	sep = "\t",
	quote = F
	);
write.table(
	tmp,
	"pathway_based_sub_networks_all.txt",
	row.names = F,
	col.names = F,
	sep = "\t",
	quote = F
	);