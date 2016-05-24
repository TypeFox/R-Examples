subs.delete <- function(part1.list, subj, group) {
	if(length(subj) != length(group)) stop("subj and group should be the same length")
	if(any(group != 1 & group != 2)) stop("all entries in group should be a 1 or 2")

	subs.1 <- subj[group == 1]
	subs.2 <- subj[group == 2]
	
	group.1 <- as.character(part1.list$groups[1])
	group.2 <- as.character(part1.list$groups[2])
	
	if(length(subs.1) != 0) {
		part1.list$data <- subset(part1.list$data, !(part1.list$data$Group == group.1 & part1.list$data$Subject %in% subs.1))
		part1.list$N.sub1 <- part1.list$N.sub1 - length(subs.1)
		part1.list$N.g1 <- part1.list$N.g1 - length(subs.1)
		part1.list$coef.id1 <- part1.list$coef.id1[-subs.1,]
		part1.list$coef.id3 <- part1.list$coef.id3[-subs.1,]
		part1.list$sdev.id1 <- part1.list$sdev.id1[-subs.1,]
		part1.list$sdev.id3 <- part1.list$sdev.id3[-subs.1,]
		part1.list$sigma.id1 <- part1.list$sigma.id1[-subs.1,]
		part1.list$sigma.id3 <- part1.list$sigma.id3[-subs.1,]
		part1.list$id.nums.g1 <- part1.list$id.nums.g1[-subs.1]
		part1.list$cor.1 <- part1.list$cor.1[-subs.1]
		part1.list$cor.3 <- part1.list$cor.3[-subs.1]
		part1.list$R2.g1.1 <- part1.list$R2.g1.1[-subs.1]
		part1.list$R2.g1.2 <- part1.list$R2.g1.2[-subs.1]
	}
	if(length(subs.2) != 0) {
		part1.list$data <- subset(part1.list$data, !(part1.list$data$Group == group.2 & part1.list$data$Subject %in% subs.2))
		part1.list$N.sub2 <- part1.list$N.sub2 - length(subs.2)
		part1.list$N.g2 <- part1.list$N.g2 - length(subs.2)
		part1.list$coef.id2 <- part1.list$coef.id2[-subs.2,]
		part1.list$coef.id4 <- part1.list$coef.id4[-subs.2,]
		part1.list$sdev.id2 <- part1.list$sdev.id2[-subs.2,]
		part1.list$sdev.id4 <- part1.list$sdev.id4[-subs.2,]
		part1.list$sigma.id2 <- part1.list$sigma.id2[-subs.2,]
		part1.list$sigma.id4 <- part1.list$sigma.id4[-subs.2,]
		part1.list$id.nums.g2 <- part1.list$id.nums.g2[-subs.2]
		part1.list$cor.2 <- part1.list$cor.2[-subs.2]
		part1.list$cor.4 <- part1.list$cor.4[-subs.2]
		part1.list$R2.g2.1 <- part1.list$R2.g2.1[-subs.2]
		part1.list$R2.g2.2 <- part1.list$R2.g2.2[-subs.2]
	}
	
	part1.list
}
