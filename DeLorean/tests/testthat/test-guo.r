library(DeLorean)
data(GuoDeLorean)
dl <- de.lorean(
    guo.expr,
    guo.gene.meta,
    guo.cell.meta)
dl <- estimate.hyper(dl, sigma.tau=0.5)
set.seed(1)
num.genes <- nrow(dl$gene.meta)  # All genes
sampled.genes <- sample_n(dl$gene.meta, num.genes)$gene
gene.filter <- function(genes) genes %in% sampled.genes
dl <- filter.genes(dl, gene.filter)
num.at.each.stage <- 9
te.sampled.cells <- (
    dl$cell.meta
    %>% filter(capture < "32C" | "TE" == cell.type)
    %>% group_by(capture)
    %>% do(sample_n(., num.at.each.stage))
)
pe.sampled.cells <- (
    dl$cell.meta
    %>% filter(capture < "32C" | "PE" == cell.type | "ICM" == cell.type)
    %>% group_by(capture)
    %>% do(sample_n(., num.at.each.stage))
)
epi.sampled.cells <- (
    dl$cell.meta
    %>% filter(capture < "32C" | "EPI" == cell.type | "ICM" == cell.type)
    %>% group_by(capture)
    %>% do(sample_n(., num.at.each.stage))
)
run.model <- function(dl, cells.sampled) {
    cell.filter <- function(cells) cells %in% cells.sampled
    dl <- filter.cells(dl, cell.filter)
    dl <- prepare.for.stan(dl)
    dl <- compile.model(dl)
    dl <- find.best.tau(dl)
    system.time(dl <- fit.model(dl, num.cores=20))
    dl <- examine.convergence(dl)
    dl <- process.posterior(dl)
    dl <- analyse.noise.levels(dl)
    dl <- make.predictions(dl)
    dl
}
#dl.te  <- run.model(dl,  te.sampled.cells$cell)
#dl.pe  <- run.model(dl,  pe.sampled.cells$cell)
#dl.epi <- run.model(dl, epi.sampled.cells$cell)
#gp <- cmp.profiles.plot(TE=dl.te, PE=dl.pe, EPI=dl.epi,
                        #genes=dl.te$gene.map$gene)
#print(gp)
sessionInfo()
