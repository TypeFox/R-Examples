context("phyndr_taxonomy")

## This is how the data were generated:
if (FALSE) {
  families <- c("Araucariaceae", "Cephalotaxaceae", "Cupressaceae",
                "Pinaceae", "Podocarpaceae", "Taxaceae")
  taxize::tpl_get("tpl", families)

  files <- dir("tpl", full.names=TRUE)
  dat <- do.call("rbind", lapply(files, read.csv, stringsAsFactors=FALSE))
  dat <- dat[!duplicated(paste(dat$Genus, dat$Species)),
             c("Species", "Genus", "Family")]
  write.csv(dat[c("Species", "Genus", "Family")], "pinales.csv",
            row.names=FALSE)

  ## A set of species to use:
  set.seed(1)
  i <- sort(sample(nrow(dat), 60))
  extra <- match(c("Larix_griffithii", "Tsuga_jeffreyi"), dat$gs)
  i <- sort(c(i, extra))
  writeLines(dat$gs[i], "pinales_sub.txt")
}

## Not really a test, but it works!
test_that("regression", {
  phy <- read.tree("pinales.tre")
  data_species <- readLines("pinales_sub.txt")

  phy2 <- phyndr_genus(phy, data_species)

  expect_that(phy2, is_a("phyndr"))
  expect_that(phy2, is_a("phylo"))

  expect_that(length(phy2$clades[["genus::Tsuga"]]), equals(1))
  expect_that(length(phy2$clades[["genus::Larix"]]), equals(1))
  expect_that(length(phy2$clades[["genus::Cedrus"]]), equals(0))

  expect_that(all(names(phy2$clades) %in% phy2$tip.label), is_true())

  expect_that(length(phy2$tip.label), equals(20))

  expect_that(phyndr_n_distinct(phy), equals(1))

  if (FALSE) {
    col <- setNames(rep("black", length(phy2$tip.label)), phy2$tip.label)
    col[names(phy2$clades)] <-
      ifelse(viapply(phy2$clades, length) > 0L, "blue", "red")
    plot(phy2, type="fan", no.margin=TRUE, cex=.5, tip.color=col)
  }
})

test_that("taxonomy", {
  phy <- read.tree("pinales.tre")
  data_species <- readLines("pinales_sub.txt")

  dat <- read.csv("pinales.csv", stringsAsFactors=FALSE)
  rownames(dat) <- paste(dat$Genus, dat$Species, sep="_")
  dat <- dat[c("Genus", "Family")]

  ## Some species that aren't in the lookup but are in the tree:
  extra <- setdiff(phy$tip.label, rownames(dat))
  extra_genus <- split_genus(extra)
  tmp <- dat[match(extra_genus, dat$Genus), ]
  rownames(tmp) <- extra
  dat2 <- rbind(dat, tmp)

  phy2 <- phyndr_taxonomy(phy, data_species, dat2)
  expect_that(length(phy2$tip.label), equals(20))
  expect_that(length(phy2$clades), equals(2))
})

test_that("corner case", {
  ## https://github.com/richfitz/phyndr/issues/5

  t_str <- "(C_a:2.509256702,((B_a:1.246000031,(A_a:1.217375313,B_b:1.217375313)nd5:0.02862471784)nd3:0.7612144472,(A_b:0.703935084,A_c:0.703935084)nd4:1.303279394)nd2:0.5020422246)nd1;"
  t <- ape::read.tree(text=t_str)

  ## If only have data for one species of B and C, then phyndr_genus
  ## will correctly return a two taxon tree containing a
  ## representative from B and C
  d_1 <- c("C_x", "B_x")
  res_1 <- phyndr_genus(t, d_1)
  expect_that(length(res_1$tip.label), equals(2L))
  expect_that(sort(res_1$tip.label), equals(c("genus::B", "genus::C")))

  ## If we don't have any data for A or B, then it should return a
  ## tree with only 1 tip, but it returns the original tree
  d_2 <- c("C_x")
  expect_that(res_2 <- phyndr_genus(t, d_2),
              throws_error("Only one species/clade in tree: genus::C"))

  t2_str <- "(D_a:1.893702895,(C_a:1.492379389,(B_a:1.34277546,(A_a:0.7803371973,(B_b:0.24427762,(A_b:0.06128815295,A_c:0.06128815295)nd6:0.182989467)nd5:0.5360595774)nd4:0.5624382631)nd3:0.1496039284)nd2:0.4013235061)nd1;"
  t2 <- ape::read.tree(text=t2_str)

  res_3 <- phyndr_genus(t2, c("D_a", "C_x"))
  expect_that(length(res_3$tip.label), equals(2))
  expect_that(sort(res_3$tip.label),
              equals(sort(c("D_a", "genus::C"))))
})
