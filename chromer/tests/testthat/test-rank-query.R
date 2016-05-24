source("helper-chromer.R")

context("Texting rank-based query")

## Query for genus Lachemilla for testing purposes
cp <- chrom_counts(taxa="Lachemilla", rank="genus", full=FALSE)
cf <- chrom_counts(taxa="Lachemilla", rank="genus", full=TRUE)

## Call function for making short species name
short_species_name <- chromer:::short_species_name

test_that("Main function returns the correct format", {

    expect_that(cp, is_a("data.frame"))
    expect_that(cf, is_a("data.frame"))
    expect_equal(ncol(cp), 4)
    expect_equal(ncol(cf), 15)
    expect_equal(nrow(cp), nrow(cf))

})

test_that("No match returns empty data frame", {

    dum <- chrom_counts(taxa="Notagenus", rank="genus")
    expect_that(nrow(dum), equals(0))
    expect_that(ncol(dum), equals(0))
})


test_that("Querying multiple taxa works", {

    mp <- chrom_counts(taxa=c("Lachemilla", "Lachemilla"), rank="genus",
                       full=FALSE)
    expect_that(mp, is_a("data.frame"))
    expect_that(nrow(mp), equals(2*nrow(cp)))
    mp <- chrom_counts(taxa=list("Lachemilla", "Lachemilla"), rank="genus",
                       full=FALSE)
    expect_that(mp, is_a("data.frame"))
    expect_that(nrow(mp), equals(2*nrow(cp)))
})
    

test_that("Query worked properly", {

    ## Full records
    gen <- unique(cf$genus)
    fam <- unique(cf$family)
    expect_that(gen, equals("Lachemilla"))
    expect_that(fam, equals("Rosaceae"))

    ## Partial records
    gen <- unique(sapply(cp$resolved_binomial, function(x)
                         {strsplit(x, split="_")[[1]][1]}))
    expect_that(gen, equals("Lachemilla"))

})

test_that("Species querying works properly for both possible queries", {

    ## With space
    min_spc <- data.frame(chrom_counts("Castilleja minata", "species"))
    min_und <- data.frame(chrom_counts("Castilleja_minata", "species"))

    expect_that(min_spc, equals(min_und))
})

test_that("Building species name worked properly",{

    spf <- cf$resolved_binomial
    spp <- cp$resolved_binomial
    expect_that(spf, equals(spp))

    sp_tmp <- cf$resolved_name[1]
    expect_that(short_species_name(sp_tmp), equals(spf[1]))

    ## Make sure strsplit is working properly
    sp_dum <- "Dummy species (Tax. auth.) Name"
    sh_dum <- short_species_name(sp_dum)
    expect_that(sh_dum, equals("Dummy_species"))

    ## Varieties
    sp_var <- "Dummy species var. x (Tax. auth.) Name"
    sh_var <- short_species_name(sp_var)
    expect_that(sh_var, equals("Dummy_species_var._x"))

    ## Subspecies
    sp_sub <- "Dummy species subsp. x (Tax. auth.) Name"
    sh_sub <- short_species_name(sp_sub)
    expect_that(sh_sub, equals("Dummy_species_subsp._x"))
})


test_that("Bad input throws error", {

    expect_that(chrom_counts("foo", c("genus", "family")),
                throws_error())
    expect_that(chrom_counts("foo", "foo"), throws_error())

})
    
    

    
    
    

