context("publish")

if(0){  # skip publishing tests.  These were all still passing at last check, but rfigshare configuration for testing is not ideal.  
  
  
  
# This loads the rOpenSci figshare sandbox credentials, so that the example 
# can run automatically during check and install.  Unlike normal figshare accounts,
# data loaded to this testing sandbox is periodically purged.  
library(rfigshare)
status <- try(fs_auth(token = "xdBjcKOiunwjiovwkfTF2QjGhROeLMw0y0nSCSgvg3YQxdBjcKOiunwjiovwkfTF2Q", token_secret = "4mdM3pfekNGO16X4hsvZdg"))
if(is(status, "try-error") || (is(status, "response") && status$status_code != 200)){
  warning("Could not authenticate figshare, skipping figshare tests")
} else {



## Create example file
  library(geiger)
  data(geospiza)

  geiger_nex <- add_trees(geospiza$phy)
  geiger_nex <- add_characters(geospiza$dat, geiger_nex)
  geiger_nex <- add_basic_meta(
    title = "Geospiza phylogeny with character data rendered as NeXML", 
    creator = "Carl Boettiger", 
    description = "This example NeXML file was created using the data originally provided in the geiger package for R to illustrate how this data can be stored, shared and distributed as NeXML.", 
    citation = citation("geiger"), 
    nexml = geiger_nex)


  test_that("We can publish to figshare", {

## Publish 
  id <- nexml_publish(geiger_nex, visibility="public", repo="figshare")


## Download and parse publication 
## Note that at present, only public files can be automatically downloaded from figshare
  library(rfigshare)
  test_nex <- nexml_read(fs_download(id))

## Extract and compare metadata from upload and download 
  m <- get_metadata(geiger_nex)
  test_m <- get_metadata(test_nex)
  expect_equal(m["dc:title"], test_m["dc:title"])
  expect_equal(m["dc:description"], test_m["dc:description"])

## Check that DOI resolves -- doesn't for the test account
#library(httr)
#page <- GET(test_m[["dc:identifier"]])
#expect_equal(page$status_code, 200)


# Check that we avoid repeated metadata entries

  expect_equal(sum(match(names(test_m), "dc:pubdate"), na.rm=TRUE), 1)
  expect_equal(sum(match(names(test_m), "cc:license"), na.rm=TRUE), 1)


  })



}

}

