
# Test efetch() -----------------------------------------------------------

context("Testing 'efetch()'")

if (getOption('reutils.test.remote')) {
  
  test_that("Fetch PMIDs 17284678 and 9997 as text abstracts", {
    a <- efetch(uid = c(17284678, 9997), db = 'pubmed', rettype = 'abstract')
    expect_equal(retmode(a), "text")
    expect_equal(rettype(a), "abstract")
    ## Content
    expect_is(content(a), 'character')
    expect_is(content(a, "text"), 'character')
    expect_error(content(a, 'xml'), "Cannot return data of retmode.+")
    expect_error(content(a, 'json'), "Cannot return data of retmode.+")
    expect_is(content(a, 'parsed'), 'character')
    expect_is(content(a, 'textConnection'), 'textConnection')
  })
  
  test_that("Fetch PMIDs 17284678 and 9997 as XML", {
    a <- efetch(c(17284678, 9997), 'pubmed', retmode = 'xml')
    expect_equal(retmode(a), "xml")
    expect_equal(rettype(a), "")
    ## Content
    expect_is(content(a), 'XMLInternalDocument')
    expect_is(content(a, "text"), 'character')
    expect_is(content(a, 'xml'), "XMLInternalDocument")
    expect_error(content(a, 'json'), "Cannot return data of retmode.+")
    expect_is(content(a, 'parsed'), 'XMLInternalDocument')
    expect_error(content(a, 'textConnection'), "Cannot return data of retmode.+")
  })
  
  test_that("Fetch 100 bases of the minus strand of GI 21614549", {
    a <- efetch('21614549', 'nuccore', strand = 1, seqstart = 1, seqstop = 100,
                rettype = "fasta", retmode = "text")
    expect_equal(retmode(a), "text")
    expect_equal(rettype(a), "fasta")
    expect_equal(nchar(paste0(strsplit(content(a), '\n')[[1]][-1], collapse = "")), 100)    
  })

  test_that("efetch works with UIDs provided as character or numeric vectors", {
    a <- efetch(c(28800982, 28628843), "protein", "fasta")
    ## sometimes NCBI returns an empty TSeqSet for whatever reason. This completely
    ## breaks the test so we better check for that scenario
    if (!all(is.na(a$xmlName("//TSeqSet/*")))) {
      expect_equal(a$xmlName("//TSeqSet/*"), c("TSeq", "TSeq"))
    }
  })
  
  test_that("efetch works with an esearch object as input", {
    x <- esearch("erythroid associated factor AND Homo sapiens", 'protein', retmax = 2)
    a <- efetch(x, rettype="fasta")
    if (!all(is.na(a$xmlValue("//TSeq_taxid")))) {
      expect_equal(a$xmlValue("//TSeq_taxid"), c("9606", "9606"))
    }
  })
  
}
