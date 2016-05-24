require(testthat)
require(compare)
require(emuR)

context("testing requeries")

.aeSampleRate=20000

.test_emu_ae_db=NULL
# .test_emu_ae_db_uuid='3f627b7b-4fb5-4b4a-8c79-b5f49df4df25'
.test_emu_ae_db_uuid = "0fc618dc-8980-414d-8c7a-144a649ce199"
.test_emu_ae_db_dir=NULL

path2demoData = file.path(tempdir(),"emuR_demoData")
path2testhatFolder = file.path(tempdir(),"emuR_testthat")

# extract internalVars from environment .emuR_pkgEnv
internalVars = get("internalVars", envir = .emuR_pkgEnv)

test_that("Copy example database ae",{
  legacyDbEmuAeTpl <- file.path(path2demoData, "legacy_ae", "ae.tpl")
  .test_emu_ae_db_dir<<-file.path(path2testhatFolder, 'test_emu_ae')
  unlink(.test_emu_ae_db_dir, recursive = T)
  
  # copy 4 faster tests
  dir.create(.test_emu_ae_db_dir)
  file.copy(file.path(path2demoData, paste0('ae', emuDB.suffix)), .test_emu_ae_db_dir, recursive = T)
})

test_that("requeries work on ae",{
  ae = load_emuDB(file.path(.test_emu_ae_db_dir, paste0('ae', emuDB.suffix)), inMemoryCache = internalVars$testingVars$inMemoryCache, verbose=FALSE)
  expect_that(ae$dbName,is_equivalent_to('ae'))
  
  test_that("Requery sequential",{
    
    # Phoneme sequences n->t
    sl1=query(ae, "[Phoneme=n -> Phoneme=t]")
    # requery two elemnts before and one after sequence
    rsl1=requery_seq(ae, sl1, offset=-2, length=5)
    rsl2=requery_seq(ae, sl1, offset=-3, length=5,offsetRef = 'END')
    
    # equivalent requery results should be equal
    expect_equal(rsl1,rsl2)
    
    expect_that(class(rsl1),is_identical_to(c('emuRsegs','emusegs','data.frame')))
    expect_that(nrow(sl1),equals(2))
    expect_that(nrow(rsl1),equals(2))
    expect_that('[.data.frame'(rsl1,1,'labels'),is_equivalent_to('l->@->n->t->l'))
    expect_that('[.data.frame'(rsl1,1,'startItemID'),equals(144))
    expect_that('[.data.frame'(rsl1,1,'endItemID'),equals(148))
    
    expect_that('[.data.frame'(rsl1,2,'labels'),is_equivalent_to('s->@->n->t->ei'))
    expect_that('[.data.frame'(rsl1,2,'startItemID'),equals(101))
    expect_that('[.data.frame'(rsl1,2,'endItemID'),equals(105))
    
    # Bug ID 42
    sl1=query(ae, "[[Phonetic=k -> Phonetic=~.*]->Phonetic=~.*]")
    sl1w=requery_hier(ae, sl1,level='Word')
    # sl1w has sequence length 1
    sl1w2=requery_seq(ae, sl1w[1,])
    # Bug startItemID != endItemID, and label is not a sequence !!
    expect_that('[.data.frame'(sl1w2,1,'startItemID'),equals(61))
    expect_that('[.data.frame'(sl1w2,1,'endItemID'),equals(61))
    
  })
  
  test_that("Requery hierarchical",{
    
    # Text beginning with 'a'
    sl1=query(ae, "Text=~'a[mn].*'")
    # requery to level Phoneme
    rsl1=requery_hier(ae, sl1,level='Phoneme')
    expect_that(class(rsl1),is_identical_to(c('emuRsegs','emusegs','data.frame')))
    expect_that(nrow(sl1),equals(3))
    expect_that(nrow(rsl1),equals(3))
    expect_that('[.data.frame'(rsl1,1,'labels'),is_equivalent_to('V->m->V->N->s->t'))
    expect_that('[.data.frame'(rsl1,1,'startItemID'),equals(114))
    expect_that('[.data.frame'(rsl1,1,'endItemID'),equals(119))
    
    expect_that('[.data.frame'(rsl1,2,'labels'),is_equivalent_to('E->n->i:'))
    expect_that('[.data.frame'(rsl1,2,'startItemID'),equals(135))
    expect_that('[.data.frame'(rsl1,2,'endItemID'),equals(137))
    
    expect_that('[.data.frame'(rsl1,3,'labels'),is_equivalent_to('@->n'))
    expect_that('[.data.frame'(rsl1,3,'startItemID'),equals(102))
    expect_that('[.data.frame'(rsl1,3,'endItemID'),equals(103))
    
  })
  
})