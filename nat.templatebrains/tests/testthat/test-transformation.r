context("Find Bridging Registrations")

test_that("we can find registrations",{

  op=options(nat.templatebrains.regdirs=NULL)
  on.exit(options(op))
  expect_error(find_reg("rhubarb.list"), "registration directories")

  td=tempfile(pattern = 'regdir1')
  dir.create(td)
  options('nat.templatebrains.regdirs'=td)

  expect_equal(find_reg("rhubarb.list"), "")
  expect_error(find_reg("rhubarb.list", mustWork = TRUE), "Unable to find")
  dir.create(file.path(td,"rhubarb.list"), recursive = T)
  expect_equal(find_reg("rhubarb.list"), file.path(td,"rhubarb.list"))

  td2=tempfile(pattern = 'regdir2')
  dir.create(file.path(td2, 'rhubarb.list'), recursive = T)
  dir.create(file.path(td2, 'crumble.list'), recursive = T)

  options(nat.templatebrains.regdirs=c(td,td2))
  expect_equal(find_reg("rhubarb.list"), file.path(td,"rhubarb.list"))
  expect_equal(find_reg("crumble.list"), file.path(td2,"crumble.list"))
  unlink(td, recursive = TRUE)
  unlink(td2, recursive = TRUE)
})


test_that("we can find bridging registrations",{
  td=tempfile(pattern = 'extrabridge')
  dir.create(td)
  op=options(nat.templatebrains.regdirs=td)
  on.exit(options(op))

  rcreg=file.path(td,"rhubarb_crumble.list")
  dir.create(rcreg, recursive = T)
  expect_true(nzchar(bridging_reg(ref='rhubarb', sample='crumble')))
  expect_error(bridging_reg(ref='crumble', sample='rhubarb', mustWork = T))
  expect_false(nzchar(bridging_reg(ref='crumble', sample='rhubarb')))

  br<-bridging_reg(reference ='crumble', sample='rhubarb', checkboth = TRUE)
  expect_true(nzchar(br))
  expect_true(attr(br,'swapped'))

  # now try equivalence of bridging_sequence and bridging_reg
  expect_equivalent(bridging_sequence(ref="rhubarb", sample = "crumble"),
                    bridging_reg(ref='rhubarb', sample='crumble'))
  expect_equivalent(bridging_sequence(ref="crumble", sample = "rhubarb", checkboth = F),
                    bridging_reg(ref='crumble', sample='rhubarb', checkboth = F))
  expect_equivalent(bridging_sequence(ref="crumble", sample = "rhubarb", checkboth=T),
                    bridging_reg(ref='crumble', sample='rhubarb', checkboth=T))
  expect_error(bridging_sequence(ref='crumble', sample='rhubarb', checkboth = F, mustWork=T))
})

context("Bridging Graph")

test_that("bridging graph and friends work",{

  expect_is(emptydf<-allreg_dataframe(NULL), 'data.frame')
  expect_equal(nrow(emptydf), 0L)

  # make a new registration directory for testing purposes
  dir.create(td<-tempfile())
  on.exit(unlink(td, recursive = TRUE))

  # set up fake registrations under that directory
  # originally made as:
  # df=subset(allreg_dataframe(), !dup)
  # saveRDS(df, 'testdata/allreg_dataframe.rds', compress='xz')
  df=readRDS("testdata/allreg_dataframe.rds")
  df$path=file.path(td, df$path)
  sapply(df$path, dir.create, recursive = TRUE, showWarnings = F)

  # note that dirs are not searched recursively - only the top level is listed
  op=options(nat.templatebrains.regdirs=unique(dirname(df$path)))
  on.exit(options(op), add = TRUE)

  expect_is(df2<-allreg_dataframe(), 'data.frame')
  # maybe dangerous to assume sort order ...
  df=df[match(df$path, df2$path),]
  expect_equal(df2, df)

  baseline=structure(file.path(td, "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/nat.flybrains/extdata/bridgingregistrations/FCWB_IS2.list"),
    swap = TRUE)
  expect_equal(shortest_bridging_seq('FCWB', 'IS2'), baseline)
  expect_is(shortest_bridging_seq('FCWB', 'IBN'), 'character')
  expect_warning(fcwb_ibn<-shortest_bridging_seq('FCWB', 'IBN', imagedata = T),
                 'very slow for image data')
  bl=structure(c(file.path(td,"/Users/jefferis/Library/Application Support/rpkg-nat.flybrains/regrepos/BridgingRegistrations/JFRC2_FCWB.list"),
                 file.path(td,"/Library/Frameworks/R.framework/Versions/3.2/Resources/library/nat.flybrains/extdata/bridgingregistrations/JFRC2_IBNWB.list"),
                 file.path(td,"/Library/Frameworks/R.framework/Versions/3.2/Resources/library/nat.flybrains/extdata/bridgingregistrations/IBNWB_IBN.list")),
               swap = c(F, T, T))
  expect_equal(fcwb_ibn, bl)
  expect_error(shortest_bridging_seq('FCWB', 'crumble'))
})

context("Transformation")

if(is.null(cmtk.bindir())){
  message("skipping transformation tests since CMTK is not installed")
} else {
test_that("can use a bridging registration in regdirs",{
  td=tempfile(pattern = 'extrabridge')
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE))
  op=options('nat.templatebrains.regdirs'=td)
  on.exit(options(op), add = TRUE)

  library(rgl)
  library(nat)
  rcreg=file.path(td,"rhubarb_crumble.list")
  m=matrix(c(1,1,1),ncol=3)
  expect_error(xform_brain(m, ref='rhubarb',sample='crumble'),
               'No bridging registrations')
  i=identityMatrix()
  cmtk.mat2dof(identityMatrix(), rcreg)
  expect_equal(xform_brain(m, ref='rhubarb',sample='crumble'), m)

  # round trip test based on affine component of JFRC2_IS2
  # x=read.cmtkreg(nat.templatebrains:::bridging_reg(ref=JFRC2, sample=IS2))
  # m=cmtkparams2affmat(x$registration$affine_xform, legacy = T)

  jfrc_is2_reg=file.path(td,"JFRC2_IS2.list")
  aff=structure(c(0.823766999979893, 0.0126432000004598, 0.0110175000032261,
                0, -0.0127625146453663, 0.862393114320825, 0.0848381903456082,
                0, -0.0314435105710272, 0.135913761768978, 1.46958245308841,
                0, -105.309000027786, 15.6708999995602, -7.7248100042806, 1),
              .Dim = c(4L, 4L))
  cmtk.mat2dof(aff, jfrc_is2_reg)
  points=matrix(c(100, 100, 50), ncol=3)
  expect_equal(xform_brain(points, ref='JFRC2',sample='JFRC2', via='IS2', checkboth = T), points)

  # round trip test with neuron and image data representing the neuron
  kc1=kcs20[[1]]
  kc1.rt=xform_brain(kc1, ref='JFRC2',sample='JFRC2', via='IS2')
  kcim=as.im3d(xyzmatrix(kc1), voxdims=c(1, 1, 1),
               BoundingBox=c(250, 410, 0, 130, 0, 120))
  dir.create(td2<-tempfile())
  owd=setwd(td2)
  on.exit(setwd(owd), add = TRUE)
  on.exit(unlink(td2, recursive = TRUE), add = TRUE)

  imfile=file.path(td2, 'kcim.nrrd')
  outfile=file.path(td2, 'kcim_roundtrip.nrrd')
  write.im3d(kcim, imfile, dtype='byte')
  xform_brain(imfile, ref='JFRC2',sample='JFRC2', via='IS2', output=outfile,
              target=imfile, checkboth = TRUE, Verbose=F)
  kc2=dotprops(outfile)
  # write an inverted registraion
  cmtk.mat2dof(solve(aff), file.path(td,"IS2_JFRC2.list"))
  outfile2=file.path(td2, 'kcim_roundtrip2.nrrd')
  # check we can still xform image
  xform_brain(imfile, ref='JFRC2',sample='JFRC2', via='IS2', output=outfile2,
              target=imfile, checkboth = TRUE, Verbose=F)
  kc3=dotprops(outfile2)
  # check mean distance between points
  library(nabor)
  expect_true(mean(knn(xyzmatrix(kc2), xyzmatrix(kc1), k = 1)$nn.dists)<1.0,
              "round trip with cmtk reformatx inversion")
  expect_true(mean(knn(xyzmatrix(kc3), xyzmatrix(kc1), k = 1)$nn.dists)<1.0,
              "round trip with pre-inverted registration")
})

}
