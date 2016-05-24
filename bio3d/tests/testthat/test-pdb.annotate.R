context("Testing pdb.annotate()")


test_that("PDB annotation works", {
  skip_on_cran()

  expected <- c('3R1C_X', '3R1C_B', '3R1C_C', '3R1C_D',
                '3R1C_E', '3R1C_F', '3R1C_G', '3R1C_H',
                '3R1C_I', '3R1C_J', '3R1C_K', '3R1C_L',
                '3R1C_M', '3R1C_N', '3R1C_O', '3R1C_P',
                '3R1C_Q', '3R1C_R', '3R1C_S', '3R1C_Y',
                '3R1C_T', '3R1C_U', '3R1C_W', '3R1C_A',
                '3R1C_V', '3R1C_Z', '3R1C_a', '3R1C_b',
                '3R1C_c', '3R1C_d', '3R1C_e', '3R1C_f',
                '3R1C_g', '3R1C_h', '3R1C_i', '3R1C_j')

  invisible(capture.output(anno <- pdb.annotate(expected)))
  expect_identical(rownames(anno), expected)

  expected <- c('3R1C_A', '3R1C_B', '3R1C_C', '3R1C_D',
                '3R1C_E', '3R1C_F', '3R1C_G', '3R1C_H',
                '3R1C_I', '3R1C_J', '3R1C_K', '3R1C_L',
                '3R1C_M', '3R1C_N', '3R1C_O', '3R1C_P',
                '3R1C_Q', '3R1C_R', '3R1C_S', '3R1C_T',
                '3R1C_U', '3R1C_V', '3R1C_W', '3R1C_X',
                '3R1C_Y', '3R1C_Z', '3R1C_a', '3R1C_b',
                '3R1C_c', '3R1C_d', '3R1C_e', '3R1C_f',
                '3R1C_g', '3R1C_h', '3R1C_i', '3R1C_j')
  
  invisible(capture.output(anno <- pdb.annotate('3R1C')))
  expect_identical(rownames(anno), expected)

  expected <- c('3R1C_A', '3R1C_B', '3R1C_C', '3R1C_D',
                '3R1C_E', '3R1C_F', '3R1C_G', '3R1C_H',
                '3R1C_I', '3R1C_J', '3R1C_K', '3R1C_L',
                '3R1C_M', '3R1C_N', '3R1C_O', '3R1C_P',
                '3R1C_Q', '3R1C_R', '3R1C_S', '3R1C_T',
                '3R1C_U', '3R1C_V', '3R1C_W', '3R1C_X',
                '3R1C_Y', '3R1C_Z', '3R1C_a', '3R1C_b',
                '3R1C_c', '3R1C_d', '3R1C_e', '3R1C_f',
                '3R1C_g', '3R1C_h', '3R1C_i', '3R1C_j',
                '1CDK_A', '1CDK_B', '1CDK_I', '1CDK_J')
  
  invisible(capture.output(anno <- pdb.annotate(c('3R1C', '1CDK'))))
  expect_identical(rownames(anno), expected)
  
  invisible(capture.output(anno <- pdb.annotate(c('3R1C_A', '3r1c_a', '3r1c_q'))))
  expect_identical(rownames(anno), expected[c(1, 27)])

  invisible(capture.output(anno <- pdb.annotate(c('3R1C_A', '3r1c_a', '3r1c_q'), unique=TRUE)))
  expect_identical(rownames(anno), "3R1C")
  expect_identical(anno$chainId, "A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,a,b,c,d,e,f,g,h,i,j")
  
  
  expected <- rep("ANP,MN,MYR,TPO", 2)
  
  invisible(capture.output(anno <- pdb.annotate(c('1cdk_A', '1cdk_B'), anno.terms="ligandId")))
  expect_identical(anno$ligandId, expected)

  
})
