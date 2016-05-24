context("makeRBPObj")

test_that("makeRBPObj", {
  # target as factor and as numeric
  ylist = list(y, as.numeric(y == positive[1]))

  for (j in 1:length(ylist)) {
    if (is.factor(ylist[[j]])) {
      pos = sonar.task$task.desc$class.levels
      for (i in pos) {
        obj = makeRBPObj(pred, ylist[[j]], i)
        expect_equal(obj$positive, i)
      }
      # by default the first class is the positive class
      obj = makeRBPObj(pred,  ylist[[j]])
      expect_equal(obj$positive, pos[1])
    } else {
      # when target is numeric, there is no positive class
      obj = makeRBPObj(pred, ylist[[j]])
      expect_equal(obj$positive, NULL)
      # when target is numeric, no positive class can be chosen
      expect_error(makeRBPObj(pred, ylist[[j]], i), 
        "Assertion on 'positive' failed: Must be NULL")
      # when target has numeric values that are not 0 or 1
      expect_error(makeRBPObj(pred, ylist[[j]] + 1), 
        "Assertion on 'y' failed: Must be a subset of \\{'0','1'\\}")
    }
  }
})