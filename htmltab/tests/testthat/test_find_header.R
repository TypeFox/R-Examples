context("Correctly identified header elements")

tab1_code <- '
<table>
<tr>
<th rowspan="2">a</th>
<th>b</th>
<th>c</th>
</tr>
<tr>
<th></th>
<th></th>
</tr>
<tr>
<td>1</td>
<td>2</td>
<td>3</td>
</tr>
</table>'

test_that("Correctly identified header 1", {

  tab1 <- XML::htmlParse(tab1_code)
  suppressMessages(tab1 <- htmltab(tab1))

  expect_that(tab1[,1], equals("1"))
  expect_that(tab1[,2], equals("2"))
  expect_that(tab1[,3], equals("3"))

  expect_that(colnames(tab1)[1], equals("a"))
  expect_that(colnames(tab1)[2], equals("b"))
  expect_that(colnames(tab1)[3], equals("c"))
})

tab2_code <- '
<table>
<thead>
<tr>
<th rowspan="2">a</th>
<th>b</th>
<th colspan="2" rowspan="2">c</th>
</tr>
<tr>
<td></td>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>2</td>
<td>3</td>
<td></td>
</tr>
</tbody>
</table>'

test_that("Correctly identified header 2", {

  tab2 <- XML::htmlParse(tab2_code)
  suppressMessages(tab2 <- htmltab(tab2))

  expect_that(tab2[,1], equals("1"))
  expect_that(tab2[,2], equals("2"))
  expect_that(tab2[,3], equals("3"))
    expect_that(is.na(tab2[,4]), is_true())

  expect_that(colnames(tab2)[1], equals("a"))
  expect_that(colnames(tab2)[2], equals("b"))
  expect_that(colnames(tab2)[3], equals("c"))
  expect_that(colnames(tab2)[4], equals("c"))
})


tab3_code <- '<table>
<thead>
<tr>
<th>a</th>
<th>b</th>
<th>c</th>
<th colspan="3">d</th>
</tr>
<tr>
<td colspan="4">e</td>
</tr>
</thead>

<tr>
<td>1</td>
<td>2</td>
<td>3</td>
<td></td>
</tr>
</table>'

test_that("Correctly identified header 3", {

  tab3 <- XML::htmlParse(tab3_code)
  suppressMessages(tab3 <- htmltab(tab3, fillNA = ""))

  expect_that(tab3[,1], equals("1"))
  expect_that(tab3[,2], equals("2"))
  expect_that(tab3[,3], equals("3"))
  expect_that(tab3[,4], equals(""))

  expect_that(colnames(tab3)[1], equals("a >> e"))
  expect_that(colnames(tab3)[2], equals("b >> e"))
  expect_that(colnames(tab3)[3], equals("c >> e"))
  expect_that(colnames(tab3)[4], equals("d >> e"))
})


tab4_code <- '
<table>
<tr>
<th>a</th>
<th>b</th>
<th>c</th>
</tr>
<tr>
<th>a2</th>
<th>b2</th>
<th>c2</th>
</tr>
<tr>
<td>1</td>
<td>2</td>
<td>3</td>
</tr>
</table>'

test_that("Correctly identified header 4", {

  tab4 <- XML::htmlParse(tab4_code)
  suppressMessages(tab4 <- htmltab(tab4))

  expect_that(tab4[,1], equals("1"))
  expect_that(tab4[,2], equals("2"))
  expect_that(tab4[,3], equals("3"))

  expect_that(colnames(tab4)[1], equals("a >> a2"))
  expect_that(colnames(tab4)[2], equals("b >> b2"))
  expect_that(colnames(tab4)[3], equals("c >> c2"))
})


test_that("Correctly identified header 5", {

  tab5 <- XML::htmlParse(tab4_code)
  suppressMessages(tab5 <- htmltab(tab5, header = 1:2))

  expect_that(tab5[,1], equals("1"))
  expect_that(tab5[,2], equals("2"))
  expect_that(tab5[,3], equals("3"))

  expect_that(colnames(tab5)[1], equals("a >> a2"))
  expect_that(colnames(tab5)[2], equals("b >> b2"))
  expect_that(colnames(tab5)[3], equals("c >> c2"))
})

test_that("Correctly identified header 6", {

  tab6 <- XML::htmlParse(tab4_code)
  suppressMessages(tab6 <- htmltab(tab6, body = 3))

  expect_that(tab6[,1], equals("1"))
  expect_that(tab6[,2], equals("2"))
  expect_that(tab6[,3], equals("3"))

  expect_that(colnames(tab6)[1], equals("a >> a2"))
  expect_that(colnames(tab6)[2], equals("b >> b2"))
  expect_that(colnames(tab6)[3], equals("c >> c2"))
})


tab5_code <- '
<table>
<thead>
<tr>
<th>a</th>
<th>b</th>
<th>c</th>
</tr>
</thead>

<tbody>
<tr>
<th>a2</th>
<th>b2</th>
<th>c2</th>
</tr>
<tr>
<td>1</td>
<td>2</td>
<td>3</td>
</tr>
</tbody>
</table>'

test_that("Header info in thead and tbody", {

  tab7 <- XML::htmlParse(tab5_code)
  suppressMessages(tab7 <- htmltab(tab7, header = 1:2))

  expect_that(tab7[,1], equals("1"))
  expect_that(tab7[,2], equals("2"))
  expect_that(tab7[,3], equals("3"))

  expect_that(colnames(tab7)[1], equals("a >> a2"))
  expect_that(colnames(tab7)[2], equals("b >> b2"))
  expect_that(colnames(tab7)[3], equals("c >> c2"))
})

test_that("Correctly identified header 8", {

  tab8 <- XML::htmlParse(tab5_code)
  suppressMessages(tab8 <- htmltab(tab8, header = 1:2, body = 3))

  expect_that(tab8[,1], equals("1"))
  expect_that(tab8[,2], equals("2"))
  expect_that(tab8[,3], equals("3"))

  expect_that(colnames(tab8)[1], equals("a >> a2"))
  expect_that(colnames(tab8)[2], equals("b >> b2"))
  expect_that(colnames(tab8)[3], equals("c >> c2"))
})
