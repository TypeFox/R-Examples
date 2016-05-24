## Asian option tests

library(testthat)
## we will generally use these parameters:
s <- 100; k <- 97; v <- 0.30; r <- 0.08;  tt <- 1.5; d <- 0.015
m <- 8; numsim <- 1e04;

x <- geomavgprice(s, k, v, r, tt, d, m)
test_that('geomavgprice works', {
    expect_equivalent(x["Call"],  11.979449562844800)
    expect_equivalent(x["Put"],  5.350469024489170000)
}
)
print('geomavgprice OK')

x <- geomavgprice(s, k, v, r, tt, d, m, cont=TRUE)
test_that('geomavgprice works, cont=TRUE', {
          expect_equivalent(x["Call"], 10.928656909353800)
          expect_equivalent(x["Put"], 4.878795485237670000)
}
)
print('geomavgprice OK with cont=TRUE')

          
x <- geomavgstrike(s, km=k, v, r, tt, d, m)
test_that('geomavgstrike works', {
    expect_equivalent(x["Call"], 11.804071467006300)
    expect_equivalent(x["Put"], 3.909402760594920000)
}
)
print('geomavgstrike OK')

x <- geomavgstrike(s, km=k, v, r, tt, d, m, cont=TRUE)
test_that('geomavgstrike works, cont=TRUE', {
          expect_equivalent(x["Call"],  12.834623645733200)
          expect_equivalent(x["Put"], 4.378209398509600000)
}
)
print('geomavgstrike OK with cont=TRUE')

set.seed(1)
x <- arithasianmc(s, k, v, r, tt, d, m, numsim=1e04)
test_that('arithasianmc works', {
    expect_equivalent(x['Avg Price', 'Call'], 12.5674068074862450572)
    expect_equivalent(x['Avg Price', 'Put'], 5.07404769559155521818)
    expect_equivalent(x['Avg Strike', 'Call'], 9.30977750803833181692)
    expect_equivalent(x['Avg Strike', 'Put'], 5.39340666762663722977)
    expect_equivalent(x['Vanilla', 'Call'], 19.7327118854260668001)
    expect_equivalent(x['Vanilla', 'Put'], 8.322981933119683262134)
})
print('arithasianmc OK')

set.seed(1); x <- geomasianmc(s, k, v, r, tt, d, m, 1e04)
test_that('geomasianmc works', {
    expect_equivalent(x['Avg Price', 'CallMC'], 11.83382442169169657120)
    expect_equivalent(x['Avg Price', 'PutMC'], 5.396934864526733655055)
    expect_equivalent(x['Avg Price', 'CallExact'], 11.97944956284476347719)
    expect_equivalent(x['Avg Price', 'PutExact'], 5.350469024489171943060)
    expect_equivalent(x['Avg Strike', 'CallMC'], 11.64281636034952960301)
    expect_equivalent(x['Avg Strike', 'PutMC'], 3.895930807646228899443)
    expect_equivalent(x['Avg Strike', 'CallExact'], 11.80407146700625276026)
    expect_equivalent(x['Avg Strike', 'PutExact'], 3.909402760594900883007)
    expect_equivalent(x['Vanilla', 'CallMC'], 19.73271188542606680016)
    expect_equivalent(x['Vanilla', 'PutMC'], 8.322981933119683262134)
    expect_equivalent(x['Vanilla', 'CallExact'], 20.06165784288828035642)
    expect_equivalent(x['Vanilla', 'PutExact'], 8.317816485118925129427)
})
print('geomasianmc OK')

set.seed(1); x <- arithavgpricecv(s, k, v, r, tt, d, m, 1e04)
test_that('arithavgpricecv works', {
    expect_equivalent(x["price"], 12.72802675462704158349)
    expect_equivalent(x["beta"], 1.053529789518180548313)
})
print('arithavgpricecv works')
