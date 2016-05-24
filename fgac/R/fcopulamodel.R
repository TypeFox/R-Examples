"fcopulamodel" <-
function(theta,delta,x,y, model = c("pCBB1", "pCBB2","pCBB3","pCBB4","pCBB5","pCBB6","pCBB7","pCMax","pCMin")) {
  model <- match.arg(model)
  switch(model,
pCBB1 = pCBB1(theta,delta,x,y),
pCBB2 = pCBB2(theta,delta,x,y),
pCBB3 = pCBB3(theta,delta,x,y),
pCBB4 = pCBB4(theta,delta,x,y),
pCBB5 = pCBB5(theta,delta,x,y),
pCBB6 = pCBB6(theta,delta,x,y),
pCBB7 = pCBB7(theta,delta,x,y),
pCMax = pCMax(theta,delta,x,y),
pCMin = pCMin(theta,delta,x,y))
}

