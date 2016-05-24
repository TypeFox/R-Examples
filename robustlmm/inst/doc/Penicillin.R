## load Penicillin data and create contaminated data
data(Penicillin, package="lme4")
Penicillin <- within(Penicillin, plate <- reorder(plate, diameter))
attr(Penicillin$plate, "scores") <- NULL
PenicillinC <- within(Penicillin, {
  contaminated <- rep.int("original", nrow(Penicillin))
  diameter[plate == "k" & sample == "F"] <- min(diameter)
  contaminated[plate == "k" & sample == "F"] <- "changed"
  diameter[plate == "g"] <- diameter[plate == "g"] / 5 * 4
  contaminated[plate == "g"] <- "changed"
  contaminated <- factor(contaminated)
})
PenicillinCE <- within(PenicillinC, {
  diameter[plate == "g"] <- diameter[plate == "g"] / 4
})
