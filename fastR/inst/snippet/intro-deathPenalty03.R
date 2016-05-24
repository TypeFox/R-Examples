require(vcd)
mosaic(~Victim+Defendant+DeathPenalty,data=deathPen)
structable(~Victim+Defendant+DeathPenalty,data=deathPen)
