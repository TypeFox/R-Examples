# NIST conversion factors

NISTabampereTOampere<-function(abampere) abampere * 10
# Input: abampere
# Output: ampere (A)

NISTampereTOabampere<-function(ampere) ampere / 10
# Input: ampere (A)
# Output: abampere

NISTabcoulombTOcoulomb<-function(abcoulomb) abcoulomb * 10
# Input: abcoulomb
# Output: coulomb (C)

NISTcoulombTOabcoulomb<-function(coulomb) coulomb / 10
# Input: coulomb (C)
# Output: abcoulomb

NISTabfaradTOfarad<-function(abfarad) abfarad * 1e+09
# Input: abfarad
# Output: farad (F)

NISTfaradTOabfarad<-function(farad) farad / 1e+09
# Input: farad (F)
# Output: abfarad

NISTabhenryTOhenry<-function(abhenry) abhenry * 1e-09
# Input: abhenry
# Output: henry (H)

NISThenryTOabhenry<-function(henry) henry / 1e-09
# Input: henry (H)
# Output: abhenry

NISTabmhoTOsiemens<-function(abmho) abmho * 1e+09
# Input: abmho
# Output: siemens (S)

NISTsiemensTOabmho<-function(siemens) siemens / 1e+09
# Input: siemens (S)
# Output: abmho

NISTabohmTOohm<-function(abohm) abohm * 1e-09
# Input: abohm
# Output: ohm (ohm)

NISTohmTOabohm<-function(ohm) ohm / 1e-09
# Input: ohm (ohm)
# Output: abohm

NISTabvoltTOvolt<-function(abvolt) abvolt * 1e-08
# Input: abvolt
# Output: volt (V)

NISTvoltTOabvolt<-function(volt) volt / 1e-08
# Input: volt (V)
# Output: abvolt

NISTaccelOfFreeFallStdTOmeterPerSecSqrd<-function(accelOfFreeFallStd) accelOfFreeFallStd * 9.80665
# Input: acceleration of free fall, standard (gn)
# Output: meter per second squared (m / s2)

NISTmeterPerSecSqrdTOaccelOfFreeFallStd<-function(meterPerSecSqrd) meterPerSecSqrd / 9.80665
# Input: meter per second squared (m / s2)
# Output: acceleration of free fall, standard (gn)

NISTacreTOsqrMeter<-function(acre) acre * 4046.873
# Input: acre (based on U.S. survey foot) 7
# Output: square meter (m2)

NISTsqrMeterTOacre<-function(sqrMeter) sqrMeter / 4046.873
# Input: square meter (m2)
# Output: acre (based on U.S. survey foot) 7

NISTacreFtTOcubMeter<-function(acreFt) acreFt * 1233.489
# Input: acre foot (based on U.S. survey foot) 7
# Output: cubic meter (m3)

NISTcubMeterTOacreFt<-function(cubMeter) cubMeter / 1233.489
# Input: cubic meter (m3)
# Output: acre foot (based on U.S. survey foot) 7

NISTampereHourTOcoulomb<-function(ampereHour) ampereHour * 3600
# Input: ampere hour (A * h)
# Output: coulomb (C)

NISTcoulombTOampereHour<-function(coulomb) coulomb / 3600
# Input: coulomb (C)
# Output: ampere hour (A * h)

NISTangstromTOmeter<-function(angstrom) angstrom * 1e-10
# Input: angstrom (a)
# Output: meter (m)

NISTmeterTOangstrom<-function(meter) meter / 1e-10
# Input: meter (m)
# Output: angstrom (a)

NISTangstromTOnanometer<-function(angstrom) angstrom * 0.1
# Input: angstrom (a)
# Output: nanometer (nm)

NISTnanometerTOangstrom<-function(nanometer) nanometer / 0.1
# Input: nanometer (nm)
# Output: angstrom (a)

NISTareTOsqrMeter<-function(are) are * 100
# Input: are (a)
# Output: square meter (m2)

NISTsqrMeterTOare<-function(sqrMeter) sqrMeter / 100
# Input: square meter (m2)
# Output: are (a)

NISTastronomicalUnTOmeter<-function(astronomicalUn) astronomicalUn * 149597900000
# Input: astronomical unit (ua)
# Output: meter (m)

NISTmeterTOastronomicalUn<-function(meter) meter / 149597900000
# Input: meter (m)
# Output: astronomical unit (ua)

NISTatmosphereStdTOpascal<-function(atmosphereStd) atmosphereStd * 101325
# Input: atmosphere, standard (atm)
# Output: pascal (Pa)

NISTpascalTOatmosphereStd<-function(pascal) pascal / 101325
# Input: pascal (Pa)
# Output: atmosphere, standard (atm)

NISTatmosphereStdTOkpascal<-function(atmosphereStd) atmosphereStd * 101.325
# Input: atmosphere, standard (atm)
# Output: kilopascal (kPa)

NISTkpascalTOatmosphereStd<-function(kpascal) kpascal / 101.325
# Input: kilopascal (kPa)
# Output: atmosphere, standard (atm)

NISTatmosphereTechnicalTOpascal<-function(atmosphereTechnical) atmosphereTechnical * 98066.5
# Input: atmosphere, technical (at) 8
# Output: pascal (Pa)

NISTpascalTOatmosphereTechnical<-function(pascal) pascal / 98066.5
# Input: pascal (Pa)
# Output: atmosphere, technical (at) 8

NISTatmosphereTechnicalTOkpascal<-function(atmosphereTechnical) atmosphereTechnical * 98.0665
# Input: atmosphere, technical (at) 8
# Output: kilopascal (kPa)

NISTkpascalTOatmosphereTechnical<-function(kpascal) kpascal / 98.0665
# Input: kilopascal (kPa)
# Output: atmosphere, technical (at) 8

NISTbarTOpascal<-function(bar) bar * 1e+05
# Input: bar (bar)
# Output: pascal (Pa)

NISTpascalTObar<-function(pascal) pascal / 1e+05
# Input: pascal (Pa)
# Output: bar (bar)

NISTbarTOkpascal<-function(bar) bar * 100
# Input: bar (bar)
# Output: kilopascal (kPa)

NISTkpascalTObar<-function(kpascal) kpascal / 100
# Input: kilopascal (kPa)
# Output: bar (bar)

NISTbarnTOsqrMeter<-function(barn) barn * 1e-28
# Input: barn (b)
# Output: square meter (m2)

NISTsqrMeterTObarn<-function(sqrMeter) sqrMeter / 1e-28
# Input: square meter (m2)
# Output: barn (b)

NISTbarrelTOcubMeter<-function(barrel) barrel * 0.1589873
# Input: barrel [for petroleum, 42 gallons (U.S.)] (bbl)
# Output: cubic meter (m3)

NISTcubMeterTObarrel<-function(cubMeter) cubMeter / 0.1589873
# Input: cubic meter (m3)
# Output: barrel [for petroleum, 42 gallons (U.S.)] (bbl)

NISTbarrelTOliter<-function(barrel) barrel * 158.9873
# Input: barrel [for petroleum, 42 gallons (U.S.)] (bbl)
# Output: liter (L)

NISTliterTObarrel<-function(liter) liter / 158.9873
# Input: liter (L)
# Output: barrel [for petroleum, 42 gallons (U.S.)] (bbl)

NISTbiotTOampere<-function(biot) biot * 10
# Input: biot (Bi)
# Output: ampere (A)

NISTampereTObiot<-function(ampere) ampere / 10
# Input: ampere (A)
# Output: biot (Bi)

NISTukThUnITtOjoule<-function(ukThUnIT) ukThUnIT * 1055.056
# Input: British thermal unitIT (BtuIT) 9
# Output: joule (J)

NISTjouleTOukThUnIT<-function(joule) joule / 1055.056
# Input: joule (J)
# Output: British thermal unitIT (BtuIT) 9

NISTukThUnthTOjoule<-function(ukThUnth) ukThUnth * 1054.35
# Input: British thermal unitth (Btuth) 9
# Output: joule (J)

NISTjouleTOukThUnth<-function(joule) joule / 1054.35
# Input: joule (J)
# Output: British thermal unitth (Btuth) 9

NISTukThUnTOjoule<-function(ukThUn) ukThUn * 1055.87
# Input: British thermal unit (mean) (Btu)
# Output: joule (J)

NISTjouleTOukThUn<-function(joule) joule / 1055.87
# Input: joule (J)
# Output: British thermal unit (mean) (Btu)

NISTukThUn39FtOjoule<-function(ukThUn39F) ukThUn39F * 1059.67
# Input: British thermal unit (39 F) (Btu)
# Output: joule (J)

NISTjouleTOukThUn39F<-function(joule) joule / 1059.67
# Input: joule (J)
# Output: British thermal unit (39 F) (Btu)

NISTukThUn59FtOjoule<-function(ukThUn59F) ukThUn59F * 1054.8
# Input: British thermal unit (59 F) (Btu)
# Output: joule (J)

NISTjouleTOukThUn59F<-function(joule) joule / 1054.8
# Input: joule (J)
# Output: British thermal unit (59 F) (Btu)

NISTukThUn60FtOjoule<-function(ukThUn60F) ukThUn60F * 1054.68
# Input: British thermal unit (60 F) (Btu)
# Output: joule (J)

NISTjouleTOukThUn60F<-function(joule) joule / 1054.68
# Input: joule (J)
# Output: British thermal unit (60 F) (Btu)

NISTukThUnITFtPerHourSqrFtDegFtOwattPerMeterK<-function(ukThUnITFtPerHourSqrFtDegF) ukThUnITFtPerHourSqrFtDegF * 1.730735
# Input: British thermal unitIT foot per hour square foot degree Fahrenheit [BtuIT * ft/(h * ft2 * F)]
# Output: watt per meter kelvin [W/(m * K)]

NISTwattPerMeterKtOukThUnITFtPerHourSqrFtDegF<-function(wattPerMeterK) wattPerMeterK / 1.730735
# Input: watt per meter kelvin [W/(m * K)]
# Output: British thermal unitIT foot per hour square foot degree Fahrenheit [BtuIT * ft/(h * ft2 * F)]

NISTukThUnthFtPerHourSqrFtDegFtOwattPerMeterK<-function(ukThUnthFtPerHourSqrFtDegF) ukThUnthFtPerHourSqrFtDegF * 1.729577
# Input: British thermal unitth foot per hour square foot degree Fahrenheit [Btuth * ft/(h * ft2 * F)]
# Output: watt per meter kelvin [W/(m * K)]

NISTwattPerMeterKtOukThUnthFtPerHourSqrFtDegF<-function(wattPerMeterK) wattPerMeterK / 1.729577
# Input: watt per meter kelvin [W/(m * K)]
# Output: British thermal unitth foot per hour square foot degree Fahrenheit [Btuth * ft/(h * ft2 * F)]

NISTukThUnITInchPerHourSqrFtDegFtOwattPerMeterK<-function(ukThUnITInchPerHourSqrFtDegF) ukThUnITInchPerHourSqrFtDegF * 0.1442279
# Input: British thermal unitIT inch per hour square foot degree Fahrenheit [BtuIT * in/(h * ft2 * F)]
# Output: watt per meter kelvin [W/(m * K)]

NISTwattPerMeterKtOukThUnITInchPerHourSqrFtDegF<-function(wattPerMeterK) wattPerMeterK / 0.1442279
# Input: watt per meter kelvin [W/(m * K)]
# Output: British thermal unitIT inch per hour square foot degree Fahrenheit [BtuIT * in/(h * ft2 * F)]

NISTukThUnthInchPerHourSqrFtDegFtOwattPerMeterK<-function(ukThUnthInchPerHourSqrFtDegF) ukThUnthInchPerHourSqrFtDegF * 0.1441314
# Input: British thermal unitth inch per hour square foot degree Fahrenheit [Btuth * in/(h * ft2 * F)]
# Output: watt per meter kelvin [W / (m * K)]

NISTwattPerMeterKtOukThUnthInchPerHourSqrFtDegF<-function(wattPerMeterK) wattPerMeterK / 0.1441314
# Input: watt per meter kelvin [W / (m * K)]
# Output: British thermal unitth inch per hour square foot degree Fahrenheit [Btuth * in/(h * ft2 * F)]

NISTukThUnITInchPerSecSqrFtDegFtOwattPerMeterK<-function(ukThUnITInchPerSecSqrFtDegF) ukThUnITInchPerSecSqrFtDegF * 519.2204
# Input: British thermal unitIT inch per second square foot degree Fahrenheit [BtuIT * in /(s * ft2 * F)]
# Output: watt per meter kelvin [W/(m * K)]

NISTwattPerMeterKtOukThUnITInchPerSecSqrFtDegF<-function(wattPerMeterK) wattPerMeterK / 519.2204
# Input: watt per meter kelvin [W/(m * K)]
# Output: British thermal unitIT inch per second square foot degree Fahrenheit [BtuIT * in /(s * ft2 * F)]

NISTukThUnthInchPerSecSqrFtDegFtOwattPerMeterK<-function(ukThUnthInchPerSecSqrFtDegF) ukThUnthInchPerSecSqrFtDegF * 518.8732
# Input: British thermal unitth inch per second square foot degree Fahrenheit [Btuth * in /(s * ft2 * F)]
# Output: watt per meter kelvin [W/(m * K)]

NISTwattPerMeterKtOukThUnthInchPerSecSqrFtDegF<-function(wattPerMeterK) wattPerMeterK / 518.8732
# Input: watt per meter kelvin [W/(m * K)]
# Output: British thermal unitth inch per second square foot degree Fahrenheit [Btuth * in /(s * ft2 * F)]

NISTukThUnITPerCubFtTOjoulePerCubMeter<-function(ukThUnITPerCubFt) ukThUnITPerCubFt * 37258.95
# Input: British thermal unitIT per cubic foot (BtuIT/ft3)
# Output: joule per cubic meter (J/m3)

NISTjoulePerCubMeterTOukThUnITPerCubFt<-function(joulePerCubMeter) joulePerCubMeter / 37258.95
# Input: joule per cubic meter (J/m3)
# Output: British thermal unitIT per cubic foot (BtuIT/ft3)

NISTukThUnthPerCubFtTOjoulePerCubMeter<-function(ukThUnthPerCubFt) ukThUnthPerCubFt * 37234.03
# Input: British thermal unitth per cubic foot (Btuth/ft3)
# Output: joule per cubic meter (J/m3)

NISTjoulePerCubMeterTOukThUnthPerCubFt<-function(joulePerCubMeter) joulePerCubMeter / 37234.03
# Input: joule per cubic meter (J/m3)
# Output: British thermal unitth per cubic foot (Btuth/ft3)

NISTukThUnITPerDegFtOjoulePerK<-function(ukThUnITPerDegF) ukThUnITPerDegF * 1899.101
# Input: British thermal unitIT per degree Fahrenheit (BtuIT/F)
# Output: joule per kelvin (J/K)

NISTjoulePerKtOukThUnITPerDegF<-function(joulePerK) joulePerK / 1899.101
# Input: joule per kelvin (J/K)
# Output: British thermal unitIT per degree Fahrenheit (BtuIT/F)

NISTukThUnthPerDegFtOjoulePerK<-function(ukThUnthPerDegF) ukThUnthPerDegF * 1897.83
# Input: British thermal unitth per degree Fahrenheit (Btuth/F)
# Output: joule per kelvin (J/K)

NISTjoulePerKtOukThUnthPerDegF<-function(joulePerK) joulePerK / 1897.83
# Input: joule per kelvin (J/K)
# Output: British thermal unitth per degree Fahrenheit (Btuth/F)

NISTukThUnITPerDegRankineTOjoulePerK<-function(ukThUnITPerDegRankine) ukThUnITPerDegRankine * 1899.101
# Input: British thermal unitIT per degree Rankine (BtuIT/R)
# Output: joule per kelvin (J/K)

NISTjoulePerKtOukThUnITPerDegRankine<-function(joulePerK) joulePerK / 1899.101
# Input: joule per kelvin (J/K)
# Output: British thermal unitIT per degree Rankine (BtuIT/R)

NISTukThUnthPerDegRankineTOjoulePerK<-function(ukThUnthPerDegRankine) ukThUnthPerDegRankine * 1897.83
# Input: British thermal unitth per degree Rankine (Btuth/R)
# Output: joule per kelvin (J/K)

NISTjoulePerKtOukThUnthPerDegRankine<-function(joulePerK) joulePerK / 1897.83
# Input: joule per kelvin (J/K)
# Output: British thermal unitth per degree Rankine (Btuth/R)

NISTukThUnITPerHourTOwatt<-function(ukThUnITPerHour) ukThUnITPerHour * 0.2930711
# Input: British thermal unitIT per hour (BtuIT/h)
# Output: watt (W)

NISTwattTOukThUnITPerHour<-function(watt) watt / 0.2930711
# Input: watt (W)
# Output: British thermal unitIT per hour (BtuIT/h)

NISTukThUnthPerHourTOwatt<-function(ukThUnthPerHour) ukThUnthPerHour * 0.2928751
# Input: British thermal unitth per hour (Btuth/h)
# Output: watt (W)

NISTwattTOukThUnthPerHour<-function(watt) watt / 0.2928751
# Input: watt (W)
# Output: British thermal unitth per hour (Btuth/h)

NISTukThUnITPerHourSqrFtDegFtOwattPerSqrMeterK<-function(ukThUnITPerHourSqrFtDegF) ukThUnITPerHourSqrFtDegF * 5.678263
# Input: British thermal unitIT per hour square foot degree Fahrenheit [BtuIT/(h * ft2 * F)]
# Output: watt per square meter kelvin [W/(m2 * K)]

NISTwattPerSqrMeterKtOukThUnITPerHourSqrFtDegF<-function(wattPerSqrMeterK) wattPerSqrMeterK / 5.678263
# Input: watt per square meter kelvin [W/(m2 * K)]
# Output: British thermal unitIT per hour square foot degree Fahrenheit [BtuIT/(h * ft2 * F)]

NISTukThUnthPerHourSqrFtDegFtOwattPerSqrMeterK<-function(ukThUnthPerHourSqrFtDegF) ukThUnthPerHourSqrFtDegF * 5.674466
# Input: British thermal unitth per hour square foot degree Fahrenheit [Btuth/(h * ft2 * F)]
# Output: watt per square meter kelvin [W/(m2 * K)]

NISTwattPerSqrMeterKtOukThUnthPerHourSqrFtDegF<-function(wattPerSqrMeterK) wattPerSqrMeterK / 5.674466
# Input: watt per square meter kelvin [W/(m2 * K)]
# Output: British thermal unitth per hour square foot degree Fahrenheit [Btuth/(h * ft2 * F)]

NISTukThUnthPerMinTOwatt<-function(ukThUnthPerMin) ukThUnthPerMin * 17.5725
# Input: British thermal unitth per minute (Btuth/min)
# Output: watt (W)

NISTwattTOukThUnthPerMin<-function(watt) watt / 17.5725
# Input: watt (W)
# Output: British thermal unitth per minute (Btuth/min)

NISTukThUnITPerPoundTOjoulePerKg<-function(ukThUnITPerPound) ukThUnITPerPound * 2326
# Input: British thermal unitIT per pound (BtuIT/lb)
# Output: joule per kilogram (J/kg)

NISTjoulePerKgTOukThUnITPerPound<-function(joulePerKg) joulePerKg / 2326
# Input: joule per kilogram (J/kg)
# Output: British thermal unitIT per pound (BtuIT/lb)

NISTukThUnthPerPoundTOjoulePerKg<-function(ukThUnthPerPound) ukThUnthPerPound * 2324.444
# Input: British thermal unitth per pound (Btuth/lb)
# Output: joule per kilogram (J/kg)

NISTjoulePerKgTOukThUnthPerPound<-function(joulePerKg) joulePerKg / 2324.444
# Input: joule per kilogram (J/kg)
# Output: British thermal unitth per pound (Btuth/lb)

NISTukThUnITPerPoundDegFtOjoulePerKgK<-function(ukThUnITPerPoundDegF) ukThUnITPerPoundDegF * 4186.8
# Input: British thermal unitIT per pound degree Fahrenheit [BtuIT/(lb * F)]
# Output: joule per kilogram kelvin [J/(kg * K)]

NISTjoulePerKgKtOukThUnITPerPoundDegF<-function(joulePerKgK) joulePerKgK / 4186.8
# Input: joule per kilogram kelvin [J/(kg * K)]
# Output: British thermal unitIT per pound degree Fahrenheit [BtuIT/(lb * F)]

NISTukThUnthPerPoundDegFtOjoulePerKgK<-function(ukThUnthPerPoundDegF) ukThUnthPerPoundDegF * 4184
# Input: British thermal unitth per pound degree Fahrenheit [Btuth/(lb * F)]
# Output: joule per kilogram kelvin [J/(kg * K)]

NISTjoulePerKgKtOukThUnthPerPoundDegF<-function(joulePerKgK) joulePerKgK / 4184
# Input: joule per kilogram kelvin [J/(kg * K)]
# Output: British thermal unitth per pound degree Fahrenheit [Btuth/(lb * F)]

NISTukThUnITPerPoundDegRankineTOjoulePerKgK<-function(ukThUnITPerPoundDegRankine) ukThUnITPerPoundDegRankine * 4186.8
# Input: British thermal unitIT per pound degree Rankine [BtuIT/(lb * R)]
# Output: joule per kilogram kelvin [J/(kg * K)]

NISTjoulePerKgKtOukThUnITPerPoundDegRankine<-function(joulePerKgK) joulePerKgK / 4186.8
# Input: joule per kilogram kelvin [J/(kg * K)]
# Output: British thermal unitIT per pound degree Rankine [BtuIT/(lb * R)]

NISTukThUnthPerPoundDegRankineTOjoulePerKgK<-function(ukThUnthPerPoundDegRankine) ukThUnthPerPoundDegRankine * 4184
# Input: British thermal unitth per pound degree Rankine [Btuth/(lb * R)]
# Output: joule per kilogram kelvin [J/(kg * K)]

NISTjoulePerKgKtOukThUnthPerPoundDegRankine<-function(joulePerKgK) joulePerKgK / 4184
# Input: joule per kilogram kelvin [J/(kg * K)]
# Output: British thermal unitth per pound degree Rankine [Btuth/(lb * R)]

NISTukThUnITPerSecTOwatt<-function(ukThUnITPerSec) ukThUnITPerSec * 1055.056
# Input: British thermal unitIT per second (BtuIT/s)
# Output: watt (W)

NISTwattTOukThUnITPerSec<-function(watt) watt / 1055.056
# Input: watt (W)
# Output: British thermal unitIT per second (BtuIT/s)

NISTukThUnthPerSecTOwatt<-function(ukThUnthPerSec) ukThUnthPerSec * 1054.35
# Input: British thermal unitth per second (Btuth/s)
# Output: watt (W)

NISTwattTOukThUnthPerSec<-function(watt) watt / 1054.35
# Input: watt (W)
# Output: British thermal unitth per second (Btuth/s)

NISTukThUnITPerSecSqrFtDegFtOwattPerSqrMeterK<-function(ukThUnITPerSecSqrFtDegF) ukThUnITPerSecSqrFtDegF * 20441.75
# Input: British thermal unitIT per second square foot degree Fahrenheit [BtuIT/(s * ft2 * F)]
# Output: watt per square meter kelvin [W/(m2 * K)]

NISTwattPerSqrMeterKtOukThUnITPerSecSqrFtDegF<-function(wattPerSqrMeterK) wattPerSqrMeterK / 20441.75
# Input: watt per square meter kelvin [W/(m2 * K)]
# Output: British thermal unitIT per second square foot degree Fahrenheit [BtuIT/(s * ft2 * F)]

NISTukThUnthPerSecSqrFtDegFtOwattPerSqrMeterK<-function(ukThUnthPerSecSqrFtDegF) ukThUnthPerSecSqrFtDegF * 20428.08
# Input: British thermal unitth per second square foot degree Fahrenheit [Btuth/(s * ft2 * F)]
# Output: watt per square meter kelvin [W/(m2 * K)]

NISTwattPerSqrMeterKtOukThUnthPerSecSqrFtDegF<-function(wattPerSqrMeterK) wattPerSqrMeterK / 20428.08
# Input: watt per square meter kelvin [W/(m2 * K)]
# Output: British thermal unitth per second square foot degree Fahrenheit [Btuth/(s * ft2 * F)]

NISTukThUnITPerSqrFtTOjoulePerSqrMeter<-function(ukThUnITPerSqrFt) ukThUnITPerSqrFt * 11356.53
# Input: British thermal unitIT per square foot (BtuIT/ft2)
# Output: joule per square meter (J/m2)

NISTjoulePerSqrMeterTOukThUnITPerSqrFt<-function(joulePerSqrMeter) joulePerSqrMeter / 11356.53
# Input: joule per square meter (J/m2)
# Output: British thermal unitIT per square foot (BtuIT/ft2)

NISTukThUnthPerSqrFtTOjoulePerSqrMeter<-function(ukThUnthPerSqrFt) ukThUnthPerSqrFt * 11348.93
# Input: British thermal unitth per square foot (Btuth/ft2)
# Output: joule per square meter (J/m2)

NISTjoulePerSqrMeterTOukThUnthPerSqrFt<-function(joulePerSqrMeter) joulePerSqrMeter / 11348.93
# Input: joule per square meter (J/m2)
# Output: British thermal unitth per square foot (Btuth/ft2)

NISTukThUnITPerSqrFtHourTOwattPerSqrMeter<-function(ukThUnITPerSqrFtHour) ukThUnITPerSqrFtHour * 3.154591
# Input: British thermal unitIT per square foot hour [BtuIT/(ft2 * h)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOukThUnITPerSqrFtHour<-function(wattPerSqrMeter) wattPerSqrMeter / 3.154591
# Input: watt per square meter (W/m2)
# Output: British thermal unitIT per square foot hour [BtuIT/(ft2 * h)]

NISTukThUnthPerSqrFtHourTOwattPerSqrMeter<-function(ukThUnthPerSqrFtHour) ukThUnthPerSqrFtHour * 3.152481
# Input: British thermal unitth per square foot hour [Btuth/(ft2 * h)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOukThUnthPerSqrFtHour<-function(wattPerSqrMeter) wattPerSqrMeter / 3.152481
# Input: watt per square meter (W/m2)
# Output: British thermal unitth per square foot hour [Btuth/(ft2 * h)]

NISTukThUnthPerSqrFtMinTOwattPerSqrMeter<-function(ukThUnthPerSqrFtMin) ukThUnthPerSqrFtMin * 189.1489
# Input: British thermal unitth per square foot minute [Btuth/(ft2 * min)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOukThUnthPerSqrFtMin<-function(wattPerSqrMeter) wattPerSqrMeter / 189.1489
# Input: watt per square meter (W/m2)
# Output: British thermal unitth per square foot minute [Btuth/(ft2 * min)]

NISTukThUnITPerSqrFtSecTOwattPerSqrMeter<-function(ukThUnITPerSqrFtSec) ukThUnITPerSqrFtSec * 11356.53
# Input: British thermal unitIT per square foot second [BtuIT/(ft2 * s)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOukThUnITPerSqrFtSec<-function(wattPerSqrMeter) wattPerSqrMeter / 11356.53
# Input: watt per square meter (W/m2)
# Output: British thermal unitIT per square foot second [BtuIT/(ft2 * s)]

NISTukThUnthPerSqrFtSecTOwattPerSqrMeter<-function(ukThUnthPerSqrFtSec) ukThUnthPerSqrFtSec * 11348.93
# Input: British thermal unitth per square foot second [Btuth/(ft2 * s)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOukThUnthPerSqrFtSec<-function(wattPerSqrMeter) wattPerSqrMeter / 11348.93
# Input: watt per square meter (W/m2)
# Output: British thermal unitth per square foot second [Btuth/(ft2 * s)]

NISTukThUnthPerSqrInchSecTOwattPerSqrMeter<-function(ukThUnthPerSqrInchSec) ukThUnthPerSqrInchSec * 1634246
# Input: British thermal unitth per square inch second [Btuth/(in2 * s)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOukThUnthPerSqrInchSec<-function(wattPerSqrMeter) wattPerSqrMeter / 1634246
# Input: watt per square meter (W/m2)
# Output: British thermal unitth per square inch second [Btuth/(in2 * s)]

NISTbushelTOcubMeter<-function(bushel) bushel * 0.03523907
# Input: bushel (U.S.) (bu)
# Output: cubic meter (m3)

NISTcubMeterTObushel<-function(cubMeter) cubMeter / 0.03523907
# Input: cubic meter (m3)
# Output: bushel (U.S.) (bu)

NISTbushelTOliter<-function(bushel) bushel * 35.23907
# Input: bushel (U.S.) (bu)
# Output: liter (L)

NISTliterTObushel<-function(liter) liter / 35.23907
# Input: liter (L)
# Output: bushel (U.S.) (bu)

NISTcalITtOjoule<-function(calIT) calIT * 4.1868
# Input: calorieIT (calIT) 10
# Output: joule (J)

NISTjouleTOcalIT<-function(joule) joule / 4.1868
# Input: joule (J)
# Output: calorieIT (calIT) 10

NISTcalthTOjoule<-function(calth) calth * 4.184
# Input: calorieth (calth) 10
# Output: joule (J)

NISTjouleTOcalth<-function(joule) joule / 4.184
# Input: joule (J)
# Output: calorieth (calth) 10

NISTcalTOjoule<-function(cal) cal * 4.19002
# Input: calorie (cal) (mean)
# Output: joule (J)

NISTjouleTOcal<-function(joule) joule / 4.19002
# Input: joule (J)
# Output: calorie (cal) (mean)

NISTcal15CtOjoule<-function(cal15C) cal15C * 4.1858
# Input: calorie (15 C) (cal15)
# Output: joule (J)

NISTjouleTOcal15C<-function(joule) joule / 4.1858
# Input: joule (J)
# Output: calorie (15 C) (cal15)

NISTcal20CtOjoule<-function(cal20C) cal20C * 4.1819
# Input: calorie (20 C) (cal20)
# Output: joule (J)

NISTjouleTOcal20C<-function(joule) joule / 4.1819
# Input: joule (J)
# Output: calorie (20 C) (cal20)

NISTcalITKgTOjoule<-function(calITKg) calITKg * 4186.8
# Input: calorieIT, kilogram (nutrition) 11
# Output: joule (J)

NISTjouleTOcalITKg<-function(joule) joule / 4186.8
# Input: joule (J)
# Output: calorieIT, kilogram (nutrition) 11

NISTcalthKgTOjoule<-function(calthKg) calthKg * 4184
# Input: calorieth , kilogram (nutrition) 11
# Output: joule (J)

NISTjouleTOcalthKg<-function(joule) joule / 4184
# Input: joule (J)
# Output: calorieth , kilogram (nutrition) 11

NISTcalKgTOjoule<-function(calKg) calKg * 4190.02
# Input: calorie (mean), kilogram (nutrition) 11
# Output: joule (J)

NISTjouleTOcalKg<-function(joule) joule / 4190.02
# Input: joule (J)
# Output: calorie (mean), kilogram (nutrition) 11

NISTcalthPerCmSecDegCtOwattPerMeterK<-function(calthPerCmSecDegC) calthPerCmSecDegC * 418.4
# Input: calorieth per centimeter second degree Celsius [calth/(cm * s * C)]
# Output: watt per meter kelvin [W/(m * K)]

NISTwattPerMeterKtOcalthPerCmSecDegC<-function(wattPerMeterK) wattPerMeterK / 418.4
# Input: watt per meter kelvin [W/(m * K)]
# Output: calorieth per centimeter second degree Celsius [calth/(cm * s * C)]

NISTcalITPerGramTOjoulePerKg<-function(calITPerGram) calITPerGram * 4186.8
# Input: calorieIT per gram (calIT/g)
# Output: joule per kilogram (J/kg)

NISTjoulePerKgTOcalITPerGram<-function(joulePerKg) joulePerKg / 4186.8
# Input: joule per kilogram (J/kg)
# Output: calorieIT per gram (calIT/g)

NISTcalthPerGramTOjoulePerKg<-function(calthPerGram) calthPerGram * 4184
# Input: calorieth per gram (calth/g)
# Output: joule per kilogram (J/kg)

NISTjoulePerKgTOcalthPerGram<-function(joulePerKg) joulePerKg / 4184
# Input: joule per kilogram (J/kg)
# Output: calorieth per gram (calth/g)

NISTcalITPerGramDegCtOjoulePerKgK<-function(calITPerGramDegC) calITPerGramDegC * 4186.8
# Input: calorieIT per gram degree Celsius [calIT/(g * C)]
# Output: joule per kilogram kelvin [J/(kg * K)]

NISTjoulePerKgKtOcalITPerGramDegC<-function(joulePerKgK) joulePerKgK / 4186.8
# Input: joule per kilogram kelvin [J/(kg * K)]
# Output: calorieIT per gram degree Celsius [calIT/(g * C)]

NISTcalthPerGramDegCtOjoulePerKgK<-function(calthPerGramDegC) calthPerGramDegC * 4184
# Input: calorieth per gram degree Celsius [calth/(g * C)]
# Output: joule per kilogram kelvin [J/(kg * K)]

NISTjoulePerKgKtOcalthPerGramDegC<-function(joulePerKgK) joulePerKgK / 4184
# Input: joule per kilogram kelvin [J/(kg * K)]
# Output: calorieth per gram degree Celsius [calth/(g * C)]

NISTcalITPerGramKtOjoulePerKgK<-function(calITPerGramK) calITPerGramK * 4186.8
# Input: calorieIT per gram kelvin [calIT/(g * K)]
# Output: joule per kilogram kelvin [J/(kg * K)]

NISTjoulePerKgKtOcalITPerGramK<-function(joulePerKgK) joulePerKgK / 4186.8
# Input: joule per kilogram kelvin [J/(kg * K)]
# Output: calorieIT per gram kelvin [calIT/(g * K)]

NISTcalthPerGramKtOjoulePerKgK<-function(calthPerGramK) calthPerGramK * 4184
# Input: calorieth per gram kelvin [calth/(g* K)]
# Output: joule per kilogram kelvin [J/(kg * K)]

NISTjoulePerKgKtOcalthPerGramK<-function(joulePerKgK) joulePerKgK / 4184
# Input: joule per kilogram kelvin [J/(kg * K)]
# Output: calorieth per gram kelvin [calth/(g* K)]

NISTcalthPerMinTOwatt<-function(calthPerMin) calthPerMin * 0.06973333
# Input: calorieth per minute (calth/min)
# Output: watt (W)

NISTwattTOcalthPerMin<-function(watt) watt / 0.06973333
# Input: watt (W)
# Output: calorieth per minute (calth/min)

NISTcalthPerSecTOwatt<-function(calthPerSec) calthPerSec * 4.184
# Input: calorieth per second (calth/s)
# Output: watt (W)

NISTwattTOcalthPerSec<-function(watt) watt / 4.184
# Input: watt (W)
# Output: calorieth per second (calth/s)

NISTcalthPerSqrCmTOjoulePerSqrMeter<-function(calthPerSqrCm) calthPerSqrCm * 41840
# Input: calorieth per square centimeter (calth/cm2)
# Output: joule per square meter (J/m2)

NISTjoulePerSqrMeterTOcalthPerSqrCm<-function(joulePerSqrMeter) joulePerSqrMeter / 41840
# Input: joule per square meter (J/m2)
# Output: calorieth per square centimeter (calth/cm2)

NISTcalthPerSqrCmMinTOwattPerSqrMeter<-function(calthPerSqrCmMin) calthPerSqrCmMin * 697.3333
# Input: calorieth per square centimeter minute [calth/(cm2 * min)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOcalthPerSqrCmMin<-function(wattPerSqrMeter) wattPerSqrMeter / 697.3333
# Input: watt per square meter (W/m2)
# Output: calorieth per square centimeter minute [calth/(cm2 * min)]

NISTcalthPerSqrCmSecTOwattPerSqrMeter<-function(calthPerSqrCmSec) calthPerSqrCmSec * 41840
# Input: calorieth per square centimeter second [calth/(cm2 * s)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOcalthPerSqrCmSec<-function(wattPerSqrMeter) wattPerSqrMeter / 41840
# Input: watt per square meter (W/m2)
# Output: calorieth per square centimeter second [calth/(cm2 * s)]

NISTcandelaPerSqrInchTOcandelaPerSqrMeter<-function(candelaPerSqrInch) candelaPerSqrInch * 1550.003
# Input: candela per square inch (cd/in2)
# Output: candela per square meter (cd/m2)

NISTcandelaPerSqrMeterTOcandelaPerSqrInch<-function(candelaPerSqrMeter) candelaPerSqrMeter / 1550.003
# Input: candela per square meter (cd/m2)
# Output: candela per square inch (cd/in2)

NISTcaratSItOkg<-function(caratSI) caratSI * 2e-04
# Input: carat, metric
# Output: kilogram (kg)

NISTkgTOcaratSI<-function(kg) kg / 2e-04
# Input: kilogram (kg)
# Output: carat, metric

NISTcaratSItOgram<-function(caratSI) caratSI * 0.2
# Input: carat, metric
# Output: gram (g)

NISTgramTOcaratSI<-function(gram) gram / 0.2
# Input: gram (g)
# Output: carat, metric

NISTcmOfMercuryTOpascal<-function(cmOfMercury) cmOfMercury * 1333.22
# Input: centimeter of mercury (0 C) 12
# Output: pascal (Pa)

NISTpascalTOcmOfMercury<-function(pascal) pascal / 1333.22
# Input: pascal (Pa)
# Output: centimeter of mercury (0 C) 12

NISTcmOfMercuryTOkpascal<-function(cmOfMercury) cmOfMercury * 1.33322
# Input: centimeter of mercury (0 C) 12
# Output: kilopascal (kPa)

NISTkpascalTOcmOfMercury<-function(kpascal) kpascal / 1.33322
# Input: kilopascal (kPa)
# Output: centimeter of mercury (0 C) 12

NISTcmOfMercuryConvtnlTOpascal<-function(cmOfMercuryConvtnl) cmOfMercuryConvtnl * 1333.224
# Input: centimeter of mercury, conventional (cmHg) 12
# Output: pascal (Pa)

NISTpascalTOcmOfMercuryConvtnl<-function(pascal) pascal / 1333.224
# Input: pascal (Pa)
# Output: centimeter of mercury, conventional (cmHg) 12

NISTcmOfMercuryConvtnlTOkpascal<-function(cmOfMercuryConvtnl) cmOfMercuryConvtnl * 1.333224
# Input: centimeter of mercury, conventional (cmHg) 12
# Output: kilopascal (kPa)

NISTkpascalTOcmOfMercuryConvtnl<-function(kpascal) kpascal / 1.333224
# Input: kilopascal (kPa)
# Output: centimeter of mercury, conventional (cmHg) 12

NISTcmOfWaterTOpascal<-function(cmOfWater) cmOfWater * 98.0638
# Input: centimeter of water (4 C) 12
# Output: pascal (Pa)

NISTpascalTOcmOfWater<-function(pascal) pascal / 98.0638
# Input: pascal (Pa)
# Output: centimeter of water (4 C) 12

NISTcmOfWaterConvtnlTOpascal<-function(cmOfWaterConvtnl) cmOfWaterConvtnl * 98.0665
# Input: centimeter of water, conventional (cmH2O) 12
# Output: pascal (Pa)

NISTpascalTOcmOfWaterConvtnl<-function(pascal) pascal / 98.0665
# Input: pascal (Pa)
# Output: centimeter of water, conventional (cmH2O) 12

NISTcentipoiseTOpascalSec<-function(centipoise) centipoise * 0.001
# Input: centipoise (cP)
# Output: pascal second (Pa * s)

NISTpascalSecTOcentipoise<-function(pascalSec) pascalSec / 0.001
# Input: pascal second (Pa * s)
# Output: centipoise (cP)

NISTcentistokesTOmeterSqrdPerSec<-function(centistokes) centistokes * 1e-06
# Input: centistokes (cSt)
# Output: meter squared per second (m2/s)

NISTmeterSqrdPerSecTOcentistokes<-function(meterSqrdPerSec) meterSqrdPerSec / 1e-06
# Input: meter squared per second (m2/s)
# Output: centistokes (cSt)

NISTchainTOmeter<-function(chain) chain * 20.11684
# Input: chain (based on U.S. survey foot) (ch) 7
# Output: meter (m)

NISTmeterTOchain<-function(meter) meter / 20.11684
# Input: meter (m)
# Output: chain (based on U.S. survey foot) (ch) 7

NISTcircularMilTOsqrMeter<-function(circularMil) circularMil * 5.067075e-10
# Input: circular mil
# Output: square meter (m2)

NISTsqrMeterTOcircularMil<-function(sqrMeter) sqrMeter / 5.067075e-10
# Input: square meter (m2)
# Output: circular mil

NISTcircularMilTOsqrMillimeter<-function(circularMil) circularMil * 0.0005067075
# Input: circular mil
# Output: square millimeter (mm2)

NISTsqrMillimeterTOcircularMil<-function(sqrMillimeter) sqrMillimeter / 0.0005067075
# Input: square millimeter (mm2)
# Output: circular mil

NISTcloTOsqrMeterKPerWatt<-function(clo) clo * 0.155
# Input: clo
# Output: square meter kelvin per watt (m2 * K/W)

NISTsqrMeterKPerWattTOclo<-function(sqrMeterKPerWatt) sqrMeterKPerWatt / 0.155
# Input: square meter kelvin per watt (m2 * K/W)
# Output: clo

NISTcordTOcubMeter<-function(cord) cord * 3.624556
# Input: cord (128 ft3)
# Output: cubic meter (m3)

NISTcubMeterTOcord<-function(cubMeter) cubMeter / 3.624556
# Input: cubic meter (m3)
# Output: cord (128 ft3)

NISTcubFtTOcubMeter<-function(cubFt) cubFt * 0.02831685
# Input: cubic foot (ft3)
# Output: cubic meter (m3)

NISTcubMeterTOcubFt<-function(cubMeter) cubMeter / 0.02831685
# Input: cubic meter (m3)
# Output: cubic foot (ft3)

NISTcubFtPerMinTOcubMeterPerSec<-function(cubFtPerMin) cubFtPerMin * 0.0004719474
# Input: cubic foot per minute (ft3/min)
# Output: cubic meter per second (m3/s)

NISTcubMeterPerSecTOcubFtPerMin<-function(cubMeterPerSec) cubMeterPerSec / 0.0004719474
# Input: cubic meter per second (m3/s)
# Output: cubic foot per minute (ft3/min)

NISTcubFtPerMinTOliterPerSec<-function(cubFtPerMin) cubFtPerMin * 0.4719474
# Input: cubic foot per minute (ft3/min)
# Output: liter per second (L/s)

NISTliterPerSecTOcubFtPerMin<-function(literPerSec) literPerSec / 0.4719474
# Input: liter per second (L/s)
# Output: cubic foot per minute (ft3/min)

NISTcubFtPerSecTOcubMeterPerSec<-function(cubFtPerSec) cubFtPerSec * 0.02831685
# Input: cubic foot per second (ft3/s)
# Output: cubic meter per second (m3/s)

NISTcubMeterPerSecTOcubFtPerSec<-function(cubMeterPerSec) cubMeterPerSec / 0.02831685
# Input: cubic meter per second (m3/s)
# Output: cubic foot per second (ft3/s)

NISTcubInchTOcubMeter<-function(cubInch) cubInch * 1.638706e-05
# Input: cubic inch (in3) 13
# Output: cubic meter (m3)

NISTcubMeterTOcubInch<-function(cubMeter) cubMeter / 1.638706e-05
# Input: cubic meter (m3)
# Output: cubic inch (in3) 13

NISTcubInchPerMinTOcubMeterPerSec<-function(cubInchPerMin) cubInchPerMin * 2.731177e-07
# Input: cubic inch per minute (in3/min)
# Output: cubic meter per second (m3/s)

NISTcubMeterPerSecTOcubInchPerMin<-function(cubMeterPerSec) cubMeterPerSec / 2.731177e-07
# Input: cubic meter per second (m3/s)
# Output: cubic inch per minute (in3/min)

NISTcubMileTOcubMeter<-function(cubMile) cubMile * 4168182000
# Input: cubic mile (mi3)
# Output: cubic meter (m3)

NISTcubMeterTOcubMile<-function(cubMeter) cubMeter / 4168182000
# Input: cubic meter (m3)
# Output: cubic mile (mi3)

NISTcubYardTOcubMeter<-function(cubYard) cubYard * 0.7645549
# Input: cubic yard (yd3)
# Output: cubic meter (m3)

NISTcubMeterTOcubYard<-function(cubMeter) cubMeter / 0.7645549
# Input: cubic meter (m3)
# Output: cubic yard (yd3)

NISTcubYardPerMinTOcubMeterPerSec<-function(cubYardPerMin) cubYardPerMin * 0.01274258
# Input: cubic yard per minute (yd3/min)
# Output: cubic meter per second (m3/s)

NISTcubMeterPerSecTOcubYardPerMin<-function(cubMeterPerSec) cubMeterPerSec / 0.01274258
# Input: cubic meter per second (m3/s)
# Output: cubic yard per minute (yd3/min)

NISTcupTOcubMeter<-function(cup) cup * 0.0002365882
# Input: cup (U.S.)
# Output: cubic meter (m3)

NISTcubMeterTOcup<-function(cubMeter) cubMeter / 0.0002365882
# Input: cubic meter (m3)
# Output: cup (U.S.)

NISTcupTOliter<-function(cup) cup * 0.2365882
# Input: cup (U.S.)
# Output: liter (L)

NISTliterTOcup<-function(liter) liter / 0.2365882
# Input: liter (L)
# Output: cup (U.S.)

NISTcupTOml<-function(cup) cup * 236.5882
# Input: cup (U.S.)
# Output: milliliter (mL)

NISTmlTOcup<-function(ml) ml / 236.5882
# Input: milliliter (mL)
# Output: cup (U.S.)

NISTcurieTObecquerel<-function(curie) curie * 3.7e+10
# Input: curie (Ci)
# Output: becquerel (Bq)

NISTbecquerelTOcurie<-function(becquerel) becquerel / 3.7e+10
# Input: becquerel (Bq)
# Output: curie (Ci)

NISTdarcyTOmeterSqrd<-function(darcy) darcy * 9.869233e-13
# Input: darcy 14
# Output: meter squared (m2)

NISTmeterSqrdTOdarcy<-function(meterSqrd) meterSqrd / 9.869233e-13
# Input: meter squared (m2)
# Output: darcy 14

NISTdayTOsec<-function(day) day * 86400
# Input: day (d)
# Output: second (s)

NISTsecTOday<-function(sec) sec / 86400
# Input: second (s)
# Output: day (d)

NISTdaySiderealTOsec<-function(daySidereal) daySidereal * 86164.09
# Input: day (sidereal)
# Output: second (s)

NISTsecTOdaySidereal<-function(sec) sec / 86164.09
# Input: second (s)
# Output: day (sidereal)

NISTdebyeTOcoulombMeter<-function(debye) debye * 3.335641e-30
# Input: debye (D)
# Output: coulomb meter (C * m)

NISTcoulombMeterTOdebye<-function(coulombMeter) coulombMeter / 3.335641e-30
# Input: coulomb meter (C * m)
# Output: debye (D)

NISTdegTOradian<-function(deg) deg * 0.01745329
# Input: degree (angle) ()
# Output: radian (rad)

NISTradianTOdeg<-function(radian) radian / 0.01745329
# Input: radian (rad)
# Output: degree (angle) ()

NISTdegCtOk<-function(degC) degC + 273.15 
# Input: degree Celsius (temperature) (C)
# Output: kelvin (K)

NISTkTOdegC<-function(k) k - 273.15 
# Input: kelvin (K)
# Output: degree Celsius (temperature) (C)

NISTdegCtOk<-function(degC) degC * 1
# Input: degree Celsius (temperature interval) (C)
# Output: kelvin (K)

NISTkTOdegC<-function(k) k / 1
# Input: kelvin (K)
# Output: degree Celsius (temperature interval) (C)

NISTdegCentigradeTOdegC<-function(degCentigrade) degCentigrade
# Input: degree centigrade (temperature) 15
# Output: degree Celsius (C)

NISTdegCtOdegCentigrade<-function(degC) degC
# Input: degree Celsius (C)
# Output: degree centigrade (temperature) 15

NISTdegCentigradeTOdegC<-function(degCentigrade) degCentigrade * 1
# Input: degree centigrade (temperature interval) 15
# Output: degree Celsius (C)

NISTdegCtOdegCentigrade<-function(degC) degC / 1
# Input: degree Celsius (C)
# Output: degree centigrade (temperature interval) 15

NISTdegFtOdegC<-function(degF) (degF - 32)/1.8  
# Input: degree Fahrenheit (temperature) (F)
# Output: degree Celsius (C)

NISTdegCtOdegF<-function(degC) (degC*1.8 + 32)
# Input: degree Celsius (C)
# Output: degree Fahrenheit (temperature) (F)

NISTdegFtOk<-function(degF) (degF + 459.67)/1.8  
# Input: degree Fahrenheit (temperature) (F)
# Output: kelvin (K)

NISTkTOdegF<-function(k) (k*1.8 - 459.67)
# Input: kelvin (K)
# Output: degree Fahrenheit (temperature) (F)

NISTdegFtOdegC<-function(degF) degF * 0.5555556
# Input: degree Fahrenheit (temperature interval) (F)
# Output: degree Celsius (C)

NISTdegCtOdegF<-function(degC) degC / 0.5555556
# Input: degree Celsius (C)
# Output: degree Fahrenheit (temperature interval) (F)

NISTdegFtOk<-function(degF) degF * 0.5555556
# Input: degree Fahrenheit (temperature interval) (F)
# Output: kelvin (K)

NISTkTOdegF<-function(k) k / 0.5555556
# Input: kelvin (K)
# Output: degree Fahrenheit (temperature interval) (F)

NISTdegFHourPerUkThUnITtOkPerWatt<-function(degFHourPerUkThUnIT) degFHourPerUkThUnIT * 1.895634
# Input: degree Fahrenheit hour per British thermal unitIT (F * h/BtuIT)
# Output: kelvin per watt (K/W)

NISTkPerWattTOdegFHourPerUkThUnIT<-function(kPerWatt) kPerWatt / 1.895634
# Input: kelvin per watt (K/W)
# Output: degree Fahrenheit hour per British thermal unitIT (F * h/BtuIT)

NISTdegFHourPerUkThUnthTOkPerWatt<-function(degFHourPerUkThUnth) degFHourPerUkThUnth * 1.896903
# Input: degree Fahrenheit hour per British thermal unitth (F * h/Btuth)
# Output: kelvin per watt (K/W)

NISTkPerWattTOdegFHourPerUkThUnth<-function(kPerWatt) kPerWatt / 1.896903
# Input: kelvin per watt (K/W)
# Output: degree Fahrenheit hour per British thermal unitth (F * h/Btuth)

NISTdegFHourSqrFtPerUkThUnITtOsqrMeterKPerWatt<-function(degFHourSqrFtPerUkThUnIT) degFHourSqrFtPerUkThUnIT * 0.1761102
# Input: degree Fahrenheit hour square foot per British thermal unitIT (F * h * ft2/BtuIT)
# Output: square meter kelvin per watt (m2 * K/W)

NISTsqrMeterKPerWattTOdegFHourSqrFtPerUkThUnIT<-function(sqrMeterKPerWatt) sqrMeterKPerWatt / 0.1761102
# Input: square meter kelvin per watt (m2 * K/W)
# Output: degree Fahrenheit hour square foot per British thermal unitIT (F * h * ft2/BtuIT)

NISTdegFHourSqrFtPerUkThUnthTOsqrMeterKPerWatt<-function(degFHourSqrFtPerUkThUnth) degFHourSqrFtPerUkThUnth * 0.176228
# Input: degree Fahrenheit hour square foot per British thermal unitth (F * h * ft2/Btuth)
# Output: square meter kelvin per watt (m2 * K/W)

NISTsqrMeterKPerWattTOdegFHourSqrFtPerUkThUnth<-function(sqrMeterKPerWatt) sqrMeterKPerWatt / 0.176228
# Input: square meter kelvin per watt (m2 * K/W)
# Output: degree Fahrenheit hour square foot per British thermal unitth (F * h * ft2/Btuth)

NISTdegFHourSqrFtPerUkThUnITInchTOmeterKPerWatt<-function(degFHourSqrFtPerUkThUnITInch) degFHourSqrFtPerUkThUnITInch * 6.933472
# Input: degree Fahrenheit hour square foot per British thermal unitIT inch [F * h * ft2/(BtuIT * in)]
# Output: meter kelvin per watt (m * K/W)

NISTmeterKPerWattTOdegFHourSqrFtPerUkThUnITInch<-function(meterKPerWatt) meterKPerWatt / 6.933472
# Input: meter kelvin per watt (m * K/W)
# Output: degree Fahrenheit hour square foot per British thermal unitIT inch [F * h * ft2/(BtuIT * in)]

NISTdegFHourSqrFtPerUkThUnthInchTOmeterKPerWatt<-function(degFHourSqrFtPerUkThUnthInch) degFHourSqrFtPerUkThUnthInch * 6.938112
# Input: degree Fahrenheit hour square foot per British thermal unitth inch [F * h * ft2/(Btuth * in)]
# Output: meter kelvin per watt (m * K/W)

NISTmeterKPerWattTOdegFHourSqrFtPerUkThUnthInch<-function(meterKPerWatt) meterKPerWatt / 6.938112
# Input: meter kelvin per watt (m * K/W)
# Output: degree Fahrenheit hour square foot per British thermal unitth inch [F * h * ft2/(Btuth * in)]

NISTdegFSecPerUkThUnITtOkPerWatt<-function(degFSecPerUkThUnIT) degFSecPerUkThUnIT * 0.0005265651
# Input: degree Fahrenheit second per British thermal unitIT (F * s/BtuIT)
# Output: kelvin per watt (K/W)

NISTkPerWattTOdegFSecPerUkThUnIT<-function(kPerWatt) kPerWatt / 0.0005265651
# Input: kelvin per watt (K/W)
# Output: degree Fahrenheit second per British thermal unitIT (F * s/BtuIT)

NISTdegFSecPerUkThUnthTOkPerWatt<-function(degFSecPerUkThUnth) degFSecPerUkThUnth * 0.0005269175
# Input: degree Fahrenheit second per British thermal unitth (F * s/Btuth)
# Output: kelvin per watt (K/W)

NISTkPerWattTOdegFSecPerUkThUnth<-function(kPerWatt) kPerWatt / 0.0005269175
# Input: kelvin per watt (K/W)
# Output: degree Fahrenheit second per British thermal unitth (F * s/Btuth)

NISTdegRankineTOk<-function(degRankine) degRankine / 1.8
# Input: degree Rankine (R)
# Output: kelvin (K)

NISTkTOdegRankine<-function(k) k * 1.8
# Input: kelvin (K)
# Output: degree Rankine (R)

NISTdegRankineTOk<-function(degRankine) degRankine * 0.5555556
# Input: degree Rankine (temperature interval) (R)
# Output: kelvin (K)

NISTkTOdegRankine<-function(k) k / 0.5555556
# Input: kelvin (K)
# Output: degree Rankine (temperature interval) (R)

NISTdenierTOkgPerMeter<-function(denier) denier * 1.111111e-07
# Input: denier
# Output: kilogram per meter (kg/m)

NISTkgPerMeterTOdenier<-function(kgPerMeter) kgPerMeter / 1.111111e-07
# Input: kilogram per meter (kg/m)
# Output: denier

NISTdenierTOgramPerMeter<-function(denier) denier * 0.0001111111
# Input: denier
# Output: gram per meter (g/m)

NISTgramPerMeterTOdenier<-function(gramPerMeter) gramPerMeter / 0.0001111111
# Input: gram per meter (g/m)
# Output: denier

NISTdyneTOnewton<-function(dyne) dyne * 1e-05
# Input: dyne (dyn)
# Output: newton (N)

NISTnewtonTOdyne<-function(newton) newton / 1e-05
# Input: newton (N)
# Output: dyne (dyn)

NISTdyneCmTOnewtonMeter<-function(dyneCm) dyneCm * 1e-07
# Input: dyne centimeter (dyn * cm)
# Output: newton meter (N * m)

NISTnewtonMeterTOdyneCm<-function(newtonMeter) newtonMeter / 1e-07
# Input: newton meter (N * m)
# Output: dyne centimeter (dyn * cm)

NISTdynePerSqrCmTOpascal<-function(dynePerSqrCm) dynePerSqrCm * 0.1
# Input: dyne per square centimeter (dyn/cm2)
# Output: pascal (Pa)

NISTpascalTOdynePerSqrCm<-function(pascal) pascal / 0.1
# Input: pascal (Pa)
# Output: dyne per square centimeter (dyn/cm2)

NISTelectronvoltTOjoule<-function(electronvolt) electronvolt * 1.602177e-19
# Input: electronvolt (eV)
# Output: joule (J)

NISTjouleTOelectronvolt<-function(joule) joule / 1.602177e-19
# Input: joule (J)
# Output: electronvolt (eV)

NISTEMUOfCapacitanceTOfarad<-function(EMUOfCapacitance) EMUOfCapacitance * 1e+09
# Input: EMU of capacitance (abfarad)
# Output: farad (F)

NISTfaradToEMUOfCapacitance<-function(farad) farad / 1e+09
# Input: farad (F)
# Output: EMU of capacitance (abfarad)

NISTEMUOfCurrentTOampere<-function(EMUOfCurrent) EMUOfCurrent * 10
# Input: EMU of current (abampere)
# Output: ampere (A)

NISTampereToEMUOfCurrent<-function(ampere) ampere / 10
# Input: ampere (A)
# Output: EMU of current (abampere)

NISTEMUOfElectricPotentialTOvolt<-function(EMUOfElectricPotential) EMUOfElectricPotential * 1e-08
# Input: EMU of electric potential (abvolt)
# Output: volt (V)

NISTvoltToEMUOfElectricPotential<-function(volt) volt / 1e-08
# Input: volt (V)
# Output: EMU of electric potential (abvolt)

NISTEMUOfInductanceTOhenry<-function(EMUOfInductance) EMUOfInductance * 1e-09
# Input: EMU of inductance (abhenry)
# Output: henry (H)

NISThenryToEMUOfInductance<-function(henry) henry / 1e-09
# Input: henry (H)
# Output: EMU of inductance (abhenry)

NISTEMUOfResistanceTOohm<-function(EMUOfResistance) EMUOfResistance * 1e-09
# Input: EMU of resistance (abohm)
# Output: ohm (ohm)

NISTohmToEMUOfResistance<-function(ohm) ohm / 1e-09
# Input: ohm (ohm)
# Output: EMU of resistance (abohm)

NISTergTOjoule<-function(erg) erg * 1e-07
# Input: erg (erg)
# Output: joule (J)

NISTjouleTOerg<-function(joule) joule / 1e-07
# Input: joule (J)
# Output: erg (erg)

NISTergPerSecTOwatt<-function(ergPerSec) ergPerSec * 1e-07
# Input: erg per second (erg/s)
# Output: watt (W)

NISTwattTOergPerSec<-function(watt) watt / 1e-07
# Input: watt (W)
# Output: erg per second (erg/s)

NISTergPerSqrCmSecTOwattPerSqrMeter<-function(ergPerSqrCmSec) ergPerSqrCmSec * 0.001
# Input: erg per square centimeter second [erg/(cm2 * s)]
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOergPerSqrCmSec<-function(wattPerSqrMeter) wattPerSqrMeter / 0.001
# Input: watt per square meter (W/m2)
# Output: erg per square centimeter second [erg/(cm2 * s)]

NISTESUOfCapacitanceTOfarad<-function(ESUOfCapacitance) ESUOfCapacitance * 1.11265e-12
# Input: ESU of capacitance (statfarad)
# Output: farad (F)

NISTfaradToESUOfCapacitance<-function(farad) farad / 1.11265e-12
# Input: farad (F)
# Output: ESU of capacitance (statfarad)

NISTESUOfCurrentTOampere<-function(ESUOfCurrent) ESUOfCurrent * 3.335641e-10
# Input: ESU of current (statampere)
# Output: ampere (A)

NISTampereToESUOfCurrent<-function(ampere) ampere / 3.335641e-10
# Input: ampere (A)
# Output: ESU of current (statampere)

NISTESUOfElectricPotentialTOvolt<-function(ESUOfElectricPotential) ESUOfElectricPotential * 299.7925
# Input: ESU of electric potential (statvolt)
# Output: volt (V)

NISTvoltToESUOfElectricPotential<-function(volt) volt / 299.7925
# Input: volt (V)
# Output: ESU of electric potential (statvolt)

NISTESUOfInductanceTOhenry<-function(ESUOfInductance) ESUOfInductance * 898755200000
# Input: ESU of inductance (stathenry)
# Output: henry (H)

NISThenryToESUOfInductance<-function(henry) henry / 898755200000
# Input: henry (H)
# Output: ESU of inductance (stathenry)

NISTESUOfResistanceTOohm<-function(ESUOfResistance) ESUOfResistance * 898755200000
# Input: ESU of resistance (statohm)
# Output: ohm (ohm)

NISTohmToESUOfResistance<-function(ohm) ohm / 898755200000
# Input: ohm (ohm)
# Output: ESU of resistance (statohm)

NISTfaradayTOcoulomb<-function(faraday) faraday * 96485.31
# Input: faraday (based on carbon 12)
# Output: coulomb (C)

NISTcoulombTOfaraday<-function(coulomb) coulomb / 96485.31
# Input: coulomb (C)
# Output: faraday (based on carbon 12)

NISTfathomTOmeter<-function(fathom) fathom * 1.828804
# Input: fathom (based on U.S. survey foot) 7
# Output: meter (m)

NISTmeterTOfathom<-function(meter) meter / 1.828804
# Input: meter (m)
# Output: fathom (based on U.S. survey foot) 7

NISTfermiTOmeter<-function(fermi) fermi * 1e-15
# Input: fermi
# Output: meter (m)

NISTmeterTOfermi<-function(meter) meter / 1e-15
# Input: meter (m)
# Output: fermi

NISTfermiTOfemtometer<-function(fermi) fermi * 1
# Input: fermi
# Output: femtometer (fm)

NISTfemtometerTOfermi<-function(femtometer) femtometer / 1
# Input: femtometer (fm)
# Output: fermi

NISTfluidOunceTOcubMeter<-function(fluidOunce) fluidOunce * 2.957353e-05
# Input: fluid ounce (U.S.) (fl oz)
# Output: cubic meter (m3)

NISTcubMeterTOfluidOunce<-function(cubMeter) cubMeter / 2.957353e-05
# Input: cubic meter (m3)
# Output: fluid ounce (U.S.) (fl oz)

NISTfluidOunceTOml<-function(fluidOunce) fluidOunce * 29.57353
# Input: fluid ounce (U.S.) (fl oz)
# Output: milliliter (mL)

NISTmlTOfluidOunce<-function(ml) ml / 29.57353
# Input: milliliter (mL)
# Output: fluid ounce (U.S.) (fl oz)

NISTftTOmeter<-function(ft) ft * 0.3048
# Input: foot (ft)
# Output: meter (m)

NISTmeterTOft<-function(meter) meter / 0.3048
# Input: meter (m)
# Output: foot (ft)

NISTftUSsurveyTOmeter<-function(ftUSsurvey) ftUSsurvey * 0.3048006
# Input: foot (U.S. survey) (ft) 7
# Output: meter (m)

NISTmeterTOftUSsurvey<-function(meter) meter / 0.3048006
# Input: meter (m)
# Output: foot (U.S. survey) (ft) 7

NISTftcandleTOlux<-function(ftcandle) ftcandle * 10.76391
# Input: footcandle
# Output: lux (lx)

NISTluxTOftcandle<-function(lux) lux / 10.76391
# Input: lux (lx)
# Output: footcandle

NISTftlambertTOcandelaPerSqrMeter<-function(ftlambert) ftlambert * 3.426259
# Input: footlambert
# Output: candela per square meter (cd/m2)

NISTcandelaPerSqrMeterTOftlambert<-function(candelaPerSqrMeter) candelaPerSqrMeter / 3.426259
# Input: candela per square meter (cd/m2)
# Output: footlambert

NISTftOfMercuryConvtnlTOpascal<-function(ftOfMercuryConvtnl) ftOfMercuryConvtnl * 40636.66
# Input: foot of mercury, conventional (ftHg) 12
# Output: pascal (Pa)

NISTpascalTOftOfMercuryConvtnl<-function(pascal) pascal / 40636.66
# Input: pascal (Pa)
# Output: foot of mercury, conventional (ftHg) 12

NISTftOfMercuryConvtnlTOkpascal<-function(ftOfMercuryConvtnl) ftOfMercuryConvtnl * 40.63666
# Input: foot of mercury, conventional (ftHg) 12
# Output: kilopascal (kPa)

NISTkpascalTOftOfMercuryConvtnl<-function(kpascal) kpascal / 40.63666
# Input: kilopascal (kPa)
# Output: foot of mercury, conventional (ftHg) 12

NISTftOfWaterTOpascal<-function(ftOfWater) ftOfWater * 2988.98
# Input: foot of water (39.2 F) 12
# Output: pascal (Pa)

NISTpascalTOftOfWater<-function(pascal) pascal / 2988.98
# Input: pascal (Pa)
# Output: foot of water (39.2 F) 12

NISTftOfWaterTOkpascal<-function(ftOfWater) ftOfWater * 2.98898
# Input: foot of water (39.2 F) 12
# Output: kilopascal (kPa)

NISTkpascalTOftOfWater<-function(kpascal) kpascal / 2.98898
# Input: kilopascal (kPa)
# Output: foot of water (39.2 F) 12

NISTftOfWaterConvtnlTOpascal<-function(ftOfWaterConvtnl) ftOfWaterConvtnl * 2989.067
# Input: foot of water, conventional (ftH2O) 12
# Output: pascal (Pa)

NISTpascalTOftOfWaterConvtnl<-function(pascal) pascal / 2989.067
# Input: pascal (Pa)
# Output: foot of water, conventional (ftH2O) 12

NISTftOfWaterConvtnlTOkpascal<-function(ftOfWaterConvtnl) ftOfWaterConvtnl * 2.989067
# Input: foot of water, conventional (ftH2O) 12
# Output: kilopascal (kPa)

NISTkpascalTOftOfWaterConvtnl<-function(kpascal) kpascal / 2.989067
# Input: kilopascal (kPa)
# Output: foot of water, conventional (ftH2O) 12

NISTftPerHourTOmeterPerSec<-function(ftPerHour) ftPerHour * 8.466667e-05
# Input: foot per hour (ft/h)
# Output: meter per second (m/s)

NISTmeterPerSecTOftPerHour<-function(meterPerSec) meterPerSec / 8.466667e-05
# Input: meter per second (m/s)
# Output: foot per hour (ft/h)

NISTftPerMinTOmeterPerSec<-function(ftPerMin) ftPerMin * 0.00508
# Input: foot per minute (ft/min)
# Output: meter per second (m/s)

NISTmeterPerSecTOftPerMin<-function(meterPerSec) meterPerSec / 0.00508
# Input: meter per second (m/s)
# Output: foot per minute (ft/min)

NISTftPerSecTOmeterPerSec<-function(ftPerSec) ftPerSec * 0.3048
# Input: foot per second (ft/s)
# Output: meter per second (m/s)

NISTmeterPerSecTOftPerSec<-function(meterPerSec) meterPerSec / 0.3048
# Input: meter per second (m/s)
# Output: foot per second (ft/s)

NISTftPerSecSqrdTOmeterPerSecSqrd<-function(ftPerSecSqrd) ftPerSecSqrd * 0.3048
# Input: foot per second squared (ft/s2)
# Output: meter per second squared (m/s2)

NISTmeterPerSecSqrdTOftPerSecSqrd<-function(meterPerSecSqrd) meterPerSecSqrd / 0.3048
# Input: meter per second squared (m/s2)
# Output: foot per second squared (ft/s2)

NISTftPoundalTOjoule<-function(ftPoundal) ftPoundal * 0.04214011
# Input: foot poundal
# Output: joule (J)

NISTjouleTOftPoundal<-function(joule) joule / 0.04214011
# Input: joule (J)
# Output: foot poundal

NISTftPoundForceTOjoule<-function(ftPoundForce) ftPoundForce * 1.355818
# Input: foot pound-force (ft * lbf)
# Output: joule (J)

NISTjouleTOftPoundForce<-function(joule) joule / 1.355818
# Input: joule (J)
# Output: foot pound-force (ft * lbf)

NISTftPoundForcePerHourTOwatt<-function(ftPoundForcePerHour) ftPoundForcePerHour * 0.0003766161
# Input: foot pound-force per hour (ft * lbf/h)
# Output: watt (W)

NISTwattTOftPoundForcePerHour<-function(watt) watt / 0.0003766161
# Input: watt (W)
# Output: foot pound-force per hour (ft * lbf/h)

NISTftPoundForcePerMinTOwatt<-function(ftPoundForcePerMin) ftPoundForcePerMin * 0.02259697
# Input: foot pound-force per minute (ft * lbf/min)
# Output: watt (W)

NISTwattTOftPoundForcePerMin<-function(watt) watt / 0.02259697
# Input: watt (W)
# Output: foot pound-force per minute (ft * lbf/min)

NISTftPoundForcePerSecTOwatt<-function(ftPoundForcePerSec) ftPoundForcePerSec * 1.355818
# Input: foot pound-force per second (ft * lbf/s)
# Output: watt (W)

NISTwattTOftPoundForcePerSec<-function(watt) watt / 1.355818
# Input: watt (W)
# Output: foot pound-force per second (ft * lbf/s)

NISTft4thPowerTOmeter4thPower<-function(ft4thPower) ft4thPower * 0.008630975
# Input: foot to the fourth power (ft4) 16
# Output: meter to the fourth power (m4)

NISTmeter4thPowerTOft4thPower<-function(meter4thPower) meter4thPower / 0.008630975
# Input: meter to the fourth power (m4)
# Output: foot to the fourth power (ft4) 16

NISTfranklinTOcoulomb<-function(franklin) franklin * 3.335641e-10
# Input: franklin (Fr)
# Output: coulomb (C)

NISTcoulombTOfranklin<-function(coulomb) coulomb / 3.335641e-10
# Input: coulomb (C)
# Output: franklin (Fr)

NISTgalTOmeterPerSecSqrd<-function(gal) gal * 0.01
# Input: gal (Gal)
# Output: meter per second squared (m/s2)

NISTmeterPerSecSqrdTOgal<-function(meterPerSecSqrd) meterPerSecSqrd / 0.01
# Input: meter per second squared (m/s2)
# Output: gal (Gal)

NISTgallonImperialTOcubMeter<-function(gallonImperial) gallonImperial * 0.00454609
# Input: gallon [Canadian and U.K. (Imperial)] (gal)
# Output: cubic meter (m3)

NISTcubMeterTOgallonImperial<-function(cubMeter) cubMeter / 0.00454609
# Input: cubic meter (m3)
# Output: gallon [Canadian and U.K. (Imperial)] (gal)

NISTgallonImperialTOliter<-function(gallonImperial) gallonImperial * 4.54609
# Input: gallon [Canadian and U.K. (Imperial)] (gal)
# Output: liter (L)

NISTliterTOgallonImperial<-function(liter) liter / 4.54609
# Input: liter (L)
# Output: gallon [Canadian and U.K. (Imperial)] (gal)

NISTgallonTOcubMeter<-function(gallon) gallon * 0.003785412
# Input: gallon (U.S.) (gal)
# Output: cubic meter (m3)

NISTcubMeterTOgallon<-function(cubMeter) cubMeter / 0.003785412
# Input: cubic meter (m3)
# Output: gallon (U.S.) (gal)

NISTgallonTOliter<-function(gallon) gallon * 3.785412
# Input: gallon (U.S.) (gal)
# Output: liter (L)

NISTliterTOgallon<-function(liter) liter / 3.785412
# Input: liter (L)
# Output: gallon (U.S.) (gal)

NISTgallonPerDayTOcubMeterPerSec<-function(gallonPerDay) gallonPerDay * 4.381264e-08
# Input: gallon (U.S.) per day (gal/d)
# Output: cubic meter per second (m3/s)

NISTcubMeterPerSecTOgallonPerDay<-function(cubMeterPerSec) cubMeterPerSec / 4.381264e-08
# Input: cubic meter per second (m3/s)
# Output: gallon (U.S.) per day (gal/d)

NISTgallonPerDayTOliterPerSec<-function(gallonPerDay) gallonPerDay * 4.381264e-05
# Input: gallon (U.S.) per day (gal/d)
# Output: liter per second (L/s)

NISTliterPerSecTOgallonPerDay<-function(literPerSec) literPerSec / 4.381264e-05
# Input: liter per second (L/s)
# Output: gallon (U.S.) per day (gal/d)

NISTgallonPerHorsepowerHourTOcubMeterPerJoule<-function(gallonPerHorsepowerHour) gallonPerHorsepowerHour * 1.410089e-09
# Input: gallon (U.S.) per horsepower hour [gal/(hp * h)]
# Output: cubic meter per joule (m3/J)

NISTcubMeterPerJouleTOgallonPerHorsepowerHour<-function(cubMeterPerJoule) cubMeterPerJoule / 1.410089e-09
# Input: cubic meter per joule (m3/J)
# Output: gallon (U.S.) per horsepower hour [gal/(hp * h)]

NISTgallonPerHorsepowerHourTOliterPerJoule<-function(gallonPerHorsepowerHour) gallonPerHorsepowerHour * 1.410089e-06
# Input: gallon (U.S.) per horsepower hour [gal/(hp * h)]
# Output: liter per joule (L/J)

NISTliterPerJouleTOgallonPerHorsepowerHour<-function(literPerJoule) literPerJoule / 1.410089e-06
# Input: liter per joule (L/J)
# Output: gallon (U.S.) per horsepower hour [gal/(hp * h)]

NISTgallonPerMinTOcubMeterPerSec<-function(gallonPerMin) gallonPerMin * 6.30902e-05
# Input: gallon (U.S.) per minute (gpm) (gal/min)
# Output: cubic meter per second (m3/s)

NISTcubMeterPerSecTOgallonPerMin<-function(cubMeterPerSec) cubMeterPerSec / 6.30902e-05
# Input: cubic meter per second (m3/s)
# Output: gallon (U.S.) per minute (gpm) (gal/min)

NISTgallonPerMinTOliterPerSec<-function(gallonPerMin) gallonPerMin * 0.0630902
# Input: gallon (U.S.) per minute (gpm) (gal/min)
# Output: liter per second (L/s)

NISTliterPerSecTOgallonPerMin<-function(literPerSec) literPerSec / 0.0630902
# Input: liter per second (L/s)
# Output: gallon (U.S.) per minute (gpm) (gal/min)

NISTgammaTOtesla<-function(gamma) gamma * 1e-09
# Input: gamma (y)
# Output: tesla (T)

NISTteslaTOgamma<-function(tesla) tesla / 1e-09
# Input: tesla (T)
# Output: gamma (y)

NISTgaussTOtesla<-function(gauss) gauss * 1e-04
# Input: gauss (Gs, G)
# Output: tesla (T)

NISTteslaTOgauss<-function(tesla) tesla / 1e-04
# Input: tesla (T)
# Output: gauss (Gs, G)

NISTgilbertTOampere<-function(gilbert) gilbert * 0.7957747
# Input: gilbert (Gi)
# Output: ampere (A)

NISTampereTOgilbert<-function(ampere) ampere / 0.7957747
# Input: ampere (A)
# Output: gilbert (Gi)

NISTgillImperialTOcubMeter<-function(gillImperial) gillImperial * 0.0001420653
# Input: gill [Canadian and U.K. (Imperial)] (gi)
# Output: cubic meter (m3)

NISTcubMeterTOgillImperial<-function(cubMeter) cubMeter / 0.0001420653
# Input: cubic meter (m3)
# Output: gill [Canadian and U.K. (Imperial)] (gi)

NISTgillImperialTOliter<-function(gillImperial) gillImperial * 0.1420653
# Input: gill [Canadian and U.K. (Imperial)] (gi)
# Output: liter (L)

NISTliterTOgillImperial<-function(liter) liter / 0.1420653
# Input: liter (L)
# Output: gill [Canadian and U.K. (Imperial)] (gi)

NISTgillTOcubMeter<-function(gill) gill * 0.0001182941
# Input: gill (U.S.) (gi)
# Output: cubic meter (m3)

NISTcubMeterTOgill<-function(cubMeter) cubMeter / 0.0001182941
# Input: cubic meter (m3)
# Output: gill (U.S.) (gi)

NISTgillTOliter<-function(gill) gill * 0.1182941
# Input: gill (U.S.) (gi)
# Output: liter (L)

NISTliterTOgill<-function(liter) liter / 0.1182941
# Input: liter (L)
# Output: gill (U.S.) (gi)

NISTgonTOradian<-function(gon) gon * 0.01570796
# Input: gon (also called grade) (gon)
# Output: radian (rad)

NISTradianTOgon<-function(radian) radian / 0.01570796
# Input: radian (rad)
# Output: gon (also called grade) (gon)

NISTgonTOdeg<-function(gon) gon * 0.9
# Input: gon (also called grade) (gon)
# Output: degree (angle) ()

NISTdegTOgon<-function(deg) deg / 0.9
# Input: degree (angle) ()
# Output: gon (also called grade) (gon)

NISTgrainTOkg<-function(grain) grain * 6.479891e-05
# Input: grain (gr)
# Output: kilogram (kg)

NISTkgTOgrain<-function(kg) kg / 6.479891e-05
# Input: kilogram (kg)
# Output: grain (gr)

NISTgrainTOmilligram<-function(grain) grain * 64.79891
# Input: grain (gr)
# Output: milligram (mg)

NISTmilligramTOgrain<-function(milligram) milligram / 64.79891
# Input: milligram (mg)
# Output: grain (gr)

NISTgrainPerGallonTOkgPerCubMeter<-function(grainPerGallon) grainPerGallon * 0.01711806
# Input: grain per gallon (U.S.) (gr/gal)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOgrainPerGallon<-function(kgPerCubMeter) kgPerCubMeter / 0.01711806
# Input: kilogram per cubic meter (kg/m3)
# Output: grain per gallon (U.S.) (gr/gal)

NISTgrainPerGallonTOmilligramPerLiter<-function(grainPerGallon) grainPerGallon * 17.11806
# Input: grain per gallon (U.S.) (gr/gal)
# Output: milligram per liter (mg/L)

NISTmilligramPerLiterTOgrainPerGallon<-function(milligramPerLiter) milligramPerLiter / 17.11806
# Input: milligram per liter (mg/L)
# Output: grain per gallon (U.S.) (gr/gal)

NISTgramForcePerSqrCmTOpascal<-function(gramForcePerSqrCm) gramForcePerSqrCm * 98.0665
# Input: gram-force per square centimeter (gf/cm2)
# Output: pascal (Pa)

NISTpascalTOgramForcePerSqrCm<-function(pascal) pascal / 98.0665
# Input: pascal (Pa)
# Output: gram-force per square centimeter (gf/cm2)

NISTgramPerCubCmTOkgPerCubMeter<-function(gramPerCubCm) gramPerCubCm * 1000
# Input: gram per cubic centimeter (g/cm3)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOgramPerCubCm<-function(kgPerCubMeter) kgPerCubMeter / 1000
# Input: kilogram per cubic meter (kg/m3)
# Output: gram per cubic centimeter (g/cm3)

NISThectareTOsqrMeter<-function(hectare) hectare * 10000
# Input: hectare (ha)
# Output: square meter (m2)

NISTsqrMeterTOhectare<-function(sqrMeter) sqrMeter / 10000
# Input: square meter (m2)
# Output: hectare (ha)

NISThorsepowerTOwatt<-function(horsepower) horsepower * 745.6999
# Input: horsepower (550 ft * lbf/s) (hp)
# Output: watt (W)

NISTwattTOhorsepower<-function(watt) watt / 745.6999
# Input: watt (W)
# Output: horsepower (550 ft * lbf/s) (hp)

NISThorsepowerBoilerTOwatt<-function(horsepowerBoiler) horsepowerBoiler * 9809.5
# Input: horsepower (boiler)
# Output: watt (W)

NISTwattTOhorsepowerBoiler<-function(watt) watt / 9809.5
# Input: watt (W)
# Output: horsepower (boiler)

NISThorsepowerElectricTOwatt<-function(horsepowerElectric) horsepowerElectric * 746
# Input: horsepower (electric)
# Output: watt (W)

NISTwattTOhorsepowerElectric<-function(watt) watt / 746
# Input: watt (W)
# Output: horsepower (electric)

NISThorsepowerSItOwatt<-function(horsepowerSI) horsepowerSI * 735.4988
# Input: horsepower (metric)
# Output: watt (W)

NISTwattTOhorsepowerSI<-function(watt) watt / 735.4988
# Input: watt (W)
# Output: horsepower (metric)

NISThorsepowerImperialTOwatt<-function(horsepowerImperial) horsepowerImperial * 745.7
# Input: horsepower (U.K.)
# Output: watt (W)

NISTwattTOhorsepowerImperial<-function(watt) watt / 745.7
# Input: watt (W)
# Output: horsepower (U.K.)

NISThorsepowerWaterTOwatt<-function(horsepowerWater) horsepowerWater * 746.043
# Input: horsepower (water)
# Output: watt (W)

NISTwattTOhorsepowerWater<-function(watt) watt / 746.043
# Input: watt (W)
# Output: horsepower (water)

NISThourTOsec<-function(hour) hour * 3600
# Input: hour (h)
# Output: second (s)

NISTsecTOhour<-function(sec) sec / 3600
# Input: second (s)
# Output: hour (h)

NISThourSiderealTOsec<-function(hourSidereal) hourSidereal * 3590.17
# Input: hour (sidereal)
# Output: second (s)

NISTsecTOhourSidereal<-function(sec) sec / 3590.17
# Input: second (s)
# Output: hour (sidereal)

NISThundredweightLongTOkg<-function(hundredweightLong) hundredweightLong * 50.80235
# Input: hundredweight (long, 112 lb)
# Output: kilogram (kg)

NISTkgTOhundredweightLong<-function(kg) kg / 50.80235
# Input: kilogram (kg)
# Output: hundredweight (long, 112 lb)

NISThundredweightShortTOkg<-function(hundredweightShort) hundredweightShort * 45.35924
# Input: hundredweight (short, 100 lb)
# Output: kilogram (kg)

NISTkgTOhundredweightShort<-function(kg) kg / 45.35924
# Input: kilogram (kg)
# Output: hundredweight (short, 100 lb)

NISTinchTOmeter<-function(inch) inch * 0.0254
# Input: inch (in)
# Output: meter (m)

NISTmeterTOinch<-function(meter) meter / 0.0254
# Input: meter (m)
# Output: inch (in)

NISTinchTOcm<-function(inch) inch * 2.54
# Input: inch (in)
# Output: centimeter (cm)

NISTcmTOinch<-function(cm) cm / 2.54
# Input: centimeter (cm)
# Output: inch (in)

NISTinchOfMercury32FtOpascal<-function(inchOfMercury32F) inchOfMercury32F * 3386.38
# Input: inch of mercury (32 F) 12
# Output: pascal (Pa)

NISTpascalTOinchOfMercury32F<-function(pascal) pascal / 3386.38
# Input: pascal (Pa)
# Output: inch of mercury (32 F) 12

NISTinchOfMercury32FtOkpascal<-function(inchOfMercury32F) inchOfMercury32F * 3.38638
# Input: inch of mercury (32 F) 12
# Output: kilopascal (kPa)

NISTkpascalTOinchOfMercury32F<-function(kpascal) kpascal / 3.38638
# Input: kilopascal (kPa)
# Output: inch of mercury (32 F) 12

NISTinchOfMercury60FtOpascal<-function(inchOfMercury60F) inchOfMercury60F * 3376.85
# Input: inch of mercury (60 F) 12
# Output: pascal (Pa)

NISTpascalTOinchOfMercury60F<-function(pascal) pascal / 3376.85
# Input: pascal (Pa)
# Output: inch of mercury (60 F) 12

NISTinchOfMercury60FtOkpascal<-function(inchOfMercury60F) inchOfMercury60F * 3.37685
# Input: inch of mercury (60 F) 12
# Output: kilopascal (kPa)

NISTkpascalTOinchOfMercury60F<-function(kpascal) kpascal / 3.37685
# Input: kilopascal (kPa)
# Output: inch of mercury (60 F) 12

NISTinchOfMercuryConvtnlTOpascal<-function(inchOfMercuryConvtnl) inchOfMercuryConvtnl * 3386.389
# Input: inch of mercury, conventional (inHg) 12
# Output: pascal (Pa)

NISTpascalTOinchOfMercuryConvtnl<-function(pascal) pascal / 3386.389
# Input: pascal (Pa)
# Output: inch of mercury, conventional (inHg) 12

NISTinchOfMercuryConvtnlTOkpascal<-function(inchOfMercuryConvtnl) inchOfMercuryConvtnl * 3.386389
# Input: inch of mercury, conventional (inHg) 12
# Output: kilopascal (kPa)

NISTkpascalTOinchOfMercuryConvtnl<-function(kpascal) kpascal / 3.386389
# Input: kilopascal (kPa)
# Output: inch of mercury, conventional (inHg) 12

NISTinchOfWater39.2FtOpascal<-function(inchOfWater39.2F) inchOfWater39.2F * 249.082
# Input: inch of water (39.2 F) 12
# Output: pascal (Pa)

NISTpascalTOinchOfWater39.2F<-function(pascal) pascal / 249.082
# Input: pascal (Pa)
# Output: inch of water (39.2 F) 12

NISTinchOfWater60FtOpascal<-function(inchOfWater60F) inchOfWater60F * 248.84
# Input: inch of water (60 F) 12
# Output: pascal (Pa)

NISTpascalTOinchOfWater60F<-function(pascal) pascal / 248.84
# Input: pascal (Pa)
# Output: inch of water (60 F) 12

NISTinchOfWaterConvtnlTOpascal<-function(inchOfWaterConvtnl) inchOfWaterConvtnl * 249.0889
# Input: inch of water, conventional (inH2O) 12
# Output: pascal (Pa)

NISTpascalTOinchOfWaterConvtnl<-function(pascal) pascal / 249.0889
# Input: pascal (Pa)
# Output: inch of water, conventional (inH2O) 12

NISTinchPerSecTOmeterPerSec<-function(inchPerSec) inchPerSec * 0.0254
# Input: inch per second (in/s)
# Output: meter per second (m/s)

NISTmeterPerSecTOinchPerSec<-function(meterPerSec) meterPerSec / 0.0254
# Input: meter per second (m/s)
# Output: inch per second (in/s)

NISTinchPerSecSqrdTOmeterPerSecSqrd<-function(inchPerSecSqrd) inchPerSecSqrd * 0.0254
# Input: inch per second squared (in/s2)
# Output: meter per second squared (m/s2)

NISTmeterPerSecSqrdTOinchPerSecSqrd<-function(meterPerSecSqrd) meterPerSecSqrd / 0.0254
# Input: meter per second squared (m/s2)
# Output: inch per second squared (in/s2)

NISTinch4thPowerTOmeter4thPower<-function(inch4thPower) inch4thPower * 4.162314e-07
# Input: inch to the fourth power (in4) 16
# Output: meter to the fourth power (m4)

NISTmeter4thPowerTOinch4thPower<-function(meter4thPower) meter4thPower / 4.162314e-07
# Input: meter to the fourth power (m4)
# Output: inch to the fourth power (in4) 16

NISTkayserTOreciprocalMeter<-function(kayser) kayser * 100
# Input: kayser (K)
# Output: reciprocal meter (m-1)

NISTreciprocalMeterTOkayser<-function(reciprocalMeter) reciprocalMeter / 100
# Input: reciprocal meter (m-1)
# Output: kayser (K)

NISTkTOdegC<-function(k) k - 273.15
# Input: kelvin (K)
# Output: degree Celsius (C)

NISTdegCtOk<-function(degC) degC + 273.15
# Input: degree Celsius (C)
# Output: kelvin (K)

NISTkilocalITtOjoule<-function(kilocalIT) kilocalIT * 4186.8
# Input: kilocalorieIT (kcalIT)
# Output: joule (J)

NISTjouleTOkilocalIT<-function(joule) joule / 4186.8
# Input: joule (J)
# Output: kilocalorieIT (kcalIT)

NISTkilocalthTOjoule<-function(kilocalth) kilocalth * 4184
# Input: kilocalorieth (kcalth)
# Output: joule (J)

NISTjouleTOkilocalth<-function(joule) joule / 4184
# Input: joule (J)
# Output: kilocalorieth (kcalth)

NISTkilocalTOjoule<-function(kilocal) kilocal * 4190.02
# Input: kilocalorie (mean) (kcal)
# Output: joule (J)

NISTjouleTOkilocal<-function(joule) joule / 4190.02
# Input: joule (J)
# Output: kilocalorie (mean) (kcal)

NISTkilocalthPerMinTOwatt<-function(kilocalthPerMin) kilocalthPerMin * 69.73333
# Input: kilocalorieth per minute (kcalth/min)
# Output: watt (W)

NISTwattTOkilocalthPerMin<-function(watt) watt / 69.73333
# Input: watt (W)
# Output: kilocalorieth per minute (kcalth/min)

NISTkilocalthPerSecTOwatt<-function(kilocalthPerSec) kilocalthPerSec * 4184
# Input: kilocalorieth per second (kcalth/s)
# Output: watt (W)

NISTwattTOkilocalthPerSec<-function(watt) watt / 4184
# Input: watt (W)
# Output: kilocalorieth per second (kcalth/s)

NISTkgForceTOnewton<-function(kgForce) kgForce * 9.80665
# Input: kilogram-force (kgf)
# Output: newton (N)

NISTnewtonTOkgForce<-function(newton) newton / 9.80665
# Input: newton (N)
# Output: kilogram-force (kgf)

NISTkgForceMeterTOnewtonMeter<-function(kgForceMeter) kgForceMeter * 9.80665
# Input: kilogram-force meter (kgf * m)
# Output: newton meter (N * m)

NISTnewtonMeterTOkgForceMeter<-function(newtonMeter) newtonMeter / 9.80665
# Input: newton meter (N * m)
# Output: kilogram-force meter (kgf * m)

NISTkgForcePerSqrCmTOpascal<-function(kgForcePerSqrCm) kgForcePerSqrCm * 98066.5
# Input: kilogram-force per square centimeter (kgf/cm2)
# Output: pascal (Pa)

NISTpascalTOkgForcePerSqrCm<-function(pascal) pascal / 98066.5
# Input: pascal (Pa)
# Output: kilogram-force per square centimeter (kgf/cm2)

NISTkgForcePerSqrCmTOkpascal<-function(kgForcePerSqrCm) kgForcePerSqrCm * 98.0665
# Input: kilogram-force per square centimeter (kgf/cm2)
# Output: kilopascal (kPa)

NISTkpascalTOkgForcePerSqrCm<-function(kpascal) kpascal / 98.0665
# Input: kilopascal (kPa)
# Output: kilogram-force per square centimeter (kgf/cm2)

NISTkgForcePerSqrMeterTOpascal<-function(kgForcePerSqrMeter) kgForcePerSqrMeter * 9.80665
# Input: kilogram-force per square meter (kgf/m2)
# Output: pascal (Pa)

NISTpascalTOkgForcePerSqrMeter<-function(pascal) pascal / 9.80665
# Input: pascal (Pa)
# Output: kilogram-force per square meter (kgf/m2)

NISTkgForcePerSqrMillimeterTOpascal<-function(kgForcePerSqrMillimeter) kgForcePerSqrMillimeter * 9806650
# Input: kilogram-force per square millimeter (kgf/mm2)
# Output: pascal (Pa)

NISTpascalTOkgForcePerSqrMillimeter<-function(pascal) pascal / 9806650
# Input: pascal (Pa)
# Output: kilogram-force per square millimeter (kgf/mm2)

NISTkgForcePerSqrMillimeterTOmpascal<-function(kgForcePerSqrMillimeter) kgForcePerSqrMillimeter * 9.80665
# Input: kilogram-force per square millimeter (kgf/mm2)
# Output: megapascal (MPa)

NISTmpascalTOkgForcePerSqrMillimeter<-function(mpascal) mpascal / 9.80665
# Input: megapascal (MPa)
# Output: kilogram-force per square millimeter (kgf/mm2)

NISTkgForceSecSqrdPerMeterTOkg<-function(kgForceSecSqrdPerMeter) kgForceSecSqrdPerMeter * 9.80665
# Input: kilogram-force second squared per meter (kgf * s2/m)
# Output: kilogram (kg)

NISTkgTOkgForceSecSqrdPerMeter<-function(kg) kg / 9.80665
# Input: kilogram (kg)
# Output: kilogram-force second squared per meter (kgf * s2/m)

NISTkmPerHourTOmeterPerSec<-function(kmPerHour) kmPerHour * 0.2777778
# Input: kilometer per hour (km/h)
# Output: meter per second (m/s)

NISTmeterPerSecTOkmPerHour<-function(meterPerSec) meterPerSec / 0.2777778
# Input: meter per second (m/s)
# Output: kilometer per hour (km/h)

NISTkilopondTOnewton<-function(kilopond) kilopond * 9.80665
# Input: kilopond (kilogram-force) (kp)
# Output: newton (N)

NISTnewtonTOkilopond<-function(newton) newton / 9.80665
# Input: newton (N)
# Output: kilopond (kilogram-force) (kp)

NISTkilowattHourTOjoule<-function(kilowattHour) kilowattHour * 3600000
# Input: kilowatt hour (kW * h)
# Output: joule (J)

NISTjouleTOkilowattHour<-function(joule) joule / 3600000
# Input: joule (J)
# Output: kilowatt hour (kW * h)

NISTkilowattHourTOmegajoule<-function(kilowattHour) kilowattHour * 3.6
# Input: kilowatt hour (kW * h)
# Output: megajoule (MJ)

NISTmegajouleTOkilowattHour<-function(megajoule) megajoule / 3.6
# Input: megajoule (MJ)
# Output: kilowatt hour (kW * h)

NISTkipTOnewton<-function(kip) kip * 4448.222
# Input: kip (1 kip= 1000 lbf)
# Output: newton (N)

NISTnewtonTOkip<-function(newton) newton / 4448.222
# Input: newton (N)
# Output: kip (1 kip= 1000 lbf)

NISTkipTOkilonewton<-function(kip) kip * 4.448222
# Input: kip (1 kip= 1000 lbf)
# Output: kilonewton (kN)

NISTkilonewtonTOkip<-function(kilonewton) kilonewton / 4.448222
# Input: kilonewton (kN)
# Output: kip (1 kip= 1000 lbf)

NISTkipPerSqrInchTOpascal<-function(kipPerSqrInch) kipPerSqrInch * 6894757
# Input: kip per square inch (ksi) (kip/in2)
# Output: pascal (Pa)

NISTpascalTOkipPerSqrInch<-function(pascal) pascal / 6894757
# Input: pascal (Pa)
# Output: kip per square inch (ksi) (kip/in2)

NISTkipPerSqrInchTOkpascal<-function(kipPerSqrInch) kipPerSqrInch * 6894.757
# Input: kip per square inch (ksi) (kip/in2)
# Output: kilopascal (kPa)

NISTkpascalTOkipPerSqrInch<-function(kpascal) kpascal / 6894.757
# Input: kilopascal (kPa)
# Output: kip per square inch (ksi) (kip/in2)

NISTknotTOmeterPerSec<-function(knot) knot * 0.5144444
# Input: knot (nautical mile per hour)
# Output: meter per second (m/s)

NISTmeterPerSecTOknot<-function(meterPerSec) meterPerSec / 0.5144444
# Input: meter per second (m/s)
# Output: knot (nautical mile per hour)

NISTlambertTOcandelaPerSqrMeter<-function(lambert) lambert * 3183.099
# Input: lambert 17
# Output: candela per square meter (cd/m2)

NISTcandelaPerSqrMeterTOlambert<-function(candelaPerSqrMeter) candelaPerSqrMeter / 3183.099
# Input: candela per square meter (cd/m2)
# Output: lambert 17

NISTlangleyTOjoulePerSqrMeter<-function(langley) langley * 41840
# Input: langley (calth/cm2)
# Output: joule per square meter (J/m2)

NISTjoulePerSqrMeterTOlangley<-function(joulePerSqrMeter) joulePerSqrMeter / 41840
# Input: joule per square meter (J/m2)
# Output: langley (calth/cm2)

NISTlightYearTOmeter<-function(lightYear) lightYear * 9.46073e+15
# Input: light year (l.y.) 18
# Output: meter (m)

NISTmeterTOlightYear<-function(meter) meter / 9.46073e+15
# Input: meter (m)
# Output: light year (l.y.) 18

NISTliterTOcubMeter<-function(liter) liter * 0.001
# Input: liter (L) 19
# Output: cubic meter (m3)

NISTcubMeterTOliter<-function(cubMeter) cubMeter / 0.001
# Input: cubic meter (m3)
# Output: liter (L) 19

NISTlumenPerSqrFtTOlux<-function(lumenPerSqrFt) lumenPerSqrFt * 10.76391
# Input: lumen per square foot (lm/ft2)
# Output: lux (lx)

NISTluxTOlumenPerSqrFt<-function(lux) lux / 10.76391
# Input: lux (lx)
# Output: lumen per square foot (lm/ft2)

NISTmaxwellTOweber<-function(maxwell) maxwell * 1e-08
# Input: maxwell (Mx)
# Output: weber (Wb)

NISTweberTOmaxwell<-function(weber) weber / 1e-08
# Input: weber (Wb)
# Output: maxwell (Mx)

NISTmhoTOsiemens<-function(mho) mho * 1
# Input: mho
# Output: siemens (S)

NISTsiemensTOmho<-function(siemens) siemens / 1
# Input: siemens (S)
# Output: mho

NISTmicroinchTOmeter<-function(microinch) microinch * 2.54e-08
# Input: microinch
# Output: meter (m)

NISTmeterTOmicroinch<-function(meter) meter / 2.54e-08
# Input: meter (m)
# Output: microinch

NISTmicroinchTOmicrometer<-function(microinch) microinch * 0.0254
# Input: microinch
# Output: micrometer (um)

NISTmicrometerTOmicroinch<-function(micrometer) micrometer / 0.0254
# Input: micrometer (um)
# Output: microinch

NISTmicronTOmeter<-function(micron) micron * 1e-06
# Input: micron (u)
# Output: meter (m)

NISTmeterTOmicron<-function(meter) meter / 1e-06
# Input: meter (m)
# Output: micron (u)

NISTmicronTOmicrometer<-function(micron) micron * 1
# Input: micron (u)
# Output: micrometer (um)

NISTmicrometerTOmicron<-function(micrometer) micrometer / 1
# Input: micrometer (um)
# Output: micron (u)

NISTmilTOmeter<-function(mil) mil * 2.54e-05
# Input: mil (0.001 in)
# Output: meter (m)

NISTmeterTOmil<-function(meter) meter / 2.54e-05
# Input: meter (m)
# Output: mil (0.001 in)

NISTmilTOmillimeter<-function(mil) mil * 0.0254
# Input: mil (0.001 in)
# Output: millimeter (mm)

NISTmillimeterTOmil<-function(millimeter) millimeter / 0.0254
# Input: millimeter (mm)
# Output: mil (0.001 in)

NISTmilTOradian<-function(mil) mil * 0.0009817477
# Input: mil (angle)
# Output: radian (rad)

NISTradianTOmil<-function(radian) radian / 0.0009817477
# Input: radian (rad)
# Output: mil (angle)

NISTmilTOdeg<-function(mil) mil * 0.05625
# Input: mil (angle)
# Output: degree ()

NISTdegTOmil<-function(deg) deg / 0.05625
# Input: degree ()
# Output: mil (angle)

NISTmileTOmeter<-function(mile) mile * 1609.344
# Input: mile (mi)
# Output: meter (m)

NISTmeterTOmile<-function(meter) meter / 1609.344
# Input: meter (m)
# Output: mile (mi)

NISTmileTOkm<-function(mile) mile * 1.609344
# Input: mile (mi)
# Output: kilometer (km)

NISTkmTOmile<-function(km) km / 1.609344
# Input: kilometer (km)
# Output: mile (mi)

NISTmileUSsurveyTOmeter<-function(mileUSsurvey) mileUSsurvey * 1609.347
# Input: mile (based on U.S. survey foot) (mi) 7
# Output: meter (m)

NISTmeterTOmileUSsurvey<-function(meter) meter / 1609.347
# Input: meter (m)
# Output: mile (based on U.S. survey foot) (mi) 7

NISTmileUSsurveyTOkm<-function(mileUSsurvey) mileUSsurvey * 1.609347
# Input: mile (based on U.S. survey foot) (mi) 7
# Output: kilometer (km)

NISTkmTOmileUSsurvey<-function(km) km / 1.609347
# Input: kilometer (km)
# Output: mile (based on U.S. survey foot) (mi) 7

NISTmileNauticalTOmeter<-function(mileNautical) mileNautical * 1852
# Input: mile, nautical 20
# Output: meter (m)

NISTmeterTOmileNautical<-function(meter) meter / 1852
# Input: meter (m)
# Output: mile, nautical 20

NISTmilePerGallonTOmeterPerCubMeter<-function(milePerGallon) milePerGallon * 425143.7
# Input: mile per gallon (U.S.) (mpg) (mi/gal)
# Output: meter per cubic meter (m/m3)

NISTmeterPerCubMeterTOmilePerGallon<-function(meterPerCubMeter) meterPerCubMeter / 425143.7
# Input: meter per cubic meter (m/m3)
# Output: mile per gallon (U.S.) (mpg) (mi/gal)

NISTmilePerGallonTOkmPerLiter<-function(milePerGallon) milePerGallon * 0.4251437
# Input: mile per gallon (U.S.) (mpg) (mi/gal)
# Output: kilometer per liter (km/L)

NISTkmPerLiterTOmilePerGallon<-function(kmPerLiter) kmPerLiter / 0.4251437
# Input: kilometer per liter (km/L)
# Output: mile per gallon (U.S.) (mpg) (mi/gal)

NISTmilePerGallonTOliterPer100Km<-function(milePerGallon) 235.215 / milePerGallon
# Input: mile per gallon (U.S.) (mpg) (mi/gal) 21
# Output: liter per 100 kilometer (L/100 km)

NISTliterPer100KmTOmilePerGallon<-function(literPer100Km) 235.215 / literPer100Km
# Input: liter per 100 kilometer (L/100 km)
# Output: mile per gallon (U.S.) (mpg) (mi/gal) 21

NISTmilePerHourTOmeterPerSec<-function(milePerHour) milePerHour * 0.44704
# Input: mile per hour (mi/h)
# Output: meter per second (m/s)

NISTmeterPerSecTOmilePerHour<-function(meterPerSec) meterPerSec / 0.44704
# Input: meter per second (m/s)
# Output: mile per hour (mi/h)

NISTmilePerHourTOkmPerHour<-function(milePerHour) milePerHour * 1.609344
# Input: mile per hour (mi/h)
# Output: kilometer per hour (km/h)

NISTkmPerHourTOmilePerHour<-function(kmPerHour) kmPerHour / 1.609344
# Input: kilometer per hour (km/h)
# Output: mile per hour (mi/h)

NISTmilePerMinTOmeterPerSec<-function(milePerMin) milePerMin * 26.8224
# Input: mile per minute (mi/min)
# Output: meter per second (m/s)

NISTmeterPerSecTOmilePerMin<-function(meterPerSec) meterPerSec / 26.8224
# Input: meter per second (m/s)
# Output: mile per minute (mi/min)

NISTmilePerSecTOmeterPerSec<-function(milePerSec) milePerSec * 1609.344
# Input: mile per second (mi/s)
# Output: meter per second (m/s)

NISTmeterPerSecTOmilePerSec<-function(meterPerSec) meterPerSec / 1609.344
# Input: meter per second (m/s)
# Output: mile per second (mi/s)

NISTmillibarTOpascal<-function(millibar) millibar * 100
# Input: millibar (mbar)
# Output: pascal (Pa)

NISTpascalTOmillibar<-function(pascal) pascal / 100
# Input: pascal (Pa)
# Output: millibar (mbar)

NISTmillibarTOkpascal<-function(millibar) millibar * 0.1
# Input: millibar (mbar)
# Output: kilopascal (kPa)

NISTkpascalTOmillibar<-function(kpascal) kpascal / 0.1
# Input: kilopascal (kPa)
# Output: millibar (mbar)

NISTmillimeterOfMercuryConvtnlTOpascal<-function(millimeterOfMercuryConvtnl) millimeterOfMercuryConvtnl * 133.3224
# Input: millimeter of mercury, conventional (mmHg) 12
# Output: pascal (Pa)

NISTpascalTOmillimeterOfMercuryConvtnl<-function(pascal) pascal / 133.3224
# Input: pascal (Pa)
# Output: millimeter of mercury, conventional (mmHg) 12

NISTmillimeterOfWaterConvtnlTOpascal<-function(millimeterOfWaterConvtnl) millimeterOfWaterConvtnl * 9.80665
# Input: millimeter of water, conventional (mm H2O) 12
# Output: pascal (Pa)

NISTpascalTOmillimeterOfWaterConvtnl<-function(pascal) pascal / 9.80665
# Input: pascal (Pa)
# Output: millimeter of water, conventional (mm H2O) 12

NISTminTOradian<-function(min) min * 0.0002908882
# Input: minute (angle) (')
# Output: radian (rad)

NISTradianTOmin<-function(radian) radian / 0.0002908882
# Input: radian (rad)
# Output: minute (angle) (')

NISTminTOsec<-function(min) min * 60
# Input: minute (min)
# Output: second (s)

NISTsecTOmin<-function(sec) sec / 60
# Input: second (s)
# Output: minute (min)

NISTminSiderealTOsec<-function(minSidereal) minSidereal * 59.83617
# Input: minute (sidereal)
# Output: second (s)

NISTsecTOminSidereal<-function(sec) sec / 59.83617
# Input: second (s)
# Output: minute (sidereal)

NISToerstedTOamperePerMeter<-function(oersted) oersted * 79.57747
# Input: oersted (Oe)
# Output: ampere per meter (A/m)

NISTamperePerMeterTOoersted<-function(amperePerMeter) amperePerMeter / 79.57747
# Input: ampere per meter (A/m)
# Output: oersted (Oe)

NISTohmCmTOohmMeter<-function(ohmCm) ohmCm * 0.01
# Input: ohm centimeter (ohm * cm)
# Output: ohm meter (ohm * m)

NISTohmMeterTOohmCm<-function(ohmMeter) ohmMeter / 0.01
# Input: ohm meter (ohm * m)
# Output: ohm centimeter (ohm * cm)

NISTohmCircularMilPerFtTOohmMeter<-function(ohmCircularMilPerFt) ohmCircularMilPerFt * 1.662426e-09
# Input: ohm circular-mil per foot
# Output: ohm meter (ohm * m)

NISTohmMeterTOohmCircularMilPerFt<-function(ohmMeter) ohmMeter / 1.662426e-09
# Input: ohm meter (ohm * m)
# Output: ohm circular-mil per foot

NISTohmCircularMilPerFtTOohmSqrMillimeterPerMeter<-function(ohmCircularMilPerFt) ohmCircularMilPerFt * 0.001662426
# Input: ohm circular-mil per foot
# Output: ohm square millimeter per meter (ohm * mm2/m)

NISTohmSqrMillimeterPerMeterTOohmCircularMilPerFt<-function(ohmSqrMillimeterPerMeter) ohmSqrMillimeterPerMeter / 0.001662426
# Input: ohm square millimeter per meter (ohm * mm2/m)
# Output: ohm circular-mil per foot

NISTounceAvoirdupoisTOkg<-function(ounceAvoirdupois) ounceAvoirdupois * 0.02834952
# Input: ounce (avoirdupois) (oz)
# Output: kilogram (kg)

NISTkgTOounceAvoirdupois<-function(kg) kg / 0.02834952
# Input: kilogram (kg)
# Output: ounce (avoirdupois) (oz)

NISTounceAvoirdupoisTOgram<-function(ounceAvoirdupois) ounceAvoirdupois * 28.34952
# Input: ounce (avoirdupois) (oz)
# Output: gram (g)

NISTgramTOounceAvoirdupois<-function(gram) gram / 28.34952
# Input: gram (g)
# Output: ounce (avoirdupois) (oz)

NISTounceTroyTOkg<-function(ounceTroy) ounceTroy * 0.03110348
# Input: ounce (troy or apothecary) (oz)
# Output: kilogram (kg)

NISTkgTOounceTroy<-function(kg) kg / 0.03110348
# Input: kilogram (kg)
# Output: ounce (troy or apothecary) (oz)

NISTounceTroyTOgram<-function(ounceTroy) ounceTroy * 31.10348
# Input: ounce (troy or apothecary) (oz)
# Output: gram (g)

NISTgramTOounceTroy<-function(gram) gram / 31.10348
# Input: gram (g)
# Output: ounce (troy or apothecary) (oz)

NISTounceImperialTOcubMeter<-function(ounceImperial) ounceImperial * 2.841306e-05
# Input: ounce [Canadian and U.K. fluid (Imperial)] (fl oz)
# Output: cubic meter (m3)

NISTcubMeterTOounceImperial<-function(cubMeter) cubMeter / 2.841306e-05
# Input: cubic meter (m3)
# Output: ounce [Canadian and U.K. fluid (Imperial)] (fl oz)

NISTounceImperialTOml<-function(ounceImperial) ounceImperial * 28.41306
# Input: ounce [Canadian and U.K. fluid (Imperial)] (fl oz)
# Output: milliliter (mL)

NISTmlTOounceImperial<-function(ml) ml / 28.41306
# Input: milliliter (mL)
# Output: ounce [Canadian and U.K. fluid (Imperial)] (fl oz)

NISTounceTOcubMeter<-function(ounce) ounce * 2.957353e-05
# Input: ounce (U.S. fluid) (fl oz)
# Output: cubic meter (m3)

NISTcubMeterTOounce<-function(cubMeter) cubMeter / 2.957353e-05
# Input: cubic meter (m3)
# Output: ounce (U.S. fluid) (fl oz)

NISTounceTOml<-function(ounce) ounce * 29.57353
# Input: ounce (U.S. fluid) (fl oz)
# Output: milliliter (mL)

NISTmlTOounce<-function(ml) ml / 29.57353
# Input: milliliter (mL)
# Output: ounce (U.S. fluid) (fl oz)

NISTounceForceTOnewton<-function(ounceForce) ounceForce * 0.2780139
# Input: ounce (avoirdupois)-force (ozf)
# Output: newton (N)

NISTnewtonTOounceForce<-function(newton) newton / 0.2780139
# Input: newton (N)
# Output: ounce (avoirdupois)-force (ozf)

NISTounceForceInchTOnewtonMeter<-function(ounceForceInch) ounceForceInch * 0.007061552
# Input: ounce (avoirdupois)-force inch (ozf * in)
# Output: newton meter (N * m)

NISTnewtonMeterTOounceForceInch<-function(newtonMeter) newtonMeter / 0.007061552
# Input: newton meter (N * m)
# Output: ounce (avoirdupois)-force inch (ozf * in)

NISTounceForceInchTOmillinewtonMeter<-function(ounceForceInch) ounceForceInch * 7.061552
# Input: ounce (avoirdupois)-force inch (ozf * in)
# Output: millinewton meter (mN * m)

NISTmillinewtonMeterTOounceForceInch<-function(millinewtonMeter) millinewtonMeter / 7.061552
# Input: millinewton meter (mN * m)
# Output: ounce (avoirdupois)-force inch (ozf * in)

NISTouncePerCubInchTOkgPerCubMeter<-function(ouncePerCubInch) ouncePerCubInch * 1729.994
# Input: ounce (avoirdupois) per cubic inch (oz/in3)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOouncePerCubInch<-function(kgPerCubMeter) kgPerCubMeter / 1729.994
# Input: kilogram per cubic meter (kg/m3)
# Output: ounce (avoirdupois) per cubic inch (oz/in3)

NISTounceImperialPerGallonTOkgPerCubMeter<-function(ounceImperialPerGallon) ounceImperialPerGallon * 6.236023
# Input: ounce (avoirdupois) per gallon [Canadian and U.K. (Imperial)] (oz/gal)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOounceImperialPerGallon<-function(kgPerCubMeter) kgPerCubMeter / 6.236023
# Input: kilogram per cubic meter (kg/m3)
# Output: ounce (avoirdupois) per gallon [Canadian and U.K. (Imperial)] (oz/gal)

NISTounceImperialPerGallonTOgramPerLiter<-function(ounceImperialPerGallon) ounceImperialPerGallon * 6.236023
# Input: ounce (avoirdupois) per gallon [Canadian and U.K. (Imperial)] (oz/gal)
# Output: gram per liter (g/L)

NISTgramPerLiterTOounceImperialPerGallon<-function(gramPerLiter) gramPerLiter / 6.236023
# Input: gram per liter (g/L)
# Output: ounce (avoirdupois) per gallon [Canadian and U.K. (Imperial)] (oz/gal)

NISTouncePerGallonTOkgPerCubMeter<-function(ouncePerGallon) ouncePerGallon * 7.489152
# Input: ounce (avoirdupois) per gallon (U.S.) (oz/gal)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOouncePerGallon<-function(kgPerCubMeter) kgPerCubMeter / 7.489152
# Input: kilogram per cubic meter (kg/m3)
# Output: ounce (avoirdupois) per gallon (U.S.) (oz/gal)

NISTouncePerGallonTOgramPerLiter<-function(ouncePerGallon) ouncePerGallon * 7.489152
# Input: ounce (avoirdupois) per gallon (U.S.) (oz/gal)
# Output: gram per liter (g/L)

NISTgramPerLiterTOouncePerGallon<-function(gramPerLiter) gramPerLiter / 7.489152
# Input: gram per liter (g/L)
# Output: ounce (avoirdupois) per gallon (U.S.) (oz/gal)

NISTouncePerSqrFtTOkgPerSqrMeter<-function(ouncePerSqrFt) ouncePerSqrFt * 0.3051517
# Input: ounce (avoirdupois) per square foot (oz/ft2)
# Output: kilogram per square meter (kg/m2)

NISTkgPerSqrMeterTOouncePerSqrFt<-function(kgPerSqrMeter) kgPerSqrMeter / 0.3051517
# Input: kilogram per square meter (kg/m2)
# Output: ounce (avoirdupois) per square foot (oz/ft2)

NISTouncePerSqrInchTOkgPerSqrMeter<-function(ouncePerSqrInch) ouncePerSqrInch * 43.94185
# Input: ounce (avoirdupois) per square inch (oz/in2)
# Output: kilogram per square meter (kg/m2)

NISTkgPerSqrMeterTOouncePerSqrInch<-function(kgPerSqrMeter) kgPerSqrMeter / 43.94185
# Input: kilogram per square meter (kg/m2)
# Output: ounce (avoirdupois) per square inch (oz/in2)

NISTouncePerSqrYardTOkgPerSqrMeter<-function(ouncePerSqrYard) ouncePerSqrYard * 0.03390575
# Input: ounce (avoirdupois) per square yard (oz/yd2)
# Output: kilogram per square meter (kg/m2)

NISTkgPerSqrMeterTOouncePerSqrYard<-function(kgPerSqrMeter) kgPerSqrMeter / 0.03390575
# Input: kilogram per square meter (kg/m2)
# Output: ounce (avoirdupois) per square yard (oz/yd2)

NISTparsecTOmeter<-function(parsec) parsec * 3.085678e+16
# Input: parsec (pc)
# Output: meter (m)

NISTmeterTOparsec<-function(meter) meter / 3.085678e+16
# Input: meter (m)
# Output: parsec (pc)

NISTpeckTOcubMeter<-function(peck) peck * 0.008809768
# Input: peck (U.S.) (pk)
# Output: cubic meter (m3)

NISTcubMeterTOpeck<-function(cubMeter) cubMeter / 0.008809768
# Input: cubic meter (m3)
# Output: peck (U.S.) (pk)

NISTpeckTOliter<-function(peck) peck * 8.809768
# Input: peck (U.S.) (pk)
# Output: liter (L)

NISTliterTOpeck<-function(liter) liter / 8.809768
# Input: liter (L)
# Output: peck (U.S.) (pk)

NISTpennyweightTOkg<-function(pennyweight) pennyweight * 0.001555174
# Input: pennyweight (dwt)
# Output: kilogram (kg)

NISTkgTOpennyweight<-function(kg) kg / 0.001555174
# Input: kilogram (kg)
# Output: pennyweight (dwt)

NISTpennyweightTOgram<-function(pennyweight) pennyweight * 1.555174
# Input: pennyweight (dwt)
# Output: gram (g)

NISTgramTOpennyweight<-function(gram) gram / 1.555174
# Input: gram (g)
# Output: pennyweight (dwt)

NISTperm0CtOkgPerPascalSecSqrMeter<-function(perm0C) perm0C * 5.72135e-11
# Input: perm (0 C)
# Output: kilogram per pascal second square meter [kg/(Pa * s * m2)]

NISTkgPerPascalSecSqrMeterTOperm0C<-function(kgPerPascalSecSqrMeter) kgPerPascalSecSqrMeter / 5.72135e-11
# Input: kilogram per pascal second square meter [kg/(Pa * s * m2)]
# Output: perm (0 C)

NISTperm23CtOkgPerPascalSecSqrMeter<-function(perm23C) perm23C * 5.74525e-11
# Input: perm (23 C)
# Output: kilogram per pascal second square meter [kg/(Pa * s * m2)]

NISTkgPerPascalSecSqrMeterTOperm23C<-function(kgPerPascalSecSqrMeter) kgPerPascalSecSqrMeter / 5.74525e-11
# Input: kilogram per pascal second square meter [kg/(Pa * s * m2)]
# Output: perm (23 C)

NISTpermInch0CtOkgPerPascalSecMeter<-function(permInch0C) permInch0C * 1.45322e-12
# Input: perm inch (0 C)
# Output: kilogram per pascal second meter [kg/(Pa*s*m)]

NISTkgPerPascalSecMeterTOpermInch0C<-function(kgPerPascalSecMeter) kgPerPascalSecMeter / 1.45322e-12
# Input: kilogram per pascal second meter [kg/(Pa*s*m)]
# Output: perm inch (0 C)

NISTpermInch23CtOkgPerPascalSecMeter<-function(permInch23C) permInch23C * 1.45929e-12
# Input: perm inch (23 C)
# Output: kilogram per pascal second meter [kg/(Pa*s*m)]

NISTkgPerPascalSecMeterTOpermInch23C<-function(kgPerPascalSecMeter) kgPerPascalSecMeter / 1.45929e-12
# Input: kilogram per pascal second meter [kg/(Pa*s*m)]
# Output: perm inch (23 C)

NISTphotTOlux<-function(phot) phot * 10000
# Input: phot (ph)
# Output: lux (lx)

NISTluxTOphot<-function(lux) lux / 10000
# Input: lux (lx)
# Output: phot (ph)

NISTpicaComputerTOmeter<-function(picaComputer) picaComputer * 0.004233333
# Input: pica (computer) (1/6 in)
# Output: meter (m)

NISTmeterTOpicaComputer<-function(meter) meter / 0.004233333
# Input: meter (m)
# Output: pica (computer) (1/6 in)

NISTpicaComputerTOmillimeter<-function(picaComputer) picaComputer * 4.233333
# Input: pica (computer) (1/6 in)
# Output: millimeter (mm)

NISTmillimeterTOpicaComputer<-function(millimeter) millimeter / 4.233333
# Input: millimeter (mm)
# Output: pica (computer) (1/6 in)

NISTpicaPrinterTOmeter<-function(picaPrinter) picaPrinter * 0.004217518
# Input: pica (printer's)
# Output: meter (m)

NISTmeterTOpicaPrinter<-function(meter) meter / 0.004217518
# Input: meter (m)
# Output: pica (printer's)

NISTpicaPrinterTOmillimeter<-function(picaPrinter) picaPrinter * 4.217518
# Input: pica (printer's)
# Output: millimeter (mm)

NISTmillimeterTOpicaPrinter<-function(millimeter) millimeter / 4.217518
# Input: millimeter (mm)
# Output: pica (printer's)

NISTpintUSdryTOcubMeter<-function(pintUSdry) pintUSdry * 0.0005506105
# Input: pint (U.S. dry) (dry pt)
# Output: cubic meter (m3)

NISTcubMeterTOpintUSdry<-function(cubMeter) cubMeter / 0.0005506105
# Input: cubic meter (m3)
# Output: pint (U.S. dry) (dry pt)

NISTpintUSdryTOliter<-function(pintUSdry) pintUSdry * 0.5506105
# Input: pint (U.S. dry) (dry pt)
# Output: liter (L)

NISTliterTOpintUSdry<-function(liter) liter / 0.5506105
# Input: liter (L)
# Output: pint (U.S. dry) (dry pt)

NISTpintUSliquidTOcubMeter<-function(pintUSliquid) pintUSliquid * 0.0004731765
# Input: pint (U.S. liquid) (liq pt)
# Output: cubic meter (m3)

NISTcubMeterTOpintUSliquid<-function(cubMeter) cubMeter / 0.0004731765
# Input: cubic meter (m3)
# Output: pint (U.S. liquid) (liq pt)

NISTpintUSliquidTOliter<-function(pintUSliquid) pintUSliquid * 0.4731765
# Input: pint (U.S. liquid) (liq pt)
# Output: liter (L)

NISTliterTOpintUSliquid<-function(liter) liter / 0.4731765
# Input: liter (L)
# Output: pint (U.S. liquid) (liq pt)

NISTpointComputerTOmeter<-function(pointComputer) pointComputer * 0.0003527778
# Input: point (computer) (1/72 in)
# Output: meter (m)

NISTmeterTOpointComputer<-function(meter) meter / 0.0003527778
# Input: meter (m)
# Output: point (computer) (1/72 in)

NISTpointComputerTOmillimeter<-function(pointComputer) pointComputer * 0.3527778
# Input: point (computer) (1/72 in)
# Output: millimeter (mm)

NISTmillimeterTOpointComputer<-function(millimeter) millimeter / 0.3527778
# Input: millimeter (mm)
# Output: point (computer) (1/72 in)

NISTpointPrinterTOmeter<-function(pointPrinter) pointPrinter * 0.0003514598
# Input: point (printer's)
# Output: meter (m)

NISTmeterTOpointPrinter<-function(meter) meter / 0.0003514598
# Input: meter (m)
# Output: point (printer's)

NISTpointPrinterTOmillimeter<-function(pointPrinter) pointPrinter * 0.3514598
# Input: point (printer's)
# Output: millimeter (mm)

NISTmillimeterTOpointPrinter<-function(millimeter) millimeter / 0.3514598
# Input: millimeter (mm)
# Output: point (printer's)

NISTpoiseTOpascalSec<-function(poise) poise * 0.1
# Input: poise (P)
# Output: pascal second (Pa * s)

NISTpascalSecTOpoise<-function(pascalSec) pascalSec / 0.1
# Input: pascal second (Pa * s)
# Output: poise (P)

NISTpoundAvoirdupoisTOkg<-function(poundAvoirdupois) poundAvoirdupois * 0.4535924
# Input: pound (avoirdupois) (lb) 22
# Output: kilogram (kg)

NISTkgTOpoundAvoirdupois<-function(kg) kg / 0.4535924
# Input: kilogram (kg)
# Output: pound (avoirdupois) (lb) 22

NISTpoundTroyTOkg<-function(poundTroy) poundTroy * 0.3732417
# Input: pound (troy or apothecary) (lb)
# Output: kilogram (kg)

NISTkgTOpoundTroy<-function(kg) kg / 0.3732417
# Input: kilogram (kg)
# Output: pound (troy or apothecary) (lb)

NISTpoundalTOnewton<-function(poundal) poundal * 0.138255
# Input: poundal
# Output: newton (N)

NISTnewtonTOpoundal<-function(newton) newton / 0.138255
# Input: newton (N)
# Output: poundal

NISTpoundalPerSqrFtTOpascal<-function(poundalPerSqrFt) poundalPerSqrFt * 1.488164
# Input: poundal per square foot
# Output: pascal (Pa)

NISTpascalTOpoundalPerSqrFt<-function(pascal) pascal / 1.488164
# Input: pascal (Pa)
# Output: poundal per square foot

NISTpoundalSecPerSqrFtTOpascalSec<-function(poundalSecPerSqrFt) poundalSecPerSqrFt * 1.488164
# Input: poundal second per square foot
# Output: pascal second (Pa * s)

NISTpascalSecTOpoundalSecPerSqrFt<-function(pascalSec) pascalSec / 1.488164
# Input: pascal second (Pa * s)
# Output: poundal second per square foot

NISTpoundFtSqrdTOkgMeterSqrd<-function(poundFtSqrd) poundFtSqrd * 0.04214011
# Input: pound foot squared (lb * ft2)
# Output: kilogram meter squared (kg * m2)

NISTkgMeterSqrdTOpoundFtSqrd<-function(kgMeterSqrd) kgMeterSqrd / 0.04214011
# Input: kilogram meter squared (kg * m2)
# Output: pound foot squared (lb * ft2)

NISTpoundForceTOnewton<-function(poundForce) poundForce * 4.448222
# Input: pound-force (lbf) 23
# Output: newton (N)

NISTnewtonTOpoundForce<-function(newton) newton / 4.448222
# Input: newton (N)
# Output: pound-force (lbf) 23

NISTpoundForceFtTOnewtonMeter<-function(poundForceFt) poundForceFt * 1.355818
# Input: pound-force foot (lbf * ft)
# Output: newton meter (N * m)

NISTnewtonMeterTOpoundForceFt<-function(newtonMeter) newtonMeter / 1.355818
# Input: newton meter (N * m)
# Output: pound-force foot (lbf * ft)

NISTpoundForceFtPerInchTOnewtonMeterPerMeter<-function(poundForceFtPerInch) poundForceFtPerInch * 53.37866
# Input: pound-force foot per inch (lbf * ft/in)
# Output: newton meter per meter (N * m/m)

NISTnewtonMeterPerMeterTOpoundForceFtPerInch<-function(newtonMeterPerMeter) newtonMeterPerMeter / 53.37866
# Input: newton meter per meter (N * m/m)
# Output: pound-force foot per inch (lbf * ft/in)

NISTpoundForceInchTOnewtonMeter<-function(poundForceInch) poundForceInch * 0.1129848
# Input: pound-force inch (lbf * in)
# Output: newton meter (N * m)

NISTnewtonMeterTOpoundForceInch<-function(newtonMeter) newtonMeter / 0.1129848
# Input: newton meter (N * m)
# Output: pound-force inch (lbf * in)

NISTpoundForceInchPerInchTOnewtonMeterPerMeter<-function(poundForceInchPerInch) poundForceInchPerInch * 4.448222
# Input: pound-force inch per inch (lbf * in/in)
# Output: newton meter per meter (N * m/m)

NISTnewtonMeterPerMeterTOpoundForceInchPerInch<-function(newtonMeterPerMeter) newtonMeterPerMeter / 4.448222
# Input: newton meter per meter (N * m/m)
# Output: pound-force inch per inch (lbf * in/in)

NISTpoundForcePerFtTOnewtonPerMeter<-function(poundForcePerFt) poundForcePerFt * 14.5939
# Input: pound-force per foot (lbf/ft)
# Output: newton per meter (N/m)

NISTnewtonPerMeterTOpoundForcePerFt<-function(newtonPerMeter) newtonPerMeter / 14.5939
# Input: newton per meter (N/m)
# Output: pound-force per foot (lbf/ft)

NISTpoundForcePerInchTOnewtonPerMeter<-function(poundForcePerInch) poundForcePerInch * 175.1268
# Input: pound-force per inch (lbf/in)
# Output: newton per meter (N/m)

NISTnewtonPerMeterTOpoundForcePerInch<-function(newtonPerMeter) newtonPerMeter / 175.1268
# Input: newton per meter (N/m)
# Output: pound-force per inch (lbf/in)

NISTpoundForcePerPoundTOnewtonPerKg<-function(poundForcePerPound) poundForcePerPound * 9.80665
# Input: pound-force per pound (lbf/lb) (thrust to mass ratio)
# Output: newton per kilogram (N/kg)

NISTnewtonPerKgTOpoundForcePerPound<-function(newtonPerKg) newtonPerKg / 9.80665
# Input: newton per kilogram (N/kg)
# Output: pound-force per pound (lbf/lb) (thrust to mass ratio)

NISTpoundForcePerSqrFtTOpascal<-function(poundForcePerSqrFt) poundForcePerSqrFt * 47.88026
# Input: pound-force per square foot (lbf/ft2)
# Output: pascal (Pa)

NISTpascalTOpoundForcePerSqrFt<-function(pascal) pascal / 47.88026
# Input: pascal (Pa)
# Output: pound-force per square foot (lbf/ft2)

NISTpoundForcePerSqrInchTOpascal<-function(poundForcePerSqrInch) poundForcePerSqrInch * 6894.757
# Input: pound-force per square inch (psi) (lbf/in2)
# Output: pascal (Pa)

NISTpascalTOpoundForcePerSqrInch<-function(pascal) pascal / 6894.757
# Input: pascal (Pa)
# Output: pound-force per square inch (psi) (lbf/in2)

NISTpoundForcePerSqrInchTOkpascal<-function(poundForcePerSqrInch) poundForcePerSqrInch * 6.894757
# Input: pound-force per square inch (psi) (lbf/in2)
# Output: kilopascal (kPa)

NISTkpascalTOpoundForcePerSqrInch<-function(kpascal) kpascal / 6.894757
# Input: kilopascal (kPa)
# Output: pound-force per square inch (psi) (lbf/in2)

NISTpoundForceSecPerSqrFtTOpascalSec<-function(poundForceSecPerSqrFt) poundForceSecPerSqrFt * 47.88026
# Input: pound-force second per square foot (lbf * s/ft2)
# Output: pascal second (Pa * s)

NISTpascalSecTOpoundForceSecPerSqrFt<-function(pascalSec) pascalSec / 47.88026
# Input: pascal second (Pa * s)
# Output: pound-force second per square foot (lbf * s/ft2)

NISTpoundForceSecPerSqrInchTOpascalSec<-function(poundForceSecPerSqrInch) poundForceSecPerSqrInch * 6894.757
# Input: pound-force second per square inch (lbf * s/in2)
# Output: pascal second (Pa * s)

NISTpascalSecTOpoundForceSecPerSqrInch<-function(pascalSec) pascalSec / 6894.757
# Input: pascal second (Pa * s)
# Output: pound-force second per square inch (lbf * s/in2)

NISTpoundInchSqrdTOkgMeterSqrd<-function(poundInchSqrd) poundInchSqrd * 0.0002926397
# Input: pound inch squared (lb * in2)
# Output: kilogram meter squared (kg *m2)

NISTkgMeterSqrdTOpoundInchSqrd<-function(kgMeterSqrd) kgMeterSqrd / 0.0002926397
# Input: kilogram meter squared (kg *m2)
# Output: pound inch squared (lb * in2)

NISTpoundPerCubFtTOkgPerCubMeter<-function(poundPerCubFt) poundPerCubFt * 16.01846
# Input: pound per cubic foot (lb/ft3)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOpoundPerCubFt<-function(kgPerCubMeter) kgPerCubMeter / 16.01846
# Input: kilogram per cubic meter (kg/m3)
# Output: pound per cubic foot (lb/ft3)

NISTpoundPerCubInchTOkgPerCubMeter<-function(poundPerCubInch) poundPerCubInch * 27679.9
# Input: pound per cubic inch (lb/in3)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOpoundPerCubInch<-function(kgPerCubMeter) kgPerCubMeter / 27679.9
# Input: kilogram per cubic meter (kg/m3)
# Output: pound per cubic inch (lb/in3)

NISTpoundPerCubYardTOkgPerCubMeter<-function(poundPerCubYard) poundPerCubYard * 0.5932764
# Input: pound per cubic yard (lb/yd3)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOpoundPerCubYard<-function(kgPerCubMeter) kgPerCubMeter / 0.5932764
# Input: kilogram per cubic meter (kg/m3)
# Output: pound per cubic yard (lb/yd3)

NISTpoundPerFtTOkgPerMeter<-function(poundPerFt) poundPerFt * 1.488164
# Input: pound per foot (lb/ft)
# Output: kilogram per meter (kg/m)

NISTkgPerMeterTOpoundPerFt<-function(kgPerMeter) kgPerMeter / 1.488164
# Input: kilogram per meter (kg/m)
# Output: pound per foot (lb/ft)

NISTpoundPerFtHourTOpascalSec<-function(poundPerFtHour) poundPerFtHour * 0.0004133789
# Input: pound per foot hour [lb/(ft * h)]
# Output: pascal second (Pa * s)

NISTpascalSecTOpoundPerFtHour<-function(pascalSec) pascalSec / 0.0004133789
# Input: pascal second (Pa * s)
# Output: pound per foot hour [lb/(ft * h)]

NISTpoundPerFtSecTOpascalSec<-function(poundPerFtSec) poundPerFtSec * 1.488164
# Input: pound per foot second [lb/(ft * s)]
# Output: pascal second (Pa * s)

NISTpascalSecTOpoundPerFtSec<-function(pascalSec) pascalSec / 1.488164
# Input: pascal second (Pa * s)
# Output: pound per foot second [lb/(ft * s)]

NISTpoundPerGallonImperialTOkgPerCubMeter<-function(poundPerGallonImperial) poundPerGallonImperial * 99.77637
# Input: pound per gallon [Canadian and U.K. (Imperial)] (lb/gal)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOpoundPerGallonImperial<-function(kgPerCubMeter) kgPerCubMeter / 99.77637
# Input: kilogram per cubic meter (kg/m3)
# Output: pound per gallon [Canadian and U.K. (Imperial)] (lb/gal)

NISTpoundPerGallonImperialTOkgPerLiter<-function(poundPerGallonImperial) poundPerGallonImperial * 0.09977637
# Input: pound per gallon [Canadian and U.K. (Imperial)] (lb/gal)
# Output: kilogram per liter (kg/L)

NISTkgPerLiterTOpoundPerGallonImperial<-function(kgPerLiter) kgPerLiter / 0.09977637
# Input: kilogram per liter (kg/L)
# Output: pound per gallon [Canadian and U.K. (Imperial)] (lb/gal)

NISTpoundPerGallonTOkgPerCubMeter<-function(poundPerGallon) poundPerGallon * 119.8264
# Input: pound per gallon (U.S.) (lb/gal)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOpoundPerGallon<-function(kgPerCubMeter) kgPerCubMeter / 119.8264
# Input: kilogram per cubic meter (kg/m3)
# Output: pound per gallon (U.S.) (lb/gal)

NISTpoundPerGallonTOkgPerLiter<-function(poundPerGallon) poundPerGallon * 0.1198264
# Input: pound per gallon (U.S.) (lb/gal)
# Output: kilogram per liter (kg/L)

NISTkgPerLiterTOpoundPerGallon<-function(kgPerLiter) kgPerLiter / 0.1198264
# Input: kilogram per liter (kg/L)
# Output: pound per gallon (U.S.) (lb/gal)

NISTpoundPerHorsepowerHourTOkgPerJoule<-function(poundPerHorsepowerHour) poundPerHorsepowerHour * 1.689659e-07
# Input: pound per horsepower hour [lb/(hp * h)]
# Output: kilogram per joule (kg/J)

NISTkgPerJouleTOpoundPerHorsepowerHour<-function(kgPerJoule) kgPerJoule / 1.689659e-07
# Input: kilogram per joule (kg/J)
# Output: pound per horsepower hour [lb/(hp * h)]

NISTpoundPerHourTOkgPerSec<-function(poundPerHour) poundPerHour * 0.0001259979
# Input: pound per hour (lb/h)
# Output: kilogram per second (kg/s)

NISTkgPerSecTOpoundPerHour<-function(kgPerSec) kgPerSec / 0.0001259979
# Input: kilogram per second (kg/s)
# Output: pound per hour (lb/h)

NISTpoundPerInchTOkgPerMeter<-function(poundPerInch) poundPerInch * 17.85797
# Input: pound per inch (lb/in)
# Output: kilogram per meter (kg/m)

NISTkgPerMeterTOpoundPerInch<-function(kgPerMeter) kgPerMeter / 17.85797
# Input: kilogram per meter (kg/m)
# Output: pound per inch (lb/in)

NISTpoundPerMinTOkgPerSec<-function(poundPerMin) poundPerMin * 0.007559873
# Input: pound per minute (lb/min)
# Output: kilogram per second (kg/s)

NISTkgPerSecTOpoundPerMin<-function(kgPerSec) kgPerSec / 0.007559873
# Input: kilogram per second (kg/s)
# Output: pound per minute (lb/min)

NISTpoundPerSecTOkgPerSec<-function(poundPerSec) poundPerSec * 0.4535924
# Input: pound per second (lb/s)
# Output: kilogram per second (kg/s)

NISTkgPerSecTOpoundPerSec<-function(kgPerSec) kgPerSec / 0.4535924
# Input: kilogram per second (kg/s)
# Output: pound per second (lb/s)

NISTpoundPerSqrFtTOkgPerSqrMeter<-function(poundPerSqrFt) poundPerSqrFt * 4.882428
# Input: pound per square foot (lb/ft2)
# Output: kilogram per square meter (kg/m2)

NISTkgPerSqrMeterTOpoundPerSqrFt<-function(kgPerSqrMeter) kgPerSqrMeter / 4.882428
# Input: kilogram per square meter (kg/m2)
# Output: pound per square foot (lb/ft2)

NISTpoundPerSqrInchTOkgPerSqrMeter<-function(poundPerSqrInch) poundPerSqrInch * 703.0696
# Input: pound per square inch (not pound-force) (lb/in2)
# Output: kilogram per square meter (kg/m2)

NISTkgPerSqrMeterTOpoundPerSqrInch<-function(kgPerSqrMeter) kgPerSqrMeter / 703.0696
# Input: kilogram per square meter (kg/m2)
# Output: pound per square inch (not pound-force) (lb/in2)

NISTpoundPerYardTOkgPerMeter<-function(poundPerYard) poundPerYard * 0.4960546
# Input: pound per yard (lb/yd)
# Output: kilogram per meter (kg/m)

NISTkgPerMeterTOpoundPerYard<-function(kgPerMeter) kgPerMeter / 0.4960546
# Input: kilogram per meter (kg/m)
# Output: pound per yard (lb/yd)

NISTpsiTOpascal<-function(psi) psi * 6894.757
# Input: psi (pound-force per square inch) (lbf/in2)
# Output: pascal (Pa)

NISTpascalTOpsi<-function(pascal) pascal / 6894.757
# Input: pascal (Pa)
# Output: psi (pound-force per square inch) (lbf/in2)

NISTpsiTOkpascal<-function(psi) psi * 6.894757
# Input: psi (pound-force per square inch) (lbf/in2)
# Output: kilopascal (kPa)

NISTkpascalTOpsi<-function(kpascal) kpascal / 6.894757
# Input: kilopascal (kPa)
# Output: psi (pound-force per square inch) (lbf/in2)

NISTquadTOjoule<-function(quad) quad * 1.055056e+18
# Input: quad (1015 BtuIT) 11
# Output: joule (J)

NISTjouleTOquad<-function(joule) joule / 1.055056e+18
# Input: joule (J)
# Output: quad (1015 BtuIT) 11

NISTquartUSdryTOcubMeter<-function(quartUSdry) quartUSdry * 0.001101221
# Input: quart (U.S. dry) (dry qt)
# Output: cubic meter (m3)

NISTcubMeterTOquartUSdry<-function(cubMeter) cubMeter / 0.001101221
# Input: cubic meter (m3)
# Output: quart (U.S. dry) (dry qt)

NISTquartUSdryTOliter<-function(quartUSdry) quartUSdry * 1.101221
# Input: quart (U.S. dry) (dry qt)
# Output: liter (L)

NISTliterTOquartUSdry<-function(liter) liter / 1.101221
# Input: liter (L)
# Output: quart (U.S. dry) (dry qt)

NISTquartUSliquidTOcubMeter<-function(quartUSliquid) quartUSliquid * 0.0009463529
# Input: quart (U.S. liquid) (liq qt)
# Output: cubic meter (m3)

NISTcubMeterTOquartUSliquid<-function(cubMeter) cubMeter / 0.0009463529
# Input: cubic meter (m3)
# Output: quart (U.S. liquid) (liq qt)

NISTquartUSliquidTOliter<-function(quartUSliquid) quartUSliquid * 0.9463529
# Input: quart (U.S. liquid) (liq qt)
# Output: liter (L)

NISTliterTOquartUSliquid<-function(liter) liter / 0.9463529
# Input: liter (L)
# Output: quart (U.S. liquid) (liq qt)

NISTradTOgray<-function(rad) rad * 0.01
# Input: rad (absorbed dose) (rad)
# Output: gray (Gy)

NISTgrayTOrad<-function(gray) gray / 0.01
# Input: gray (Gy)
# Output: rad (absorbed dose) (rad)

NISTremTOsievert<-function(rem) rem * 0.01
# Input: rem (rem)
# Output: sievert (Sv)

NISTsievertTOrem<-function(sievert) sievert / 0.01
# Input: sievert (Sv)
# Output: rem (rem)

NISTrevolutionTOradian<-function(revolution) revolution * 6.283185
# Input: revolution (r)
# Output: radian (rad)

NISTradianTOrevolution<-function(radian) radian / 6.283185
# Input: radian (rad)
# Output: revolution (r)

NISTrevolutionPerMinTOradianPerSec<-function(revolutionPerMin) revolutionPerMin * 0.1047198
# Input: revolution per minute (rpm) (r/min)
# Output: radian per second (rad/s)

NISTradianPerSecTOrevolutionPerMin<-function(radianPerSec) radianPerSec / 0.1047198
# Input: radian per second (rad/s)
# Output: revolution per minute (rpm) (r/min)

NISTrheTOreciprocalPascalSec1<-function(rhe) rhe * 10
# Input: rhe
# Output: reciprocal pascal second (Pa * s)-1

NISTreciprocalPascalSec1TOrhe<-function(reciprocalPascalSec1) reciprocalPascalSec1 / 10
# Input: reciprocal pascal second (Pa * s)-1
# Output: rhe

NISTrodTOmeter<-function(rod) rod * 5.02921
# Input: rod (based on U.S. survey foot) (rd) 7
# Output: meter (m)

NISTmeterTOrod<-function(meter) meter / 5.02921
# Input: meter (m)
# Output: rod (based on U.S. survey foot) (rd) 7

NISTroentgenTOcoulombPerKg<-function(roentgen) roentgen * 0.000258
# Input: roentgen (R)
# Output: coulomb per kilogram (C/kg)

NISTcoulombPerKgTOroentgen<-function(coulombPerKg) coulombPerKg / 0.000258
# Input: coulomb per kilogram (C/kg)
# Output: roentgen (R)

NISTrpmTOradianPerSec<-function(rpm) rpm * 0.1047198
# Input: rpm (revolution per minute) (r/min)
# Output: radian per second (rad/s)

NISTradianPerSecTOrpm<-function(radianPerSec) radianPerSec / 0.1047198
# Input: radian per second (rad/s)
# Output: rpm (revolution per minute) (r/min)

NISTsecTOradian<-function(sec) sec * 4.848137e-06
# Input: second (angle) ('')
# Output: radian (rad)

NISTradianTOsec<-function(radian) radian / 4.848137e-06
# Input: radian (rad)
# Output: second (angle) ('')

NISTsecTOsec<-function(sec) sec * 0.9972696
# Input: second (sidereal)
# Output: second (s)

NISTsecTOsec<-function(sec) sec / 0.9972696
# Input: second (s)
# Output: second (sidereal)

NISTshakeTOsec<-function(shake) shake * 1e-08
# Input: shake
# Output: second (s)

NISTsecTOshake<-function(sec) sec / 1e-08
# Input: second (s)
# Output: shake

NISTshakeTOnanosec<-function(shake) shake * 10
# Input: shake
# Output: nanosecond (ns)

NISTnanosecTOshake<-function(nanosec) nanosec / 10
# Input: nanosecond (ns)
# Output: shake

NISTslugTOkg<-function(slug) slug * 14.5939
# Input: slug (slug)
# Output: kilogram (kg)

NISTkgTOslug<-function(kg) kg / 14.5939
# Input: kilogram (kg)
# Output: slug (slug)

NISTslugPerCubFtTOkgPerCubMeter<-function(slugPerCubFt) slugPerCubFt * 515.3788
# Input: slug per cubic foot (slug/ft3)
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOslugPerCubFt<-function(kgPerCubMeter) kgPerCubMeter / 515.3788
# Input: kilogram per cubic meter (kg/m3)
# Output: slug per cubic foot (slug/ft3)

NISTslugPerFtSecTOpascalSec<-function(slugPerFtSec) slugPerFtSec * 47.88026
# Input: slug per foot second [slug/(ft * s)]
# Output: pascal second (Pa * s)

NISTpascalSecTOslugPerFtSec<-function(pascalSec) pascalSec / 47.88026
# Input: pascal second (Pa * s)
# Output: slug per foot second [slug/(ft * s)]

NISTsqrFtTOsqrMeter<-function(sqrFt) sqrFt * 0.09290304
# Input: square foot (ft2)
# Output: square meter (m2)

NISTsqrMeterTOsqrFt<-function(sqrMeter) sqrMeter / 0.09290304
# Input: square meter (m2)
# Output: square foot (ft2)

NISTsqrFtPerHourTOsqrMeterPerSec<-function(sqrFtPerHour) sqrFtPerHour * 2.58064e-05
# Input: square foot per hour (ft2/h)
# Output: square meter per second (m2/s)

NISTsqrMeterPerSecTOsqrFtPerHour<-function(sqrMeterPerSec) sqrMeterPerSec / 2.58064e-05
# Input: square meter per second (m2/s)
# Output: square foot per hour (ft2/h)

NISTsqrFtPerSecTOsqrMeterPerSec<-function(sqrFtPerSec) sqrFtPerSec * 0.09290304
# Input: square foot per second (ft2/s)
# Output: square meter per second (m2/s)

NISTsqrMeterPerSecTOsqrFtPerSec<-function(sqrMeterPerSec) sqrMeterPerSec / 0.09290304
# Input: square meter per second (m2/s)
# Output: square foot per second (ft2/s)

NISTsqrInchTOsqrMeter<-function(sqrInch) sqrInch * 0.00064516
# Input: square inch (in2)
# Output: square meter (m2)

NISTsqrMeterTOsqrInch<-function(sqrMeter) sqrMeter / 0.00064516
# Input: square meter (m2)
# Output: square inch (in2)

NISTsqrInchTOsqrCm<-function(sqrInch) sqrInch * 6.4516
# Input: square inch (in2)
# Output: square centimeter (cm2)

NISTsqrCmTOsqrInch<-function(sqrCm) sqrCm / 6.4516
# Input: square centimeter (cm2)
# Output: square inch (in2)

NISTsqrMileTOsqrMeter<-function(sqrMile) sqrMile * 2589988
# Input: square mile (mi2)
# Output: square meter (m2)

NISTsqrMeterTOsqrMile<-function(sqrMeter) sqrMeter / 2589988
# Input: square meter (m2)
# Output: square mile (mi2)

NISTsqrMileTOsqrKm<-function(sqrMile) sqrMile * 2.589988
# Input: square mile (mi2)
# Output: square kilometer (km2)

NISTsqrKmTOsqrMile<-function(sqrKm) sqrKm / 2.589988
# Input: square kilometer (km2)
# Output: square mile (mi2)

NISTsqrMileUSsurveyTOsqrMeter<-function(sqrMileUSsurvey) sqrMileUSsurvey * 2589998
# Input: square mile (based on U.S. survey foot) (mi2) 7
# Output: square meter (m2)

NISTsqrMeterTOsqrMileUSsurvey<-function(sqrMeter) sqrMeter / 2589998
# Input: square meter (m2)
# Output: square mile (based on U.S. survey foot) (mi2) 7

NISTsqrMileUSsurveyTOsqrKm<-function(sqrMileUSsurvey) sqrMileUSsurvey * 2.589998
# Input: square mile (based on U.S. survey foot) (mi2) 7
# Output: square kilometer (km2)

NISTsqrKmTOsqrMileUSsurvey<-function(sqrKm) sqrKm / 2.589998
# Input: square kilometer (km2)
# Output: square mile (based on U.S. survey foot) (mi2) 7

NISTsqrYardTOsqrMeter<-function(sqrYard) sqrYard * 0.8361274
# Input: square yard (yd2)
# Output: square meter (m2)

NISTsqrMeterTOsqrYard<-function(sqrMeter) sqrMeter / 0.8361274
# Input: square meter (m2)
# Output: square yard (yd2)

NISTstatampereTOampere<-function(statampere) statampere * 3.335641e-10
# Input: statampere
# Output: ampere (A)

NISTampereTOstatampere<-function(ampere) ampere / 3.335641e-10
# Input: ampere (A)
# Output: statampere

NISTstatcoulombTOcoulomb<-function(statcoulomb) statcoulomb * 3.335641e-10
# Input: statcoulomb
# Output: coulomb (C)

NISTcoulombTOstatcoulomb<-function(coulomb) coulomb / 3.335641e-10
# Input: coulomb (C)
# Output: statcoulomb

NISTstatfaradTOfarad<-function(statfarad) statfarad * 1.11265e-12
# Input: statfarad
# Output: farad (F)

NISTfaradTOstatfarad<-function(farad) farad / 1.11265e-12
# Input: farad (F)
# Output: statfarad

NISTstathenryTOhenry<-function(stathenry) stathenry * 898755200000
# Input: stathenry
# Output: henry (H)

NISThenryTOstathenry<-function(henry) henry / 898755200000
# Input: henry (H)
# Output: stathenry

NISTstatmhoTOsiemens<-function(statmho) statmho * 1.11265e-12
# Input: statmho
# Output: siemens (S)

NISTsiemensTOstatmho<-function(siemens) siemens / 1.11265e-12
# Input: siemens (S)
# Output: statmho

NISTstatohmTOohm<-function(statohm) statohm * 898755200000
# Input: statohm
# Output: ohm (ohm)

NISTohmTOstatohm<-function(ohm) ohm / 898755200000
# Input: ohm (ohm)
# Output: statohm

NISTstatvoltTOvolt<-function(statvolt) statvolt * 299.7925
# Input: statvolt
# Output: volt (V)

NISTvoltTOstatvolt<-function(volt) volt / 299.7925
# Input: volt (V)
# Output: statvolt

NISTstereTOcubMeter<-function(stere) stere * 1
# Input: stere (st)
# Output: cubic meter (m3)

NISTcubMeterTOstere<-function(cubMeter) cubMeter / 1
# Input: cubic meter (m3)
# Output: stere (st)

NISTstilbTOcandelaPerSqrMeter<-function(stilb) stilb * 10000
# Input: stilb (sb)
# Output: candela per square meter (cd/m2)

NISTcandelaPerSqrMeterTOstilb<-function(candelaPerSqrMeter) candelaPerSqrMeter / 10000
# Input: candela per square meter (cd/m2)
# Output: stilb (sb)

NISTstokesTOmeterSqrdPerSec<-function(stokes) stokes * 1e-04
# Input: stokes (St)
# Output: meter squared per second (m2/s)

NISTmeterSqrdPerSecTOstokes<-function(meterSqrdPerSec) meterSqrdPerSec / 1e-04
# Input: meter squared per second (m2/s)
# Output: stokes (St)

NISTtablespoonTOcubMeter<-function(tablespoon) tablespoon * 1.478676e-05
# Input: tablespoon
# Output: cubic meter (m3)

NISTcubMeterTOtablespoon<-function(cubMeter) cubMeter / 1.478676e-05
# Input: cubic meter (m3)
# Output: tablespoon

NISTtablespoonTOml<-function(tablespoon) tablespoon * 14.78676
# Input: tablespoon
# Output: milliliter (mL)

NISTmlTOtablespoon<-function(ml) ml / 14.78676
# Input: milliliter (mL)
# Output: tablespoon

NISTteaspoonTOcubMeter<-function(teaspoon) teaspoon * 4.928922e-06
# Input: teaspoon
# Output: cubic meter (m3)

NISTcubMeterTOteaspoon<-function(cubMeter) cubMeter / 4.928922e-06
# Input: cubic meter (m3)
# Output: teaspoon

NISTteaspoonTOml<-function(teaspoon) teaspoon * 4.928922
# Input: teaspoon
# Output: milliliter (mL)

NISTmlTOteaspoon<-function(ml) ml / 4.928922
# Input: milliliter (mL)
# Output: teaspoon

NISTtexTOkgPerMeter<-function(tex) tex * 1e-06
# Input: tex
# Output: kilogram per meter (kg/m)

NISTkgPerMeterTOtex<-function(kgPerMeter) kgPerMeter / 1e-06
# Input: kilogram per meter (kg/m)
# Output: tex

NISTthermECtOjoule<-function(thermEC) thermEC * 105506000
# Input: therm (EC) 24
# Output: joule (J)

NISTjouleTOthermEC<-function(joule) joule / 105506000
# Input: joule (J)
# Output: therm (EC) 24

NISTthermUStOjoule<-function(thermUS) thermUS * 105480400
# Input: therm (U.S.) 24
# Output: joule (J)

NISTjouleTOthermUS<-function(joule) joule / 105480400
# Input: joule (J)
# Output: therm (U.S.) 24

NISTtonAssayTOkg<-function(tonAssay) tonAssay * 0.02916667
# Input: ton, assay (AT)
# Output: kilogram (kg)

NISTkgTOtonAssay<-function(kg) kg / 0.02916667
# Input: kilogram (kg)
# Output: ton, assay (AT)

NISTtonAssayTOgram<-function(tonAssay) tonAssay * 29.16667
# Input: ton, assay (AT)
# Output: gram (g)

NISTgramTOtonAssay<-function(gram) gram / 29.16667
# Input: gram (g)
# Output: ton, assay (AT)

NISTtonForceTOnewton<-function(tonForce) tonForce * 8896.443
# Input: ton-force (2000 lbf)
# Output: newton (N)

NISTnewtonTOtonForce<-function(newton) newton / 8896.443
# Input: newton (N)
# Output: ton-force (2000 lbf)

NISTtonForceTOkilonewton<-function(tonForce) tonForce * 8.896443
# Input: ton-force (2000 lbf)
# Output: kilonewton (kN)

NISTkilonewtonTOtonForce<-function(kilonewton) kilonewton / 8.896443
# Input: kilonewton (kN)
# Output: ton-force (2000 lbf)

NISTtonLongTOkg<-function(tonLong) tonLong * 1016.047
# Input: ton, long (2240 lb)
# Output: kilogram (kg)

NISTkgTOtonLong<-function(kg) kg / 1016.047
# Input: kilogram (kg)
# Output: ton, long (2240 lb)

NISTtonLongPerCubYardTOkgPerCubMeter<-function(tonLongPerCubYard) tonLongPerCubYard * 1328.939
# Input: ton, long, per cubic yard
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOtonLongPerCubYard<-function(kgPerCubMeter) kgPerCubMeter / 1328.939
# Input: kilogram per cubic meter (kg/m3)
# Output: ton, long, per cubic yard

NISTtonSItOkg<-function(tonSI) tonSI * 1000
# Input: ton, metric (t)
# Output: kilogram (kg)

NISTkgTOtonSI<-function(kg) kg / 1000
# Input: kilogram (kg)
# Output: ton, metric (t)

NISTtonneTOkg<-function(tonne) tonne * 1000
# Input: tonne (called metric ton in U.S.) (t)
# Output: kilogram (kg)

NISTkgTOtonne<-function(kg) kg / 1000
# Input: kilogram (kg)
# Output: tonne (called metric ton in U.S.) (t)

NISTtonOfRefrigerationTOwatt<-function(tonOfRefrigeration) tonOfRefrigeration * 3516.853
# Input: ton of refrigeration (12 000 BtuIT/h)
# Output: watt (W)

NISTwattTOtonOfRefrigeration<-function(watt) watt / 3516.853
# Input: watt (W)
# Output: ton of refrigeration (12 000 BtuIT/h)

NISTtonOfTNTtOjoule<-function(tonOfTNT) tonOfTNT * 4.184e+09
# Input: ton of TNT (energy equivalent) 25
# Output: joule (J)

NISTjouleTOtonOfTNT<-function(joule) joule / 4.184e+09
# Input: joule (J)
# Output: ton of TNT (energy equivalent) 25

NISTtonRegisterTOcubMeter<-function(tonRegister) tonRegister * 2.831685
# Input: ton, register
# Output: cubic meter (m3)

NISTcubMeterTOtonRegister<-function(cubMeter) cubMeter / 2.831685
# Input: cubic meter (m3)
# Output: ton, register

NISTtonShortTOkg<-function(tonShort) tonShort * 907.1847
# Input: ton, short (2000 lb)
# Output: kilogram (kg)

NISTkgTOtonShort<-function(kg) kg / 907.1847
# Input: kilogram (kg)
# Output: ton, short (2000 lb)

NISTtonShortPerCubYardTOkgPerCubMeter<-function(tonShortPerCubYard) tonShortPerCubYard * 1186.553
# Input: ton, short, per cubic yard
# Output: kilogram per cubic meter (kg/m3)

NISTkgPerCubMeterTOtonShortPerCubYard<-function(kgPerCubMeter) kgPerCubMeter / 1186.553
# Input: kilogram per cubic meter (kg/m3)
# Output: ton, short, per cubic yard

NISTtonShortPerHourTOkgPerSec<-function(tonShortPerHour) tonShortPerHour * 0.2519958
# Input: ton, short, per hour
# Output: kilogram per second (kg/s)

NISTkgPerSecTOtonShortPerHour<-function(kgPerSec) kgPerSec / 0.2519958
# Input: kilogram per second (kg/s)
# Output: ton, short, per hour

NISTtorrTOpascal<-function(torr) torr * 133.3224
# Input: torr (Torr)
# Output: pascal (Pa)

NISTpascalTOtorr<-function(pascal) pascal / 133.3224
# Input: pascal (Pa)
# Output: torr (Torr)

NISTunPoleTOweber<-function(unPole) unPole * 1.256637e-07
# Input: unit pole
# Output: weber (Wb)

NISTweberTOunPole<-function(weber) weber / 1.256637e-07
# Input: weber (Wb)
# Output: unit pole

NISTwattHourTOjoule<-function(wattHour) wattHour * 3600
# Input: watt hour (W * h)
# Output: joule (J)

NISTjouleTOwattHour<-function(joule) joule / 3600
# Input: joule (J)
# Output: watt hour (W * h)

NISTwattPerSqrCmTOwattPerSqrMeter<-function(wattPerSqrCm) wattPerSqrCm * 10000
# Input: watt per square centimeter (W/cm2)
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOwattPerSqrCm<-function(wattPerSqrMeter) wattPerSqrMeter / 10000
# Input: watt per square meter (W/m2)
# Output: watt per square centimeter (W/cm2)

NISTwattPerSqrInchTOwattPerSqrMeter<-function(wattPerSqrInch) wattPerSqrInch * 1550.003
# Input: watt per square inch (W/in2)
# Output: watt per square meter (W/m2)

NISTwattPerSqrMeterTOwattPerSqrInch<-function(wattPerSqrMeter) wattPerSqrMeter / 1550.003
# Input: watt per square meter (W/m2)
# Output: watt per square inch (W/in2)

NISTwattSecTOjoule<-function(wattSec) wattSec * 1
# Input: watt second (W * s)
# Output: joule (J)

NISTjouleTOwattSec<-function(joule) joule / 1
# Input: joule (J)
# Output: watt second (W * s)

NISTyardTOmeter<-function(yard) yard * 0.9144
# Input: yard (yd)
# Output: meter (m)

NISTmeterTOyard<-function(meter) meter / 0.9144
# Input: meter (m)
# Output: yard (yd)

NISTyearTOsec<-function(year) year * 31536000
# Input: year (365 days)
# Output: second (s)

NISTsecTOyear<-function(sec) sec / 31536000
# Input: second (s)
# Output: year (365 days)

NISTyearSiderealTOsec<-function(yearSidereal) yearSidereal * 31558150
# Input: year (sidereal)
# Output: second (s)

NISTsecTOyearSidereal<-function(sec) sec / 31558150
# Input: second (s)
# Output: year (sidereal)

NISTyearTropicalTOsec<-function(yearTropical) yearTropical * 31556930
# Input: year (tropical)
# Output: second (s)

NISTsecTOyearTropical<-function(sec) sec / 31556930
# Input: second (s)
# Output: year (tropical)
