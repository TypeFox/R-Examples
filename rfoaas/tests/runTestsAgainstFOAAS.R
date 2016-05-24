##  runTestsAgainstFOOAS --- simple Continuous Integration tests
##
##  Copyright (C) 2015  Dirk Eddelbuettel <edd@debian.org>
##
##  This file is part of rfoaas
##
##  rfoaas is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 2 of the License, or
##  (at your option) any later version.
##
##  rfoaas is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with rfoaas.  If not, see <http://www.gnu.org/licenses/>.

## By default, do not run the tests
## which also means do not run on CRAN
runTests <- FALSE

## Use the Travis / GitHub integrations as we set this
## environment variable to "yes" in .travis.yml
##
## Set this variable manually if you want to run the tests
##
if (Sys.getenv("RunFOAASTests=yes") == "yes") runTests <- TRUE

## Also run the tests when building on Dirk's box, even whem
## the environment variable is not set
if (isTRUE(unname(Sys.info()["user"])=="edd")) runTests <- TRUE

if (runTests) {

    library(rfoaas)

    ## placeholders
    name      <- "Someone"
    from      <- "Me"
    reference <- "Something somewhere"
    company   <- "XYZ Corp"
    tool      <- "Some magic thing"
    do        <- "Get"
    something <- "Something"

    ## basic operations
    stopifnot(off         (name=name, from=from)  == "Fuck off, Someone. - Me") 
    stopifnot(you         (name=name, from=from)  == "Fuck you, Someone. - Me")
    stopifnot(this        (from=from)             == "Fuck this. - Me")
    stopifnot(that        (from=from)             == "Fuck that. - Me")
    stopifnot(everything  (from=from)             == "Fuck everything. - Me")
    stopifnot(everyone    (from=from)             == "Everyone can go and fuck off. - Me")
    stopifnot(donut       (name=name, from=from)  == "Someone, go and take a flying fuck at a rolling donut. - Me")
    stopifnot(shakespeare (name=name, from=from)  == "Someone, Thou clay-brained guts, thou knotty-pated fool, thou whoreson obscene greasy tallow-catch! - Me")
    stopifnot(linus       (name=name, from=from)  == "Someone, there aren't enough swear-words in the English language, so now I'll have to call you perkeleen vittupää just to express my disgust and frustration with this crap. - Me" )
    stopifnot(king        (name=name, from=from)  == "Oh fuck off, just really fuck off you total dickface. Christ Someone, you are fucking thick. - Me" )
    stopifnot(pink        (name=name)             == "Well, Fuck me pink. - Someone")
    stopifnot(life        (name=name)             == "Fuck my life. - Someone")
    stopifnot(chainsaw    (name=name, from=from)  == "Fuck me gently with a chainsaw, Someone. Do I look like Mother Teresa? - Me")
    stopifnot(outside     (name=name, from=from)  == "Someone, why don't you go outside and play hide-and-go-fuck-yourself? - Me")
    stopifnot(thanks      (from=from)             == "Fuck you very much. - Me")
    stopifnot(flying      (from=from)             == "I don't give a flying fuck. - Me")
    stopifnot(fascinating (from=from)             == "Fascinating story, in what chapter do you shut the fuck up? - Me")
    stopifnot(madison     (name=name, from=from)  == "What you've just said is one of the most insanely idiotic things I have ever heard, Someone. At no point in your rambling, incoherent response were you even close to anything that could be considered a rational thought. Everyone in this room is now dumber for having listened to it. I award you no points Someone, and may God have mercy on your soul. - Me")
    stopifnot(cool        (from=from)             == "Cool story, bro. - Me")
    stopifnot(field       (name=name, from=from,
                           reference=reference)   == "And Someone said unto Me, 'Verily, cast thine eyes upon the field in which I grow my fucks', and Me gave witness unto the field, and saw that it was barren. - Something somewhere")
    
    stopifnot(nugget      (name=name, from=from)  == "Well Someone, aren't you a shining example of a rancid fuck-nugget. - Me")
    stopifnot(yoda        (name=name, from=from)  == "Fuck off, you must, Someone. - Me")
    stopifnot(ballmer     (name=name, company=company,
                           from=from)             == "Fucking Someone is a fucking pussy. I'm going to fucking bury that guy, I have done it before, and I will do it again. I'm going to fucking kill XYZ Corp. - Me")
    
    stopifnot(what        (from=from)             == "What the fuck‽ - Me")
    stopifnot(because     (from=from)             == "Why? Because Fuck you, that's why. - Me")
    stopifnot(caniuse     (tool=tool, from=from)  == "Can you use Some magic thing? Fuck no! - Me")
    stopifnot(bye         (from=from)             == "Fuckity bye! - Me")
    stopifnot(diabetes    (from=from)             == "I'd love to stop and chat to you but I'd rather have type 2 diabetes. - Me")
    stopifnot(bus         (from=from)             == "Fuck bus. - Me")
    stopifnot(xmas        (name=name, from=from)  == "Merry Fucking Christmas, Someone. - Me")
    stopifnot(awesome     (from=from)             == "This is Fucking Awesome. - Me")
    
    stopifnot(tucker      (from=from)             == "Come the fuck in or fuck the fuck off. - Me")
    stopifnot(bucket      (from=from)             == "Please choke on a bucket of cocks. - Me")
    stopifnot(bday        (name=name, from=from)  == "Happy Fucking Birthday, Someone. - Me")
    stopifnot(family_     (from=from)             == "Fuck you, your whole family, your pets, and your feces. - Me")
    stopifnot(shutup      (name=name, from=from)  == "Someone, shut the fuck up. - Me")
    stopifnot(zayn        (from=from)             == "Ask me if I give a motherfuck ?!! - Me")
    stopifnot(dalton      (name=name, from=from)  == "Someone: A fucking problem solving super-hero. - Me")
    stopifnot(dosomething (do=do, something=something,
                           from=from)             == "Get the fucking Something! - Me")
    ##off_with
    stopifnot(retard      (from=from)             == "You Fucktard! - Me")
    stopifnot(thumbs      (name=name, from=from)  == "Who has two thumbs and doesn't give a fuck? Someone. - Me")
    stopifnot(thing       (name=name, from=from)  == "Fuck Someone. - Me")
    
    ## shoutcloud
    stopifnot(off         (name=name, from=from, filter="shoutcloud")  == "FUCK OFF, SOMEONE. - ME") 

    ## language
    #stopifnot(off         (name=name, from=from, language="de")  == "Verpiss dich, jemand. - Me")
    
    ## shoutcloud and language
    #stopifnot(off         (name=name, from=from, filter="shoutcloud", language="de")  == "VERPISS DICH, JEMAND. - ME")
    
}

