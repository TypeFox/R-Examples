#Function written by Robert Pearce for NCCT, December, 2013
#SBML version of vLiver PBPK model by James Sluka
export_pbtk_sbml <- function(chem.cas=NULL,
                             chem.name=NULL,
                             species="Human",
                             initial.amounts=list(Agutlumen=0),
                             filename="default.xml", 
                             digits = 4)
{
  Agutlumen <- Aart <- Aven <- Alung <- Agut <- Aliver <- Akidney <- Arest <- Atubules <- Ametabolized <- NULL
  for (this.compartment in c("Agutlumen","Aart","Aven","Alung","Agut","Aliver","Akidney","Arest","Atubules","Ametabolized"))
  {
    if (this.compartment %in% names(initial.amounts)) 
    {
      eval(parse(text=paste(this.compartment,"<-",initial.amounts[[this.compartment]])))
    }
    else eval(parse(text=paste(this.compartment,"<- 0")))
  }
  
  inlist <- parameterize_pbtk(chem.cas=chem.cas,chem.name=chem.name,species=species)
  inlist[["Qcardiac"]] <- inlist[["Qcardiacc"]] * 24 * inlist[["BW"]]^0.75
  out <- get_chem_id(chem.cas=chem.cas,chem.name=chem.name)
  chem.cas <- out$chem.cas
  chem.name <- out$chem.name
  
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<sbml xmlns = \"http://www.sbml.org/sbml/level2/version4\" level = \"2\" version = \"4\">
   <model id = \"cell\">
      <listOfCompartments>
         <compartment id = \"compartment\" size = \"1\"/>
      </listOfCompartments>
      <listOfSpecies>
         <species id = \"Aart\" boundaryCondition = \"false\" initialConcentration = \"",Aart,"\" compartment = \"compartment\"/>
         <species id = \"Agut\" boundaryCondition = \"false\" initialConcentration = \"",Agut,"\" compartment = \"compartment\"/>
         <species id = \"Agutlumen\" boundaryCondition = \"false\" initialConcentration = \"",Agutlumen,"\" compartment = \"compartment\"/>
         <species id = \"Alung\" boundaryCondition = \"false\" initialConcentration = \"",Alung,"\" compartment = \"compartment\"/>
         <species id = \"Aven\" boundaryCondition = \"false\" initialConcentration = \"",Aven,"\" compartment = \"compartment\"/>
         <species id = \"Arest\" boundaryCondition = \"false\" initialConcentration = \"",Arest,"\" compartment = \"compartment\"/>
         <species id = \"Aliver\" boundaryCondition = \"false\" initialConcentration = \"",Aliver,"\" compartment = \"compartment\"/>
         <species id = \"Ametabolized\" boundaryCondition = \"false\" initialConcentration = \"",Ametabolized,"\" compartment = \"compartment\"/>
         <species id = \"Akidney\" boundaryCondition = \"false\" initialConcentration = \"",Akidney,"\" compartment = \"compartment\"/>
         <species id = \"Atubules\" boundaryCondition = \"false\" initialConcentration = \"",Atubules,"\" compartment = \"compartment\"/>
        </listOfSpecies>
        <listOfParameters>
         <parameter id = \"Qgut\" value = \"",signif(inlist$Qgutf * inlist$Qcardiac,digits),"\"/>
         <parameter id = \"Vart\" value = \"",signif(inlist$Vartc * inlist$BW,digits),"\"/>
         <parameter id = \"kgutabs\" value = \"",signif(24*inlist$kgutabs,digits),"\"/>
         <parameter id = \"Qcardiac\" value = \"",signif(inlist$Qcardiac,digits),"\"/>
         <parameter id = \"Vlung\" value = \"",signif(inlist$Vlungc * inlist$BW,digits),"\"/>
         <parameter id = \"Vven\" value = \"",signif(inlist$Vvenc * inlist$BW,digits),"\"/>
         <parameter id = \"Qrest\" value = \"",signif(inlist$Qcardiac*(1-inlist$Qgutf-inlist$Qliverf-inlist$Qkidneyf),digits),"\"/>
         <parameter id = \"Vrest\" value = \"",signif(inlist$Vrestc * inlist$BW,digits),"\"/>
         <parameter id = \"Rblood2plasma\" value = \"",signif(inlist$Rblood2plasma,digits),"\"/>
         <parameter id = \"Krest2pu\" value = \"",signif(inlist$Krest2pu,digits),"\"/>
         <parameter id = \"Fraction_unbound_plasma\" value = \"",signif(inlist$Funbound.plasma,digits),"\"/>
         <parameter id = \"Qliver\" value = \"",signif(inlist$Qliverf * inlist$Qcardiac,digits),"\"/>
         <parameter id = \"Clmetabolism\" value = \"",signif(24*inlist$Clmetabolism*inlist$BW,digits),"\"/>
         <parameter id = \"Vliver\" value = \"",signif(inlist$Vliverc * inlist$BW,digits),"\"/>
         <parameter id = \"Kliver2pu\" value = \"",signif(inlist$Kliver2pu,digits),"\"/>
         <parameter id = \"Vgut\" value = \"",signif(inlist$Vgutc * inlist$BW,digits),"\"/>
         <parameter id = \"Qkidney\" value = \"",signif(inlist$Qkidneyf * inlist$Qcardiac,digits),"\"/>
         <parameter id = \"Qgfr\" value = \"",signif(24* inlist$Qgfrc * inlist$BW^0.75 ,digits),"\"/>
         <parameter id = \"Vkidney\" value = \"",signif(inlist$Vkidneyc * inlist$BW,digits),"\"/>
         <parameter id = \"Kkidney2pu\" value = \"",signif(inlist$Kkidney2pu,digits),"\"/>
         <parameter id = \"Klung2pu\" value = \"",signif(inlist$Klung2pu,digits),"\"/>
         <parameter id = \"Kgut2pu\" value = \"",signif(inlist$Kgut2pu,digits),"\"/>
      </listOfParameters>
      <listOfReactions>
         <reaction id = \"J1\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Aart\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Agut\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <times/>
                        <ci>
                              Qgut
                        </ci>
                        <ci>
                              Aart
                        </ci>
                     </apply>
                     <ci>
                           Vart
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J2\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Agutlumen\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Agut\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <times/>
                     <ci>
                           kgutabs
                     </ci>
                     <ci>
                           Agutlumen
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         
         <reaction id = \"J3\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Alung\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Aart\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <divide/>
                        <apply>
                           <times/>
                           <apply>
                              <divide/>
                              <apply>
                                 <times/>
                                 <ci>
                                       Qcardiac
                                 </ci>
                                 <ci>
                                       Alung
                                 </ci>
                              </apply>
                              <ci>
                                    Vlung
                              </ci>
                           </apply>
                           <ci>
                                 Rblood2plasma
                           </ci>
                        </apply>
                        <ci>
                              Klung2pu
                        </ci>
                     </apply>
                     <ci>
                           Fraction_unbound_plasma
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J4\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Aven\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Alung\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <times/>
                        <ci>
                              Qcardiac
                        </ci>
                        <ci>
                              Aven
                        </ci>
                     </apply>
                     <ci>
                           Vven
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J5\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Aart\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Arest\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <times/>
                        <ci>
                              Qrest
                        </ci>
                        <ci>
                              Aart
                        </ci>
                     </apply>
                     <ci>
                           Vart
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J6\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Arest\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Aven\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <divide/>
                        <apply>
                           <times/>
                           <apply>
                              <divide/>
                              <apply>
                                 <times/>
                                 <ci>
                                       Qrest
                                 </ci>
                                 <ci>
                                       Arest
                                 </ci>
                              </apply>
                              <ci>
                                    Vrest
                              </ci>
                           </apply>
                           <ci>
                                 Rblood2plasma
                           </ci>
                        </apply>
                        <ci>
                              Krest2pu
                        </ci>
                     </apply>
                     <ci>
                           Fraction_unbound_plasma
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J7\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Aart\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Aliver\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <times/>
                        <ci>
                              Qliver
                        </ci>
                        <ci>
                              Aart
                        </ci>
                     </apply>
                     <ci>
                           Vart
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J13\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Aliver\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Ametabolized\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <divide/>
                           <apply>
                              <times/>
                              <ci>
                                    Clmetabolism
                              </ci>
                              <ci>
                                    Aliver
                              </ci>
                           </apply>
                           <ci>
                                 Vliver
                           </ci>
                        </apply>
                        <ci>
                              Kliver2pu
                        </ci>
                     </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J8\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Agut\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Aliver\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <divide/>
                        <apply>
                           <times/>
                           <apply>
                              <divide/>
                              <apply>
                                 <times/>
                                 <ci>
                                       Qgut
                                 </ci>
                                 <ci>
                                       Agut
                                 </ci>
                              </apply>
                              <ci>
                                    Vgut
                              </ci>
                           </apply>
                           <ci>
                                 Rblood2plasma
                           </ci>
                        </apply>
                        <ci>
                              Kgut2pu
                        </ci>
                     </apply>
                     <ci>
                           Fraction_unbound_plasma
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J9\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Aliver\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Aven\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <divide/>
                        <apply>
                           <times/>
                           <apply>
                              <divide/>
                              <apply>
                                 <times/>
                                 <apply>
                                    <plus/>
                                    <ci>
                                          Qliver
                                    </ci>
                                    <ci>
                                          Qgut
                                    </ci>
                                 </apply>
                                 <ci>
                                       Aliver
                                 </ci>
                              </apply>
                              <ci>
                                    Vliver
                              </ci>
                           </apply>
                           <ci>
                                 Rblood2plasma
                           </ci>
                        </apply>
                        <ci>
                              Kliver2pu
                        </ci>
                     </apply>
                     <ci>
                           Fraction_unbound_plasma
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J10\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Aart\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Akidney\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <times/>
                        <ci>
                              Qkidney
                        </ci>
                        <ci>
                              Aart
                        </ci>
                     </apply>
                     <ci>
                           Vart
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J11\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Akidney\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Atubules\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <divide/>
                        <apply>
                           <times/>
                           <ci>
                                 Qgfr
                           </ci>
                           <ci>
                                 Akidney
                           </ci>
                        </apply>
                        <ci>
                              Vkidney
                        </ci>
                     </apply>
                     <ci>
                           Kkidney2pu
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
         <reaction id = \"J12\" reversible = \"false\">
            <listOfReactants>
               <speciesReference species = \"Akidney\" stoichiometry = \"1\"/>
            </listOfReactants>
            <listOfProducts>
               <speciesReference species = \"Aven\" stoichiometry = \"1\"/>
            </listOfProducts>
            <kineticLaw>
               <math xmlns = \"http://www.w3.org/1998/Math/MathML\">
                  <apply>
                     <divide/>
                     <apply>
                        <divide/>
                        <apply>
                           <times/>
                           <apply>
                              <divide/>
                              <apply>
                                 <times/>
                                 <ci>
                                       Qkidney
                                 </ci>
                                 <ci>
                                       Akidney
                                 </ci>
                              </apply>
                              <ci>
                                    Vkidney
                              </ci>
                           </apply>
                           <ci>
                                 Rblood2plasma
                           </ci>
                        </apply>
                        <ci>
                              Kkidney2pu
                        </ci>
                     </apply>
                     <ci>
                           Fraction_unbound_plasma
                     </ci>
                  </apply>
               </math>
            </kineticLaw>
         </reaction>
      </listOfReactions>
   </model>
</sbml>",file=filename)
}
