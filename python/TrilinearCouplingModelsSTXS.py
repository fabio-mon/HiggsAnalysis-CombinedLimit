from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from HiggsAnalysis.CombinedLimit.LHCHCGModels import LHCHCGBaseModel
import ROOT, os, json

def getGenProdDecMode(bin,process,options):
    """Return a triple of (production, decay, energy)"""
    #assuming that process names have a form like 'STXSprocname_year_decay', e.g. 'ggH_0J_PTH_GT10_2016_hgg' 
    print "hi I am getGenProdDecMode"
    print "doing bin ", bin, "process ", process
    processSource = process.rsplit('_',2)[0] 
    decaySource = process.rsplit('_',2)[2] 
    foundEnergy = None
    for D in [ '7TeV', '8TeV', '13TeV', '14TeV' ]:
	if D in decaySource:
	    if foundEnergy: raise RuntimeError, "Validation Error: decay string %s contains multiple known energies" % decaySource
	    foundEnergy = D
    if not foundEnergy:
	for D in [ '7TeV', '8TeV', '13TeV', '14TeV' ]:
	    if D in options.fileName+":"+bin:
		if foundEnergy: raise RuntimeError, "Validation Error: decay string %s contains multiple known energies" % decaySource
		foundEnergy = D
    if foundEnergy:
        if foundEnergy != '13TeV':
            raise RuntimeError, "Validation Error: only 13TeV energy supported"
    else:
        print "[WARNING]: no energy found -> Assuming 13 TeV"
        foundEnergy = '13TeV'

    print processSource, decaySource, foundEnergy
    return (processSource, decaySource, foundEnergy)

class TrilinearHiggsKappaVKappaFSTXS12(LHCHCGBaseModel):
    "assume the SM coupling but let the Higgs mass to float"
    def __init__(self,BRU=True):
        print "hi I am TrilinearHiggsKappaVKappaFSTXS12"
        LHCHCGBaseModel.__init__(self) 
        self.doBRU = BRU
        self.STXSScalingFunctions = {}
    def setPhysicsOptions(self,physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = (po.replace("BRU=","") in [ "yes", "1", "Yes", "True", "true" ])
        print "BR uncertainties in partial widths: %s " % self.doBRU
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.modelBuilder.doVar("kappa_V[1,0.0,2.0]")
        self.modelBuilder.doVar("kappa_F[1,-2.0,2.0]")
        self.modelBuilder.doVar("kappa_lambda[1,-20,20]")
        pois = 'kappa_V,kappa_F,kappa_lambda'
        self.doMH()
        self.modelBuilder.doSet("POI",pois)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        
        print "setup"
	#self.dobbH()
        # SM BR

        for d in SM_HIGG_DECAYS + [ "hss" ]:
            self.SMH.makeBR(d)

        # BR uncertainties
        if self.doBRU:
            self.SMH.makePartialWidthUncertainties()
        else:
            for d in SM_HIGG_DECAYS:
                self.modelBuilder.factory_('HiggsDecayWidth_UncertaintyScaling_%s[1.0]' % d)

        # fix to have all BRs add up to unity
        self.modelBuilder.factory_("sum::c7_SMBRs(%s)" %  (",".join("SM_BR_"+X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())))
        self.modelBuilder.out.function("c7_SMBRs").Print("")

        # get VBF, tHq, tHW, ggZH cross section and resolved loops
        self.SMH.makeScaling('qqH', CW='kappa_V', CZ='kappa_V')
        self.SMH.makeScaling("tHq", CW='kappa_V', Ctop="kappa_F")
        self.SMH.makeScaling("tHW", CW='kappa_V', Ctop="kappa_F")
        self.SMH.makeScaling("ggZH",CZ='kappa_V', Ctop="kappa_F",Cb="kappa_F")
        self.SMH.makeScaling('ggH', Cb='kappa_F', Ctop='kappa_F', Cc="kappa_F")
        self.SMH.makeScaling('hgluglu', Cb='kappa_F', Ctop='kappa_F')
        self.SMH.makeScaling('hgg', Cb='kappa_F', Ctop='kappa_F', CW='kappa_V', Ctau='kappa_F')
        self.SMH.makeScaling('hzg', Cb='kappa_F', Ctop='kappa_F', CW='kappa_V', Ctau='kappa_F')


	cGammap = {"hgg":0.49e-2,"hzz":0.83e-2,"hww":0.73e-2,"hgluglu":0.66e-2,"htt":0,"hbb":0,"hcc":0,"hmm":0}
	
	# First we need to create the terms that account for the self-coupling --> Just scale partial width first - https://arxiv.org/abs/1709.08649 Eq 22.
	# probably a better way to code this since the partial width expressions are being repeated when we write the BR 
        for dec in cGammap.keys(): 
	   valC1 = cGammap[dec]
	   self.modelBuilder.factory_('expr::kl_scalBR_%s("(@0-1)*%g",kappa_lambda)' % (dec,valC1))

	# next make the partial widths, also including the kappas -> we want to include the term from the normal kappas and the one from the self-coupling 
        self.modelBuilder.factory_('expr::kVkFkl_Gscal_Z("(@0*@0+@3)*@1*@2", kappa_V, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz, kl_scalBR_hzz)')
        self.modelBuilder.factory_('expr::kVkFkl_Gscal_W("(@0*@0+@3)*@1*@2", kappa_V, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww, kl_scalBR_hww)')
        self.modelBuilder.factory_('expr::kVkFkl_Gscal_tau("(@0*@0+@6)*@1*@4 + (@2*@2+@7)*@3*@5", kappa_F, SM_BR_htt, kappa_F, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm,kl_scalBR_htt, kl_scalBR_hmm)')
        self.modelBuilder.factory_('expr::kVkFkl_Gscal_top("(@0*@0+@3)*@1*@2", kappa_F, SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc, kl_scalBR_hcc)')
        self.modelBuilder.factory_('expr::kVkFkl_Gscal_bottom("(@0*@0+@4) * (@1*@3+@2)", kappa_F, SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb, kl_scalBR_hbb)')
        self.modelBuilder.factory_('expr::kVkFkl_Gscal_gluon("  (@0+@3)  * @1 * @2", Scaling_hgluglu, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu, kl_scalBR_hgluglu)')
        self.modelBuilder.factory_('expr::kVkFkl_Gscal_gamma("(@0+@6)*@1*@4 + @2*@3*@5",  Scaling_hgg, SM_BR_hgg, Scaling_hzg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg, kl_scalBR_hgg)') # no kappa_lambda dependance on H->zg known yet ?
        # fix to have all BRs add up to unity
        self.modelBuilder.factory_("sum::kVkFkl_SMBRs(%s)" %  (",".join("SM_BR_"+X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())))
        self.modelBuilder.out.function("kVkFkl_SMBRs").Print("")        

        ## total witdh, normalized to the SM one (just the sum over the partial widths/SM total BR)
        self.modelBuilder.factory_('expr::kVkFkl_Gscal_tot("(@0+@1+@2+@3+@4+@5+@6)/@7", kVkFkl_Gscal_Z, kVkFkl_Gscal_W, kVkFkl_Gscal_tau, kVkFkl_Gscal_top, kVkFkl_Gscal_bottom, kVkFkl_Gscal_gluon, kVkFkl_Gscal_gamma, kVkFkl_SMBRs)') 

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM) 
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_hww("(@0*@0+@3)*@2/@1", kappa_V, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hww, kl_scalBR_hww)')
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_hzz("(@0*@0+@3)*@2/@1", kappa_V, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzz, kl_scalBR_hzz)')
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_htt("(@0*@0+@3)*@2/@1", kappa_F, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt, kl_scalBR_htt)')
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_hmm("(@0*@0+@3)*@2/@1", kappa_F, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hmm, kl_scalBR_hmm)')
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_hbb("(@0*@0+@3)*@2/@1", kappa_F, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hbb, kl_scalBR_hbb)')
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_hcc("(@0*@0+@3)*@2/@1", kappa_F, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc, kl_scalBR_hcc)')
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_hgg("(@0+@3)*@2/@1", Scaling_hgg, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgg,kl_scalBR_hgg)')
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_hzg("@0*@2/@1", Scaling_hzg, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzg)')
        self.modelBuilder.factory_('expr::kVkFkl_BRscal_hgluglu("(@0+@3)*@2/@1", Scaling_hgluglu, kVkFkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgluglu, kl_scalBR_hgluglu)')
	
        #Build the production XS scalings for the STXS bins 
        trilinearcoeffs = {}
        jsonfile = open(os.path.join(self.SMH.datadir, "../trilinearHiggsModel/TrilinearCoeffSTXS.json"),"r")
        trilinearcoeffs = json.load(jsonfile) # trilinear_coeff[STXSproc]["C1"] and trilinear_coeff[STXSproc]["EWK"]  
        dZH = -1.536e-3
        energy="13TeV"
        jsonfile.close()

        for production in SM_HIGG_PROD:
            print "building scaling for ",production
            if production in ["WPlusH","WMinusH","VH"]: continue
            
            elif production in  [ "ggZH", "tHW"]: #trilinear scaling is not available --> use only scaling from SMH
                self.STXSScalingFunctions[production]="Scaling_%s_%s"%(production,energy)

            elif production in [ "ggH", "qqH", "tHq"]: #the scaling built by SMH combined with inclusive trilinear scaling
                C1 =  trilinearcoeffs[production]["C1"]
                EWK = trilinearcoeffs[production]["EWK"]
                self.modelBuilder.factory_("expr::kVkFkl_XSscal_%s_%s(\"(@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,Scaling_%s_%s)"\
	       				%(production,energy,C1,EWK,dZH,production,energy))
                self.STXSScalingFunctions[production]="kVkFkl_XSscal_%s_%s"%(production,energy)

	    elif production in [ "ZH", "WH"]: #k-scaling combined with trilinear scaling specific for each stxs bin 
                for STXSprocname in trilinearcoeffs.keys():
                    if not STXSprocname.startswith(production): continue
                    C1  =  trilinearcoeffs[str(STXSprocname)]["C1"]
                    EWK =  trilinearcoeffs[str(STXSprocname)]["EWK"]
                    self.modelBuilder.factory_("expr::kVkFkl_XSscal_%s_%s(\"(@1*@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,kappa_V)"\
                                               %(str(STXSprocname),energy,C1,EWK,dZH))
                    self.STXSScalingFunctions[str(STXSprocname)]="kVkFkl_XSscal_%s_%s"%(str(STXSprocname),energy)

            elif production == "ttH":  #k-scaling combined with trilinear scaling specific for each stxs bin 
                for STXSprocname in trilinearcoeffs.keys():
                    if not STXSprocname.startswith(production): continue
                    C1  =  trilinearcoeffs[str(STXSprocname)]["C1"]
                    EWK =  trilinearcoeffs[str(STXSprocname)]["EWK"]
                    self.modelBuilder.factory_("expr::kVkFkl_XSscal_%s_%s(\"(@1*@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,kappa_F)"\
                                               %(str(STXSprocname),energy,C1,EWK,dZH))
                    self.STXSScalingFunctions[str(STXSprocname)]="kVkFkl_XSscal_%s_%s"%(str(STXSprocname),energy)

            elif production == "bbH":  #k-scaling
                self.modelBuilder.factory_("expr::kVkFkl_XSscal_%s_%s(\"@0*@0\",kappa_F)"\
                                           %(production,energy))
                self.STXSScalingFunctions[production]="kVkFkl_XSscal_%s_%s"%(production,energy)

            else: raise RuntimeError, "Production %s not supported" % production
            

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        print "Hi, this is getYieldScale"
        if not self.DC.isSignal[process]: return 1
        (processSource, foundDecay, foundEnergy) = getGenProdDecMode(bin,process,self.options)
        return self.getHiggsSignalYieldScale(processSource, foundDecay, foundEnergy)


    def getHiggsSignalYieldScale(self,processSource,decay,energy):

        print "Hi, this is getHiggsSignalYieldScale"
        print "retriving scaling for ", processSource, decay, energy
        production = processSource.split('_')[0]

        if production in [ "ggZH", "tHq", "tHW", "ggH", "qqH", "bbH"]: # scale the inclusive XS
            XSscal = self.STXSScalingFunctions[production]
        elif production in [ "ZH", "WH", "ttH"]: # scale the specific STXS bin
            XSscal = self.STXSScalingFunctions[processSource]
        else:
            raise RuntimeError, "Production %s not supported" % production


        if decay == "hss": 
            decay = "hbb"
        BRscal = "kVkFkl_BRscal_"+decay
        if not self.modelBuilder.out.function(BRscal):
            raise RuntimeError, "Decay mode %s not supported" % decay
        
        XSBRscaling = "%s_%s"%(XSscal,BRscal)
        if self.modelBuilder.out.function(XSBRscaling) == None:
	    self.modelBuilder.factory_('expr::%s("@0*@1", %s, %s)'%(XSBRscaling, XSscal, BRscal))
            self.modelBuilder.out.function(XSBRscaling).Print("")
        print "XSscal = ", XSscal, type(XSscal)
        print "BRscal = ", BRscal, type(BRscal)
        print "XSBRscal = ", XSBRscaling, type(XSBRscaling)
        print "------------"
        return XSBRscaling

TrilinearHiggsKappaVKappaFSTXS12 = TrilinearHiggsKappaVKappaFSTXS12()
