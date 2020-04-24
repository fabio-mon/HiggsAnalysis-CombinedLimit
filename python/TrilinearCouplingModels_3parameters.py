from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from HiggsAnalysis.CombinedLimit.LHCHCGModels import LHCHCGBaseModel
import ROOT, os
import json

## Naming conventions
CMS_to_LHCHCG_Dec = { 
    'hww': 'WW',
    'hzz': 'ZZ',
    'hgg': 'gamgam',
    'hbb': 'bb',
    'hcc': 'cc',
    'htt': 'tautau',
    'hmm': 'mumu',
    'hzg': 'Zgam',
    'hgluglu': 'gluglu',
    'hinv': 'inv',
}
CMS_to_LHCHCG_DecSimple = { 
    'hww': 'WW',
    'hzz': 'ZZ',
    'hgg': 'gamgam',
    'hbb': 'bb',
    'hcc': 'bb',
    'htt': 'tautau',
    'hmm': 'tautau',
    'hzg': 'gamgam',
    'hgluglu': 'bb',
    'hinv': 'inv',
}
CMS_to_LHCHCG_Prod = { 
    'ggH': 'ggF',
    'qqH': 'VBF',
    'WH': 'WH',
    'ZH': 'qqZH',
    'ggZH': 'ggZH',
    'ttH': 'ttH',
    'tHq': 'tHjb',
    'tHW': 'WtH',
    'bbH': 'bbH',
 } 

datadir = os.environ['CMSSW_BASE']+"/src/HiggsAnalysis/CombinedLimit/data/"
fit_dict={}
with open(datadir+"/trilinearHiggsModel/kl_fit_parameters_MAR1620.json","r") as fit_json: 
    fit_dict = json.load(fit_json)


class TrilinearHiggsKappaVKappaF(LHCHCGBaseModel):
    "assume the SM coupling but let the Higgs mass to float"
    def __init__(self,BRU=True):
        LHCHCGBaseModel.__init__(self) 
        self.doBRU = BRU
        self.fit_dict = fit_dict

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
        
 
    def getHiggsSignalYieldScale(self,production,decay,energy):

        name = "kVkFkl_XSBRscal_%s_%s_%s" % (production,decay,energy)
        if self.modelBuilder.out.function(name) == None:
                
            # now make production scaling --> taken from Tab. 2 of https://arxiv.org/pdf/1607.04251v1.pdf, using formula from https://arxiv.org/pdf/1709.08649.pdf (eqn 18)
	    cXSmap_7   = {"ggH":0.66e-2,"qqH":0.65e-2,"WH":1.06e-2,"ZH":1.23e-2,"ttH":3.87e-2,"VH":(0.5*(1.06e-2+1.23e-2))}
	    cXSmap_8   = {"ggH":0.66e-2,"qqH":0.65e-2,"WH":1.05e-2,"ZH":1.22e-2,"ttH":3.78e-2,"VH":(0.5*(1.05e-2+1.22e-2))}
	    cXSmap_13  = {"ggH":0.66e-2,"qqH":0.64e-2,"WH":1.03e-2,"ZH":1.19e-2,"ttH":3.51e-2,"VH":(0.5*(1.03e-2+1.19e-2))}
	    EWKmap_13  = {"ggH":1.049,"qqH":0.932,"WH":0.93,"ZH":0.947,"ttH":1.014,"VH":(0.5*(0.93+0.947))}
            cXSmaps = {"7TeV":cXSmap_7, "8TeV":cXSmap_8, "13TeV":cXSmap_13}
            dZH = -1.536e-3
             
            if production in [ "ggZH", "tHq", "tHW"]: 
                XSscal = ("@0", "Scaling_%s_%s" % (production,energy) )

            elif production in [ "ggH", "qqH" ]:
                C1_map = cXSmaps[energy]
                EWK = EWKmap_13[production]
                self.modelBuilder.factory_("expr::kVkFkl_XSscal_%s_%s(\"(@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,Scaling_%s_%s)"\
                                           %(production,energy,C1_map[production],EWK,dZH,production,energy))
                XSscal = ("@0", "kVkFkl_XSscal_%s_%s, " % (production,energy) )

            elif production in [ "ZH", "WH","VH"]: 
                C1_map = cXSmaps[energy]
                EWK = EWKmap_13[production]
                self.modelBuilder.factory_("expr::kVkFkl_XSscal_%s_%s(\"(@1*@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,kappa_V)"\
                                           %(production,energy,C1_map[production],EWK,dZH))
                XSscal = ("@0", "kVkFkl_XSscal_%s_%s, " % (production,energy) )

            elif production == "ttH": 
                C1_map = cXSmaps[energy]
                EWK = EWKmap_13[production]
                self.modelBuilder.factory_("expr::kVkFkl_XSscal_%s_%s(\"(@1*@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,kappa_F)"\
                                           %(production,energy,C1_map[production],EWK,dZH))
                XSscal = ("(0.+(@0>0)*@0)", "kVkFkl_XSscal_%s_%s, " % (production,energy) )

            elif production == "bbH": 
                XSscal = ("@0*@0", "kappa_F")

            elif "HH_" in production: #only for the HH case the production is HH_binname without _13TeV 
                #Get the bin name
                binname = production.replace("HH_","")
                p=[]
                for ipar in range(0,3):
                    p.append(self.fit_dict[binname]['p%d'%ipar])
                self.modelBuilder.factory_("expr::kVkFkl_XSscal_%s_%s(\"(@0*@0*@1*@1*(%g)+@0*@1*@1*@1*(%g)+@1*@1*@1*@1*(%g))\",kappa_lambda,kappa_F)"%(production,energy,p[0],p[1],p[2]))
                XSscal = ("(0.+(@0>0)*@0)", "kVkFkl_XSscal_%s_%s, " % (production,energy) )

            else: raise RuntimeError, "Production %s not supported" % production
            
            BRscal = decay
            if decay == "hss": BRscal = "hbb"
            if not self.modelBuilder.out.function("kVkFkl_BRscal_"+BRscal):
                raise RuntimeError, "Decay mode %s not supported" % decay
            
            self.modelBuilder.factory_('expr::%s("%s*@1", %s, kVkFkl_BRscal_%s)' % (name, XSscal[0], XSscal[1], BRscal))
            print '[LHC-HCG Kappas]', name, production, decay, energy,": ",
            self.modelBuilder.out.function(name).Print("")
        return name

##################################################################################################################
##################################################################################################################

class TrilinearHiggsKappaVKappat(LHCHCGBaseModel):
    "assume the SM coupling but let the Higgs mass to float"
    def __init__(self,BRU=True):
        LHCHCGBaseModel.__init__(self) 
        self.doBRU = BRU
        self.fit_dict = fit_dict
    def setPhysicsOptions(self,physOptions):
        self.setPhysicsOptionsBase(physOptions)
        for po in physOptions:
            if po.startswith("BRU="):
                self.doBRU = (po.replace("BRU=","") in [ "yes", "1", "Yes", "True", "true" ])
        print "BR uncertainties in partial widths: %s " % self.doBRU
    def doParametersOfInterest(self):
        """Create POI out of signal strength and MH"""
        self.modelBuilder.doVar("kappa_V[1,0.0,2.0]")
        self.modelBuilder.doVar("kappa_t[1,-2.0,2.0]")
        self.modelBuilder.doVar("kappa_lambda[1,-20,20]")
        pois = 'kappa_V,kappa_t,kappa_lambda'
        self.doMH()
        self.modelBuilder.doSet("POI",pois)
        self.SMH = SMHiggsBuilder(self.modelBuilder)
        self.setup()

    def setup(self):
        
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
        self.SMH.makeScaling("tHq", CW='kappa_V', Ctop="kappa_t")
        self.SMH.makeScaling("tHW", CW='kappa_V', Ctop="kappa_t")
        self.SMH.makeScaling("ggZH",CZ='kappa_V', Ctop="kappa_t",Cb="1")
        self.SMH.makeScaling('ggH', Cb='1', Ctop='kappa_t', Cc="1")
        self.SMH.makeScaling('hgluglu', Cb='1', Ctop='kappa_t')
        self.SMH.makeScaling('hgg', Cb='1', Ctop='kappa_t', CW='kappa_V', Ctau='1')
        self.SMH.makeScaling('hzg', Cb='1', Ctop='kappa_t', CW='kappa_V', Ctau='1')


        cGammap = {"hgg":0.49e-2,"hzz":0.83e-2,"hww":0.73e-2,"hgluglu":0.66e-2,"htt":0,"hbb":0,"hcc":0,"hmm":0}
        
        # First we need to create the terms that account for the self-coupling --> Just scale partial width first - https://arxiv.org/abs/1709.08649 Eq 22.
        # probably a better way to code this since the partial width expressions are being repeated when we write the BR 
        for dec in cGammap.keys(): 
            valC1 = cGammap[dec]
            self.modelBuilder.factory_('expr::kl_scalBR_%s("(@0-1)*%g",kappa_lambda)' % (dec,valC1))

        # next make the partial widths, also including the kappas -> we want to include the term from the normal kappas and the one from the self-coupling 
        self.modelBuilder.factory_('expr::kVktkl_Gscal_Z("(@0*@0+@3)*@1*@2", kappa_V, SM_BR_hzz, HiggsDecayWidth_UncertaintyScaling_hzz, kl_scalBR_hzz)')
        self.modelBuilder.factory_('expr::kVktkl_Gscal_W("(@0*@0+@3)*@1*@2", kappa_V, SM_BR_hww, HiggsDecayWidth_UncertaintyScaling_hww, kl_scalBR_hww)')
        self.modelBuilder.factory_('expr::kVktkl_Gscal_tau("(1+@4)*@0*@2 + (1+@5)*@1*@3", SM_BR_htt, SM_BR_hmm, HiggsDecayWidth_UncertaintyScaling_htt, HiggsDecayWidth_UncertaintyScaling_hmm,kl_scalBR_htt, kl_scalBR_hmm)')
        self.modelBuilder.factory_('expr::kVktkl_Gscal_top("(1+@2)*@0*@1", SM_BR_hcc, HiggsDecayWidth_UncertaintyScaling_hcc, kl_scalBR_hcc)')
        self.modelBuilder.factory_('expr::kVktkl_Gscal_bottom("(1+@3) * (@0*@2+@1)", SM_BR_hbb, SM_BR_hss, HiggsDecayWidth_UncertaintyScaling_hbb, kl_scalBR_hbb)')
        self.modelBuilder.factory_('expr::kVktkl_Gscal_gluon("  (@0+@3)  * @1 * @2", Scaling_hgluglu, SM_BR_hgluglu, HiggsDecayWidth_UncertaintyScaling_hgluglu, kl_scalBR_hgluglu)')
        self.modelBuilder.factory_('expr::kVktkl_Gscal_gamma("(@0+@6)*@1*@4 + @2*@3*@5",  Scaling_hgg, SM_BR_hgg, Scaling_hzg, SM_BR_hzg, HiggsDecayWidth_UncertaintyScaling_hgg, HiggsDecayWidth_UncertaintyScaling_hzg, kl_scalBR_hgg)') # no kappa_lambda dependance on H->zg known yet ?
        # fix to have all BRs add up to unity
        self.modelBuilder.factory_("sum::kVktkl_SMBRs(%s)" %  (",".join("SM_BR_"+X for X in "hzz hww htt hmm hcc hbb hss hgluglu hgg hzg".split())))
        self.modelBuilder.out.function("kVktkl_SMBRs").Print("")        

        ## total witdh, normalized to the SM one (just the sum over the partial widths/SM total BR)
        self.modelBuilder.factory_('expr::kVktkl_Gscal_tot("(@0+@1+@2+@3+@4+@5+@6)/@7", kVktkl_Gscal_Z, kVktkl_Gscal_W, kVktkl_Gscal_tau, kVktkl_Gscal_top, kVktkl_Gscal_bottom, kVktkl_Gscal_gluon, kVktkl_Gscal_gamma, kVktkl_SMBRs)') 

        ## BRs, normalized to the SM ones: they scale as (partial/partial_SM) / (total/total_SM) 
        self.modelBuilder.factory_('expr::kVktkl_BRscal_hww("(@0*@0+@3)*@2/@1", kappa_V, kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hww, kl_scalBR_hww)')
        self.modelBuilder.factory_('expr::kVktkl_BRscal_hzz("(@0*@0+@3)*@2/@1", kappa_V, kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzz, kl_scalBR_hzz)')
        self.modelBuilder.factory_('expr::kVktkl_BRscal_htt("(1+@2)*@1/@0", kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_htt, kl_scalBR_htt)')
        self.modelBuilder.factory_('expr::kVktkl_BRscal_hmm("(1+@2)*@1/@0", kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hmm, kl_scalBR_hmm)')
        self.modelBuilder.factory_('expr::kVktkl_BRscal_hbb("(1+@2)*@1/@0", kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hbb, kl_scalBR_hbb)')
        self.modelBuilder.factory_('expr::kVktkl_BRscal_hcc("(1+@2)*@1/@0", kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hcc, kl_scalBR_hcc)')
        self.modelBuilder.factory_('expr::kVktkl_BRscal_hgg("(@0+@3)*@2/@1", Scaling_hgg, kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgg,kl_scalBR_hgg)')
        self.modelBuilder.factory_('expr::kVktkl_BRscal_hzg("@0*@2/@1", Scaling_hzg, kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hzg)')
        self.modelBuilder.factory_('expr::kVktkl_BRscal_hgluglu("(@0+@3)*@2/@1", Scaling_hgluglu, kVktkl_Gscal_tot, HiggsDecayWidth_UncertaintyScaling_hgluglu, kl_scalBR_hgluglu)')
        
 
    def getHiggsSignalYieldScale(self,production,decay,energy):

        name = "kVktkl_XSBRscal_%s_%s_%s" % (production,decay,energy)
        if self.modelBuilder.out.function(name) == None:
                
            # now make production scaling --> taken from Tab. 2 of https://arxiv.org/pdf/1607.04251v1.pdf, using formula from https://arxiv.org/pdf/1709.08649.pdf (eqn 18)
	    cXSmap_7   = {"ggH":0.66e-2,"qqH":0.65e-2,"WH":1.06e-2,"ZH":1.23e-2,"ttH":3.87e-2,"VH":(0.5*(1.06e-2+1.23e-2))}
	    cXSmap_8   = {"ggH":0.66e-2,"qqH":0.65e-2,"WH":1.05e-2,"ZH":1.22e-2,"ttH":3.78e-2,"VH":(0.5*(1.05e-2+1.22e-2))}
	    cXSmap_13  = {"ggH":0.66e-2,"qqH":0.64e-2,"WH":1.03e-2,"ZH":1.19e-2,"ttH":3.51e-2,"VH":(0.5*(1.03e-2+1.19e-2))}
	    EWKmap_13  = {"ggH":1.049,"qqH":0.932,"WH":0.93,"ZH":0.947,"ttH":1.014,"VH":(0.5*(0.93+0.947))}
            cXSmaps = {"7TeV":cXSmap_7, "8TeV":cXSmap_8, "13TeV":cXSmap_13}
            dZH = -1.536e-3
             
            if production in [ "ggZH", "tHq", "tHW"]: 
                XSscal = ("@0", "Scaling_%s_%s" % (production,energy) )
            
            elif production in [ "ggH", "qqH" ]:
                C1_map = cXSmaps[energy]
                EWK = EWKmap_13[production]
                self.modelBuilder.factory_("expr::kVktkl_XSscal_%s_%s(\"(@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,Scaling_%s_%s)"\
                                           %(production,energy,C1_map[production],EWK,dZH,production,energy))
                XSscal = ("@0", "kVktkl_XSscal_%s_%s, " % (production,energy) )

            elif production in [ "ZH", "WH", "VH"]: 
                C1_map = cXSmaps[energy]
                EWK = EWKmap_13[production]
                self.modelBuilder.factory_("expr::kVktkl_XSscal_%s_%s(\"(@1*@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,kappa_V)"\
                                           %(production,energy,C1_map[production],EWK,dZH))
                XSscal = ("@0", "kVktkl_XSscal_%s_%s, " % (production,energy) )

            elif production == "ttH": 
                C1_map = cXSmaps[energy]
                EWK = EWKmap_13[production]
                self.modelBuilder.factory_("expr::kVktkl_XSscal_%s_%s(\"(@1*@1+(@0-1)*%g/%g)/((1-(@0*@0-1)*%g))\",kappa_lambda,kappa_t)"\
                                           %(production,energy,C1_map[production],EWK,dZH))
                XSscal = ("(0.+(@0>0)*@0)", "kVktkl_XSscal_%s_%s, " % (production,energy) )

            elif "HH_" in production: #only for the HH case the production is HH_binname without _13TeV 
                #Get the bin name
                binname = production.replace("HH_","")
                p=[]
                for ipar in range(0,3):
                    p.append(self.fit_dict[binname]['p%d'%ipar])
                self.modelBuilder.factory_("expr::kVktkl_XSscal_%s_%s(\"(@0*@0*@1*@1*(%g)+@0*@1*@1*@1*(%g)+@1*@1*@1*@1*(%g))\",kappa_lambda,kappa_t)"%(production,energy,p[0],p[1],p[2]))
                XSscal = ("(0.+(@0>0)*@0)", "kVktkl_XSscal_%s_%s, " % (production,energy) )

            else: raise RuntimeError, "Production %s not supported" % production
            
            BRscal = decay
            if decay == "hss": BRscal = "hbb"
            if not self.modelBuilder.out.function("kVktkl_BRscal_"+BRscal):
                raise RuntimeError, "Decay mode %s not supported" % decay
            
            self.modelBuilder.factory_('expr::%s("%s*@1", %s, kVktkl_BRscal_%s)' % (name, XSscal[0], XSscal[1], BRscal))
            print '[LHC-HCG Kappas]', name, production, decay, energy,": ",
            self.modelBuilder.out.function(name).Print("")
        return name

trilinearHiggskVkF = TrilinearHiggsKappaVKappaF()
trilinearHiggskVkt = TrilinearHiggsKappaVKappat()

