import os,sys,copy,math
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--outdir",  default='/afs/cern.ch/user/f/fmonti/work/NewflashggFinalFit/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/out/MAR1620_allsignals_onlyHHcats/',help="output directory" )
parser.add_option("--datacard",default='/afs/cern.ch/user/f/fmonti/work/NewflashggFinalFit/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/out/MAR1620_allsignals_onlyHHcats/SMdatacard.root',      help="datacard directory" )
parser.add_option("--klmin_scan",    type="float", default=-20., help="kl min for scan" )
parser.add_option("--klmax_scan",    type="float", default=+20., help="kl max for scan" )
parser.add_option("--ktmin_scan",    type="float", default=-2.,  help="kt min for scan" )
parser.add_option("--ktmax_scan",    type="float", default=+2.,  help="kt max for scan" )
parser.add_option("--Njobs",         type="int",   default=1000, help="Total number of jobs to submit")
parser.add_option("--Npointsperjob", type="int",   default=10,   help="Total number of points to compute per jobs")
parser.add_option("--unblind",       action="store_true", default=False, help="Unblind")
parser.add_option("--jonathontricks",action="store_true", default=False, help="Jonathon tricks to speed up the likelihood scan")
parser.add_option("--resetpdfindexes",action="store_true", default=False, help="re-set pdf indexes using post-fit result")
parser.add_option("--maskTTHcategories",action="store_true", default=False, help="maskTTHcategories")
(options,args)=parser.parse_args()

##############################################                                                                                    
#text2workspace command
text2workspace_command = 'text2workspace.py %s -P HiggsAnalysis.CombinedLimit.TrilinearCouplingModels:trilinearHiggskVkF  -m 125 --PO \'BRU=false\' --channel-masks '%options.datacard.replace(".root",".txt")
if options.maskTTHcategories:
    text2workspace_command +=" --channel-masks "
print text2workspace_command

options.datacard = options.datacard.replace(".txt",".root")

##############################################                                                                                                       
#generate an Asimov distribution for the SM expectation of sig (klambda = 1) + bk
#freezing_option_str = " --freezeNuisances kappa_V "
asimov_command =  "combine -M GenerateOnly -t -1 -m 125.00 --saveToys -n SM_toys "
asimov_command += " --setParameters kappa_lambda=1,kappa_t=1"

if options.maskTTHcategories:
    asimov_command += ",mask_TTHHadronicTag_0_13TeV=1,mask_TTHHadronicTag_1_13TeV=1,mask_TTHHadronicTag_2_13TeV=1,mask_TTHHadronicTag_3_13TeV=1,mask_TTHLeptonicTag_0_13TeV=1,mask_TTHLeptonicTag_1_13TeV=1,mask_TTHLeptonicTag_2_13TeV=1,mask_TTHLeptonicTag_3_13TeV=1"

if options.resetpdfindexes:
#    asimov_command += ",pdfindex_DoubleHTag_0_13TeV=0,pdfindex_DoubleHTag_1_13TeV=1,pdfindex_DoubleHTag_2_13TeV=0,pdfindex_DoubleHTag_3_13TeV=0,pdfindex_DoubleHTag_4_13TeV=0,pdfindex_DoubleHTag_5_13TeV=0,pdfindex_DoubleHTag_6_13TeV=1,pdfindex_DoubleHTag_7_13TeV=0,pdfindex_DoubleHTag_8_13TeV=0,pdfindex_DoubleHTag_9_13TeV=0,pdfindex_DoubleHTag_10_13TeV=0,pdfindex_DoubleHTag_11_13TeV=2,pdfindex_TTHHadronicTag_0_13TeV=2,pdfindex_TTHHadronicTag_1_13TeV=0,pdfindex_TTHHadronicTag_2_13TeV=2,pdfindex_TTHHadronicTag_3_13TeV=4,pdfindex_TTHLeptonicTag_0_13TeV=1,pdfindex_TTHLeptonicTag_1_13TeV=3,pdfindex_TTHLeptonicTag_2_13TeV=3,pdfindex_TTHLeptonicTag_3_13TeV=3 "
    asimov_command += ",pdfindex_DoubleHTag_0_13TeV=0,pdfindex_DoubleHTag_1_13TeV=4,pdfindex_DoubleHTag_2_13TeV=9,pdfindex_DoubleHTag_3_13TeV=5,pdfindex_DoubleHTag_4_13TeV=3,pdfindex_DoubleHTag_5_13TeV=9,pdfindex_DoubleHTag_6_13TeV=7,pdfindex_DoubleHTag_7_13TeV=9,pdfindex_DoubleHTag_8_13TeV=6,pdfindex_DoubleHTag_9_13TeV=12,pdfindex_DoubleHTag_10_13TeV=6,pdfindex_DoubleHTag_11_13TeV=6,pdfindex_TTHHadronicTag_0_13TeV=2,pdfindex_TTHHadronicTag_1_13TeV=0,pdfindex_TTHHadronicTag_2_13TeV=2,pdfindex_TTHHadronicTag_3_13TeV=4,pdfindex_TTHLeptonicTag_0_13TeV=1,pdfindex_TTHLeptonicTag_1_13TeV=3,pdfindex_TTHLeptonicTag_2_13TeV=3,pdfindex_TTHLeptonicTag_3_13TeV=3 "
else:
    asimov_command += " "

if options.jonathontricks:
    asimov_command += " --X-rt MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --X-rtd TMCSO_AdaptivePseudoAsimov=0 --X-rtd TMCSO_PseudoAsimov=0  --cminDefaultMinimizerStrategy 0  --cminFallbackAlgo Minuit2,Migrad,0:0.1 "

#asimov_command += freezing_option_str+" "
asimov_command += options.datacard
asimov_name = "higgsCombineSM_toys.GenerateOnly.mH125.123456.root"
if not options.unblind:
    print "asimov command: "+asimov_command+" && mv "+asimov_name+" "+options.outdir
asimov_name = os.path.abspath(options.outdir+"/"+asimov_name)
#os.system(asimov_command)

##############################################                                                                                               

#run the likelihood scan splitting in multiple jobs
#i want to run a scan with kt fixed to 1, one with kl fixed to 1, and a 2d scan
likelihoodscan2D_jobs_dir = "%s/likelihood_jobs/scan2D_exp/"%options.outdir
likelihoodscan2D_out_dir = "%s/likelihood_out/scan2D_exp/"%options.outdir
if options.unblind:
    likelihoodscan2D_jobs_dir = "%s/likelihood_jobs/scan2D_obs/"%options.outdir
    likelihoodscan2D_out_dir = "%s/likelihood_out/scan2D_obs/"%options.outdir

os.system("mkdir -p "+likelihoodscan2D_jobs_dir)
os.system("mkdir -p "+likelihoodscan2D_out_dir)

#freezing_option_str_2D = freezing_option_str.replace("kl,kt,","")

likelihood_command = "combine -M MultiDimFit --algo grid -m 125.00 "
likelihood_command += " --floatOtherPOIs 0 "
if not options.unblind:
    likelihood_command += " -t -1 --toysFile %s "%asimov_name
likelihood_command += " -P kappa_lambda -P kappa_t "
if options.maskTTHcategories:
    likelihood_command += " --setParameters mask_TTHHadronicTag_0_13TeV=1,mask_TTHHadronicTag_1_13TeV=1,mask_TTHHadronicTag_2_13TeV=1,mask_TTHHadronicTag_3_13TeV=1,mask_TTHLeptonicTag_0_13TeV=1,mask_TTHLeptonicTag_1_13TeV=1,mask_TTHLeptonicTag_2_13TeV=1,mask_TTHLeptonicTag_3_13TeV=1 "


likelihood_command += " --setParameterRanges kappa_lambda=%f,%f:kappa_t=%f,%f "%(options.klmin_scan,options.klmax_scan,options.ktmin_scan,options.ktmax_scan)
#likelihood_command += " --setParameters kappa_V=1 "
likelihood_command += " --cminDefaultMinimizerStrategy 0  --cminFallbackAlgo Minuit2,Migrad,0:0.1 "
likelihood_command += options.datacard

Njobs=options.Njobs
Npointsperjob_2Dscan=options.Npointsperjob

for ijob in range(0,Njobs):

    likelihood_command_2D = likelihood_command
#    likelihood_command_2D += freezing_option_str_2D+" "
    likelihood_command_2D += " --points %i "%(Njobs*Npointsperjob_2Dscan)
    likelihood_command_2D += " --firstPoint=%i --lastPoint=%i "%(ijob*Npointsperjob_2Dscan,(ijob+1)*Npointsperjob_2Dscan-1)
    likelihood_command_2D += " -n MultiDim_2Dscan_Job%d "%(ijob)
    if options.jonathontricks:
        likelihood_command_2D += " --X-rt MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --X-rtd TMCSO_AdaptivePseudoAsimov=0 --X-rtd TMCSO_PseudoAsimov=0 "
    mv_command = "mv *MultiDim_2Dscan_Job%d*.root %s"%(ijob,os.path.abspath(likelihoodscan2D_out_dir))

    likelihoodscan2D_scriptname = likelihoodscan2D_jobs_dir+"/job_2Dscan_%i.sh"%(ijob)
    likelihoodscan2D_script = open(likelihoodscan2D_scriptname ,'w')
    likelihoodscan2D_script.write('#!/bin/sh\ncd %s\neval `scramv1 runtime -sh`\ncd -\n'%os.getcwd())
    likelihoodscan2D_script.write('touch %s.run\n'%os.path.abspath(likelihoodscan2D_scriptname))
    likelihoodscan2D_script.write('if ( %s && %s ) then\n'%(likelihood_command_2D,mv_command))
    likelihoodscan2D_script.write('\t echo "DONE" \n')
    likelihoodscan2D_script.write('\t touch %s.done\n'%os.path.abspath(likelihoodscan2D_scriptname))
    likelihoodscan2D_script.write('else\n')
    likelihoodscan2D_script.write('\t echo "FAIL" \n')
    likelihoodscan2D_script.write('\t touch %s.fail\n'%os.path.abspath(likelihoodscan2D_scriptname))
    likelihoodscan2D_script.write('fi\n')
    likelihoodscan2D_script.write('cd -\n')
    likelihoodscan2D_script.write('echo "RM RUN "\n')
    likelihoodscan2D_script.write('rm -f %s.run\n'%os.path.abspath(likelihoodscan2D_scriptname))
    likelihoodscan2D_script.write("echo DONE \n")
    likelihoodscan2D_script.close()
    os.system("chmod +x "+likelihoodscan2D_scriptname)

print("writing .sub file")
condorsubname = likelihoodscan2D_jobs_dir+"/submit_jobs.sub"
condorsub = open(condorsubname,'w')
#condorsub.write('requirements = (OpSysAndVer =?= "SLCern6")\n')
condorsub.write("executable            = $(scriptname)\n")
condorsub.write("output                = $(scriptname).out\n")
condorsub.write("error                 = $(scriptname).err\n")
condorsub.write("log                   = $(scriptname).log\n")
condorsub.write('+JobFlavour           = "workday"\n')
condorsub.write("queue scriptname matching %s/job_2Dscan_*.sh"%(os.path.abspath(likelihoodscan2D_jobs_dir)))
condorsub.close()
print("condor_submit "+condorsubname)
