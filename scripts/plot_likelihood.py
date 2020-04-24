#!/bin/python
import os,sys,copy,math
import json
from array import array
import ROOT as r
from optparse import OptionParser
from math import sqrt

r.gStyle.SetPalette(1)
r.gStyle.SetOptStat(0)
r.gStyle.SetOptTitle(0)
r.gStyle.SetNumberContours(200)
r.gROOT.SetBatch(True)

def extra_texts():
    # print "... drawing extra texts"
    ### extra text
    cmstextfont   = 61  # font of the "CMS" label
    cmstextsize   = 0.05  # font size of the "CMS" label
    chantextsize = 18
    extratextfont = 52     # for the "preliminary"
    extratextsize = 0.76 * cmstextsize # for the "preliminary"
    lumitextfont  = 42
    cmstextinframe = False

    yoffset = -0.046

    lumibox = r.TLatex  (0.9, 0.964+yoffset, "136.8 fb^{-1} (13 TeV)")
    lumibox.SetNDC()
    lumibox.SetTextAlign(31)
    lumibox.SetTextSize(extratextsize)
    lumibox.SetTextFont(lumitextfont)
    lumibox.SetTextColor(r.kBlack)

    # xpos  = 0.177
    xpos  = 0.137
    if cmstextinframe:
        ypos  = 0.94 ## inside the frame
    else:
        ypos  = 0.995  ## ouside the frame
    CMSbox = r.TLatex  (xpos, ypos+yoffset+0.01, "CMS")       
    CMSbox.SetNDC()
    CMSbox.SetTextSize(cmstextsize)
    CMSbox.SetTextFont(cmstextfont)
    CMSbox.SetTextColor(r.kBlack)
    CMSbox.SetTextAlign(13) ## inside the frame

    # simBox = r.TLatex  (xpos, ypos - 0.05+yoffset, "Simulation")
    simBox = r.TLatex  (xpos + 0.12, ypos+yoffset, "")
    simBox.SetNDC()
    simBox.SetTextSize(extratextsize)
    simBox.SetTextFont(extratextfont)
    simBox.SetTextColor(r.kBlack)
    simBox.SetTextAlign(13)


    channelLabel = r.TLatex  (0.6-0.05, 0.8, "")#"HH #rightarrow #gamma#gammab#bar{b}")
    channelLabel.SetNDC()
    # channelLabel.SetTextAlign(31)
    channelLabel.SetTextSize(1.15*extratextsize)
    channelLabel.SetTextFont(lumitextfont)
    channelLabel.SetTextColor(r.kBlack)


    return [lumibox, CMSbox, simBox, channelLabel]
    # lumibox.Draw()
    # CMSbox.Draw()
    # simBox.Draw()
    # channelLabel.Draw()

def fixemptybins(histo2D):
    print "in function fixemptybins" 
    xbinwidth = histo2D.GetXaxis().GetBinWidth(1)
    ybinwidth = histo2D.GetYaxis().GetBinWidth(1)
    for xbin in range(1,histo2D.GetXaxis().GetNbins()+1):
        x = histo2D.GetXaxis().GetBinCenter(xbin)
        for ybin in range(1,histo2D.GetYaxis().GetNbins()+1):
            y = histo2D.GetYaxis().GetBinCenter(ybin)
            ibin = histo2D.FindBin(x,y)
            if histo2D.GetBinEntries(ibin)==0:
                print "found empty bin at %f,%f"%(x,y)
                ibinxp1 = histo2D.FindBin(x+xbinwidth,y)
                ibinxm1 = histo2D.FindBin(x-xbinwidth,y)
                ibinyp1 = histo2D.FindBin(x,y+ybinwidth)
                ibinym1 = histo2D.FindBin(x,y-ybinwidth)
#                avg=0.
#                if histo2D.GetBinEntries(ibinxm1)!=0 and histo2D.GetBinEntries(ibinxp1)!=0 and histo2D.GetBinEntries(ibinym1)!=0 and histo2D.GetBinEntries(ibinyp1)!=0:
#                    avg = (histo2D.GetBinContent(ibinxm1)+histo2D.GetBinContent(ibinxp1)+histo2D.GetBinContent(ibinym1)+histo2D.GetBinContent(ibinyp1))/4.
 
                xavg=0.
                yavg=0.
                avg=0.
                if histo2D.GetBinEntries(ibinxm1)!=0 and histo2D.GetBinEntries(ibinxp1)!=0:
                    xavg = 0.5*(histo2D.GetBinContent(ibinxm1) + histo2D.GetBinContent(ibinxp1))
                if histo2D.GetBinEntries(ibinym1)!=0 and histo2D.GetBinEntries(ibinyp1)!=0:
                    yavg = 0.5*(histo2D.GetBinContent(ibinym1) + histo2D.GetBinContent(ibinyp1))
                if xavg!=0 and yavg!=0:
                    avg = 0.5*(xavg+yavg)
                elif xavg!=0:
                    avg = xavg
                elif yavg!=0:
                    avg = yavg
 
                if avg!=0:
                    histo2D.SetBinContent(ibin,avg)
                    histo2D.SetBinEntries(ibin,1)
                    print "fixed with value %f!"%avg
                else: 
                    print "can't fix this bin !!!"


parser = OptionParser()
parser.add_option("--infilename_2Dscan",               default="",         help="Input file(s) for 2D scan")
parser.add_option("--infilename_klscan",               default="",         help="Input file(s) for kl scan")
parser.add_option("--infilename_ktscan",               default="",         help="Input file(s) for kt scan")
parser.add_option("--infilename_obs",                  default="",         help="Input file(s) with observed result")
parser.add_option("--infilename_exp",                  default="",         help="Input file(s) with expected result")
parser.add_option("--outdir",                          default="./plots/", help="Output dir for plots")
parser.add_option("--Npoints_2Dscan",    type="int",   default=10000,      help="Number of points for 2D scan")
parser.add_option("--klmin",             type="float", default=1.,         help="kl min")
parser.add_option("--klmax",             type="float", default=1.,         help="kl max")
parser.add_option("--ktmin",             type="float", default=1.,         help="kt min")
parser.add_option("--ktmax",             type="float", default=1.,         help="kt max")
(options,args)=parser.parse_args()

klmin = options.klmin
klmax = options.klmax
ktmin = options.ktmin
ktmax = options.ktmax

outfile = r.TFile("%s/likelihood_scans.root"%options.outdir,"RECREATE")
outfile.cd()

if options.infilename_obs != "" and options.infilename_exp != "":
    ''' compare observed with expected result '''
    c = r.TCanvas('c', 'c', 600, 600)
    c.SetFrameLineWidth(3)
    c.SetBottomMargin(0.13)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.13)
    c.SetFrameLineWidth(3)
    c.SetBottomMargin(0.13)
    c.SetLeftMargin(0.13)    
    c.SetGridx()    
    c.SetGridy()    

    obsfile = r.TFile(options.infilename_obs,"READ")
    contour68_obs = obsfile.Get("contour68")
    contour95_obs = obsfile.Get("contour95")
    Bestfit_obs   = obsfile.Get("Bestfit")
    klscan_obs    = obsfile.Get("klscan") 
    ktscan_obs    = obsfile.Get("ktscan") 

    expfile = r.TFile(options.infilename_exp,"READ")
    contour68_exp = expfile.Get("contour68")
    contour95_exp = expfile.Get("contour95")
    Bestfit_exp   = expfile.Get("Bestfit")
    klscan_exp    = expfile.Get("klscan") 
    ktscan_exp    = expfile.Get("ktscan") 

    contour68_obs.SetLineColor(1)
    contour95_obs.SetLineColor(1)
    Bestfit_obs.SetMarkerColor(1)
    contour68_exp.SetLineColor(4)
    contour95_exp.SetLineColor(4)
    Bestfit_exp.SetMarkerColor(4)

    legend_2D = r.TLegend()
    legend_2D.AddEntry(contour68_obs,"obs. 68%C.L.","L")
    legend_2D.AddEntry(contour95_obs,"obs. 95%C.L.","L")
    legend_2D.AddEntry(contour68_exp,"exp. 68%C.L.","L")
    legend_2D.AddEntry(contour95_exp,"exp. 95%C.L.","L")

    contour68_obs.Draw("cont3")
    contour95_obs.Draw("cont3 SAME")
    Bestfit_obs.Draw("P SAME")
    contour68_exp.Draw("cont3 SAME")
    contour95_exp.Draw("cont3 SAME")
    Bestfit_exp.Draw("P SAME")
    legend_2D.Draw()
    contour68_obs.GetXaxis().SetTitle("#kappa_{#lambda}")
    contour68_obs.GetYaxis().SetTitle("#kappa_{t}")
    et = extra_texts()
    for t in et: t.Draw()
    c.Print("%s/compare_2Dscan.pdf"%options.outdir)
    c.Print("%s/compare_2Dscan.png"%options.outdir)

    klscan_obs.SetLineColor(1)    
    klscan_exp.SetLineColor(4)    
    legend_1D = r.TLegend()
    legend_1D.AddEntry(klscan_obs,"obs.","L")
    legend_1D.AddEntry(klscan_exp,"exp.","L")

    klscan_obs.Draw("APC")
    klscan_exp.Draw("PC same")
    legend_1D.Draw()
    klscan_obs.GetXaxis().SetTitle("#kappa_{#lambda}")
    klscan_obs.GetYaxis().SetTitle("-2#Deltaln(L)")
    klscan_obs.GetXaxis().SetLimits(-7,13)
    klscan_obs.GetYaxis().SetRangeUser(0,5)
    for t in et: t.Draw()

    xmin=klscan_obs.GetXaxis().GetXmin()
    xmax=klscan_obs.GetXaxis().GetXmax()
    ### lines
    sigmas = [1,1.96]
    CL = [68,95]
    lines = []
    for s in sigmas:
        l = r.TLine(xmin, s*s, xmax, s*s)
        l.SetLineStyle(7)
        l.SetLineWidth(1)
        l.SetLineColor(r.kGray+2)
        lines.append(l)
    for l in lines: l.Draw()
    c.SetGridy(0)
    c.Print("%s/compare_klscan.pdf"%options.outdir)
    c.Print("%s/compare_klscan.png"%options.outdir)


    ktscan_obs.SetLineColor(1)    
    ktscan_exp.SetLineColor(4)    
    legend_1D = r.TLegend()
    legend_1D.AddEntry(ktscan_obs,"obs.","L")
    legend_1D.AddEntry(ktscan_exp,"exp.","L")

    ktscan_obs.Draw("APC")
    ktscan_exp.Draw("PC same")
    legend_1D.Draw()
    ktscan_obs.GetXaxis().SetTitle("#kappa_{#lambda}")
    ktscan_obs.GetYaxis().SetTitle("-2#Deltaln(L)")
    ktscan_obs.GetXaxis().SetLimits(-1.,2.2)
    ktscan_obs.GetYaxis().SetRangeUser(0,5)
    for t in et: t.Draw()

    xmin=ktscan_obs.GetXaxis().GetXmin()
    xmax=ktscan_obs.GetXaxis().GetXmax()
    ### lines
    sigmas = [1,1.96]
    CL = [68,95]
    lines = []
    for s in sigmas:
        l = r.TLine(xmin, s*s, xmax, s*s)
        l.SetLineStyle(7)
        l.SetLineWidth(1)
        l.SetLineColor(r.kGray+2)
        lines.append(l)
    for l in lines: l.Draw()
    c.SetGridy(0)
    c.Print("%s/compare_ktscan.pdf"%options.outdir)
    c.Print("%s/compare_ktscan.png"%options.outdir)


if options.infilename_2Dscan != "":
    c0 = r.TCanvas('c0', 'c0', 600, 600)
    c0.SetFrameLineWidth(3)
    c0.SetBottomMargin(0.13)
    c0.SetLeftMargin(0.13)
    c0.SetRightMargin(0.13)
    c0.SetFrameLineWidth(3)
    c0.SetBottomMargin(0.13)
    c0.SetLeftMargin(0.13)    
    c0.SetGridx()    
    c0.SetGridy()    
    inchain = r.TChain("limit","limit")
    for infilename in options.infilename_2Dscan.split(","):
        inchain.Add(infilename)
    Npoints = int(sqrt(options.Npoints_2Dscan))
    inchain.Draw("kappa_t:kappa_lambda","quantileExpected == -1","P same")
    best_fit = r.TGraph(r.gROOT.FindObject("Graph"))
    best_fit.SetName("Bestfit")
    best_fit.SetTitle("Bestfit")
    best_fit.SetMarkerSize(1)
    best_fit.SetMarkerStyle(33)
    inchain.Draw("2*deltaNLL:kappa_t:kappa_lambda>>likelihood_vs_kl_kt(%i,%f,%f,%i,%f,%f)"%(Npoints,klmin,klmax,Npoints,ktmin,ktmax),
                 "2*deltaNLL<50",
                 "prof colz goff")
    inchain.Draw("2*deltaNLL:kappa_t:kappa_lambda>>likelihood_vs_kl_kt_contour(%i,%f,%f,%i,%f,%f)"%(Npoints,klmin,klmax,Npoints,ktmin,ktmax),
                 "",
                 "prof colz goff")
    
    likelihood_vs_kl_kt = r.TProfile2D(r.gDirectory.Get("likelihood_vs_kl_kt"))
    fixemptybins(likelihood_vs_kl_kt)
    likelihood_vs_kl_kt.GetXaxis().SetTitle("k_{#lambda}")
    likelihood_vs_kl_kt.GetYaxis().SetTitle("k_{t}")
    likelihood_vs_kl_kt.GetZaxis().SetTitle("-2#Deltaln(L)'")
    likelihood_vs_kl_kt.GetZaxis().SetRangeUser(0.,50.)
    ###################################################
    #likelihood_vs_kl_kt.GetYaxis().SetRangeUser(0.65,1.35)
    #likelihood_vs_kl_kt.GetXaxis().SetRangeUser(-10,15)
    ###################################################
    likelihood_vs_kl_kt.Draw("COLZ")

    likelihood_vs_kl_kt_contour = r.TProfile2D(r.gDirectory.Get("likelihood_vs_kl_kt_contour"))
    fixemptybins(likelihood_vs_kl_kt_contour)
    #draw 1 and 2 sigma contours
    contours =  array('d',[2.3,5.99])
    contour68 = array('d',[2.3])
    contour95 = array('d',[5.99])

    outfile.cd()
    best_fit.Write()

    outfile.cd()
    likelihood_vs_kl_kt_contour68 = r.TProfile2D(likelihood_vs_kl_kt_contour)
    likelihood_vs_kl_kt_contour68.SetName("contour68")
    likelihood_vs_kl_kt_contour68.SetTitle("contour68")
    likelihood_vs_kl_kt_contour68.SetContour(1,contour68);
    likelihood_vs_kl_kt_contour68.SetLineWidth(2);
    likelihood_vs_kl_kt_contour68.SetLineStyle(1);
    likelihood_vs_kl_kt_contour68.Write()

    outfile.cd()
    likelihood_vs_kl_kt_contour95 = r.TProfile2D(likelihood_vs_kl_kt_contour)
    likelihood_vs_kl_kt_contour95.SetName("contour95")
    likelihood_vs_kl_kt_contour95.SetTitle("contour95")
    likelihood_vs_kl_kt_contour95.SetContour(1,contour95);
    likelihood_vs_kl_kt_contour95.SetLineWidth(2);
    likelihood_vs_kl_kt_contour95.SetLineStyle(2);
    likelihood_vs_kl_kt_contour95.Write()

    likelihood_vs_kl_kt_contour.SetContour(2,contours);
    likelihood_vs_kl_kt_contour.SetLineColor(2);
    likelihood_vs_kl_kt_contour.SetLineWidth(2);
    likelihood_vs_kl_kt_contour.Draw("cont3 same");
    best_fit.Draw("p same")
    et = extra_texts()
    for t in et: t.Draw()

    box = r.TBox()
    box.SetFillStyle(1001)
    box.SetFillColor(16)
    box.SetLineColor(16)
    #box.DrawBox(klmin,ktmin,klmax,-1.35) 
    #box.DrawBox(klmin,-0.65,klmax,0.65)
    #box.DrawBox(klmin,1.35,klmax,ktmax)  
    #box.DrawBox(klmin+0.1,-2.1,klmax-0.1,-1.35) 
    #box.DrawBox(klmin+0.1,-0.65,klmax-0.1,0.65)
    #box.DrawBox(klmin+0.1,1.35,klmax-0.1,2.1)  

    c0.Print("%s/likelihood_2Dscan.pdf"%options.outdir)
    c0.Print("%s/likelihood_2Dscan.png"%options.outdir)
    c0.SaveAs("%s/likelihood_2Dscan_canvas.root"%options.outdir)
    outfile.cd()
    likelihood_vs_kl_kt.Write()

    '''
    #kl likelihood scan with kt fixed to a certain value
    ktrange = 0.5*(ktmax-ktmin)/Npoints
    for ikt in range(0,Npoints):
        ikt+=5
        kt = options.ktmin + (0.5+ikt)*(ktmax-ktmin)/Npoints
        c1 = r.TCanvas()
        inchain.Draw("2*deltaNLL:kappa_lambda>>likelihood_vs_kl_ikt%i(%i,%f,%f)"%(ikt,Npoints,klmin,klmax),
                     "kappa_t>%f && kappa_t<%f && 2*deltaNLL<100"%(kt-ktrange,kt+ktrange),
                     "prof goff")
        likelihood_vs_kl = r.TProfile(r.gDirectory.Get("likelihood_vs_kl_ikt%i"%ikt))
        likelihood_vs_kl.GetXaxis().SetTitle("k_{#lambda}")
        likelihood_vs_kl.GetYaxis().SetTitle("-2#Deltaln(L)'")
        likelihood_vs_kl.Draw("histl")
        xmin=likelihood_vs_kl.GetXaxis().GetXmin()
        xmax=likelihood_vs_kl.GetXaxis().GetXmax()
        ### lines
        sigmas = [1,1.96]
        CL = [68,95]
        lines = []
        for s in sigmas:
            l = r.TLine(xmin, s*s, xmax, s*s)
            l.SetLineStyle(7)
            l.SetLineWidth(1)
            l.SetLineColor(r.kGray+2)
            lines.append(l)
        for l in lines: l.Draw()
        c1.Print("%s/likelihood_klscan_kt%.3f.png"%(options.outdir,kt))

    c2 = r.TCanvas()
    inchain.Draw("2*deltaNLL:kappa_t>>likelihood_vs_kt(%i,%f,%f)"%(Npoints,ktmin,ktmax),
                 "kappa_lambda>0.95 && kappa_lambda<1.05 && 2*deltaNLL<20",
                 "prof goff")
    likelihood_vs_kt = r.TProfile(r.gDirectory.Get("likelihood_vs_kt"))
    likelihood_vs_kt.GetXaxis().SetTitle("k_{t}")
    likelihood_vs_kt.GetYaxis().SetTitle("-2#Deltaln(L)'")
    likelihood_vs_kt.Draw("histl")
    xmin=likelihood_vs_kt.GetXaxis().GetXmin()
    xmax=likelihood_vs_kt.GetXaxis().GetXmax()
    ### lines
    sigmas = [1,1.96]
    CL = [68,95]
    lines = []
    for s in sigmas:
        l = r.TLine(xmin, s*s, xmax, s*s)
        l.SetLineStyle(7)
        l.SetLineWidth(1)
        l.SetLineColor(r.kGray+2)
        lines.append(l)
    for l in lines: l.Draw()
    c2.Print("%s/likelihood_ktscan.pdf"%options.outdir)
    c2.Print("%s/likelihood_ktscan.png"%options.outdir)
    '''

if options.infilename_klscan != "":
    c1 = r.TCanvas()
    inchain = r.TChain("limit","limit")
    for infilename in options.infilename_klscan.split(","):
        inchain.Add(infilename)
    Npoints=inchain.Draw("2*deltaNLL:kappa_lambda","2*deltaNLL<6 && deltaNLL>0","goff")
    print "draw %i points"%Npoints

    deltaNLL_values=inchain.GetV1()
    deltaNLL_values.SetSize(Npoints)
    a_deltaNLL = array('d',deltaNLL_values)  

    kl_values=inchain.GetV2()
    kl_values.SetSize(Npoints)  
    a_kl = array('d',kl_values) 
 
    likelihood_vs_kl = r.TGraph(Npoints,a_kl,a_deltaNLL)
    likelihood_vs_kl.SetName("klscan")
    likelihood_vs_kl.SetTitle("klscan")
    likelihood_vs_kl.Sort()
    outfile.cd()
    likelihood_vs_kl.Write()
    likelihood_vs_kl.Draw("AL")
    xmin=likelihood_vs_kl.GetXaxis().GetXmin()
    xmax=likelihood_vs_kl.GetXaxis().GetXmax()
    ### lines
    sigmas = [1,1.96]
    CL = [68,95]
    lines = []
    for s in sigmas:
        l = r.TLine(xmin, s*s, xmax, s*s)
        l.SetLineStyle(7)
        l.SetLineWidth(1)
        l.SetLineColor(r.kGray+2)
        lines.append(l)
    for l in lines: l.Draw()
    c1.Print("%s/likelihood_klscan.pdf"%options.outdir)
    c1.Print("%s/likelihood_klscan.png"%options.outdir)

if options.infilename_ktscan != "":
    c1 = r.TCanvas()
    inchain = r.TChain("limit","limit")
    for infilename in options.infilename_ktscan.split(","):
        inchain.Add(infilename)
    Npoints=inchain.Draw("2*deltaNLL:kappa_t","2*deltaNLL<6 && deltaNLL>0","goff")
    deltaNLL_values=inchain.GetV1()
    kt_values=inchain.GetV2()

    deltaNLL_values=inchain.GetV1()
    deltaNLL_values.SetSize(Npoints)
    a_deltaNLL = array('d',deltaNLL_values)  

    kt_values=inchain.GetV2()
    kt_values.SetSize(Npoints)  
    a_kt = array('d',kt_values) 

    likelihood_vs_kt = r.TGraph(Npoints,a_kt,a_deltaNLL)
    likelihood_vs_kt.SetName("ktscan")
    likelihood_vs_kt.SetTitle("ktscan")
    likelihood_vs_kt.Sort()
    outfile.cd()
    likelihood_vs_kt.Write()
    likelihood_vs_kt.Draw("AL")
    xmin=likelihood_vs_kt.GetXaxis().GetXmin()
    xmax=likelihood_vs_kt.GetXaxis().GetXmax()
    ### lines
    sigmas = [1,1.96]
    CL = [68,95]
    lines = []
    for s in sigmas:
        l = r.TLine(xmin, s*s, xmax, s*s)
        l.SetLineStyle(7)
        l.SetLineWidth(1)
        l.SetLineColor(r.kGray+2)
        lines.append(l)
    for l in lines: l.Draw()
    c1.Print("%s/likelihood_ktscan.pdf"%options.outdir)
    c1.Print("%s/likelihood_ktscan.png"%options.outdir)

outfile.Close()
