import ROOT
import math
import sys

if int(sys.argv[1]) == 0:
    fIn=ROOT.TFile('MinBias3D_DATA.root')
    fOut=ROOT.TFile('RadLength_DATA.root', 'RECREATE')
else:
    fIn=ROOT.TFile('MinBias3D_MC.root')
    fOut=ROOT.TFile('RadLength_MC.root', 'RECREATE')

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# loop over histograms
for key in fIn.GetListOfKeys():
    print key.GetName()
    if not 'minus' in key.GetName(): continue

    # take only the minus histograms fro the list
    th3=ROOT.TH3
    th3=fIn.Get(key.GetName())
    triplet=[p for p in th3.GetName().split('_') if not 'minus' in p][0]
    print triplet

    #find the paired plus
    th3pl=ROOT.TH3
    th3pl=fIn.Get('sag3Dplus_'+triplet)

    print th3.GetName(), th3pl.GetName()
    fOut.mkdir(triplet)
    # define TH1 to fill
    slope=ROOT.TH1D("slope"+triplet, "slope"+triplet, 41, -4, 4)
    radlength=ROOT.TH1D("radlength"+triplet, "radlength"+triplet, 41, -4, 4)

    fOut.cd(triplet)
    #loop over the x and y bins
    for i in range(1, th3.GetNbinsX()+1):
        rms=ROOT.TH1D("{i}".format(i=i), "{i}".format(i=i), 10, 0.5, 2.25);
        for j in range(1, th3.GetNbinsY()+1):

            th1=ROOT.TH1
            th1pl=ROOT.TH1

            th1=th3.ProjectionZ("binMinus_{i}_{j}".format(i=i,j=j), i, i, j, j)
            th1pl=th3pl.ProjectionZ("binPlus_{i}_{j}".format(i=i,j=j), i, i, j, j)

            th1.Write()
            th1pl.Write()

            mean=(1./2*(th1.GetRMS()+th1pl.GetRMS()))
            mean_err=mean*math.sqrt(th1.GetRMSError()*th1.GetRMSError()+th1pl.GetRMSError()*th1pl.GetRMSError())

            rms.SetBinContent(j, mean*mean)
            rms.SetBinError(j, mean_err)

        if(rms.GetEntries()<7): continue
        fit=ROOT.TF1("fit", "pol1");
        rms.Fit(fit, "Q");
        rms.Write()
        
        slope.SetBinContent(i,fit.GetParameter(1))
        slope.SetBinError(i, fit.GetParError(1))

        radlength.SetBinContent(i,fit.GetParameter(0)/0.013/0.013)
        radlength.SetBinError(i, fit.GetParError(0)/0.013/0.013)

    fOut.cd()
    slope.Write()
    radlength.GetXaxis().SetTitle("local #eta")
    radlength.GetYaxis().SetTitle("x/X_{0}")
    radlength.Write()
