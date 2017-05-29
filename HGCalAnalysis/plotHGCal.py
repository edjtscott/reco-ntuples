import ROOT as r
import os
from math import exp,cos,sin,sqrt


from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inDir",default="/afs/cern.ch/work/e/escott/public/HGCStudies/Ntuples/")
parser.add_option("-f","--fullFilePath",default="")
parser.add_option("-o","--outDir",default="HGCPlots/")
parser.add_option("-w","--webDir",default="")
parser.add_option("-p","--particleType",default="Photon")
parser.add_option("-m","--ptVal",default="35")
parser.add_option("-e","--ext",default="UpdatedNtuple_20mm")
parser.add_option("-n","--nClus",type="int",default=3,help="min number of 2D clusters for a multicluster to be included in sums. Default is three")
parser.add_option("-s","--superClus",type="int",default=30,help="set baseline superclustering radius. Default is thirty")
parser.add_option("-c","--doConverted",type="int",default=0,help="for photons: set to one to only plot unconverted, two for unconverted. Default is zero (all)")
parser.add_option("-b","--doBasics",dest="doBasics",default=False,action="store_true",help="make basic plots with no cuts applied")
parser.add_option("-v","--verbosity",type="int",default=0)
(opts,args) = parser.parse_args()


class BasicPlots:
  """Class to make the simple diagnostic plots"""
  types = ["genpart","cluster2d","multiclus","rechit"]
  quantities = ["eta","phi","pt","energy","z"]

  def __init__(self):
    self.hists = {}
    for atype in self.types:
      for quantity in self.quantities:
        histName = atype+'_'+quantity
        if   'eta'    in histName: self.hists[histName] = r.TH1F(histName,histName,100,-5.,5.) 
        elif 'phi'    in histName: self.hists[histName] = r.TH1F(histName,histName,64,-3.2,3.2)
        elif 'pt'     in histName: self.hists[histName] = r.TH1F(histName,histName,100,0.,50.)
        elif 'energy' in histName: self.hists[histName] = r.TH1F(histName,histName,100,0.,500.)
        elif 'z'      in histName and 'multiclus' in histName: self.hists[histName] = r.TH1F(histName,histName,100,-500.,500.)
   
  def fillHists(self,tree):
    print "about to loop through tree with %d entries"%tree.GetEntries()
    for i in range(tree.GetEntries()):
      if i%100==0: print "processing entry %d"%i
      tree.GetEntry(i)
      for histName in self.hists.keys():
        collection = getattr(tree,histName)
        for val in collection:
          self.hists[histName].Fill(val)

  def drawHists(self,outDir):
    canv = r.TCanvas("c","c")
    canv.SetTicks(1,1)
    for histName,hist in self.hists.iteritems():
      hist.Draw()
      canv.Print(outDir+histName+".pdf")
      canv.Print(outDir+histName+".png")


def printOpts(fileName):
  print "Running HGC plotting code for a %s with pT %s"%(opts.particleType,opts.ptVal)
  print "Infile:     %s"%fileName
  print "Outdir:     %s"%opts.outDir
  print "Conversion: %s"%opts.doConverted

def makeOutDir(thedir):
  if not thedir[-1]=='/':
    thedir+='/'
  if not os.path.exists(thedir): 
    os.makedirs(thedir)
  if not os.path.isfile(thedir+'index.php'): 
    os.system('cp /afs/cern.ch/user/e/escott/www/HGCclustering/FixClusteringTest/index.php %s'%thedir)

def initGenHists():
  pass

def initMultiHists():
  hists = {}
 
  hists['bestEnFrac'] = r.TH1F('hMulti_bestEnFrac','hMulti_bestEnFrac',100,0.,2.)
  enBins = 50
  enLow  = 0.
  enHigh = 0.
  for i in range(1,opts.superClus+1):
    hists['superEnFrac_EEfirst%d'%i] = r.TH1F('hMulti_superEnFrac_EEfirst%d'%i,'hMulti_superEnFrac_EEfirst%d'%i,enBins,enLow,enHigh)
    hists['superEnFrac_EEsecond%d'%i] = r.TH1F('hMulti_superEnFrac_EEsecond%d'%i,'hMulti_superEnFrac_EEsecond%d'%i,enBins,enLow,enHigh)
    hists['superEnFrac_FH%d'%i] = r.TH1F('hMulti_superEnFrac_FH%d'%i,'hMulti_superEnFrac_FH%d'%i,enBins,enLow,enHigh)
    hists['superEnFrac_BH%d'%i] = r.TH1F('hMulti_superEnFrac_BH%d'%i,'hMulti_superEnFrac_BH%d'%i,enBins,enLow,enHigh)

  return hists

def setXTitle(key,hist):
  names = {}
  names['pt'] = 'p_{T}'
  names['eta'] = '#eta'
  names['en'] = 'E'
  for val in names.keys():
    if val in key.lower():
      if 'frac' in key.lower():
        hist.GetXaxis().SetTitle('%s / gen %s'%(names[val],names[val]))
      else:
         hist.GetXaxis().SetTitle(names[val])
  

def printHists(canv,hists,outdir):
  for key,hist in hists.iteritems():
    canv.cd() 
    setXTitle(key,hist)
    hist.Draw()
    canv.Print(outdir+hist.GetName()+".pdf")
    canv.Print(outdir+hist.GetName()+".png")

def fitHists(canv,hists,outdir):
  theGraphs = {'EEfirst':[],'EEsecond':[],'FH':[],'BH':[]}
  for section in theGraphs.keys():
    theGraphs[section].append(r.TGraph(opts.superClus))
    theGraphs[section][0].SetName('gSuper_%s_mean'%section)
    theGraphs[section].append(r.TGraph(opts.superClus))
    theGraphs[section][1].SetName('gSuper_%s_sigma'%section)
    theGraphs[section].append(r.TGraph(opts.superClus))
    theGraphs[section][2].SetName('gSuper_%s_res'%section)
  for key,hist in hists.iteritems():
    canv.cd() 
    setXTitle(key,hist)
    hist.Draw()
    #fit here
    fit = r.TF1("fit","gaus")
    hist.Fit(fit)
    mean  = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    for section in theGraphs.keys():
      if 'super' in key and section in key:
        point = int(key.split(section)[1])
        theGraphs[section][0].SetPoint( point-1, point, mean )
        theGraphs[section][1].SetPoint( point-1, point, sigma )
        if mean > 0.:
          theGraphs[section][2].SetPoint( point-1, point, sigma / mean )
        else:
          theGraphs[section][2].SetPoint( point-1, point, 1. )
  for section,graphs in theGraphs.iteritems():
    for graph in graphs:
      canv.cd() 
      canv.Clear() 
      graph.GetXaxis().SetName('Radius / cm')
      graph.GetYaxis().SetName(graph.GetName().split('_')[2])
      graph.Draw()
      canv.Print(outdir+graph.GetName()+".pdf")
      canv.Print(outdir+graph.GetName()+".png")

def etaPhiZtoX(eta,phi,z):
  t = exp( -1. * eta );
  x = z * 2. * t * cos(phi) / ( 1. - t*t );
  return x;

def etaPhiZtoY(eta,phi,z):
  t = exp( -1. * eta );
  y = z * 2. * t * sin(phi) / ( 1. - t*t );
  return y;

def deltaX(x1,x2,y1,y2):
  dX = x1 - x2;
  dY = y1 - y2;
  return sqrt( dX*dX + dY*dY );

def etasPhisZsToDeltaX(eta1,eta2,phi1,phi2,z1,z2):
  t1 = exp( -1. * eta1 );
  x1 = z1 * 2. * t1 * cos(phi1) / ( 1. - t1*t1 );
  y1 = z1 * 2. * t1 * sin(phi1) / ( 1. - t1*t1 );
  t2 = exp( -1 * eta2 );
  x2 = z2 * 2. * t2 * cos(phi2) / ( 1. - t2*t2 );
  y2 = z2 * 2. * t2 * sin(phi2) / ( 1. - t2*t2 );
  return deltaX(x1,x2,y1,y2);


def main():
  r.gROOT.SetBatch(True)

  #get file and tree
  if opts.inDir[-1]!='/':
    opts.inDir+='/'
  fileName = opts.inDir+'partGun_'+opts.particleType+'_Pt'+opts.ptVal+'_'+opts.ext+'.root'
  if opts.fullFilePath:
    fileName = opts.fullFilePath
  printOpts(fileName)
  assert(os.path.isfile(fileName) and "Error: file does not exist")
  theFile = r.TFile(fileName)
  if(opts.verbosity>-1): print "Successfully opened file",fileName
  theTree = theFile.Get("ana/hgc")
  assert(theTree and "Error: tree does not exist")
  if(opts.verbosity>-1): print "Successfully got the tree"
  if opts.outDir[-1]!='/':
    opts.outDir+='/'
  makeOutDir(opts.outDir)

  #do the normal set of basic plots
  if opts.doBasics:
    plots = BasicPlots()
    plots.fillHists(theTree)
    plots.drawHists(opts.outDir)

  #setup hists for each type of object
  genHists   = initGenHists()
  multiHists = initMultiHists()

  #loop over entries in tree and actully do things
  for iEntry in range(theTree.GetEntries()):
    #setup gen values
    if iEntry%100==0 or iEntry==0: print "Processing entry %d"%iEntry
    theTree.GetEntry(iEntry)
    genEtas = getattr(theTree,"genpart_eta")
    genPts = getattr(theTree,"genpart_pt")
    genEnergies = getattr(theTree,"genpart_energy")
    genDvzs = getattr(theTree,"genpart_dvz")

    #loop over gen particles (assuming one per endcap)
    for iGen in range(genEtas.size()):
      #apply selections
      genPt  = genPts[iGen]
      if genPt < float(opts.ptVal): #FIXME: why are these here?!
        continue
      genEta = genEtas[iGen]
      if abs(genEta) < 1.7 or abs(genEta) > 2.7: 
        continue
      genDvz = genDvzs[iGen]
      if opts.doConverted==1 and abs(genDvz) < 300.: #FIXME is 300 exactly right?
        continue
      if opts.doConverted==2 and abs(genDvz) > 300.:
        continue
      genEnergy  = genEnergies[iGen]
      
      #setup multi values
      multiEtas = getattr(theTree,"multiclus_eta")
      multiPhis = getattr(theTree,"multiclus_phi")
      multiPts = getattr(theTree,"multiclus_pt")
      multiEnergies = getattr(theTree,"multiclus_energy")
      multiZees = getattr(theTree,"multiclus_z")
      multiN2Ds = getattr(theTree,"multiclus_cluster2d")
      bestMultiIndex = -9999
      bestMultiEnergy = -9999.
      selectedMultiIndices = []

      #loop over multis
      for iMulti in range(multiEtas.size()):
        #check same endcap
        multiEta = multiEtas[iMulti]
        if( multiEta*genEta < 0. ):
          continue

        #analysis
        multiEnergy = multiEnergies[iMulti]
        if multiEnergy > bestMultiEnergy:
          bestMultiIndex = iMulti
          bestMultiEnergy = multiEnergy
        multiN2D = len(multiN2Ds[iMulti])
        if multiN2D >= opts.nClus:
          #print "multi with %d 2D being added"%multiN2D
          selectedMultiIndices.append(iMulti)

      #end multi loop
      multiHists['bestEnFrac'].Fill(bestMultiEnergy/genEnergy)
      bestMultiZ = multiZees[bestMultiIndex]
      bestMultiPhi = multiPhis[bestMultiIndex]
      bestMultiEta = multiEtas[bestMultiIndex]

      #setup superclustering
      numSuperPoints = opts.superClus
      scanEnergies = {'EEfirst':[-9999.],'EEsecond':[-9999.],'FH':[-9999.],'BH':[-9999.]}
      for key in scanEnergies.keys():
        for i in range(numSuperPoints):
          scanEnergies[key].append(bestMultiEnergy)
     
      #loop over selected multis again for superclustering
      #idea is to add together the multis within a default distance of 30cm, 
      #then tightening that threshold on each subdetector individually
      for iSel in selectedMultiIndices:
        multiEta = multiEtas[iSel]
        if multiEta*genEta < 0.:
          continue
        if iSel == bestMultiIndex:
          continue
        multiZ = multiZees[iSel]
        multiPhi = multiPhis[iSel]
        deltaX = etasPhisZsToDeltaX(multiEta,bestMultiEta,multiPhi,bestMultiPhi,multiZ,bestMultiZ) #FIXME check this works
        #print "multiEta,multiPhi,multiZ",multiEta,multiPhi,multiZ
        #print "bestEta,bestPhi,bestZ",bestMultiEta,bestMultiPhi,bestMultiZ
        #print "deltaX",deltaX
        multiEnergy = multiEnergies[iSel]
        for ncm in range(1,numSuperPoints+1):
          if abs(multiZ) < 335.:
            if deltaX < float(ncm):
              #print deltaX,float(ncm)
              #print "debug A"
              #print deltaX,multiZ,multiEnergy
              scanEnergies['EEfirst'][ncm] += multiEnergy
              #print bestMultiEnergy,scanEnergies['EEfirst'][ncm]
            if deltaX < float(numSuperPoints):
              #print deltaX,float(ncm)
              scanEnergies['EEsecond'][ncm] += multiEnergy
              scanEnergies['FH'][ncm] += multiEnergy
              scanEnergies['BH'][ncm] += multiEnergy
          elif abs(multiZ) < 384.:
            if deltaX < float(ncm):
              scanEnergies['EEsecond'][ncm] += multiEnergy
            if deltaX < float(numSuperPoints):
              scanEnergies['EEfirst'][ncm] += multiEnergy
              scanEnergies['FH'][ncm] += multiEnergy
              scanEnergies['BH'][ncm] += multiEnergy
          elif abs(multiZ) < 425.:
            if deltaX < float(ncm):
              scanEnergies['FH'][ncm] += multiEnergy
            if deltaX < float(numSuperPoints):
              scanEnergies['EEfirst'][ncm] += multiEnergy
              scanEnergies['EEsecond'][ncm] += multiEnergy
              scanEnergies['BH'][ncm] += multiEnergy
          elif abs(multiZ) < 999.:
            if deltaX < float(ncm):
              scanEnergies['BH'][ncm] += multiEnergy
            if deltaX < float(numSuperPoints):
              scanEnergies['EEfirst'][ncm] += multiEnergy
              scanEnergies['EEsecond'][ncm] += multiEnergy
              scanEnergies['FH'][ncm] += multiEnergy

      #end superclustering loop
      for icm in range(1,numSuperPoints+1):
        #print scanEnergies['EEfirst'][icm]
        #print genEnergy
        multiHists['superEnFrac_EEfirst%d'%icm].Fill( scanEnergies['EEfirst'][icm] / genEnergy )
        multiHists['superEnFrac_EEsecond%d'%icm].Fill( scanEnergies['EEsecond'][icm] / genEnergy )
        multiHists['superEnFrac_FH%d'%icm].Fill( scanEnergies['FH'][icm] / genEnergy )
        multiHists['superEnFrac_BH%d'%icm].Fill( scanEnergies['BH'][icm] / genEnergy )

    #end gen loop

  #end loop over entries

  canv = r.TCanvas('canv','canv')
  printHists(canv,multiHists,opts.outDir)
  fitHists(canv,multiHists,opts.outDir)
  
  #copy plots across
  if opts.webDir:
    if opts.webDir[-1]!='/':
      opts.webDir+='/'
    makeOutDir(opts.webDir)
    if opts.doBasics:
      opts.webDir+="Basics/"
    elif opts.doConverted==0:
      opts.webDir+="All/"
    elif opts.doConverted==1:
      opts.webDir+="Unconverted/"
    elif opts.doConverted==2:
      opts.webDir+="Converted/"
    makeOutDir(opts.webDir)
    if opts.verbosity>-1: print "moving plots to",opts.webDir
    os.system("mv %s* %s"%(opts.outDir,opts.webDir))

  #end main
  return 0


if __name__ == '__main__':
  main()
