import ROOT as r
import os
from math import exp,cos,sin,sqrt,atan


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
parser.add_option("--endText",dest="endText",default="\n",help="text to print at end")
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

def makeOutDir(thedir,index=True):
  if not thedir[-1]=='/':
    thedir+='/'
  if not os.path.exists(thedir): 
    os.makedirs(thedir)
  if not os.path.isfile(thedir+'index.php'): 
    if index: os.system('cp /afs/cern.ch/user/e/escott/www/HGCclustering/FixClusteringTest/index.php %s'%thedir)

def initGenHists():
  pass

def initMultiHists():
  hists = {}
 
  #hists['bestEnFrac'] = r.TH1F('hMulti_bestEnFrac','hMulti_bestEnFrac',100,0.,2.)
  #hists['bestPfEnFrac'] = r.TH1F('hMulti_bestEnFrac','hMulti_bestEnFrac',100,0.,2.)
  hists['bestEnFrac'] = r.TH1F('hMulti_bestEnFrac','hMulti_bestEnFrac',100,1.,0.)
  hists['bestPfEnFrac'] = r.TH1F('hMulti_bestPfEnFrac','hMulti_bestPfEnFrac',100,1.,0.)
  hists['totalEnFrac'] = r.TH1F('hMulti_totalEnFrac','hMulti_totalEnFrac',100,1.,0.)
  hists['suggestedEnFrac'] = r.TH1F('hMulti_suggestedEnFrac','hMulti_suggestedEnFrac',100,1.,0.)
  hists['layersInBest'] = r.TH1F('hMulti_layersInBest','hMulti_layersInBest',53,-0.5,52.5) #TODO: fill these
  hists['layersInBestWeighted'] = r.TH1F('hMulti_layersInBestWeighted','hMulti_layersInBestWeighted',53,-0.5,52.5)
  hists['drGenBestAtFace'] = r.TH1F('hMulti_drGenBestAtFace','hMulti_drGenBestAtFace',100,0.,5.)
  hists['drGenBestAtBestZ'] = r.TH1F('hMulti_drGenBestAtBestZ','hMulti_drGenBestAtBestZ',100,0.,5.)
  hists['drGenBestAtFaceCorr'] = r.TH1F('hMulti_drGenBestAtFaceCorr','hMulti_drGenBestAtFaceCorr',100,0.,5.)
  hists['drGenBestAtBestZCorr'] = r.TH1F('hMulti_drGenBestAtBestZCorr','hMulti_drGenBestAtBestZCorr',100,0.,5.)
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
  

#def writeHists(hists,outdir,mode="UPDATE"):
def writeHists(hists,outdir,mode="RECREATE"):
  makeOutDir(outdir,False)
  convType = ['All','Unconverted','Converted']
  fileName = outdir+'/partGun_'+opts.particleType+'_Pt'+opts.ptVal+'_'+opts.ext+'_'+convType[opts.doConverted]+'.root'
  if opts.fullFilePath:
    fileName = outdir+'/'+opts.fullFilePath.split('/')[-1].replace('.root','_'+convType[opts.doConverted]+'.root')
  outFile = r.TFile(fileName,mode)
  for key,hist in hists.iteritems():
    hist.Write()
  outFile.Close()

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
    theGraphs[section][0].SetTitle('gSuper_%s_mean'%section)
    theGraphs[section].append(r.TGraph(opts.superClus))
    theGraphs[section][1].SetName('gSuper_%s_sigma'%section)
    theGraphs[section][1].SetTitle('gSuper_%s_sigma'%section)
    theGraphs[section].append(r.TGraph(opts.superClus))
    theGraphs[section][2].SetName('gSuper_%s_res'%section)
    theGraphs[section][2].SetTitle('gSuper_%s_res'%section)
  for key,hist in hists.iteritems():
    canv.cd() 
    setXTitle(key,hist)
    hist.Draw()
    #fit here
    fit = r.TF1("fit","gaus")
    hist.Fit(fit)
    mean  = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    if 'total' in key or 'Pf' in key or 'suggested' in key:
      opts.endText += "For hist %s the mean is %1.2f and the width is %1.2f\n"%(key,mean,sigma)
      opts.endText += "hence resolution is %1.2f\n\n"%(sigma/mean)
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
      graphType = graph.GetName().split('_')[2]
      graph.GetXaxis().SetTitle('Radius / cm')
      graph.GetYaxis().SetTitle(graphType)
      if graphType == 'mean':
        graph.GetYaxis().SetRangeUser(0.6,0.7)
        graph.SetLineColor(r.kBlue+2)
      elif graphType == 'sigma':
        graph.GetYaxis().SetRangeUser(0.13,0.17)
      elif graphType == 'res':
        graph.GetYaxis().SetRangeUser(0.18,0.28)
        graph.SetLineColor(r.kGreen+2)
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

def cot(val):
  return cos(val)/sin(val)


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
    genPhis = getattr(theTree,"genpart_phi")
    genReachedEEs= getattr(theTree,"genpart_reachedEE")

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
      genReachedEE = genReachedEEs[iGen]
      if opts.doConverted==1 and not genReachedEE:
        continue
      if opts.doConverted==2 and genReachedEE:
        continue
      genEnergy = genEnergies[iGen]
      genPhi = genPhis[iGen]
      genDz = genDvzs[iGen]

      #setup track values - not working atm
      #trackEtas = getattr(theTree,"track_eta")
      #print len(genEtas)
      #print len(trackEtas)
      #exit("testing")

      #setup pfcluster values
      pfEtas = getattr(theTree,"pfcluster_eta")
      pfPhis = getattr(theTree,"pfcluster_phi")
      pfPts = getattr(theTree,"pfcluster_pt")
      pfEnergies = getattr(theTree,"pfcluster_energy")
      bestPfIndex = -9999
      bestPfEnergy = -9999.
      for iPf in range(pfEnergies.size()):
        pfEta = pfEtas[iPf]
        if( pfEta*genEta < 0. ):
          continue
        pfEnergy = pfEnergies[iPf]
        if pfEnergy > bestPfEnergy:
          bestPfEnergy = pfEnergy
          bestPfIndex = iPf
      multiHists['bestPfEnFrac'].Fill(bestPfEnergy/genEnergy)
      
      #setup multi values
      multiEtas = getattr(theTree,"multiclus_eta")
      multiPhis = getattr(theTree,"multiclus_phi")
      multiPts = getattr(theTree,"multiclus_pt")
      multiEnergies = getattr(theTree,"multiclus_energy")
      multiZees = getattr(theTree,"multiclus_z")
      multi2Ds = getattr(theTree,"multiclus_cluster2d")
      bestMultiIndex = -9999
      bestMultiEnergy = -9999.
      totalMultiEnergy = 0.
      selectedMultiIndices = []

      #loop over multis
      for iMulti in range(multiEtas.size()):
        #check same endcap
        multiEta = multiEtas[iMulti]
        if( multiEta*genEta < 0. ):
          continue

        #analysis
        multiEnergy = multiEnergies[iMulti]
        totalMultiEnergy += multiEnergy
        if multiEnergy > bestMultiEnergy:
          bestMultiIndex = iMulti
          bestMultiEnergy = multiEnergy
        multiN2D = len(multi2Ds[iMulti])
        #multiN2D = multiN2Ds[iMulti]
        if multiN2D >= opts.nClus:
          #print "multi with %d 2D being added"%multiN2D
          selectedMultiIndices.append(iMulti)

      #end multi loop
      multiHists['bestEnFrac'].Fill(bestMultiEnergy/genEnergy)
      multiHists['totalEnFrac'].Fill(totalMultiEnergy/genEnergy)
      bestMultiZ = multiZees[bestMultiIndex]
      bestMultiPhi = multiPhis[bestMultiIndex]
      bestMultiEta = multiEtas[bestMultiIndex]
      multiHists['drGenBestAtFace'].Fill( etasPhisZsToDeltaX(genEta,bestMultiEta,genPhi,bestMultiPhi,320.,320.) )
      multiHists['drGenBestAtBestZ'].Fill( etasPhisZsToDeltaX(genEta,bestMultiEta,genPhi,bestMultiPhi,bestMultiZ,bestMultiZ) )
      #dPhi = q*B*dZ/pZ; prefactor is the conversion of GeV to normal units
      genDeltaPhiFace = 0.2998 * 1. * 3.8 * 0.01*(320.-genDz) * (1./(genPt*cot(2*atan(exp(-1*genEta)))))
      genDeltaPhiBestZ = 0.2998 * 1. * 3.8 * 0.01*(bestMultiZ-genDz) * (1./(genPt*cot(2*atan(exp(-1*genEta)))))
      multiHists['drGenBestAtFaceCorr'].Fill( etasPhisZsToDeltaX(genEta,bestMultiEta,genPhi+genDeltaPhiFace,bestMultiPhi,320.,320.) )
      multiHists['drGenBestAtBestZCorr'].Fill( etasPhisZsToDeltaX(genEta,bestMultiEta,genPhi+genDeltaPhiBestZ,bestMultiPhi,bestMultiZ,bestMultiZ) )

      #setup superclustering
      numSuperPoints = opts.superClus
      scanEnergies = {'EEfirst':[-9999.],'EEsecond':[-9999.],'FH':[-9999.],'BH':[-9999.]}
      for key in scanEnergies.keys():
        for i in range(numSuperPoints):
          scanEnergies[key].append(bestMultiEnergy)
     
      #loop over selected multis again for superclustering
      #idea is to add together the multis within a default distance of 30cm, 
      #then tightening that threshold on each subdetector individually
      suggestedEnergy = bestMultiEnergy
      suggestedRadii = {'EEfirst':2.,'EEsecond':10.,'FH':15.,'BH':15.}
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
            if ncm==1 and deltaX < suggestedRadii['EEfirst']: suggestedEnergy += multiEnergy
          elif abs(multiZ) < 384.:
            if deltaX < float(ncm):
              scanEnergies['EEsecond'][ncm] += multiEnergy
            if deltaX < float(numSuperPoints):
              scanEnergies['EEfirst'][ncm] += multiEnergy
              scanEnergies['FH'][ncm] += multiEnergy
              scanEnergies['BH'][ncm] += multiEnergy
            if ncm==1 and deltaX < suggestedRadii['EEsecond']: suggestedEnergy += multiEnergy
          elif abs(multiZ) < 425.:
            if deltaX < float(ncm):
              scanEnergies['FH'][ncm] += multiEnergy
            if deltaX < float(numSuperPoints):
              scanEnergies['EEfirst'][ncm] += multiEnergy
              scanEnergies['EEsecond'][ncm] += multiEnergy
              scanEnergies['BH'][ncm] += multiEnergy
            if ncm==1 and deltaX < suggestedRadii['FH']: suggestedEnergy += multiEnergy
          elif abs(multiZ) < 999.:
            if deltaX < float(ncm):
              scanEnergies['BH'][ncm] += multiEnergy
            if deltaX < float(numSuperPoints):
              scanEnergies['EEfirst'][ncm] += multiEnergy
              scanEnergies['EEsecond'][ncm] += multiEnergy
              scanEnergies['FH'][ncm] += multiEnergy
            if ncm==1 and deltaX < suggestedRadii['BH']: suggestedEnergy += multiEnergy

      #end superclustering loop
      for icm in range(1,numSuperPoints+1):
        #print scanEnergies['EEfirst'][icm]
        #print genEnergy
        multiHists['superEnFrac_EEfirst%d'%icm].Fill( scanEnergies['EEfirst'][icm] / genEnergy )
        multiHists['superEnFrac_EEsecond%d'%icm].Fill( scanEnergies['EEsecond'][icm] / genEnergy )
        multiHists['superEnFrac_FH%d'%icm].Fill( scanEnergies['FH'][icm] / genEnergy )
        multiHists['superEnFrac_BH%d'%icm].Fill( scanEnergies['BH'][icm] / genEnergy )
      multiHists['suggestedEnFrac'].Fill( suggestedEnergy / genEnergy )

    #end gen loop

  #end loop over entries

  canv = r.TCanvas('canv','canv')
  printHists(canv,multiHists,opts.outDir)
  fitHists(canv,multiHists,opts.outDir)
  writeHists(multiHists,'Output/')
  
  #copy plots across
  if opts.webDir:
    if opts.webDir[-1]!='/':
      opts.webDir+='/'
    makeOutDir(opts.webDir)
    convType = ['All/','Unconverted/','Converted/']
    if opts.doBasics:
      opts.webDir+="Basics/"
    else:
      opts.webDir+=convType[opts.doConverted]
    makeOutDir(opts.webDir)
    if opts.verbosity>-1: print "moving plots to",opts.webDir
    os.system("rm %s*"%opts.webDir)
    os.system("mv %s* %s"%(opts.outDir,opts.webDir))

  #end main
  print opts.endText
  return 0


if __name__ == '__main__':
  main()
