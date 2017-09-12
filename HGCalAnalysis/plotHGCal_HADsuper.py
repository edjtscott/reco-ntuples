import ROOT as r
import os
import sys
from math import exp,cos,sin,tan,sqrt,atan,log


from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inDir",default="/afs/cern.ch/work/e/escott/public/HGCStudies/Ntuples/")
parser.add_option("-f","--fullFilePath",default="")
parser.add_option("-o","--outDir",default="HGCPlots/")
parser.add_option("-w","--webDir",default="")
parser.add_option("-p","--particleType",default="Photon")
parser.add_option("-m","--ptVal",default="35")
parser.add_option("-e","--ext",default="UpdatedNtuple_20mm")
parser.add_option("--enRadii",default="2,4,6,8")
parser.add_option("--scRadii",default="4,6,8,10")
parser.add_option("--etaWindow",type="float",default=0.05,help="if less than one, just the eta window. If greater than one, the number of cm the window should correspond to")
parser.add_option("-n","--nClus",type="int",default=3,help="min number of 2D clusters for a multicluster to be included in sums. Default is three")
parser.add_option("-s","--superClus",type="int",default=30,help="set number of points to scan for superclustering. Default is thirty")
parser.add_option("-c","--doConverted",type="int",default=0,help="for photons: set to one to only plot unconverted, two for unconverted. Default is zero (all)")
parser.add_option("-b","--doBasics",dest="doBasics",default=False,action="store_true",help="make basic plots with no cuts applied")
parser.add_option("--endText",dest="endText",default="Done! \n",help="text to print at end")
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
      hist.Draw("hist")
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
    if index: os.system('cp /afs/cern.ch/user/e/escott/www/HGCclustering/Pass23/index.php %s'%thedir)

def initGenHists():
  pass

def initMultiHists(enRadii,scRadii):
  hists = {}
 
  hists['nMultis'] = r.TH1F('hMulti_nMultis','hMulti_nMultis',51,-0.5,50.5)
  hists['nSeedMultis'] = r.TH1F('hMulti_nSeedMultis','hMulti_nSeedMultis',51,-0.5,50.5)
  hists['twodEnFrac'] = r.TH1F('hMulti_twodEnFrac','hMulti_twodEnFrac',100,1.,0.)
  hists['bestEnFrac'] = r.TH1F('hMulti_bestEnFrac','hMulti_bestEnFrac',100,1.,0.)
  hists['bestPfEnFrac'] = r.TH1F('hMulti_bestPfEnFrac','hMulti_bestPfEnFrac',100,1.,0.)
  hists['twoBestPfEnFrac'] = r.TH1F('hMulti_twoBestPfEnFrac','hMulti_twoBestPfEnFrac',100,1.,0.)
  hists['totalEnFrac'] = r.TH1F('hMulti_totalEnFrac','hMulti_totalEnFrac',100,1.,0.)
  hists['total2DEnFrac'] = r.TH1F('hMulti_total2DEnFrac','hMulti_total2DEnFrac',100,1.,0.)
  hists['totalRechitsInMultisEnFrac'] = r.TH1F('hMulti_totalRechitsInMultisEnFrac','hMulti_totalRechitsInMultisEnFrac',100,1.,0.)
  hists['suggestedEnFrac'] = r.TH1F('hMulti_suggestedEnFrac','hMulti_suggestedEnFrac',100,1.,0.)
  hists['layersInBest'] = r.TH1F('hMulti_layersInBest','hMulti_layersInBest',53,-0.5,52.5)
  hists['layersInBestWeighted'] = r.TH1F('hMulti_layersInBestWeighted','hMulti_layersInBestWeighted',53,-0.5,52.5)
  hists['num2DsInBest'] = r.TH1F('hMulti_num2DsInBest','hMulti_num2DsInBest',53,-0.5,52.5)
  hists['num2DsInBestWeighted'] = r.TH1F('hMulti_num2DsInBestWeighted','hMulti_num2DsInBestWeighted',53,-0.5,52.5)
  hists['hOverEInBest'] = r.TH1F('hMulti_hOverEInBest','hMulti_hOverEInBest',100,0.,0.1)
  hists['hOverEInBestWeighted'] = r.TH1F('hMulti_hOverEInBestWeighted','hMulti_hOverEInBestWeighted',100,0.,0.1)
  hists['maxClusterInBestX'] = r.TH1F('hMulti_maxClusterInBestX','hMulti_hMulti_maxClusterInBestX',100,1.,0.)
  hists['maxClusterInBestY'] = r.TH1F('hMulti_maxClusterInBestY','hMulti_hMulti_maxClusterInBestY',100,1.,0.)
  hists['maxClusterInBestZ'] = r.TH1F('hMulti_maxClusterInBestZ','hMulti_hMulti_maxClusterInBestZ',100,1.,0.)
  hists['maxClusterInBestLayer'] = r.TH1F('hMulti_maxClusterInBestLayer','hMulti_hMulti_maxClusterInBestZ',53,-0.5,52.5)
  hists['maxClusterInBestGenDrho'] = r.TH1F('hMulti_maxClusterInBestGenDrho','hMulti_hMulti_maxClusterInBestGenDrho',100,0.,5.)
  hists['maxClusterInBestGenDeta'] = r.TH1F('hMulti_maxClusterInBestGenDeta','hMulti_hMulti_maxClusterInBestGenDeta',100,-0.01,0.01)
  hists['maxClusterInBestGenDphi'] = r.TH1F('hMulti_maxClusterInBestGenDphi','hMulti_hMulti_maxClusterInBestGenDphi',100,-0.1,0.1)
  hists['maxClusterInBestGenDetaCm'] = r.TH1F('hMulti_maxClusterInBestGenDetaCm','hMulti_hMulti_maxClusterInBestGenDetaCm',100,-0.5,0.5)
  hists['maxClusterInBestGenDphiScaled'] = r.TH1F('hMulti_maxClusterInBestGenDphiScaled','hMulti_hMulti_maxClusterInBestGenDphiScaled',100,-0.1,0.1)
  hists['drGenBestAtFace'] = r.TH1F('hMulti_drGenBestAtFace','hMulti_drGenBestAtFace',100,0.,5.)
  hists['drGenBestAtBestZ'] = r.TH1F('hMulti_drGenBestAtBestZ','hMulti_drGenBestAtBestZ',100,0.,5.)
  hists['drGenBestAtFaceCorr'] = r.TH1F('hMulti_drGenBestAtFaceCorr','hMulti_drGenBestAtFaceCorr',100,0.,5.)
  hists['drGenBestAtBestZCorr'] = r.TH1F('hMulti_drGenBestAtBestZCorr','hMulti_drGenBestAtBestZCorr',100,0.,5.)

  for enRadius in enRadii:
    hists['bestEnFrac_Alt%s'%str(enRadius)] = r.TH1F('hMulti_bestEnFrac_Alt%s'%str(enRadius),'hMulti_bestEnFrac_Alt%s'%str(enRadius),100,1.,0.)
    hists['totalEnFrac_Alt%s'%str(enRadius)] = r.TH1F('hMulti_totalEnFrac_Alt%s'%str(enRadius),'hMulti_totalEnFrac_Alt%s'%str(enRadius),100,1.,0.)
    hists['total2DEnFrac_Alt%s'%str(enRadius)] = r.TH1F('hMulti_total2DEnFrac_Alt%s'%str(enRadius),'hMulti_total2DEnFrac_Alt%s'%str(enRadius),100,1.,0.)
    for scRadius in scRadii:
      nameStr = 'superEnFrac_Alt%s_Sc%s'%(str(enRadius),str(scRadius))
      hists[nameStr] = r.TH1F('hMulti_%s'%nameStr,'hMulti_%s'%nameStr,100,1.,0.)

  hists['num2DsInSelected'] = r.TH1F('hMulti_num2DsInSelected','hMulti_num2DsInSelected',53,-0.5,52.5)
  hists['num2DsInSelectedWeighted'] = r.TH1F('hMulti_num2DsInSelectedWeighted','hMulti_num2DsInSelectedWeighted',53,-0.5,52.5)
  hists['hOverEInSelected'] = r.TH1F('hMulti_hOverEInSelected','hMulti_hOverEInSelected',100,0.,0.1)
  hists['hOverEInSelectedWeighted'] = r.TH1F('hMulti_hOverEInSelectedWeighted','hMulti_hOverEInSelectedWeighted',100,0.,0.1)

  hists['TotalAlt2VsBestFrac'] = r.TH2F('hMulti_TotalAlt2VsBestFrac','hMulti_TotalAlt2VsBestFrac',20,1.,0.,20,1.,0.)
  hists['TotalAlt2VsNumMultis'] = r.TH2F('hMulti_TotalAlt2VsNumMultis','hMulti_TotalAlt2VsNumMultis',20,1.,0.,20,1.,0.)
  hists['TotalAlt2VsLargestScaledDPhi'] = r.TH2F('hMulti_TotalAlt2VsLargestScaledDPhi','hMulti_TotalAlt2VsLargestScaledDPhi',20,1.,0.,20,1.,0.)
  hists['TotalAlt2VsRmsPhi'] = r.TH2F('hMulti_TotalAlt2VsRmsPhi','hMulti_TotalAlt2VsRmsPhi',20,1.,0.,20,1.,0.)
  hists['BestFracVsNumMultis'] = r.TH2F('hMulti_BestFracVsNumMultis','hMulti_BestFracVsNumMultis',20,1.,0.,20,1.,0.)
  hists['BestFracVsLargestScaledDPhi'] = r.TH2F('hMulti_BestFracVsLargestScaledDPhi','hMulti_BestFracVsLargestScaledDPhi',20,1.,0.,20,1.,0.)
  hists['BestFracVsRmsPhi'] = r.TH2F('hMulti_BestFracVsRmsPhi','hMulti_BestFracVsRmsPhi',20,1.,0.,20,1.,0.)

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
    if str(type(hist)) == '<class \'ROOT.TH2F\'>':
      hist.Draw("colz")
    else:
      hist.Draw("hist")
    canv.Print(outdir+hist.GetName()+".pdf")
    canv.Print(outdir+hist.GetName()+".png")

def fitHists(canv,hists,outdir):
  fitFile = open('Output/fits_%s_Pt%s_Conv%s.txt'%(opts.particleType,opts.ptVal,opts.doConverted),'w')
  for key,hist in hists.iteritems():
    canv.cd() 
    setXTitle(key,hist)
    hist.Draw("hist")
    #fit here
    fit = r.TF1("fit","gaus")
    hist.Fit(fit)
    mean  = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    res = 0.
    if mean > 0.: res = sigma/mean
    effSigma = getEffSigma(hist)
    effRes = 0.
    if mean > 0.: effRes = effSigma/mean
    fitFile.write('hist %s: mean %1.4f, sigma %1.4f, effSigma %1.4f, res %1.4f, effRes %1.4f\n'%(key,mean,sigma,effSigma,res,effRes))
  fitFile.close()

def getEffSigma( theHist, wmin=0.2, wmax=1.8, step=0.001, epsilon=0.007 ):
#def getEffSigma( theHist, wmin=0.1, wmax=2.8, step=0.0002, epsilon=0.005 ):
  point = wmin
  weight = 0.
  points = [] #vector<pair<double,double> > 
  thesum = theHist.Integral()
  for i in range(theHist.GetNbinsX()):
    weight += theHist.GetBinContent(i)
    if weight > epsilon:
      points.append( [theHist.GetBinCenter(i),weight/thesum] )
  low = wmin
  high = wmax
  width = wmax-wmin
  for i in range(len(points)):
    for j in range(i,len(points)):
      wy = points[j][1] - points[i][1]
      if abs(wy-0.683) < epsilon:
        wx = points[j][0] - points[i][0]
        if wx < width:
          low = points[i][0]
          high = points[j][0]
          width=wx
  return 0.5*(high-low)

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

def etaToTheta(val):
  return 2 * atan( exp(-1*val) )

def thetaToEta(val):
  return -1 * log( tan(val/2.) )

def deltaPhi( phi1, phi2):
  dPhi = phi1 - phi2;
  pi = 3.14159265;
  if     ( dPhi <=-pi): dPhi += 2.0*pi;
  elif( dPhi >  pi): dPhi -= 2.0*pi;
  return dPhi;


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
  #convert comma-separated list of radii into list of floats
  enRadii = [int(val) for val in opts.enRadii.split(',')]
  scRadii = [int(val) for val in opts.scRadii.split(',')]

  #do the normal set of basic plots
  if opts.doBasics:
    plots = BasicPlots()
    plots.fillHists(theTree)
    plots.drawHists(opts.outDir)

  #setup hists for each type of object
  genHists   = initGenHists()
  multiHists = initMultiHists(enRadii,scRadii)

  #loop over entries in tree and actully do things
  for iEntry in range(theTree.GetEntries()):
  #for iEntry in range(1):
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
      #if opts.doConverted==1 and not abs(genDvz)>220.:
        continue
      if opts.doConverted==2 and genReachedEE:
      #if opts.doConverted==2 and abs(genDvz)>220.:
        continue
      genEnergy = genEnergies[iGen]
      genPhi = genPhis[iGen]

      #setup track values - not working atm
      #trackEtas = getattr(theTree,"track_eta")
      #print len(genEtas)
      #print len(trackEtas)
      #exit("testing")

      #get simcluster quantities
      simEtas = getattr(theTree,"simcluster_eta")
      simPhis = getattr(theTree,"simcluster_phi")
      simPts = getattr(theTree,"simcluster_pt")
      simEnergies = getattr(theTree,"simcluster_energy")
      totSimEnergy = 0.
      for iSim in range(simEnergies.size()):
        simEta = simEtas[iSim]
        if simEta*genEta < 0.: continue
        simEnergy = simEnergies[iSim]
        totSimEnergy += simEnergy
      #print "totSimEnergy",totSimEnergy

      #get pfcluster quantities
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
      secondBestPfIndex = -9999
      secondBestPfEnergy = -9999.
      for iPf in range(pfEnergies.size()):
        pfEta = pfEtas[iPf]
        if( pfEta*genEta < 0. ):
          continue
        if iPf == bestPfIndex:
          continue
        pfEnergy = pfEnergies[iPf]
        if pfEnergy > secondBestPfEnergy:
          secondBestPfEnergy = pfEnergy
          secondBestPfIndex = iPf
      if bestPfEnergy>0.: multiHists['bestPfEnFrac'].Fill(bestPfEnergy/genEnergy)
      if bestPfEnergy+secondBestPfEnergy>0.: multiHists['twoBestPfEnFrac'].Fill((bestPfEnergy+secondBestPfEnergy)/genEnergy)
      #print "bestPfEnergy",bestPfEnergy

      #get rechit quantities
      rechitEnergies = getattr(theTree,"rechit_energy")
      rechitEtas = getattr(theTree,"rechit_eta")
      rechitPhis = getattr(theTree,"rechit_phi")
      rechitPts = getattr(theTree,"rechit_pt")
      rechitZees = getattr(theTree,"rechit_z")
      rechitLayers = getattr(theTree,"rechit_layer")
      rechit2Ds = getattr(theTree,"rechit_cluster2d")
      rechitXs = getattr(theTree,"rechit_x")
      rechitYs = getattr(theTree,"rechit_y")
      
      #get 2D quantities
      twodRechits = getattr(theTree,"cluster2d_rechits")
      twodEnergies = getattr(theTree,"cluster2d_energy")
      twodEtas = getattr(theTree,"cluster2d_eta")
      totalRechitsInMultisEnergy = 0.
      total2Denergy = 0.
      totalAlt2Denergies = {enRadius: 0. for enRadius in enRadii}
      #alt2Denergies = [{enRadius: 0. for enRadius in enRadii}] * len(twodEtas)
      alt2Denergies = [{enRadius: 0. for enRadius in enRadii} for twodEta in twodEtas]
      #then loop over to get total energy estimate, various ways
      #new version takes sum of rechits within r cm of the highest rechit in the layer cluster
      for i2D in range(len(twodEtas)):
        twodEta = twodEtas[i2D]
        if twodEta*genEta<0.: 
          continue
        total2Denergy += twodEnergies[i2D]
        recIndices = twodRechits[i2D]
        highestRecEnergy = -9999.
        highestRecIndex = -9999.
        highestRecEta = -9999.
        highestRecPhi = -9999.
        highestRecZ = -9999.
        highestRecX = -9999.
        highestRecY = -9999.
        for iRec in recIndices:
          recEnergy = rechitEnergies[iRec]
          recEta = rechitEtas[iRec]
          recPhi = rechitPhis[iRec]
          recZ = rechitZees[iRec]
          recX = rechitXs[iRec]
          recY = rechitYs[iRec]
          totalRechitsInMultisEnergy += recEnergy
          if recEnergy > highestRecEnergy:
            highestRecEnergy = recEnergy
            highestRecIndex = iRec
            highestRecEta = recEta
            highestRecPhi = recPhi
            highestRecZ = recZ
            highestRecX = recX
            highestRecY = recY
        for iRec in recIndices:
          recEnergy = rechitEnergies[iRec]
          recEta = rechitEtas[iRec]
          recPhi = rechitPhis[iRec]
          recZ = rechitZees[iRec]
          recX = rechitXs[iRec]
          recY = rechitYs[iRec]
          drRecBest = etasPhisZsToDeltaX(recEta,highestRecEta,recPhi,highestRecPhi,recZ,highestRecZ)
          actualDr = deltaX(recX,highestRecX,recY,highestRecY)
          for enRadius in enRadii:
            if drRecBest < enRadius:
              alt2Denergies[i2D][enRadius] += recEnergy
              totalAlt2Denergies[enRadius] += recEnergy
      #end loop over 2Ds
      for enRadius in enRadii: 
        multiHists['total2DEnFrac_Alt%s'%str(enRadius)].Fill( totalAlt2Denergies[enRadius] / genEnergy )
      multiHists['total2DEnFrac'].Fill( total2Denergy / genEnergy )
      multiHists['totalRechitsInMultisEnFrac'].Fill( totalRechitsInMultisEnergy / genEnergy )

      #get multi quantities
      multiEtas = getattr(theTree,"multiclus_eta")
      multiPhis = getattr(theTree,"multiclus_phi")
      multiPts = getattr(theTree,"multiclus_pt")
      multiEnergies = getattr(theTree,"multiclus_energy")
      multiZees = getattr(theTree,"multiclus_z")
      multi2Ds = getattr(theTree,"multiclus_cluster2d")

      #setup multi values, collections
      bestMultiIndex = -9999
      bestMultiEnergy = -9999.
      totalMultiEnergy = 0.
      nMultis = 0
      selectedMultiIndices = []
      seedMultiIndices = []
      #altMultiEnergies = [{enRadius: 0. for enRadius in enRadii}] * len(multiEtas)
      altMultiEnergies = [{enRadius: 0. for enRadius in enRadii} for multiEta in multiEtas]
      totalAltMultiEnergies = { enRadius: 0. for enRadius in enRadii }

      #loop over multis
      for iMulti in range(multiEtas.size()):
        #check same endcap
        multiEta = multiEtas[iMulti]
        if( multiEta*genEta < 0. ):
          continue
        nMultis += 1

        #analysis
        multiEnergy = multiEnergies[iMulti]
        multiPt = multiPts[iMulti]
        multiPhi = multiPhis[iMulti]
        multiZ = multiZees[iMulti]
        totalMultiEnergy += multiEnergy
        if multiEnergy > bestMultiEnergy:
          bestMultiIndex = iMulti
          bestMultiEnergy = multiEnergy
        multiN2D = len(multi2Ds[iMulti])
        if multiN2D >= opts.nClus:
          #print "multi with %d 2D being added"%multiN2D
          selectedMultiIndices.append(iMulti)
        if multiN2D >= 7 and multiPt > 1.: #seed multi selection
          seedMultiIndices.append(iMulti)
        #sanity check sum of all rechits in the sum-of-all-multis plot (will be SLOW)
        #and get alternative energy definition for multis
        #print "number of 2Ds for this multi is %d"%len(multi2Ds[iMulti])
        for i2D in multi2Ds[iMulti]:
          for enRadius in enRadii:
            altMultiEnergies[iMulti][enRadius] += alt2Denergies[i2D][enRadius]
            totalAltMultiEnergies[enRadius] += alt2Denergies[i2D][enRadius]

      #end multi loop
      if nMultis == 0:
        print "ALERT: zero multiclusters found in endcap"
        continue
      multiHists['nMultis'].Fill(nMultis)
      multiHists['nSeedMultis'].Fill(len(seedMultiIndices))
      multiHists['bestEnFrac'].Fill(bestMultiEnergy/genEnergy)
      multiHists['totalEnFrac'].Fill(totalMultiEnergy/genEnergy)

      bestMultiZ = multiZees[bestMultiIndex]
      bestMultiPhi = multiPhis[bestMultiIndex]
      bestMultiEta = multiEtas[bestMultiIndex]
      #following are used in superclustering
      bestMultiTheta = etaToTheta(bestMultiEta)
      bestMultiRho = deltaX( etaPhiZtoX(bestMultiEta,bestMultiPhi,bestMultiZ), 0., etaPhiZtoY(bestMultiEta,bestMultiPhi,bestMultiZ), 0. )
      bestMultiR = sqrt( bestMultiZ*bestMultiZ + bestMultiRho*bestMultiRho )
      bestMulti2Ds = multi2Ds[bestMultiIndex]

      for enRadius in enRadii:
        multiHists['bestEnFrac_Alt%s'%str(enRadius)].Fill(altMultiEnergies[bestMultiIndex][enRadius]/genEnergy)
        multiHists['totalEnFrac_Alt%s'%str(enRadius)].Fill(totalAltMultiEnergies[enRadius]/genEnergy)

      multiHists['num2DsInBest'].Fill(len(bestMulti2Ds))
      multiHists['num2DsInBestWeighted'].Fill(len(bestMulti2Ds),bestMultiEnergy)
      multiHists['drGenBestAtFace'].Fill( etasPhisZsToDeltaX(genEta,bestMultiEta,genPhi,bestMultiPhi,320.,320.) )
      multiHists['drGenBestAtBestZ'].Fill( etasPhisZsToDeltaX(genEta,bestMultiEta,genPhi,bestMultiPhi,bestMultiZ,bestMultiZ) )
      #dPhi = q*B*dZ/pZ; prefactor is the conversion of GeV to normal units NB DOES NOT WORK
      genDeltaPhiFace = 0.2998 * 1. * 3.8 * 0.01*(320.-genDvz) * (1./(genPt*cot(2*atan(exp(-1*genEta)))))
      genDeltaPhiBestZ = 0.2998 * 1. * 3.8 * 0.01*(bestMultiZ-genDvz) * (1./(genPt*cot(2*atan(exp(-1*genEta)))))
      multiHists['drGenBestAtFaceCorr'].Fill( etasPhisZsToDeltaX(genEta,bestMultiEta,genPhi+genDeltaPhiFace,bestMultiPhi,320.,320.) )
      multiHists['drGenBestAtBestZCorr'].Fill( etasPhisZsToDeltaX(genEta,bestMultiEta,genPhi+genDeltaPhiBestZ,bestMultiPhi,bestMultiZ,bestMultiZ) )

      #quick loop over 2Ds in best multi to get layer info
      twodEnergies = getattr(theTree,"cluster2d_energy")
      twodLayers = getattr(theTree,"cluster2d_layer")
      twodXs = getattr(theTree,"cluster2d_x")
      twodYs = getattr(theTree,"cluster2d_y")
      twodZs = getattr(theTree,"cluster2d_z")
      twodEtas = getattr(theTree,"cluster2d_eta")
      twodPhis = getattr(theTree,"cluster2d_phi")
      maxClusterInBestIndex = -9999
      maxClusterInBestEnergy = -9999.
      bestEmEnergy = 0.
      bestHadEnergy = 0.
      for iBest2D in range(len(bestMulti2Ds)):
        twodIndex = bestMulti2Ds[iBest2D]
        twodLayer = twodLayers[twodIndex]
        twodEnergy = twodEnergies[twodIndex]
        multiHists['layersInBest'].Fill(twodLayer)
        multiHists['layersInBestWeighted'].Fill(twodLayer,twodEnergy)
        if twodEnergy>maxClusterInBestEnergy:
          maxClusterInBestEnergy = twodEnergy
          maxClusterInBestIndex = twodIndex
        if twodLayer<=28:
          bestEmEnergy += twodEnergy
        elif twodLayer>28:
          bestHadEnergy += twodEnergy
      #end quick loop over 2Ds
      maxClusterInBestX = twodXs[maxClusterInBestIndex]
      maxClusterInBestY = twodYs[maxClusterInBestIndex]
      maxClusterInBestZ = twodZs[maxClusterInBestIndex]
      maxClusterInBestLayer = twodLayers[maxClusterInBestIndex]
      maxClusterInBestEta = twodEtas[maxClusterInBestIndex]
      maxClusterInBestPhi = twodPhis[maxClusterInBestIndex]
      maxClusterInBestGenDrho = etasPhisZsToDeltaX(genEta, maxClusterInBestEta, genPhi, maxClusterInBestPhi, maxClusterInBestZ, maxClusterInBestZ)
      maxClusterInBestGenDeta = maxClusterInBestEta - genEta
      maxClusterInBestGenDphi = deltaPhi(maxClusterInBestPhi, genPhi)
      multiHists['maxClusterInBestX'].Fill( maxClusterInBestX )
      multiHists['maxClusterInBestY'].Fill( maxClusterInBestY )
      multiHists['maxClusterInBestZ'].Fill( maxClusterInBestZ )
      multiHists['maxClusterInBestLayer'].Fill( maxClusterInBestLayer )
      multiHists['maxClusterInBestGenDrho'].Fill( maxClusterInBestGenDrho )
      multiHists['maxClusterInBestGenDeta'].Fill( maxClusterInBestGenDeta )
      multiHists['maxClusterInBestGenDphi'].Fill( maxClusterInBestGenDphi )
      
      recipRho = 1. / deltaX( etaPhiZtoX(1.5,0.,320.), 0., etaPhiZtoY(1.5,0.,320.), 0. )
      maxClusterInBestRho = deltaX( etaPhiZtoX(maxClusterInBestEta,maxClusterInBestPhi,maxClusterInBestZ), 0., etaPhiZtoY(maxClusterInBestEta,maxClusterInBestPhi,maxClusterInBestZ), 0. )
      maxClusterInBestR = sqrt( maxClusterInBestZ*maxClusterInBestZ + maxClusterInBestRho*maxClusterInBestRho )
      maxClusterInBestGenDetaCm = (etaToTheta(maxClusterInBestEta) - etaToTheta(genEta)) * maxClusterInBestR
      maxClusterInBestGenDphiScaled = maxClusterInBestGenDphi*(1./(maxClusterInBestRho*recipRho))
      multiHists['maxClusterInBestGenDetaCm'].Fill( maxClusterInBestGenDetaCm )
      multiHists['maxClusterInBestGenDphiScaled'].Fill( maxClusterInBestGenDphiScaled )
      if bestEmEnergy>0.:
        if bestHadEnergy/bestEmEnergy <0.1:
          multiHists['hOverEInBest'].Fill(bestHadEnergy/bestEmEnergy)
          multiHists['hOverEInBestWeighted'].Fill(bestHadEnergy/bestEmEnergy, bestMultiEnergy)
        else:
          multiHists['hOverEInBest'].Fill(0.099)
          multiHists['hOverEInBestWeighted'].Fill(0.099, bestMultiEnergy)
      else:
        multiHists['hOverEInBest'].Fill(0.099)
        multiHists['hOverEInBestWeighted'].Fill(0.099, bestMultiEnergy)

      #setup superclustering
      superAltEnergies = {(scRadius,enRadius): altMultiEnergies[bestMultiIndex][enRadius] for enRadius in enRadii for scRadius in scRadii}
      #print 'superAltEnergies',superAltEnergies

      #loop over selected multis again for superclustering
      #for hadrons just use 'simple' deltaX superclustering
      for iSel in selectedMultiIndices:
        multiEta = multiEtas[iSel]
        if multiEta*genEta < 0.:
          continue
        if iSel == bestMultiIndex:
          continue
        multiZ = multiZees[iSel]
        multiPhi = multiPhis[iSel]
        multiEnergy = multiEnergies[iSel]
        multiDx = etasPhisZsToDeltaX(multiEta,bestMultiEta,multiPhi,bestMultiPhi,multiZ,bestMultiZ)
        for scRadius in scRadii:
          if multiDx < scRadius: 
            for enRadius in enRadii: 
              superAltEnergies[(scRadius,enRadius)] += altMultiEnergies[iSel][enRadius]
      #end of superclustering loop
      for scRadius in scRadii:
        for enRadius in enRadii:
          multiHists['superEnFrac_Alt%s_Sc%s'%(str(enRadius),str(scRadius))].Fill( superAltEnergies[(scRadius,enRadius)] / genEnergy )

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
    os.system('mv Output/fits_%s_Pt%s_Conv%s.txt %s'%(opts.particleType,opts.ptVal,opts.doConverted,opts.webDir))

  #end main
  print opts.endText
  return 0


if __name__ == '__main__':
  main()
