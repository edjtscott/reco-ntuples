import os
#passName = "Pass25"
#passName = "Pass25a"
passName = "Pass25_Cap12"

#particles = ["Photon","Pion"]
particles = ["Pion"]
#particles = ["Photon"]
#particles = ["Electron"]
#particles = ["Photon","Electron"]

#ptvals = ["35","25","15"]
#ptvals = ["35"]
ptvals = ["25"]
#ptvals = ["15"]
#ptvals = ["25","35","15"]
#ptvals = ["25","15"]

#names = {"Photon":"D17","Pion":"D17_255_225","Electron":"D17"}
#names = {"Photon":["LogWeightingOff","LogWeightingOn","DropNoMultis"],"Pion":["D17_255_225"],"Electron":["D17"]}
#names = {"Photon":["NewRadii"],"Pion":["D17_255_225"],"Electron":["D17"]}
#names = {"Photon":["LogWeightingOn","LogWeightingOff","DropNoMultis"],"Pion":["D17_255_225"],"Electron":["D17"]}
#names = {"Photon":["93X"],"Pion":["93X"],"Electron":["93X"]}
#names = {"Photon":["93X"],"Pion":["93XDropNoMultis","93X"],"Electron":["93X"]}
names = {"Photon":["93X"],"Pion":["93X_FH4cm","93X"],"Electron":["93X"]}
#prenames = ["","PU200_"]
#prenames = ["PU200_"]
prenames = [""]

for particle in particles:
  for ptval in ptvals:
    for prename in prenames:
      for postname in names[particle]:
        name = prename + postname
        infile = "/afs/cern.ch/work/e/escott/public/HGCStudies/Ntuples/partGun_%s_Pt%s_%s.root"%(particle,ptval,name)
        web = "/afs/cern.ch/user/e/escott/www/HGCclustering/%s/%s/%s_Pt%s/"%(passName,name,particle,ptval)

        if 'Photon' in particle:
          web = "/afs/cern.ch/user/e/escott/www/HGCclustering/%s/%s/%s_Pt%s/"%(passName,name+'_2cm',particle,ptval)
          os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 2 --etaWindow 2."%(infile,web,particle,ptval))
          os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 1 --etaWindow 2."%(infile,web,particle,ptval))
          #os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 2 --etaWindow 2."%(infile,web,particle,ptval))
          web = "/afs/cern.ch/user/e/escott/www/HGCclustering/%s/%s/%s_Pt%s/"%(passName,name+'_3cm',particle,ptval)
          #os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 1 --etaWindow 3."%(infile,web,particle,ptval))
          #os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 2 --etaWindow 3."%(infile,web,particle,ptval))
          web = "/afs/cern.ch/user/e/escott/www/HGCclustering/%s/%s/%s_Pt%s/"%(passName,name+'_4cm',particle,ptval)
          #os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 1 --etaWindow 4."%(infile,web,particle,ptval))
          #os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 2 --etaWindow 4."%(infile,web,particle,ptval))
          web = "/afs/cern.ch/user/e/escott/www/HGCclustering/%s/%s/%s_Pt%s/"%(passName,name+'_constEta',particle,ptval)
          #os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 1 --etaWindow 0.05"%(infile,web,particle,ptval))
          #os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 2 --etaWindow 0.05"%(infile,web,particle,ptval))

        web = "/afs/cern.ch/user/e/escott/www/HGCclustering/%s/%s/%s_Pt%s/"%(passName,name,particle,ptval)
        if 'Pion' in particle:
          #os.system("python plotHGCal_HADsuper.py -f %s -w %s -p %s -m %s -c 0"%(infile,web,particle,ptval))
          os.system("python plotHGCal_HADsuper.py -f %s -w %s -p %s -m %s -c 1"%(infile,web,particle,ptval))
          #os.system("python plotHGCal_HADsuper.py -f %s -w %s -p %s -m %s -c 1 --nClus 0"%(infile,web,particle,ptval))
          #os.system("python plotHGCal_HADsuper.py -f %s -w %s -p %s -m %s -c 2"%(infile,web,particle,ptval))

        web = "/afs/cern.ch/user/e/escott/www/HGCclustering/%s/%s/%s_Pt%s/"%(passName,name,particle,ptval)

        if 'Electron' in particle:
          web = "/afs/cern.ch/user/e/escott/www/HGCclustering/%s/%s/%s_Pt%s/"%(passName,name+'_2cm',particle,ptval)
          os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 1 --etaWindow 2."%(infile,web,particle,ptval))
          os.system("python plotHGCal_EMsuper.py -f %s -w %s -p %s -m %s -c 2 --etaWindow 2."%(infile,web,particle,ptval))
