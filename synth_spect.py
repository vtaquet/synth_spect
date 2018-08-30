#################################################################################
# pred_line.py
# read JPL/CDMS data file, compute Eup, Aij and fluxes of transitions from
# PD analysis
#################################################################################

from scipy.interpolate import interp1d
from numpy import array, arange, sin
import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import requests
import pandas as pd
import os

##---------------------------------
##---------------------------------
## input parameters
##---------------------------------
##---------------------------------

# check if directories exist
if os.path.isdir('predictions') == False:
  os.system("mkdir predictions")
if os.path.isdir('data') == False:
  os.system("mkdir data")
if os.path.isdir('data/cdms') == False:
  os.system("mkdir data/cdms")
if os.path.isdir('data/jpl') == False:
  os.system("mkdir data/jpl")

# read input file and store input params into dictionary
inpfile = "input.in"
inpascii = open(inpfile,'r')
inpparam = [] ; inpvalue = []
for line in inpascii:
  s = line.split()
  inpparam.append(s[0])
  inpvalue.append(s[2])
inp_dict = dict(zip(inpparam,inpvalue))

# assign input variables
choice_obs = inp_dict["choice_obs"]       # Use an observed spectrum? If yes, specify fileobs
choice_data = inp_dict["choice_data"]     # Where to look for spectro data (local or online)
filemod = inp_dict["filemod"]             # species name
fileobs = inp_dict["fileobs"]             # species name
prefix   = inp_dict["prefix"]             # directory location 
choice_y = inp_dict["choice_y"]           # Unit of intensity value: Tmb (K), Fpeak (Jy)
beamsize = float(inp_dict["beamsize"])    # Beam size arcsec
choice_tau = inp_dict["choice_tau"]       # Include opacity in calculations (yes or no)
numin = float(inp_dict["numin"])          # Minimum frequency in GHz 
numax = float(inp_dict["numax"])          # Maximum frequency in GHz
dnu = float(inp_dict["dnu"])              # Spectral resolution in MHz
Vmin = float(inp_dict["Vmin"])            # Minimal value for param chosen above
rms = float(inp_dict["rms"])              # Minimal value for param chosen above
choice_plot = int(inp_dict["choice_plot"])# Output for the plot: 1) pdf, 2) python window

numin = numin*1e3 ; numax = numax*1e3

# read mod file
filemodascii = open(filemod,'r')
species = [] ; Ntot = [] ; Tex = [] ; ss = [] ; voff = [] ; dv = [] ; datab = [] ; vibf = [] ; ref = []
for line in filemodascii:
  if line[0:1] != "!" and line[0:1] != "#" and line[0:1] != "*":
    #s = line.split()
    species.append(   line[0:22].replace(" ",""))
    Ntot.append(float(line[22:28]))
    Tex.append( float(line[28:33]))
    ss.append(  float(line[33:38]))
    voff.append(float(line[38:44]))
    dv.append(  float(line[44:49]))
    datab.append(    (line[49:55].replace(" ","")))
    vibf.append(float(line[55:61]))
    ref.append(      (line[61:-1].replace(" ","")))
  elif line[0:1] == '*':
    break
modfile_dict = {"species" : species,
                "Ntot" : Ntot,
                "Tex" : Tex,
                "ss" : ss,
                "voff" : voff,
                "dv" : dv,
                "datab" : datab,
                "vibf" : vibf,
                "ref" : ref}
mod_pd = pd.DataFrame(modfile_dict)  
mod_pd.index = mod_pd["species"]
Nsp = len(species)

# read obs file
if choice_obs == 'yes':
  fileobsascii = open(fileobs,'r')
  nuobs = [] ; Iobs = []
  for line in fileobsascii:
    s = line.split()
    nuobs.append(float(s[0]))
    Iobs.append(float(s[1]))
  Nchan = len(nuobs)
  dnu = nuobs[1]-nuobs[0]
  numin = [] ; numax = []
  numin.append(nuobs[0])
  for ic in range(Nchan-1):
    if np.isclose(nuobs[ic+1],nuobs[ic]+dnu,1e-4) == False:
      numax.append(nuobs[ic])
      numin.append(nuobs[ic+1])
numax.append(nuobs[Nchan-1])
obs_pd = pd.DataFrame({"nu" : nuobs, choice_y : Iobs})
if choice_y == "Fpeak":
  obs_pd["Tpeak"] = obs_pd[choice_y]/((1.222e6)**(-1)*beamsize**2*((numin[0]+numax[-1])/2./1e3)**2)*1e-3
elif choice_y == "Tpeak":
  obs_pd["Fpeak"] = (1.222e6)**(-1)*beamsize**2*((numin[0]+numax[-1])/2./1e3)**2*obs_pd[choice_y]*1e3

# constants
c=2.99e10; kb = 1.3803e-16; hb = 6.62e-27; lun = 1
Ktocm = 0.695; JtoeV = 1.602e-19; cmtoeV = 8065.54 


##------------------------------------
##------------------------------------
## read partition functions data
##------------------------------------
##------------------------------------

# read catdir from JPL and CDMS websites or locally
catdir_jpl = "https://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.cat"
catdir_cdms = "http://www.astro.uni-koeln.de/site/vorhersagen/catalog/partition_function.html"
#
if choice_data == 'online': # read the catdir online
  resp_jpl = requests.get(catdir_jpl)
  catdir_file_jpl = resp_jpl.text.split("\n")
  del(catdir_file_jpl[-1])
  resp_cdms = requests.get(catdir_cdms)
  catdir_file_cdms = resp_cdms.text.split("\n")
  catdir_file_cdms = [line.replace("---", "0.0") for line in catdir_file_cdms]
  del(catdir_file_cdms[0:14])
  del(catdir_file_cdms[-5:-1])
  del(catdir_file_cdms[-1])
  # write the catdir files locally
  with open("data/catdir_jpl.cat", 'w') as catjpl:
    for line in catdir_file_jpl:
      catjpl.write(line+'\n')
  with open("data/catdir_cdms.cat", 'w') as catcdms:
    for line in catdir_file_cdms:
      catcdms.write(line+'\n')
elif choice_data == 'local': # read the catdir locally
  catdir_file_jpl = [] ; catdir_file_cdms = []
  with open("data/catdir_jpl.cat",'r') as catjpl:
    for line in catjpl.readlines():
      catdir_file_jpl.append(line[0:-1])
  with open("data/catdir_cdms.cat",'r') as catcdms:
    for line in catcdms.readlines():
      catdir_file_cdms.append(line[0:-1])

# create dictionary and panda frames
catdir_dic_jpl = {'tag': [int(line[0:7]) for line in catdir_file_jpl], 
              'species': [line[7:20].replace(" ","") for line in catdir_file_jpl], 
              'Nl':      [int(line[20:27]) for line in catdir_file_jpl],
              'T_300K':  [pow(10,float(line[27:33])) for line in catdir_file_jpl],
              'T_225K':  [pow(10,float(line[34:40])) for line in catdir_file_jpl],
              'T_150K':  [pow(10,float(line[41:47])) for line in catdir_file_jpl],
              'T_75K':   [pow(10,float(line[48:54])) for line in catdir_file_jpl],
              'T_37.5K': [pow(10,float(line[55:61])) for line in catdir_file_jpl],
              'T_18.75': [pow(10,float(line[62:68])) for line in catdir_file_jpl],
              'T_9.375K':[pow(10,float(line[69:75])) for line in catdir_file_jpl]}
catdir_dic_cdms = {'tag':[int(line[0:7]) for line in catdir_file_cdms], 
              'species': [line[7:33].replace(" ","") for line in catdir_file_cdms], 
              'Nl':      [int(line[33:39]) for line in catdir_file_cdms],
              'T_1000K': [pow(10,float(line[46:52])) for line in catdir_file_cdms],
              'T_500K':  [pow(10,float(line[59:65])) for line in catdir_file_cdms],
              'T_300K':  [pow(10,float(line[71:78])) for line in catdir_file_cdms],
              'T_225K':  [pow(10,float(line[84:91])) for line in catdir_file_cdms],
              'T_150K':  [pow(10,float(line[97:104])) for line in catdir_file_cdms],
              'T_75K':   [pow(10,float(line[110:117])) for line in catdir_file_cdms],
              'T_37.5K': [pow(10,float(line[123:130])) for line in catdir_file_cdms],
              'T_18.75': [pow(10,float(line[136:143])) for line in catdir_file_cdms],
              'T_9.375K':[pow(10,float(line[149:156])) for line in catdir_file_cdms],
              'T_5K':    [pow(10,float(line[162:169])) for line in catdir_file_cdms],
              'T_2.725K':[pow(10,float(line[175:183])) for line in catdir_file_cdms]}

catdir_pd_jpl = pd.DataFrame(catdir_dic_jpl)
catdir_pd_jpl.index = [catdir_pd_jpl["species"]]
catdir_pd_cdms = pd.DataFrame(catdir_dic_cdms)
catdir_pd_cdms.index = [catdir_pd_cdms["species"]]

        
##----------------------------------------
##----------------------------------------
## read spectro data from JPL/CDMS
##----------------------------------------
##----------------------------------------

for isp in range(Nsp):

  namespec = mod_pd.ix[[isp],["species"]].values.tolist()[0][0]
  datab = mod_pd.ix[[isp],["datab"]].values.tolist()[0][0]
  Tex = float(mod_pd.ix[[isp],["Tex"]].values.tolist()[0][0])
  Ntot = float(mod_pd.ix[[isp],["Ntot"]].values.tolist()[0][0])
  dv = float(mod_pd.ix[[isp],["dv"]].values.tolist()[0][0])
  ss = float(mod_pd.ix[[isp],["ss"]].values.tolist()[0][0])
  bd = (ss)**2/(beamsize**2+ss**2)
  if bd > 1.:
      bd = 1.

  print("Reading "+namespec+"...")

  # read spectro file from JPL or CDMS website
  if datab == 'jpl':
    try:
      pf = catdir_pd_jpl.loc[[namespec],['T_300K','T_225K','T_150K','T_75K','T_37.5K','T_18.75','T_9.375K']].values.tolist()[0]
    except:
      print("%15s is not in the JPL database: STOP" % (namespec))
      exit()
    Tpf = [300., 225., 150., 75., 37.5, 18.75, 9.375] ; Npf = len(Tpf)
    specurl = "https://spec.jpl.nasa.gov/ftp/pub/catalog/c%06i.cat" % (catdir_pd_jpl.loc[namespec,"tag"])
    specpath = 'data/'+datab+'/c%06i.cat' % (catdir_pd_jpl.loc[namespec,"tag"])
  elif datab == 'cdms':
    try:
      pf = catdir_pd_cdms.loc[[namespec],['T_1000K','T_500K','T_300K','T_225K','T_150K','T_75K','T_37.5K','T_18.75','T_9.375K','T_5K','T_2.725K']].values.tolist()[0]
    except:
      print("%15s is not in the CDMS database: STOP" % (namespec))
      exit()
    Tpf = [1000., 500., 300., 225., 150., 75., 37.5, 18.75, 9.375, 5., 2.725] ; Npf = len(Tpf)
    #print(str(catdir_pd_cdms.loc[namespec,"tag"]).zfill(6))
    specurl = "http://www.astro.uni-koeln.de/site/vorhersagen/catalog/c%06i.cat" % (catdir_pd_cdms.loc[namespec,"tag"])
    specpath = 'data/'+datab+'/c%06i.cat' % (catdir_pd_cdms.loc[namespec,"tag"])
    # assign pf values to closest ones when there is no data    
    for it in range(Npf):
        if pf[it] == 1.0:
            if it > 1:
                if pf[it-1] > 1.0:
                    pf[it] = pf[it-1]
                elif pf[it-1] == 1.0:
                    pf[it] = pf[it-2]
            else:
                if pf[it+1] > 1.0:
                    pf[it] = pf[it+1]
                elif pf[it+1] == 1.0:
                    pf[it] = pf[it+2]

  # compute partition function
  pf_interp = interp1d(Tpf,pf)
  Z = pf_interp(mod_pd.loc[namespec,"Tex"])*mod_pd.loc[namespec,"vibf"] # pf from specified Trot
  Zint = float(pf_interp(300.)) # pf at 300 K
  mod_pd["Z_300K"] = Zint
  mod_pd["Z_Trot"] = Z

  # read spectro datafile online or locally
  if choice_data == 'online':
    response = requests.get(specurl)
    specfile = response.text.split("\n")
    #open(specpath, 'w').write(response.content)
    with open(specpath, 'w') as catfile:
      for line in specfile:
        catfile.write(line+'\n')
  elif choice_data == 'local':
    specfile = []
    with open(specpath,'r') as catfile:
      for line in catfile.readlines():
        specfile.append(line[0:-1])
  del(specfile[-1])

  # create dictionary and panda frame from spectro file
  try:
    spec_pd2 = pd.DataFrame(
               {'species' : [namespec for line in specfile], 
                'nu':       [float(line[0:13]) for line in specfile], 
                'dnu':      [float(line[13:20]) for line in specfile],
                'intens':   [10**float(line[21:29]) for line in specfile],
                'Elow':     [float(line[31:41]) for line in specfile],
                'gup':      [float(line[41:44]) for line in specfile],
                'line':     [(line[44:-1]) for line in specfile],
                'Ntot':     [Ntot for line in specfile],
                'Tex':      [Tex for line in specfile],
                'dv':       [dv for line in specfile],
                'ss':       [ss for line in specfile],
                'bd':       [bd for line in specfile],
                'Z':        [Z for line in specfile],
                'Zint':     [Zint for line in specfile]})
  except: # gup is not a int in some CDMS files
    spec_pd2 = pd.DataFrame(
               {'species' : [namespec for line in specfile], 
                'nu':       [float(line[0:13]) for line in specfile], 
                'dnu':      [float(line[13:20]) for line in specfile],
                'intens':   [10**float(line[21:29]) for line in specfile],
                'Elow':     [float(line[31:41]) for line in specfile],
                'gup':      [float(line[42:44]) for line in specfile],
                'line':     [(line[44:-1]) for line in specfile],
                'Ntot':     [Ntot for line in specfile],
                'Tex':      [Tex for line in specfile],
                'dv':       [dv for line in specfile],
                'ss':       [ss for line in specfile],
                'bd':       [bd for line in specfile],
                'Z':        [Z for line in specfile],
                'Zint':     [Zint for line in specfile]})

  spec_pd2 = spec_pd2[(spec_pd2["nu"] >= numin[0]) & (spec_pd2["nu"] <= numax[-1])]

  # add transitions to dataframe
  if isp == 0:
    spec_pd = spec_pd2
  else:
    spec_pd = spec_pd.append(spec_pd2,ignore_index=True)


# compute intensities from spectro properties
spec_pd["Eup"] = spec_pd["Elow"] + hb*1e-7*spec_pd["nu"]*1e6*cmtoeV/JtoeV
spec_pd["Aij"] = spec_pd["intens"]*((spec_pd["nu"])**2)*spec_pd["Zint"]/spec_pd["gup"]*\
                 (np.exp(-spec_pd["Elow"]/Ktocm/3e2)-np.exp(-spec_pd["Eup"]/Ktocm/3e2))**(-1)*2.7964e-16

if choice_tau == 'yes':
    spec_pd["tau"] = c**3*spec_pd["Aij"]*spec_pd["gup"]*spec_pd["Ntot"]/(8.*np.pi*(spec_pd["nu"]*1e6)**3*spec_pd["dv"]*1e5*spec_pd["Z"])*\
                     np.exp(-spec_pd["Eup"]/Ktocm/spec_pd["Tex"])*(np.exp(hb*spec_pd["nu"]*1e6/(kb*spec_pd["Tex"]))-1) 
    spec_pd["Nup"] = 0.
    tau_0 = spec_pd["tau"] <= 0
    tau_1 = spec_pd["tau"] > 0
    spec_pd.loc[tau_1,"Nup"] = spec_pd.loc[tau_1,"gup"]*spec_pd.loc[tau_1,"Ntot"]/spec_pd.loc[tau_1,"Z"]*np.exp(-spec_pd.loc[tau_1,"Eup"]/Ktocm/spec_pd.loc[tau_1,"Tex"])*(1-np.exp(-spec_pd.loc[tau_1,"tau"]))/spec_pd.loc[tau_1,"tau"]
    spec_pd.loc[tau_0,"Nup"] = spec_pd.loc[tau_0,"gup"]*spec_pd.loc[tau_0,"Ntot"]/spec_pd.loc[tau_0,"Z"]*np.exp(-spec_pd.loc[tau_0,"Eup"]/Ktocm/spec_pd.loc[tau_0,"Tex"])
    #spec_pd["Nup"] = spec_pd["gup"]*spec_pd["Ntot"]/spec_pd["Z"]*np.exp(-spec_pd["Eup"]/Ktocm/spec_pd["Tex"])*(1-np.exp(-1e0*spec_pd["tau"]))/spec_pd["tau"]
                           
elif choice_tau == 'no':
    spec_pd["tau"] = 0.
    spec_pd["Nup"] = spec_pd["gup"]*spec_pd["Ntot"]/spec_pd["Z"]*np.exp(-spec_pd["Eup"]/Ktocm/spec_pd["Tex"])

spec_pd["Tint"] = spec_pd["Nup"]*spec_pd["bd"]*hb*c**3*spec_pd["Aij"]/(8e0*np.pi*kb*(spec_pd["nu"]*1e6)**2)/1e5
spec_pd["Tpeak"]= spec_pd["Tint"]/(1.06447*dv)
spec_pd["Fint"] = (1.222e6)**(-1)*beamsize**2*(spec_pd["nu"]/1e3)**2*spec_pd["Tint"]*1e3
spec_pd["Fpeak"]= (1.222e6)**(-1)*beamsize**2*(spec_pd["nu"]/1e3)**2*spec_pd["Tpeak"]*1e3

# sort dataframe according to frequency
spec_pd = spec_pd.sort_values(by=["nu"])


##----------------------------------------
##----------------------------------------
## write/plot output
##----------------------------------------
##----------------------------------------

# write transitions properties in ASCII file
where = spec_pd[choice_y] > Vmin
spec_pd_out = spec_pd[where]

filepf = open('predictions/'+prefix+"list_spect.cat",'w')
filepf.write('!            Species  Freq. (MHz)  Eup (K)     Aij(s-1)    Nup(cm-2)  '+\
             'Tmbdv (K km/s) Tpeak (K) Fint (mJy km/s) Fpeak (mJy)        tau\n')
for index, row in spec_pd_out.iterrows():
    filepf.write('%20s %12.3f %8.2f %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n' \
                % (row["species"], row["nu"], row["Eup"]/Ktocm, row["Aij"], row["Nup"], \
                   row["Tint"], row["Tpeak"], row["Fint"], row["Fpeak"], row["tau"]))

# compute the overall spectrum
def gaus(x,a,sigma,x0):
    #x0 = 0.
    offset = 0.
    #sigma = 0.5/2.355
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+offset

nu_plot = np.arange(numin[0],numax[-1],dnu)
T_plot = np.full_like(nu_plot, 0.)
F_plot = np.full_like(nu_plot, 0.) 
for index, row in spec_pd_out.iterrows():
    T_plot += gaus(nu_plot,row["Tpeak"],dv/2.355,row["nu"])
    F_plot += gaus(nu_plot,row["Fpeak"],dv/2.355,row["nu"])

if choice_y == 'Tpeak':
  rmst = rms
  rmsf = (1.222e6)**(-1)*beamsize**2*((numin[0]+numax[-1])/2./1e3)**2*rmst*1e3
  T_plot += np.random.rand(len(nu_plot))*rmst/2.-rmst/4.
  F_plot += np.random.rand(len(nu_plot))*rmsf/2.-rmsf/4.
elif choice_y == 'Fpeak':
  rmsf = rms
  rmst = rmsf/((1.222e6)**(-1)*beamsize**2*((numin[0]+numax[-1])/2./1e3)**2)*1e-3
  F_plot += np.random.rand(len(nu_plot))*rmsf/2.-rmsf/4.
  T_plot += np.random.rand(len(nu_plot))*rmst/2.-rmst/4.

# plot the spectrum
hfont = {'family' : 'normal',
         'weight' : 'medium',
         'size'   : 10}
haxes = {'linewidth' : 1.5}
hticks = {'major.size' : 6,
          'major.width' : 1.5,
          'minor.size' : 3,
          'minor.width' : 1}
plt.rc('font', **hfont)
plt.rc('axes',**haxes)
plt.rc('xtick',**hticks)
plt.rc('ytick',**hticks)
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

# Tmb plot
max_y = T_plot.max() ; min_y = -0.1*max_y
plt.figure(1,figsize=(20,5))
plt.xlabel('Frequency [GHz]')
plt.xlim(numin[0],numax[-1])
plt.ylabel('Tmb [K]')
plt.ylim(min_y,max_y)
plt.step(nu_plot,T_plot,'r')#,'fontname':'Helvetica')
if choice_obs == 'yes':
  plt.step(obs_pd["nu"].values.tolist(),obs_pd["Tpeak"].values.tolist(),'k')#,'fontname':'Helvetica')
if choice_plot == 2:
  for index, row in spec_pd_out.iterrows():
    plt.text(row["nu"], 0.5*min_y, row["species"], rotation=90, horizontalalignment='center',verticalalignment='center',fontsize=5)
if choice_plot == 1:
    plt.savefig('predictions/'+prefix+'spect_Tmb_all.eps',bbox_inches='tight')
    #plt.savefig('predictions/'+prefix+'spect_Tmb_all.pdf',bbox_inches='tight')
if choice_plot == 2:
    plt.show()

# Flux plot
max_y = F_plot.max() ; min_y = -0.1*max_y
plt.figure(1,figsize=(20,5))
plt.figure(2,figsize=(20,5))
plt.xlabel('Frequency [GHz]')
plt.xlim(numin[0],numax[-1])
plt.ylabel('F [mJy]')
plt.ylim(min_y,max_y)
plt.step(nu_plot,F_plot,'r')#,'fontname':'Helvetica')
if choice_obs == 'yes':
  plt.step(obs_pd["nu"].values.tolist(),obs_pd["Fpeak"].values.tolist(),'k')#,'fontname':'Helvetica')
if choice_plot == 2:
  for index, row in spec_pd_out.iterrows():
    plt.text(row["nu"], 0.5*min_y, row["species"], rotation=90, horizontalalignment='center',verticalalignment='center',fontsize=5)
if choice_plot == 1:
    plt.savefig('predictions/'+prefix+'spect_F_all.eps',bbox_inches='tight')
    #plt.savefig('predictions/'+prefix+'spect_F.pdf',bbox_inches='tight')
if choice_plot == 2:
    plt.show()
    
# 
if choice_plot == 1: 
  for iw in range(10):
    numin2 = numin[0] + iw*(numax[-1]-numin[0])/10.
    numax2 = numin[0] + (iw+1)*(numax[-1]-numin[0])/10.
    max_y = obs_pd["Tpeak"].max() ; min_y = -0.1*max_y
    plt.figure(iw,figsize=(20,5))
    plt.xlabel('Frequency [GHz]')
    plt.xlim(numin2,numax2)
    plt.ylabel('Tmb [K]')
    plt.ylim(min_y,max_y)
    plt.step(nu_plot,T_plot,'r')#,'fontname':'Helvetica')
    if choice_obs == 'yes':
      plt.step(obs_pd["nu"].values.tolist(),obs_pd["Tpeak"].values.tolist(),'k')#,'fontname':'Helvetica')
    plt.savefig('predictions/'+prefix+'spect_Tmb_zoom_'+str(iw)+'.eps',bbox_inches='tight')
    
    max_y = obs_pd["Fpeak"].max() ; min_y = -0.1*max_y
    plt.figure(iw+10,figsize=(20,5))
    plt.xlabel('Frequency [GHz]')
    plt.xlim(numin2,numax2)
    plt.ylabel('F [mJy]')
    plt.ylim(min_y,max_y)
    plt.step(nu_plot,F_plot,'r')#,'fontname':'Helvetica')
    if choice_obs == 'yes':
      plt.step(obs_pd["nu"].values.tolist(),obs_pd["Fpeak"].values.tolist(),'k')#,'fontname':'Helvetica')
    plt.savefig('predictions/'+prefix+'spect_F_zoom_'+str(iw)+'.eps',bbox_inches='tight')
    



