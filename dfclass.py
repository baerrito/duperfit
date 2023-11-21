"""
    MODULE: DFCLASS
    *************************************
    Contains a Class object definition for Duperfit
    Intended to be generally more user-intuitive than the core functions
    Feedback on ease of use is greatly appreciated
"""

# Import modules
import os
import time
import json
import extinction
import astropy.io.ascii as tasc
from astropy.table import Table
# Custom packages
from spec_handling import *
from temp_handling import *
from duperfit import *
# Plotting
import matplotlib.pyplot as plt
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('font', family='serif')
plt.rcParams['errorbar.capsize'] = 3
opts = {'mec':'k', 'mew': 0.5, 'lw': 1}
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

try:
    dfdir=os.environ['DFDIR']
except KeyError:
    raise Exception("$DFDIR not recognized")


#############################################
# Most of this module is still being tested #
#############################################

def get_cum_frac(scores, typa, typb, ord='pos'):
	if ord=='pos':
		wsort = np.argsort(scores)
	elif ord=='neg':
		wsort = np.argsort(-scores)
	cscores = scores[wsort]
	#wf = np.isfinite(cscores)
	#cscores = cscores[wf]
	tot = len(cscores)
	cfd = []; dtype = [('frac',float),('score',float)]
	#cfd.append(np.array((0, 0), dtype=dtype))
	#cscore = 0
	for score in cscores:
		if ord=='pos':
			ws = cscores <= score
		elif ord=='neg':
			ws = cscores >= score
		frac = len(cscores[ws])/tot
		cfd.append(np.array((frac, score), dtype=dtype))
	cfd = np.array(cfd)
	return cfd

class Duperfit():

    def __init__(self,pars):

        odir=pars['IO']['object_dir']
        ofile=pars['IO']['object_file']
        self.fname=odir+ofile
        self.tempchoice=pars['IO']['SN_templates']
        self.uinptemp=pars['IO']['user_SN_template']
        self.galpicks=pars['IO']['gal_templates']
        self.uinpgal=pars['IO']['user_gal_template']
        self.output_path=pars['IO']['output_path']
        self.outfile=pars['IO']['output_file']

        if self.tempchoice=='User Input':
            if self.uinptemp=='':
                raise ValueError('Provide a template file')
        if 'User Input' in self.galpicks:
            if self.uinpgal=='':
                raise ValueError('Provide a galaxy file')

        self.exactz=pars['fit_params']['use_exact_z']
        self.zstart=pars['fit_params']['z_min']
        self.zstop=pars['fit_params']['z_max']
        self.zstep=pars['fit_params']['delta_z']
        self.Avmin=pars['fit_params']['Av_min']
        self.Avmax=pars['fit_params']['Av_max']
        self.Rv=pars['fit_params']['Rv']
        self.tscale=pars['fit_params']['max_template_scale']
        self.gscale=pars['fit_params']['max_galaxy_scale']

        self.wsrc=pars['fit_weight']['weight_source']
        self.isest=pars['fit_weight']['estimate_error']

        self.sclip=pars['sigma_clipping']['sigma_clip']
        self.sigsrc=pars['sigma_clipping']['sigma_source']
        self.niter=pars['sigma_clipping']['n_iterations']
        self.grow=pars['sigma_clipping']['n_grow']
        self.nsig=pars['sigma_clipping']['n_sigma']

        self.beginw=pars['wavelength_range']['min_wavelength']
        self.endw=pars['wavelength_range']['max_wavelength']
        self.resolution=pars['wavelength_range']['wavelength_bin']
        self.sfractol=pars['wavelength_range']['minimum_wavelength_fraction']

        self.silence=pars['options']['silence_messages']
        self.save=pars['options']['save_output']
        self.optimizer=pars['options']['optimizer']

        # Header
        print("*"*60)
        print("*"*3+" "*23+"DUPERFIT"+" "*23+"*"*3)
        print("*"*3+" "*11+"(Duplicated Superfit in Python)"+" "*12+"*"*3)
        print("*"*60)
        print("*"+" "*58+"*")
        print("*"+" "*4+"Written by Michael J. Baer"+" "*28+"*")
        print("*"+" "*4+"Adapted from Superfit by D. Andrew Howell"+" "*13+"*")
        print("*"+" "*58+"*")
        print("*"+" "*4+"When using this program, please cite:"+" "*17+"*")
        print("*"+" "*9+"D.A. Howell et al (2005), ApJ, 634, 1190"+" "*9+"*")
        print("*"+" "*9+"M.J. Baer (in Preparation)"+" "*23+"*")
        print("*"+" "*58+"*")
        print("*"*60)
        # Parameters
        print()
        print(" "*3+"PARAMETERS:")
        print()
        print(" "*6+f"Begin wavelength: {self.beginw:.1f} Angstrom")
        print(" "*6+f"End wavelength: {self.endw:.1f} Angstrom")
        print(" "*6+f"Binning: {self.resolution:.1f} Angstrom")
        print(" "*6+f"Minimum template wavelength coverage: {self.sfractol*100:.0f}%")
        print()
        if self.exactz:
            print(" "*6+f"Redshift: {self.zstart:.3f}")
        else:
            print(" "*6+f"zmin: {self.zstart:.3f}, zmax: {self.zstop:.3f}, dz: {self.zstep:.3f}")
        print(" "*6+f"Avmin: {self.Avmin:.1f}, Avmax: {self.Avmax:.1f}, Rv: {self.Rv:.1f}")
        print(" "*6+f"Maximum template scale: {self.tscale:.1f}")
        print(" "*6+f"Maximum galaxy scale: {self.gscale:.1f}")
        print()
        print(" "*6+f"Sigma source: {self.sigsrc}")
        print(" "*6+f"nsigma: {self.nsig}, grow: {self.grow}, niter: {self.niter}")
        print(" "*6+f"Weight source: {self.wsrc}, Estimate error flux: {self.isest}")
        print()
        print(" "*6+f"Galaxies: {self.galpicks}")
        print(" "*6+f"Templates: {self.tempchoice}")
        print(" "*6+f"Optimizer: {self.optimizer}")
        print()

        self.ospec,self.isest=read_spec(self.fname,isest=self.isest,silence=self.silence)
        if self.tempchoice=="All SNe":
            self.tempchoice="allsne"
        self.tempdict=get_tempdict(self.tempchoice,self.uinptemp)
        self.galdict=get_galdict(self.galpicks,self.uinpgal)

    def fit(self, refit=False):

        # if not refitting and an old output exists, load it
        if not refit:
            outfile=self.output_path+self.outfile
            if os.path.isfile(outfile):
                parts=self.fname.split("/")
                self.objname=parts[-1]
                print("Existing fit found!")
                self.output=tasc.read(outfile,format='fixed_width')
                print("Top 10 matches:")
                print(self.output[:10])
                return

        parts=self.fname.split("/")
        parts2=parts[-1].split(".")
        self.objname=''
        for part in parts2[:-1]:
            self.objname+=part
        start=time.time()
        self.output=duperfit(self.zstart,self.zstop,self.zstep,self.ospec,self.objname,
                        self.tempdict,self.galdict,Avmin=self.Avmin,Avmax=self.Avmax,
                        wsrc=self.wsrc,sclip=self.sclip,sigsrc=self.sigsrc,
                        niter=self.niter,grow=self.grow,nsig=self.nsig,
                        exactz=self.exactz,beginw=self.beginw,endw=self.endw,
                        resolution=self.resolution,sfractol=self.sfractol,
                        Rv=self.Rv,tscale=self.tscale,gscale=self.gscale,
                        silence=self.silence,save=self.save,
                        output_path=self.output_path,outfile=self.outfile)
        end=time.time()
        print(f"Runtime: {end-start:.1f}s")
        print("Top 10 matches:")
        print(self.output[:10])
        return

    def plot(self,iparam,sav=False):

        try:

            trow = self.output[iparam]

            for tfile, tspec in zip(self.tempdict['tfiles'],self.tempdict['tspecs']):
                if trow['SN'] in tfile:
                    Tempspec=np.copy(tspec)
                    temname=tfile.split("/")[-1]
                    break
            for gfile, gspec in zip(self.galdict['gfiles'],self.galdict['gspecs']):
                if trow['gal'] in gfile:
                    Galspec=np.copy(gspec)
                    galname=gfile.split("/")[-1]
                    break

            Tempspec['wav']*=(1+trow['z'])
            Galspec['wav']*=(1+trow['z'])
            Tempspec['flux']/=np.median(Tempspec['flux'])
            Galspec['flux']/=np.median(Galspec['flux'])
            minw=np.max([Tempspec['wav'][0],Galspec['wav'][0]])
            maxw=np.min([Tempspec['wav'][-1],Galspec['wav'][-1]])
            rt=np.max(Tempspec['wav'][1:] - Tempspec['wav'][:-1])
            rg=np.max(Galspec['wav'][1:] - Galspec['wav'][:-1])
            ares=np.max([rt,rg])
            Glspec = bin_spec_weighted(Galspec,minw,maxw,10.)
            Tmspec = bin_spec_weighted(Tempspec,minw,maxw,10.)
            redlaw=extinction.ccm89(Tmspec['wav'],1.,self.Rv)
            cc=trow['cSN']
            ff=trow['cgal']
            aa=-0.4*trow['Av']
            params = (cc,ff,aa)
            X = (Tmspec['flux'],Glspec['flux'],redlaw)
            modflux = sf_model(X,*params)

            S=trow['S']
            C=trow['cSN'];D=trow['cgal'];A_V=trow['Av']
            norm = 1./np.median(self.ospec['flux'])
            fig=plt.figure(figsize=(20,10))
            fig.set_facecolor('k')
            fig.set_edgecolor('k')
            ax=fig.add_subplot(111)
            ax.set_facecolor('k')
            wo=(self.ospec['wav']>=self.beginw)&(self.ospec['wav']<=self.endw)
            ax.plot(self.ospec['wav'][wo],self.ospec['flux'][wo]*norm,'-w',
                     drawstyle='steps-mid',label=self.objname)
            wt=(Tmspec['wav']>=self.beginw)&(Tmspec['wav']<=self.endw)
            ax.plot(Tmspec['wav'][wt],modflux[wt],'-r',drawstyle='steps-mid',
                     label=temname+', '+galname+f', $C={C:.4f}$, $D={D:.4f}$'\
                     +f', $A_V={A_V:.4f}$')
            ax.text(self.ospec['wav'].max(),(self.ospec['flux']*norm).max(),
                     f'S={S:.3f}',fontsize=16)
            ax.spines['bottom'].set_color('w')
            ax.spines['top'].set_color('w')
            ax.spines['left'].set_color('w')
            ax.spines['right'].set_color('w')
            ax.xaxis.label.set_color('w')
            ax.tick_params(axis='x', colors='w')
            ax.yaxis.label.set_color('w')
            ax.tick_params(axis='y', colors='w')
            ax.legend(shadow=True,fontsize=14,facecolor='k',labelcolor='w')
            ax.set_xlabel('Wavelength ($\\AA$)')
            ax.set_ylabel('Normalized F$_{\\lambda}$')
            ax.grid(True, which='major', alpha=0.8, ls='-',color='w')
            ax.grid(True, which='minor', alpha=0.4, ls='--',color='w')
            if sav:
                plt.savefig(self.output_path+'df_plot_'+self.objname.replace(".dat","")+'.jpg')
            plt.show()
            return

        except AttributeError:

            print("I don't have output to plot... use the fit function")
            return

    def MIDscore(self,snphase=0.,evaltypes=['SNIa','SNIb','SNIc','SNII','SNIIn','SLSN-I','SLSN-II'],
                 plot=False):

        try:

            SNlist=self.output['SN']

            self.MIDs={}
            Xarr=[]
            iarrs={}

            Xscors=[];Yscors=[]
            for itype in evaltypes:
                Xarr.append(itype)
                self.MIDs[itype]={}
                wtyp1=wheretype(itype,SNlist,conIIn=True)
                ii1=np.mean(wtyp1[:5])
                iXscors=[];iYscors=[];iarrs[itype]=[]
                for jtype in evaltypes:
                    if jtype==itype:
                        iarrs[itype].append(0)
                        continue
                    wtyp2=wheretype(jtype,SNlist,conIIn=True)
                    ii2=np.mean(wtyp2[:5])
                    self.MIDs[itype][jtype]=ii1-ii2
                    iXscors.append(ii1-ii2)
                    iYscors.append(ii2-ii1)
                    iarrs[itype].append(ii2-ii1)
                iarrs[itype]=np.array(iarrs[itype])
                Xscors.append(np.mean(iXscors))
                Yscors.append(np.mean(iYscors))
            Xarr=np.array(Xarr)

            data=(Xarr,);names=('X',);dtype=('str',)
            for itype in evaltypes:
                data+=(iarrs[itype],)
                names+=(itype,)
                dtype+=('float',)

            self.scortab=Table(data,names=names,dtype=dtype)
            for itype in evaltypes:
                self.scortab[itype].info.format='.1f'
            fname=self.objname+".midscore"
            fname=fname.replace(".dat","")
            tasc.write(self.scortab,self.output_path+fname,format='fixed_width',overwrite=True)
            fname="MIDscore_"+self.objname+".tex"
            fname=fname.replace(".dat","")
            tasc.write(self.scortab,self.output_path+fname,format='latex',overwrite=True)


            imn=np.argmin(Xscors);imx=np.argmax(Yscors)
            if imn==imx:
                self.besttype=evaltypes[imx]
            else:
                self.besttype='ambiguous'

            print()
            print(' '*3+'MID score evaluations:')

            print(' '*6+'Best-matched type is '+self.besttype)
            print(' '*6+'However, you should check the scores to be sure.')
            print()

            if plot:
                supscore = {}; ds = []; allkeys = set()
                for type in evaltypes:
                    fname = dfdir+f"/tempbank/score_eval/scores_{type}.pickle"
                    with open(fname, "rb") as f:
                        scores = pickle.load(f)
                    f.close()
                    ds.append(scores)
                    supscore[f"for_{type}"] = scores

                type1=evaltypes[imn]
                type2=evaltypes.copy()
                type2.remove(type1)

                alltypes=['SNIa','SNIb','SNIc','SNII','SNIIn','SLSN-I','SLSN-II']
                markers = ["o","D","^","s","*","v",">"]
                colors = ["blue","orange","g","r","magenta","cyan","yellow"]
                marks=[];cols=[]
                for i,itype in enumerate(alltypes):
                    if type1==itype:
                        mark1=markers[i]
                        col1=colors[i]
                for itype in type2:
                    for j,jtype in enumerate(alltypes):
                        if itype==jtype:
                            marks.append(markers[j])
                            cols.append(colors[j])
                fig, axs = plt.subplots(len(type2),2,figsize=(10,6*len(type2)),sharex='col',\
            							gridspec_kw={'width_ratios': [2, 1]})
                fig.set_facecolor('k')
                fig.set_edgecolor('k')

                ig=0
                for ll, ty in enumerate(type2):
                    phasea = np.array(supscore[f'for_{type1}'][type1][ty]['phases'])
                    scorea = np.array(supscore[f'for_{type1}'][type1][ty]['scores'])
                    phaseb = np.array(supscore[f'for_{ty}'][type1][ty]['phases'])
                    scoreb = np.array(supscore[f'for_{ty}'][type1][ty]['scores'])
                    snscore = self.MIDs[type1][ty]
                    #quit()

                    w=(phasea>=-20)&(phasea<=200)
                    phasea=phasea[w];scorea=scorea[w]
                    w=(phaseb>=-20)&(phaseb<=200)
                    phaseb=phaseb[w];scoreb=scoreb[w]

                    cfda = get_cum_frac(scorea, type1, ty, ord='neg')
                    cfdb = get_cum_frac(scoreb, type1, ty)

                    if len(cfda['frac'])>len(cfdb['frac']):
                        finterp=interp1d(cfda['frac'],cfda['score'])
                        cfaob=finterp(cfdb['frac'])
                        cfdif=cfdb['score']-cfaob
                        wmin=np.argmin(abs(cfdif))
                        fdmin=cfdb['score'][wmin]
                    else:
                        finterp=interp1d(cfdb['frac'],cfdb['score'])
                        cfboa=finterp(cfda['frac'])
                        cfdif=cfboa-cfda['score']
                        wmin=np.argmin(abs(cfdif))
                        fdmin=cfda['score'][wmin]

                    if len(type2)>1:
                        axs[ig,0].set_facecolor('k')
                        axs[ig,1].set_facecolor('k')
                        axs[ig,0].axhline(fdmin,ls=':',color='k')
                        axs[ig,0].scatter(phasea,scorea,marker=mark1,edgecolors=col1,\
                                        facecolors='none',label=type1,linewidths=0.8,\
                                        s=85)
                        axs[ig,0].scatter(phaseb,scoreb,marker=marks[ll],\
                                        edgecolors=cols[ll],facecolors='none',\
                                        label=f'{ty}',linewidths=1.2,s=75)
                        axs[ig,0].plot(snphase,snscore,marker='d',color='w',markersize=10)
                        axs[ig,0].set_ylabel(f"{type1}$-${ty}")
                        axs[ig,0].grid(True,which='major',ls='-',alpha=0.8,color='w')
                        axs[ig,0].grid(True,which='minor',ls='--',alpha=0.4,color='w')
                        axs[ig,0].set_xlim(-60,230)
                        axs[ig,0].spines['bottom'].set_color('w')
                        axs[ig,0].spines['top'].set_color('w')
                        axs[ig,0].spines['left'].set_color('w')
                        axs[ig,0].spines['right'].set_color('w')
                        axs[ig,0].xaxis.label.set_color('w')
                        axs[ig,0].tick_params(axis='x', colors='w')
                        axs[ig,0].yaxis.label.set_color('w')
                        axs[ig,0].tick_params(axis='y', colors='w')
                        axs[ig,0].legend(shadow=True,loc='upper right',fontsize=14,facecolor='k',labelcolor='w')
                        axs[ig,1].axhline(fdmin,ls=':',color='w')
                        axs[ig,1].plot(cfda['frac'],cfda['score'],c=col1,drawstyle='steps-mid',\
                                    label=type1)
                        axs[ig,1].plot(cfdb['frac'],cfdb['score'],c=cols[ll],drawstyle='steps-mid',\
                                    label=f'{ty}',ls='--')
                        axs[ig,1].set_xlim(0,1)
                        axs[ig,1].grid(True,which='major',ls='-',alpha=0.8,color='w')
                        axs[ig,1].grid(True,which='minor',ls='--',alpha=0.4,color='w')
                        #axs[ig,1].set_xlabel("Fraction")
                        axs[ig,1].spines['bottom'].set_color('w')
                        axs[ig,1].spines['top'].set_color('w')
                        axs[ig,1].spines['left'].set_color('w')
                        axs[ig,1].spines['right'].set_color('w')
                        axs[ig,1].xaxis.label.set_color('w')
                        axs[ig,1].tick_params(axis='x', colors='w')
                        axs[ig,1].yaxis.label.set_color('w')
                        axs[ig,1].tick_params(axis='y', colors='w')
                        axs[ig,1].legend(shadow=True,loc='upper left',fontsize=14,facecolor='k',labelcolor='w')
                        ig += 1
                    else:
                        axs[0].scatter(phasea, scorea, marker=mark1, edgecolors=col1, \
                                        facecolors='none', label=type1)
                        axs[0].scatter(phaseb, scoreb, marker=marks[ll], \
                                        edgecolors=cols[ll], facecolors='none', \
                                        label=f'{ty}')
                        axs[0].axhline(fdmin,ls=':',color='k')
                        axs[0].set_ylabel(f"{type1}$-${ty}")
                        axs[0].grid(True,which='major',ls='-',alpha=0.8)
                        axs[0].grid(True,which='minor',ls='--',alpha=0.4)
                        axs[1].plot(cfda['frac'], cfda['score'], c=col1, drawstyle='steps-mid', \
                                    label=type1)
                        axs[1].plot(cfdb['frac'], cfdb['score'], c=cols[ll], drawstyle='steps-mid', \
                                    label=f'{ty}')
                        axs[1].axhline(fdmin,ls='--',color='k',alpha=0.5)
                        axs[1].set_xlim(0,1)
                        axs[1].grid(True,which='major',ls='-',alpha=0.8)
                        axs[1].grid(True,which='minor',ls='--',alpha=0.4)
                        axs[1].legend(shadow=True)
                if len(type2)>1:
                    for ig in range(len(type2)):
                        axs[ig,1].set_yticklabels([])
                    #fig.tight_layout()
                    #fig.subplots_adjust(wspace=0, hspace=0)
                    axs[ig,0].set_xlabel("Light-curve phase (days)")
                    axs[ig,1].set_xlabel("Fraction")
                else:
                    axs[1].set_yticklabels([])
                    #fig.tight_layout()
                    #fig.subplots_adjust(wspace=0, hspace=0)
                    axs[0].set_xlabel("Light-curve phase (days)")
                    axs[1].set_xlabel("Fraction")

                fig.tight_layout(pad=0.5)
                fig.subplots_adjust(wspace=0, hspace=0)
                fname="MIDscore_"+self.objname+".png"
                fname=fname.replace(".dat","")
                #print(fname)
                #quit()
                plt.savefig(self.output_path+fname,dpi=150)
                #plt.show()

            return


        except AttributeError:

            print("I don't have output to score... use the fit function")
            return
