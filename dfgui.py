"""
    SCRIPT: DFGUI
    *************************************
    A GUI run of Duperfit
"""

# Import modules
import os
import time
import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
# Custom packages
#from spec_handling import read_spec
#from temp_handling import *
#from duperfit import duperfit
from dfclass import Duperfit
# Plotting
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('font', family='serif')
plt.rcParams['errorbar.capsize'] = 3
opts = {'mec':'k', 'mew': 0.5, 'lw': 1}
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

#*************************************************************
# Classes
#*************************************************************
class FileBrowser(tk.Frame):
    """
    CLASS: FILEBROWSER
    ****************************************************************************
    File browser object, capable of optionally displaying output files
    ****************************************************************************
    Written by Michael Baer
    """
    def __init__(self,parent=None,fblbl="",fdlbl="",w=60,
                 filetypes=(("ASCII files","*.dat *.ascii *.flm"),("All files","*")),
                 ofd=False,offunc=None,btnstate=tk.NORMAL,sticky=tk.W):
        self.filename="";self.outfile=tk.StringVar()
        self.filetypes=filetypes
        tk.Frame.__init__(self,parent)
        self.initlbl=tk.Label(self,text=fblbl).grid(row=0,column=0,sticky=sticky)
        self.dispfilelabel=tk.Label(self,text="",width=w,height=2,relief="sunken",bd=4)
        self.filebrowseBtn=tk.Button(self,
                                text="Browse Files",state=btnstate,
                                command=lambda: self.browseFiles(ofd,offunc))
        self.dispfilelabel.grid(row=0,column=1,sticky=sticky)
        self.filebrowseBtn.grid(row=0,column=2,sticky=sticky)
        if ofd and offunc is not None:
            self.fdinlbl=tk.Label(self,text=fdlbl)
            self.fdinlbl.grid(row=1,column=0,sticky=sticky)
            #self.fddflbl=tk.Label(self,text="",width=60,height=2,relief="sunken",bd=4)
            self.fddflbl=ttk.Entry(self,textvariable=self.outfile,width=40,font=('Arial 12'))#,relief="sunken",bd=4)
            self.fddflbl.grid(row=1,column=1,columnspan=2,sticky=sticky)
    def browseFiles(self,ofd=False,offunc=None):
        self.filename=filedialog.askopenfilename(initialdir="./",
                                              title="Select a File",
                                              filetypes=self.filetypes)
        self.dispfilelabel.configure(text=self.filename)
        if ofd and offunc is not None:
            #self.outfile=offunc(self.filename)
            self.outfile.set(offunc(self.filename))
            self.fddflbl.configure(text=self.outfile)

class EntryGroup(tk.Frame):
    """
    CLASS: ENTRYGROUP
    ****************************************************************************
    Entry group object
    ****************************************************************************
    Written by Michael Baer
    """
    def __init__(self,parent=None,entries=[],invals=[],ew=5,eh=2,order="row",
                 sticky=tk.W):
        #check list lengths match
        if len(invals)!=len(entries):
            raise ValueError("invals and entries not same length")
        tk.Frame.__init__(self,parent)
        self.entries=[]
        ii=0
        for i,entry in enumerate(entries):
            #check width
            if isinstance(ew,int):
                wid=ew
            elif isinstance(ew,list):
                if len(ew)!=len(entries):
                    raise ValueError("w and entries not same length")
                wid=ew[i]
            else:
                raise ValueError("w not of right data type")
            #check height
            if isinstance(eh,int):
                hei=eh
            elif isinstance(eh,list):
                if len(eh)!=len(entries):
                    raise ValueError("h and entries not same length")
                hei=eh[i]
            else:
                raise ValueError("h not of right data type")
            if order=="row":
                elabeli=tk.Label(self,text=entry,height=hei).grid(row=ii,column=0,
                                                                  sticky=sticky)
                ei=tk.Entry(self,text="",width=wid)
                ei.insert(10,invals[i])
                ei.grid(row=ii,column=1,sticky=sticky)
                ii+=1
            elif order=="col":
                elabeli=tk.Label(self,text=entry,height=hei).grid(row=0,column=ii,
                                                                  sticky=sticky)
                ei=tk.Entry(self,text="",width=wid)
                ei.insert(10,invals[i])
                ei.grid(row=0,column=ii+1,sticky=sticky,padx=5)
                ii+=2
            self.entries.append(ei)
    def varstate(self):
        return map((lambda entry: entry.get()), self.entries)

class RbtnGroup(tk.Frame):
    """
    CLASS: RBTNGROUP
    ****************************************************************************
    Radial button group object, capable of adding a user input option with a
    FileBrowser
    ****************************************************************************
    Written by Michael Baer
    """
    def __init__(self,parent=None,picks=[],lblname="",order="col",dfalt=0,
                 uinpopt=False,sticky=tk.W):
        tk.Frame.__init__(self,parent)
        Rbtnlabel=tk.Label(self,text=lblname).grid(row=0,column=0,sticky=sticky)
        self.var=tk.IntVar()
        if uinpopt:
            def update_uinp():
                if self.var.get()==maxval+1:
                    self.browser.filebrowseBtn.configure(state=tk.ACTIVE)
                else:
                    self.browser.filebrowseBtn.configure(state=tk.DISABLED)
            for i,pick in enumerate(picks):
                if order=="col":
                    Rbtn=tk.Radiobutton(self,text=pick,variable=self.var,value=i+1,
                                        command=update_uinp).grid(row=0,column=i+1,
                                        sticky=sticky)
                elif order=="row":
                    Rbtn=tk.Radiobutton(self,text=pick,variable=self.var,value=i+1,
                                        command=update_uinp).grid(row=i+1,column=0,
                                        sticky=sticky)
            maxval=i+1 #capture the last option's value
            # Add the user input option
            Rbtn=tk.Radiobutton(self,text="User Input:",variable=self.var,
                                value=maxval+1,command=update_uinp)
            if order=="col":
                Rbtn.grid(row=0,column=maxval+1,sticky=sticky)
            elif order=="row":
                Rbtn.grid(row=maxval+1,column=0,sticky=sticky)
            self.browser=FileBrowser(self,btnstate=tk.DISABLED,w=30)
            if order=="col":
                self.browser.grid(row=1,column=0,columnspan=maxval+1,sticky=sticky)
            elif order=="row":
                self.browser.grid(row=maxval+2,column=0,columnspan=1,sticky=sticky)
            self.uinpfile=self.browser.filename
        else:
            for i,pick in enumerate(picks):
                if order=="col":
                    Rbtn=tk.Radiobutton(self,text=pick,variable=self.var,
                                        value=i+1).grid(row=0,column=i+1,
                                        sticky=sticky)
                elif order=="row":
                    Rbtn=tk.Radiobutton(self,text=pick,variable=self.var,
                                        value=i+1).grid(row=i+1,column=0,
                                        sticky=sticky)
    def varstate(self):
        return self.var.get()

class CbtnGroup(tk.Frame):
    """
    CLASS: FILEBROWSER
    ****************************************************************************
    Check button group object, capable of adding a user input option with a
    FileBrowser
    ****************************************************************************
    Written by Michael Baer
    """
    def __init__(self,parent=None,picks=[],lblname="",sticky=tk.W,maxrows=4,uinpopt=False):
        tk.Frame.__init__(self,parent)
        subframe=tk.Frame(self)
        subframe.grid(row=1)
        self.vars=[]
        ii=0;jj=0
        for pick in picks:
            var=tk.IntVar()
            Cbtn=tk.Checkbutton(subframe,text=pick,variable=var,bd=1)
            Cbtn.grid(row=ii,column=jj,sticky=sticky)
            self.vars.append(var)
            ii+=1
            if ii+1>maxrows:
                ii=0
                jj+=1
        maxcol=jj;lastrow=ii
        Cbtnlabel=tk.Label(self,text=lblname).grid(row=0,column=0,columnspan=maxcol,
                                                   padx=0,sticky=sticky)
        if uinpopt:
            def update_uinp():
                if uinpvar.get()==1:
                    self.browser.filebrowseBtn.configure(state=tk.ACTIVE)
                else:
                    self.browser.filebrowseBtn.configure(state=tk.DISABLED)
            uinpvar=tk.IntVar()
            Cbtn=tk.Checkbutton(subframe,text="User Input:",variable=uinpvar,
                                command=update_uinp)
            Cbtn.grid(row=lastrow,column=maxcol,sticky=sticky)
            self.browser=FileBrowser(subframe,btnstate=tk.DISABLED,w=30)
            self.browser.grid(row=maxrows+1,column=0,columnspan=maxcol+1,
                              sticky=sticky)
            self.vars.append(uinpvar)
    def varstate(self):
        return map((lambda var: var.get()), self.vars)

class Slider(tk.Frame):
    """
    CLASS: SLIDER
    ****************************************************************************
    Value slider object which actively displays values as the slider moves
    ****************************************************************************
    Written by Michael Baer
    """
    def __init__(self,parent=None,lblname="",slmin=0.0,slmax=1.0,slstep=0.01,
                 defval=0.7,slwidth=500,sticky=tk.W):
        tk.Frame.__init__(self,parent)
        var=tk.DoubleVar()
        Slabel=tk.Label(self,text=lblname).grid(row=0,column=0,sticky=sticky)
        self.Etext=ttk.Entry(self,textvariable=var)
        self.Etext.grid(row=1,column=0,sticky=sticky)
        self.slider=tk.Scale(self,from_=slmin,to=slmax,resolution=slstep,
                             orient=tk.HORIZONTAL,length=slwidth,showvalue=0,
                             variable=var)
        self.slider.set(defval)
        self.slider.grid(row=2,column=0,sticky=sticky)
    def varstate(self):
        return self.slider.get()

class DropDown(tk.Frame):
    """
    CLASS: DROPDOWN
    ****************************************************************************
    Dropdown option menu object
    ****************************************************************************
    Written by Michael Baer
    """
    def __init__(self,parent=None,lblname="",options=[],uinpopt=False,
                 filetypes=(("ASCII files","*.dat *.ascii *.flm"),("All files","*")),
                 sticky=tk.W):
        tk.Frame.__init__(self,parent)
        self.choice=tk.StringVar()
        dpdlabel=tk.Label(self,text=lblname).grid(row=0,column=0,sticky=sticky)
        if uinpopt:
            def update_uinp(event):
                if self.choice.get()=="User Input":
                    self.browser.filebrowseBtn.configure(state=tk.ACTIVE)
                else:
                    self.browser.filebrowseBtn.configure(state=tk.DISABLED)
            options.append("User Input")
            dpddrop=tk.OptionMenu(self,self.choice,*(options),command=update_uinp)
            dpddrop.grid(row=0,column=1,sticky=sticky)
            self.browser=FileBrowser(self,btnstate=tk.DISABLED,w=30,filetypes=filetypes)
            self.browser.grid(row=1,column=1,sticky=sticky)
        else:
            dpddrop=tk.OptionMenu(self,self.choice,*options).grid(row=0,column=1,
                                                                  sticky=sticky)

    def varstate(self):
        return self.choice.get()

class WindowFigure(tk.Frame):
    """
    CLASS: WINDOWFIGURE
    ****************************************************************************
    Figure object, capable of drawing spectra when user clicks a button to
    redraw
    ****************************************************************************
    Written by Michael Baer
    """
    def __init__(self,parent=None,w=5,h=4):
        tk.Frame.__init__(self,parent)
        self.fig=plt.Figure(figsize=(w,h),dpi=100)
        self.canvas=FigureCanvasTkAgg(self.fig,master=parent)
        self.canvas.draw()
        # placing the canvas on the Tkinter window
        self.canvas.get_tk_widget().pack()
        # creating the Matplotlib toolbar
        #self.toolbar = NavigationToolbar2Tk(self.canvas,
        #                               parent)
        #self.toolbar.update()
        # placing the toolbar on the Tkinter window
        #self.canvas.get_tk_widget().pack()
        # adding the subplot
        self.plot1=self.fig.add_subplot(111)
        self.plot1.grid(True,which='major',alpha=0.8,ls='-')
        self.plot1.grid(True,which='minor',alpha=0.4,ls='--')
    def redraw(self,expobj=None):# plotting the graph
        self.plot1.clear()
        if expobj is not None:
            fname=expobj.filename
        if os.path.isfile(fname):
            spec,isest=read_spec(fname,isest=True,silence=True)
            if spec is not None:
                self.plot1.plot(spec['wav'],spec['flux'],drawstyle='steps-mid',c='r',lw=0.5)
                self.plot1.grid(True,which='major',alpha=0.8,ls='-')
                self.plot1.grid(True,which='minor',alpha=0.4,ls='--')
                #self.plot1.set_xlabel('Wavelength (Angstrom)')
                self.canvas.draw()

#*************************************************************
# Functions
#*************************************************************
def getoutfile(filename):
    """
    FUNCTION: GETOUTFILE
    ****************************************************************************
    Function for making an output filename based on the input filename
    ****************************************************************************
    INPUT:
        filename    :   Name of the input file | character
    OUTPUT:
        outfile     :   Name of the output file | character
    ****************************************************************************
    Written by Michael Baer
    """
    parts=filename.split("/")
    basefile=parts[-1]
    parts=basefile.split(".")
    outfilebase=""
    for i,part in enumerate(parts[:-1]):
        if i==0:
            outfilebase+=part
        else:
            outfilebase+="."+part
    outfile="sf_result_"+outfilebase+".dfo"
    return outfile

def getentval(ent,dfalt,dtype="float"):
    """
    FUNCTION: GETENTVAL
    ****************************************************************************
    Function for obtaining values from text entries
    ****************************************************************************
    INPUT:
        ent     :   Entry value | character
        dfalt   :   Default value | arbitrary
    OPTIONAL INPUT:
        dtype   :   Data type | character, default `float`
    OUTPUT:
        ent     :   Entry value | arbitrary
    ****************************************************************************
    Written by Michael Baer
    """
    if ent=="":
        return dfalt
    else:
        if dtype=="int":
            return int(ent)
        elif dtype=="float":
            return float(ent)

def getrbtval(wvar,picks,dfalt):
    """
    FUNCTION: GETRBTNVAL
    ****************************************************************************
    Function for obtaining values from radial buttons
    ****************************************************************************
    INPUT:
        wvar    :   Location of choice | integer
        picks   :   List of choices | arbitrary
        dfalt   :   Default value | arbitrary
    ****************************************************************************
    Written by Michael Baer
    """
    if wvar==0:
        return dfalt
    else:
        return picks[wvar-1]

def runfunc():
    """
    FUNCTION: RUNFUNC
    ****************************************************************************
    Function on user clicking the Run button
    ****************************************************************************
    Written by Michael Baer
    """
    # Parameters
    ioe=0
    fname=obsBrowse.filename
    outfile=obsBrowse.outfile.get()
    print(outfile)
    if fname=="":
        messagebox.showerror('I/O Error','Provide an input file')
        ioe=1
    elif not os.path.isfile(fname):
        messagebox.showerror('I/O Error','File not found')
        ioe=2
    # Wavelength range
    wavpicks=list(wavEntry.varstate())
    beginw=getentval(wavpicks[0],2000.)
    endw=getentval(wavpicks[1],10000.)
    # Sigma-clipping
    wsig=sigRbtn.varstate()
    sigsrc=getrbtval(wsig,["incl","calc","none"],"calc")
    if sigsrc=="none":
        sclip=False
    else:
        sclip=True
    sigpicks=list(sigPars.varstate())
    nsig=getentval(sigpicks[0],2.7)
    grow=getentval(sigpicks[1],0,dtype="int")
    niter=getentval(sigpicks[2],5,dtype="int")
    # Error flux
    werr=esrcRbtn.varstate()
    isest=getrbtval(werr,[False,True],False)
    # Redshift grid
    zpicks=list(zEntry.varstate())
    zstart=getentval(zpicks[0],0.)
    zstop=getentval(zpicks[1],0.)
    zstep=getentval(zpicks[2],0.01)
    if zstop==0.:
        exactz=True
    else:
        exactz=False
    # Extinction fitting
    extpicks=list(extEntry.varstate())
    Avmin=getentval(extpicks[0],-2)
    Avmax=getentval(extpicks[1],2)
    Rv=getentval(extpicks[2],3.1)
    # Template library
    tempchoice=tempDrop.varstate()
    if tempchoice=="":
        tempchoice="All SNe"
    uinptemp=tempDrop.browser.filename
    # Wavelength binning
    resolution=getentval(list(binEntry.varstate())[0],20.)

    # Optimizer choice
    optimizer=optDrop.varstate()
    if optimizer=="":
        optimizer="L-BFGS-B"

    # Min fractional wavelength coverage
    sfractol=sfracSlide.varstate()
    # Galaxy selection
    galchoices=list(galChecks.varstate())
    allgals=gals+["User Input"]
    galpicks=[]
    for wi,g in zip(galchoices,allgals):
        if wi==1:
            galpicks.append(g)
    if len(galpicks)==0:
        galpicks=['E']
    uinpgal=galChecks.browser.filename
    if "User Input" in galpicks:
        if uinpgal=="":
            messagebox.showerror('I/O Error','Provide a galaxy file')
            ioe=1
        elif not os.path.isfile(uinpgal):
            messagebox.showerror('I/O Error','User-provided galaxy file not found')
            ioe=2
    # Weight source
    weightpick=weightRbtn.varstate()
    wsrc=getrbtval(weightpick,["none","tell","incl","usr"],"incl")
    if wsrc=="usr":
        wsrc=weightRbtn.browser.filename
        if wsrc=="":
            messagebox.showerror('I/O Error','Provide a weight file')
            ioe=1
        elif not os.path.isfile(wsrc):
            messagebox.showerror('I/O Error','User-provided weight file not found')
            ioe=2
    # Scale parameter limits
    tscale=CscaleSlide.varstate()
    gscale=DscaleSlide.varstate()
    if ioe==0:
        root.destroy()

        dfParams={}
        dfParams['IO']={}
        parts=fname.split("/")
        indir=""
        for part in parts[:-1]:
            indir+=part+"/"
        dfParams['IO']['object_dir']=indir
        dfParams['IO']['object_file']=fname.split("/")[-1]
        dfParams['IO']['output_path']=os.getcwd()+"/"
        dfParams['IO']['output_file']=outfile
        dfParams['IO']['SN_templates']=tempchoice
        dfParams['IO']['user_SN_template']=uinptemp
        dfParams['IO']['gal_templates']=galpicks
        dfParams['IO']['user_gal_template']=uinpgal

        dfParams['fit_params']={}
        dfParams['fit_params']['use_exact_z']=exactz
        dfParams['fit_params']['z_min']=zstart
        dfParams['fit_params']['z_max']=zstop
        dfParams['fit_params']['delta_z']=zstep
        dfParams['fit_params']['Av_min']=Avmin
        dfParams['fit_params']['Av_max']=Avmax
        dfParams['fit_params']['Rv']=3.1
        dfParams['fit_params']['max_template_scale']=tscale
        dfParams['fit_params']['max_galaxy_scale']=gscale

        dfParams['fit_weight']={}
        dfParams['fit_weight']['weight_source']=wsrc
        dfParams['fit_weight']['estimate_error']=isest

        dfParams['sigma_clipping']={}
        dfParams['sigma_clipping']['sigma_clip']=sclip
        dfParams['sigma_clipping']['sigma_source']=sigsrc
        dfParams['sigma_clipping']['n_iterations']=niter
        dfParams['sigma_clipping']['n_grow']=grow
        dfParams['sigma_clipping']['n_sigma']=nsig

        dfParams['wavelength_range']={}
        dfParams['wavelength_range']['min_wavelength']=beginw
        dfParams['wavelength_range']['max_wavelength']=endw
        dfParams['wavelength_range']['wavelength_bin']=resolution
        dfParams['wavelength_range']['minimum_wavelength_fraction']=sfractol

        dfParams['options']={}
        dfParams['options']['silence_messages']=True
        dfParams['options']['save_output']=True
        dfParams['options']['optimizer']=optimizer

        mySupernova = Duperfit(dfParams)
        mySupernova.fit(runRefit=True)

def cancelfunc():
    """
    FUNCTION: CANCELFUNC
    """
    root.destroy()

#*************************************************************
# Main Loop
#*************************************************************
if __name__ == '__main__':

    root=tk.Tk()
    root.title("Duperfit")
    #root.geometry("800x800")
    root.grid_columnconfigure(0,weight=1)
    root.grid_columnconfigure(1,weight=1)

    leftframe=tk.Frame(root)
    leftframe.grid(row=0,column=0)
    leftframe.rowconfigure(0,weight=1)
    leftframe.rowconfigure(1,weight=1)
    leftframe.rowconfigure(2,weight=1)
    leftframe.rowconfigure(3,weight=1)
    leftframe.rowconfigure(4,weight=1)
    leftframe.rowconfigure(5,weight=1)
    leftframe.rowconfigure(6,weight=1)
    leftframe.rowconfigure(7,weight=1)
    leftframe.rowconfigure(8,weight=1)

    ioframe=tk.Frame(leftframe,relief="groove",bd=1)
    ioframe.grid(row=0,sticky="nswe")
    obsBrowse=FileBrowser(ioframe,"Observation","Output File",ofd=True,offunc=getoutfile)
    obsBrowse.pack()

    wavsigframe=tk.Frame(leftframe)
    wavsigframe.grid(row=1,sticky="nswe")
    wavsigframe.columnconfigure(0,weight=1)
    wavsigframe.columnconfigure(1,weight=1)
    wavsigframe.columnconfigure(2,weight=1)
    wavframe=tk.Frame(wavsigframe,relief="groove",bd=1)
    wavframe.grid(row=0,column=0,sticky="nswe")
    wavEntry=EntryGroup(wavframe,entries=["Begin Wavelength:","End Wavelength:"],
                        invals=["2000","10000"])
    wavEntry.pack(padx=10)
    sigframe=tk.Frame(wavsigframe,relief="groove",bd=1)
    sigframe.grid(row=0,column=1,sticky="nswe")
    sigframe.rowconfigure(0,weight=1)
    sigframe.rowconfigure(1,weight=1)
    sigRbtn=RbtnGroup(sigframe,picks=["input","calculate","none"],lblname="Sigma:")
    sigRbtn.grid(row=0,padx=10)
    sigPars=EntryGroup(sigframe,entries=["nsigma","grow","niter"],invals=["2.7","0","5"],
                       order="col")
    sigPars.grid(row=1,padx=10)
    esrcframe=tk.Frame(wavsigframe,relief="groove",bd=1)
    esrcframe.grid(row=0,column=2,sticky="nswe")
    esrcRbtn=RbtnGroup(esrcframe,lblname="Error Flux:",picks=["included","estimate"],order="row")
    esrcRbtn.pack(padx=10)

    zextframe=tk.Frame(leftframe)
    zextframe.grid(row=2,sticky="nswe")
    zextframe.columnconfigure(0,weight=1)
    zextframe.columnconfigure(1,weight=1)
    zframe=tk.Frame(zextframe,relief="groove",bd=1)
    zframe.grid(row=0,column=0,sticky="nswe")
    zframe.columnconfigure(0,weight=1)
    zframe.columnconfigure(1,weight=1)
    zEntry=EntryGroup(zframe,entries=["Minimum z:","Maximum z:","Delta z:"],
                      invals=["0.00","0.00",""],order="col")
    zEntry.pack()
    extframe=tk.Frame(zextframe,relief="groove",bd=1)
    extframe.grid(row=0,column=1,sticky="nswe")
    extEntry=EntryGroup(extframe,entries=["Minimum Av:","Maximum Av:","Rv:"],
                        invals=["-2","2","3.1"],order="col")
    extEntry.pack()

    tempbinoptframe=tk.Frame(leftframe)
    tempbinoptframe.grid(row=3,sticky="nswe")
    tempbinoptframe.columnconfigure(0,weight=1)
    tempbinoptframe.columnconfigure(1,weight=1)
    #tempbinoptframe.columnconfigure(2,weight=1)
    tempframe=tk.Frame(tempbinoptframe,relief="groove",bd=1)
    tempframe.grid(row=0,column=0,sticky="nswe")
    tempopts=["All SNe","SNe <=10d","SNe Ia","SNe Ib","SNe Ic","SNe II","SNe IIn",
              "SLSN-I","SLSN-II"]
    tempDrop=DropDown(tempframe,lblname="Templates:",options=tempopts,
                      uinpopt=True,filetypes=(("Pickle files","*.pickle"),
                                              ("All files","*")))
    tempDrop.pack()
    binframe=tk.Frame(tempbinoptframe,relief="groove",bd=1)
    binframe.grid(row=0,column=1,sticky="nswe")
    binEntry=EntryGroup(binframe,entries=["Binning in Angstroms"],invals=["20"])
    binEntry.pack()
    optframe=tk.Frame(tempbinoptframe,relief="groove",bd=1)
    optframe.grid(row=0,column=2,sticky="nswe")
    optimizers=["L-BFGS-B","TRF"]
    optDrop=DropDown(optframe,lblname="Optimizer:",options=optimizers)
    optDrop.pack()

    sfracframe=tk.Frame(leftframe,relief="groove",bd=1)
    sfracframe.grid(row=4,sticky="nswe")
    sfracSlide=Slider(sfracframe,lblname="Minimum Template Wavelength Coverage:")
    sfracSlide.pack()

    galwgtframe=tk.Frame(leftframe)
    galwgtframe.grid(row=5,sticky="nswe")
    galwgtframe.columnconfigure(0,weight=1)
    galwgtframe.columnconfigure(1,weight=1)
    galframe=tk.Frame(galwgtframe,relief="groove",bd=1)
    galframe.grid(row=0,column=0,sticky="nswe")
    gals=['E','S0','Sa','Sb','Sc','SB1','SB2','SB3','SB4','SB5','SB6']
    galChecks=CbtnGroup(galframe,gals,lblname="Galaxies:",uinpopt=True)
    galChecks.grid(sticky="nswe")
    galChecks.rowconfigure(0,weight=1)
    galChecks.rowconfigure(1,weight=1)
    galChecks.rowconfigure(2,weight=1)
    galChecks.rowconfigure(3,weight=1)
    galChecks.rowconfigure(4,weight=1)
    galChecks.rowconfigure(5,weight=1)
    galChecks.columnconfigure(0,weight=1)
    galChecks.columnconfigure(1,weight=1)
    galChecks.columnconfigure(2,weight=1)
    wgtframe=tk.Frame(galwgtframe,relief="groove",bd=1)
    wgtframe.grid(row=0,column=1,sticky="nswe")
    weightRbtn=RbtnGroup(wgtframe,picks=["Unweighted","Telluric Deweighted","Error Flux"],
                         lblname="Weights:",uinpopt=True,order="row")
    weightRbtn.grid(sticky="nswe")
    weightRbtn.rowconfigure(0,weight=1)
    weightRbtn.rowconfigure(1,weight=1)
    weightRbtn.rowconfigure(2,weight=1)
    weightRbtn.rowconfigure(3,weight=1)
    weightRbtn.rowconfigure(4,weight=1)
    weightRbtn.rowconfigure(5,weight=1)

    scaleframe=tk.Frame(leftframe)
    scaleframe.grid(row=6,sticky="nswe")
    scaleframe.columnconfigure(0,weight=1)
    scaleframe.columnconfigure(1,weight=1)
    Csframe=tk.Frame(scaleframe,relief="groove",bd=1)
    Csframe.grid(row=0,column=0,sticky="nswe")
    CscaleSlide=Slider(Csframe,lblname="Maximum Template Scale:",slmin=0.01,slmax=3.0,
                       defval=3.0,slwidth=250)
    CscaleSlide.pack()
    Dsframe=tk.Frame(scaleframe,relief="groove",bd=1)
    Dsframe.grid(row=0,column=1,sticky="nswe")
    DscaleSlide=Slider(Dsframe,lblname="Maximum Galaxy Scale:",slmin=0.0,slmax=3.0,
                       defval=3.0,slwidth=250)
    DscaleSlide.pack()

    Btnframe=tk.Frame(leftframe)
    Btnframe.grid(row=7,sticky="nswe")
    Run=tk.Button(Btnframe,text="Run",command=runfunc)
    Run.pack(side=tk.LEFT)
    Cancel=tk.Button(Btnframe,text="Cancel",command=cancelfunc)
    Cancel.pack(side=tk.LEFT)

    rightframe=tk.Frame(root)
    rightframe.grid(row=0,column=1)

    specplot=WindowFigure(rightframe)
    specplot.pack()
    Redraw=tk.Button(Btnframe,text="Draw",command=lambda: specplot.redraw(obsBrowse))
    Redraw.pack(side=tk.LEFT)

    root.mainloop()
