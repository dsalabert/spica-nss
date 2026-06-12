#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on

@author: dsalabert
"""
# Developed in Python 3.8.12
# Last modified december 2023

# import DataCursorPlot
import astropy.io.votable
import json
import platform
import tempfile
import warnings

import matplotlib.artist as artists
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pyvo as vo

# Handle Aspro2's model representation using a2p2 Models
from a2p2.jmmc import Models  # since a2p2 V0.3.3
from a2p2.jmmc import Catalog
from a2p2.jmmc.models import (
    _model,
)  # used to convert catalog's models to aspro's XML format
from astroplan import FixedTarget, Observer
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable.tree import Param
from astropy.samp import SAMPIntegratedClient
from astropy.table import MaskedColumn
from astropy.table import vstack, Table, join, join_skycoord
from astropy.time import Time
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from datetime import date, datetime
from functools import reduce
from idlelib.tooltip import Hovertip
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from pathlib import Path
from scipy.special import jv
from tkinter import *
from tkinter import ttk
from tkinter import messagebox, simpledialog

if platform.system() == "Darwin":
    from tkmacosx import Button

Simbad.add_votable_fields("otype")
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message=".*partition.*mask.*")
plt.ioff()

# from functools import partial
# import datetime
# import time
# from astropy.io import fits
# from astroquery.vizier import Vizier
# from astropy.coordinates import Angle
# from scipy.special import jv
# from astropy.coordinates import EarthLocation
# from astropy.table import Table, vstack
# from astroplan import AltitudeConstraint, AirmassConstraint, AtNightConstraint
# from astroplan import is_observable, is_always_observable, months_observable
# from astroplan import observability_table
# from matplotlib.figure import Figure
# from itertools import compress
# from PIL import ImageTk, Image


class spica_NSS:

    def __init__(self):
        """


        Returns
        -------
        None.

        """

        print(f"*** Welcome to the SPICA-NSS tool ({nssVersion}) ***")

        self.iter = 0

        # Initialize the main Tkinter window
        self.root = Tk()
        self.root.resizable(True, True)  # Width, Height
        self.root.title(f"SPICA-NSS TOOL ({nssVersion})")
        # self.root.geometry('1000x900')
        self.root.configure(background="#dddddd")
        self.myFontLabelFrame = ("Courier", 11, "bold")
        self.myFont = (
            "Courier",
            11,
        )  # 'bold')
        self.root.eval("::msgcat::mclocale en")  # Message in english

        # self.root.rowconfigure(0, weight=1)
        # self.root.columnconfigure(0, weight=1)

        """
        frame8 = Frame(self.root)
        frame8.grid()

        width = 180
        height = 75
        img = Image.open("/home/dsalabert/Workspace/SPICA-DB/GUI/erc_eu.png")
        img = img.resize((width,height), Image.ANTIALIAS)
        photoImg =  ImageTk.PhotoImage(img)
        label = Label(frame8, image=photoImg, width=180)
        label.grid(column=0, row=7)

        width = 180
        height = 75
        img2 = Image.open("/home/dsalabert/Workspace/SPICA-DB/GUI/OCAlogoQlarge_WEB_250px.png")
        # img = img.resize((width,height), Image.ANTIALIAS)
        photoImg2 =  ImageTk.PhotoImage(img2)
        label = Label(frame8, image=photoImg2)#, width=180)
        label.grid(column=1, row=7)

        width = 180
        height = 75
        img3 = Image.open("/home/dsalabert/Workspace/SPICA-DB/GUI/téléchargement.jpeg")
        #img3 = img3.resize((width,height), Image.ANTIALIAS)
        photoImg3 =  ImageTk.PhotoImage(img3)
        label = Label(frame8, image=photoImg3)#, width=180)
        label.grid(column=2, row=7)
        """

        # Define StringVar and IntVar for various inputs
        strDate = StringVar()
        self.FinalPriorityId = IntVar()
        self.strDecMin = StringVar()
        self.strDecMax = StringVar()
        self.strVmagMin = StringVar()
        self.strVmagMax = StringVar()
        self.strDecMean = StringVar()
        self.strVmagMean = StringVar()
        self.strRaRangePrim = StringVar()
        self.strDecRangePrim = StringVar()
        self.strVmagRangePrim = StringVar()
        self.strRaRangeSec = StringVar()
        self.strDecRangeSec = StringVar()
        self.strVmagRangeSec = StringVar()
        self.strLDDChiSec = StringVar()
        self.strRelErrorSec = StringVar()
        self.strMinVisSec = StringVar()
        self.intMaxBaseline = IntVar()

        # Initialisation of the different frames
        FrameDate = LabelFrame(
            self.root, text="Date (YYYY-MM-DD)", font=self.myFontLabelFrame
        )
        FrameActions = LabelFrame(
            self.root, text="Action Buttons", font=self.myFontLabelFrame
        )
        FrameWorkPackages = LabelFrame(self.root)
        FrameInstModes = LabelFrame(self.root)
        FramePriorities = LabelFrame(
            self.root, text="Priorities", font=self.myFontLabelFrame
        )
        self.FrameObjects = LabelFrame(
            self.root, text="Science objects", font=self.myFontLabelFrame
        )
        FrameAddTarget = LabelFrame(self.root, font=self.myFontLabelFrame)
        self.FrameCalPrims = LabelFrame(
            self.root, text="Primary calibrators", font=self.myFontLabelFrame
        )
        self.FrameCalSecs = LabelFrame(
            self.root, text="Secondary calibrators", font=self.myFontLabelFrame
        )
        FrameLog = Frame(self.root)

        # Date
        # labelDate = Label(Frame1a, font=self.myFont)
        entryDate = Entry(
            FrameDate,
            textvariable=strDate,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=10,
        )
        entryDate.insert(END, str(self.date))
        entryDate.bind("<Return>", (lambda _: self.entryDateCallback(entryDate)))
        self.labelValDate = Label(
            FrameDate, text=str(self.date), font=self.myFont, fg="red"
        )
        entryDate.focus()  # Place cursor into date entry at start
        myTip = Hovertip(entryDate, "Enter the chosen date.", hover_delay=1000)

        # Target query button
        buttonQuery = Button(
            FrameActions,
            text="QUERY_CATALOG",
            font=self.myFont,
            fg="white",
            bg="green",
            command=self.onQuery,
            cursor="hand1",
        )
        myTip = Hovertip(buttonQuery, "Click to query SPICA-DB.", hover_delay=1000)

        # Apply_filters button
        # buttonFilters = Button(FrameActions, text='APPLY_FILTERS', font=self.myFont, fg='white', bg='lightgoldenrod4',
        #                     command=self.onFilters, cursor='hand1', state=DISABLED)

        # Apply Best declination button
        buttonBestDec = Button(
            FrameActions,
            text="BEST_DEC.",
            font=self.myFont,
            fg="white",
            bg="lightgoldenrod4",
            command=self.open_popupBestDec,
            cursor="hand1",
            justify="center",
        )
        myTip = Hovertip(
            buttonBestDec,
            "Click to get the number of stars by declinaison ranges.",
            hover_delay=1000,
        )

        # Reset button
        buttonReset = Button(
            FrameActions,
            text="RESET",
            font=self.myFont,
            fg="white",
            bg="orange",
            command=self.onReset,
            cursor="hand1",
        )
        myTip = Hovertip(buttonReset, "Click to reset the selection.", hover_delay=1000)

        # Apply Info targets button
        buttonInfoTargets = Button(
            FrameActions,
            text="INFO_TARGETS",
            font=self.myFont,
            fg="white",
            bg="purple",
            command=self.open_popupInfoTargets,
            cursor="hand1",
        )
        myTip = Hovertip(
            buttonInfoTargets,
            "Click to get information about the selected targets.",
            hover_delay=1000,
        )

        # Send_to_Aspro button
        buttonAspro = Button(
            FrameActions,
            text="SEND2ASPRO",
            font=self.myFont,
            fg="white",
            bg="blue",
            command=self.onAspro,
            cursor="hand1",
        )
        myTip = Hovertip(
            buttonAspro,
            "Click to send the selected targets to Aspro2.",
            hover_delay=1000,
        )

        # DRS button
        # buttonDRS = Button(FrameActions, text='FAKE_DRS', font=self.myFont, fg='white', bg='purple',
        #                   command=self.onDRS, cursor='hand1', state=DISABLED)

        # Quit button
        buttonQuit = Button(
            FrameActions,
            text="QUIT",
            font=self.myFont,
            fg="white",
            bg="red",
            command=self.onQuit,
            cursor="X_cursor",
        )
        myTip = Hovertip(
            buttonQuit, "Click to quit the SPICA-NSS tool.", hover_delay=1000
        )

        # Log button
        buttonLog = Button(
            FrameActions,
            text="LOG",
            font=self.myFont,
            fg="black",
            command=self.onLog,
            cursor="hand1",
        )
        myTip = Hovertip(
            buttonLog,
            "Click to get statistics about the ISSP survey.",
            hover_delay=1000,
        )

        """
        # Workpackages
        labelProgName = Button(Frame2, text='Workpackages', relief='groove', font=self.myFont, activebackground = 'green',
                               command=self.open_popupProgName, cursor='hand1')
        """
        """
        # Instrumental modes
        labelInstMode = Button(Frame3, text='Modes', relief='groove', font=self.myFont, activebackground = 'green',
                               command=self.open_popupInstMode, cursor='hand1')
        """
        # Workpackages (checkbuttons)
        labelProgName = Label(
            FrameWorkPackages,
            text="ProgNames",
            font=self.myFont,
            fg="black",
            bg="orange",
            width=15,
        )
        # self.ProgName = ['WP01', 'WP02', 'WP03', 'WP07', 'WP08', 'WP09', 'WP10', 'WP11', 'WP12', 'WP13', 'WP15']
        self.ProgName = ["S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08"]
        self.buttonProgName = []
        self.SelectedProgName = []
        for iProgName in self.ProgName:
            ProgNameId = StringVar()
            ProgNameId.set("0")
            self.buttonProgName.append(
                Checkbutton(
                    FrameWorkPackages,
                    text=iProgName,
                    onvalue=iProgName,
                    offvalue=0,
                    indicatoron=1,
                    variable=ProgNameId,
                    selectcolor="red",
                    activebackground="green",
                    command=self.plotSelectedProgName,
                )
            )
            self.SelectedProgName.append(ProgNameId)
        myTip = Hovertip(
            labelProgName, "Selection of the ISSP programs.", hover_delay=1000
        )

        # Instrumental Modes (checkbuttons)
        labelInstMode = Label(
            FrameInstModes,
            text="Inst. Modes",
            font=self.myFont,
            fg="black",
            bg="orange",
            width=15,
        )
        """
        self.buttonInstMode = []
        self.SelectedInstMode = []
        for iInstMode in range(len(self.InstMode)):
            InstModeId = IntVar()
            InstModeId.set(0)
            self.buttonInstMode.append(Checkbutton(Frame3, onvalue=self.InstMode[iInstMode],
                                              text=iInstMode+1, state='normal',
                                              variable=InstModeId, font=self.myFont,
                                              selectcolor='red', activebackground='green',
                                              command=self.plotSelectedInstMode, cursor='hand1',
                                              ))
            self.SelectedInstMode.append(InstModeId)
            print( self.SelectedInstMode, InstModeId)
        """

        self.InstMode = ["DIA", "DLD", "IMA", "TMP", "SPI"]
        self.buttonInstMode = []
        self.SelectedInstMode = []
        for iInstMode in self.InstMode:
            InstModeId = StringVar()
            InstModeId.set("0")
            self.buttonInstMode.append(
                Checkbutton(
                    FrameInstModes,
                    text=iInstMode,
                    onvalue=iInstMode,
                    offvalue=0,
                    indicatoron=1,
                    variable=InstModeId,
                    selectcolor="red",
                    activebackground="green",
                    command=self.plotSelectedInstMode,
                )
            )
            self.SelectedInstMode.append(InstModeId)

        # Final priority (checkbutton)
        labelFinalPriority = Label(
            FramePriorities,
            text="Priority_final",
            font=self.myFont,
            fg="black",
            bg="orange",
            width=15,
        )
        self.buttonFinalPriority = []
        self.SelectedFinalPriority = []
        for iFinalPriority in range(4):
            FinalPriorityId = IntVar()
            FinalPriorityId.set(0)
            self.buttonFinalPriority.append(
                Checkbutton(
                    FramePriorities,
                    onvalue=iFinalPriority + 1,
                    text=iFinalPriority + 1,
                    state="normal",
                    variable=FinalPriorityId,
                    font=self.myFont,
                    selectcolor="red",
                    activebackground="green",
                    command=self.plotSelectedFinalPriority,
                )
            )
            self.SelectedFinalPriority.append(FinalPriorityId)

        # Declination
        labelDeclination = Label(
            self.FrameObjects,
            text="Declination (deg)",
            font=self.myFont,
            fg="black",
            bg="light blue",
            width=20,
        )

        # Magnitude
        labelMagnitude = Label(
            self.FrameObjects,
            text="Magnitude",
            font=self.myFont,
            fg="black",
            bg="light green",
            width=20,
        )

        # Dec_Min
        labelDecMin = Label(
            self.FrameObjects, text="DEC_MIN", font=self.myFont, width=10
        )
        self.entryDecMin = Entry(
            self.FrameObjects,
            textvariable=self.strDecMin,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryDecMin.insert(END, str(self.decmin))
        self.entryDecMin.bind(
            "<Return>", (lambda _: self.entryDecMinCallback(self.entryDecMin))
        )

        # Dec_Max
        labelDecMax = Label(
            self.FrameObjects, text="DEC_MAX", font=self.myFont, width=10
        )
        self.entryDecMax = Entry(
            self.FrameObjects,
            textvariable=self.strDecMax,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryDecMax.insert(END, str(self.decmax))
        self.entryDecMax.bind(
            "<Return>", (lambda _: self.entryDecMaxCallback(self.entryDecMax))
        )

        # Vmag_Min
        labelVmagMin = Label(
            self.FrameObjects, text="VMAG_MIN", font=self.myFont, width=10
        )
        self.entryVmagMin = Entry(
            self.FrameObjects,
            textvariable=self.strVmagMin,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryVmagMin.insert(END, str(self.vmagmin))
        self.entryVmagMin.bind(
            "<Return>", (lambda _: self.entryVmagMinCallback(self.entryVmagMin))
        )

        # Vmag_Max
        labelVmagMax = Label(
            self.FrameObjects, text="VMAG_MAX", font=self.myFont, width=10
        )
        self.entryVmagMax = Entry(
            self.FrameObjects,
            text=self.strVmagMax,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryVmagMax.insert(END, str(self.vmagmax))
        self.entryVmagMax.bind(
            "<Return>", (lambda _: self.entryVmagMaxCallback(self.entryVmagMax))
        )  # self.labelValVmagMax = Label(self.frame5, text=str(self.vmagmax), fg='red', font=self.myFont)

        # Dec_Mean
        labelDecMean = Label(
            self.FrameObjects, text="<DEC> =", font=self.myFont, width=10
        )
        DecMeanValue = Label(
            self.FrameObjects,
            textvariable=self.strDecMean,
            justify="right",
            font=self.myFont,
        )

        # Vmag_Mean
        labelVmagMean = Label(
            self.FrameObjects, text="<VMAG> =", font=self.myFont, width=10
        )
        VmagMeanValue = Label(
            self.FrameObjects,
            textvariable=self.strVmagMean,
            justify="right",
            font=self.myFont,
        )

        # Add star to send to Aspro2
        buttonAddTarget = Button(
            FrameAddTarget,
            text="Add a star",
            font=self.myFont,
            command=self.open_popupAddTarget,
            cursor="hand1",
        )
        myTip = Hovertip(
            buttonAddTarget, "To add a star in selection for Aspro2.", hover_delay=1000
        )

        # Primary calibrators
        # labelCalPrim = Label(frame6, text='Primary calibrators:', font=self.myFont, fg='black', bg='light blue', width=23)
        # RA range
        labelRaRangePrim = Label(
            self.FrameCalPrims,
            text="RA range (min)",
            font=self.myFont,
            bg="light blue",
            width=20,
        )
        self.entryRaRangePrim = Entry(
            self.FrameCalPrims,
            textvariable=self.strRaRangePrim,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryRaRangePrim.insert(END, str(self.rarangeprim))
        self.entryRaRangePrim.bind(
            "<Return>", (lambda _: self.entryRaRangePrimCallback(self.entryRaRangePrim))
        )

        # DEC range
        labelDecRangePrim = Label(
            self.FrameCalPrims,
            text="DEC range (deg)",
            font=self.myFont,
            bg="light blue",
            width=20,
        )
        self.entryDecRangePrim = Entry(
            self.FrameCalPrims,
            textvariable=self.strDecRangePrim,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryDecRangePrim.insert(END, str(self.decrangeprim))
        self.entryDecRangePrim.bind(
            "<Return>",
            (lambda _: self.entryDecRangePrimCallback(self.entryDecRangePrim)),
        )

        # Vmag range
        labelVmagRangePrim = Label(
            self.FrameCalPrims,
            text="Vmag range",
            font=self.myFont,
            bg="light green",
            width=20,
        )
        self.entryVmagRangePrim = Entry(
            self.FrameCalPrims,
            textvariable=self.strVmagRangePrim,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryVmagRangePrim.insert(END, str(self.vmagrangeprim))
        self.entryVmagRangePrim.bind(
            "<Return>",
            (lambda _: self.entryVmagRangePrimCallback(self.entryVmagRangePrim)),
        )

        # Secondary calibrators
        # labelCalSec = Label(frame7, text='Secondary calibrators:', font=self.myFont, fg='black', bg='light green', width=23)
        # RA range
        labelRaRangeSec = Label(
            self.FrameCalSecs,
            text="RA range (min)",
            font=self.myFont,
            bg="light blue",
            width=20,
        )
        self.entryRaRangeSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strRaRangeSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryRaRangeSec.insert(END, str(self.rarangesec))
        self.entryRaRangeSec.bind(
            "<Return>", (lambda _: self.entryRaRangeSecCallback(self.entryRaRangeSec))
        )

        # DEC range
        labelDecRangeSec = Label(
            self.FrameCalSecs,
            text="DEC range (deg)",
            font=self.myFont,
            bg="light blue",
            width=20,
        )
        self.entryDecRangeSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strDecRangeSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryDecRangeSec.insert(END, str(self.decrangesec))
        self.entryDecRangeSec.bind(
            "<Return>", (lambda _: self.entryDecRangeSecCallback(self.entryDecRangeSec))
        )

        # Vmag range
        labelVmagRangeSec = Label(
            self.FrameCalSecs,
            text="Vmag range",
            font=self.myFont,
            bg="light green",
            width=20,
        )
        self.entryVmagRangeSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strVmagRangeSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryVmagRangeSec.insert(END, str(self.vmagrangesec))
        self.entryVmagRangeSec.bind(
            "<Return>",
            (lambda _: self.entryVmagRangeSecCallback(self.entryVmagRangeSec)),
        )

        # LDD Chi2
        labelLDDChiSec = Label(
            self.FrameCalSecs,
            text="Max. LDD Chi2",
            font=self.myFont,
            bg="wheat1",
            width=20,
        )
        self.entryLDDChiSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strLDDChiSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryLDDChiSec.insert(END, str(self.lddchisec))
        self.entryLDDChiSec.bind(
            "<Return>", (lambda _: self.entryLDDChiSecCallback(self.entryLDDChiSec))
        )

        # Relative Error on diameter
        labelRelErrorSec = Label(
            self.FrameCalSecs,
            text="Max. rel. error (%)",
            font=self.myFont,
            bg="wheat1",
            width=20,
        )
        self.entryRelErrorSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strRelErrorSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryRelErrorSec.insert(END, str(self.relerrorsec))
        self.entryRelErrorSec.bind(
            "<Return>", (lambda _: self.entryRelErrorSecCallback(self.entryRelErrorSec))
        )

        # Visibility limit
        labelMinVisSec = Label(
            self.FrameCalSecs, text="Min. vis2", font=self.myFont, bg="wheat1", width=20
        )
        self.entryMinVisSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strMinVisSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        # self.entryMinVisSec =Entry(self.FrameCalSecs, textvariable=self.intMaxBaseline, justify='right', font=self.myFont, cursor='pencil', width=5)
        self.entryMinVisSec.insert(END, str(self.minvissec))
        self.entryMinVisSec.bind(
            "<Return>", (lambda _: self.entryMinVisSecCallback(self.entryMinVisSec))
        )

        # Baseline limit
        labelMaxBaseline = Label(
            self.FrameCalSecs,
            text="Max. Baseline (m)",
            font=self.myFont,
            bg="wheat1",
            width=20,
        )
        self.entryMaxBaseline = Scale(
            self.FrameCalSecs,
            variable=self.intMaxBaseline,
            orient=HORIZONTAL,
            from_=330,
            to_=30,
            resolution=10,
            tickinterval=300,
        )  # , command= (lambda _:  self.entryMinVisSecCallback(self.entryMinVisSec)))
        self.entryMaxBaseline.set(self.maxbaseline)
        # self.entryMaxBaseline.insert(END, str(self.maxbaseline))
        self.entryMaxBaseline.bind(
            "<ButtonRelease-1>",
            (lambda _: self.entryMinVisSecCallback2(self.entryMaxBaseline)),
        )

        # Log button
        # buttonLog = Button(FrameLog, text='Log', relief='groove', bg='AntiqueWhite3', font=self.myFont, activebackground='green',
        #                       command=self.onLog, cursor='hand1', justify='center')

        # Gridding widgets
        FrameDate.rowconfigure(0, weight=1)
        FrameDate.columnconfigure(0, weight=1)
        FrameDate.grid(row=0, column=0, padx=15, pady=5, ipady=4)

        FrameActions.rowconfigure(0, weight=1)
        FrameActions.columnconfigure(0, weight=15)
        FrameActions.grid(row=0, column=1, pady=5, ipady=1, padx=5)

        FrameWorkPackages.rowconfigure(1, weight=1)
        FrameWorkPackages.columnconfigure(0, weight=1)
        FrameWorkPackages.grid(
            sticky=W, row=1, column=0, columnspan=2, padx=15, ipady=1
        )

        self.FrameObjects.rowconfigure(0, weight=1)
        self.FrameObjects.columnconfigure(0, weight=1)
        self.FrameObjects.grid(
            sticky=W, row=4, column=0, columnspan=2, padx=15, pady=1, ipady=1
        )
        # self.FrameObjects.grid(row=4, column=0, padx=15, pady=5, ipady=4)

        FrameAddTarget.rowconfigure(0, weight=1)
        FrameAddTarget.columnconfigure(0, weight=1)
        FrameAddTarget.grid(sticky=E, row=4, column=1, padx=15, pady=1, ipady=1)

        FrameInstModes.grid(sticky=W, row=2, column=0, columnspan=2, padx=15, ipady=1)
        FramePriorities.grid(
            sticky=W, row=3, column=0, columnspan=2, padx=15, pady=1, ipady=1
        )
        # self.FrameObjects.grid(sticky=W, row=4, column=0, columnspan=2, padx=15, pady=1, ipady=1)
        # self.FrameAddTarget.grid(sticky=N, row=4, column=0, columnspan=2, padx=15, pady=1, ipady=1)
        self.FrameCalPrims.grid(
            sticky=W, row=5, column=0, columnspan=2, padx=15, ipadx=3, pady=1, ipady=1
        )
        self.FrameCalSecs.grid(
            sticky=W, row=6, column=0, columnspan=2, padx=15, ipadx=3, pady=1, ipady=1
        )
        FrameLog.grid(sticky=W, row=8, column=0, columnspan=2, padx=15, pady=1, ipady=1)

        # labelDate.grid(column=0, row=0, padx=5)
        entryDate.grid(column=0, row=0, padx=5)
        self.labelValDate.grid(column=1, row=0, padx=5)

        # l0 = Label(frame1a, text='                      ')
        # l0.grid(column=3, row=iRow)

        buttonQuery.grid(column=0, row=0, padx=5)
        # buttonFilters.grid(column=1, row=0, padx=5)
        buttonBestDec.grid(column=1, row=0, padx=5)
        buttonInfoTargets.grid(column=2, row=0, padx=5)
        buttonAspro.grid(column=3, row=0, padx=5)
        buttonReset.grid(column=4, row=0, padx=5)
        # buttonCalPrim.grid(column=7, row=iRow)
        # buttonCalSec.grid(column=8, row=iRow)
        # buttonDRS.grid(column=4, row=0, padx=5)
        buttonQuit.grid(column=5, row=0, padx=5)
        buttonLog.grid(column=6, row=0, padx=5)

        labelProgName.grid(column=0, row=1, padx=5)
        iRow = 0
        iRow += 1
        labelProgName.grid(column=0, row=0, padx=5)
        for iProgName in list(range(len(self.ProgName))):
            self.buttonProgName[iProgName].grid(column=iProgName + 1, row=0, padx=5)

        labelInstMode.grid(column=0, row=0, padx=5)
        for iInstMode in list(range(len(self.InstMode))):
            self.buttonInstMode[iInstMode].grid(column=iInstMode + 1, row=0, padx=5)

        labelFinalPriority.grid(column=0, row=0, padx=5)
        for iPriority in list(range(4)):
            self.buttonFinalPriority[iPriority].grid(
                column=iPriority + 1, row=0, padx=5
            )

        labelDeclination.grid(column=0, row=0, padx=5)
        labelDecMin.grid(column=1, row=0, padx=5)
        self.entryDecMin.grid(column=2, row=0, padx=5)
        labelDecMax.grid(column=3, row=0, padx=5)
        self.entryDecMax.grid(column=4, row=0, padx=5)
        labelDecMean.grid(column=5, row=0)
        DecMeanValue.grid(column=6, row=0)

        buttonAddTarget.grid(column=7, row=0, padx=5)

        labelMagnitude.grid(column=0, row=1)
        labelVmagMin.grid(column=1, row=1)
        self.entryVmagMin.grid(column=2, row=1)
        labelVmagMax.grid(column=3, row=1)
        self.entryVmagMax.grid(column=4, row=1)
        labelVmagMean.grid(column=5, row=1)
        VmagMeanValue.grid(column=6, row=1)

        # labelCalPrim.grid(column=0, row=iRow, sticky='w', pady=15)
        labelRaRangePrim.grid(column=0, row=0, padx=5)
        self.entryRaRangePrim.grid(column=1, row=0)
        labelDecRangePrim.grid(column=2, row=0, padx=5)
        self.entryDecRangePrim.grid(column=3, row=0)
        labelVmagRangePrim.grid(column=4, row=0, padx=5)
        self.entryVmagRangePrim.grid(column=5, row=0)
        # Button(self.FrameCalPrims, text='SEARCH', font=self.myFont, fg='white', bg = 'green', command=self.goCalPrim, cursor='hand1').grid(column=6, row=0, padx=5)
        Button(
            self.FrameCalPrims,
            text="UNDO",
            font=self.myFont,
            fg="white",
            bg="red",
            command=self.delCalPrim,
            cursor="hand1",
        ).grid(column=6, row=0, padx=5)

        # labelCalSec.grid(column=0, row=iRow, sticky='w')
        labelRaRangeSec.grid(column=0, row=0, padx=5)
        self.entryRaRangeSec.grid(column=1, row=0)
        labelDecRangeSec.grid(column=2, row=0, padx=5)
        self.entryDecRangeSec.grid(column=3, row=0)
        labelVmagRangeSec.grid(column=4, row=0, padx=5)
        self.entryVmagRangeSec.grid(column=5, row=0)

        labelLDDChiSec.grid(column=0, row=1, padx=5)
        self.entryLDDChiSec.grid(column=1, row=1)
        labelRelErrorSec.grid(column=2, row=1, padx=5)
        self.entryRelErrorSec.grid(column=3, row=1)
        labelMinVisSec.grid(column=4, row=1, padx=5)
        self.entryMinVisSec.grid(column=5, row=1)
        labelMaxBaseline.grid(column=4, row=2, padx=5)
        self.entryMaxBaseline.grid(column=5, row=2, padx=5)
        # Button(self.FrameCalSecs, text='SEARCH', font=self.myFont, fg='white', bg = 'green', command=self.goCalSec, cursor='hand1').grid(column=6, row=0, rowspan=1, padx=5)
        Button(
            self.FrameCalSecs,
            text="UNDO",
            font=self.myFont,
            fg="white",
            bg="red",
            command=self.delCalSec,
            cursor="hand1",
        ).grid(column=6, row=0, rowspan=1, padx=5)

        # buttonLog.grid(column=0, row=0)

        # self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.root.mainloop()

    def open_popupAddTarget(self):
        # Check if the SPICA catalog has been queried
        if self.spica_catg is None:
            messagebox.showinfo(title="Info", message="Please query SPICA-DB first.")
            return

        # Prompt the user to input the target name
        name = simpledialog.askstring("Add a target", "Enter your target name:")
        # If no name is entered, show a message and exit
        if not name:
            messagebox.showinfo(title="Info", message="No name entered.")
            return

        join_same_target_angle = 5 * u.arcsec

        # Attempt to query Simbad for the target
        try:
            t_simbad = Simbad.query_object(name)
            t_simbad.rename_columns(
                t_simbad.colnames, [col.lower() for col in t_simbad.colnames]
            )
            self.name_simbad = t_simbad["main_id"][0]
            coo_new_target = Table(
                [SkyCoord(t_simbad["ra"], t_simbad["dec"], unit=(u.hourangle, u.deg))],
                names=["sc"],
            )
        except Exception as exc:
            messagebox.showerror(
                title="Error",
                message=f"Error querying Simbad for {name}.",
                detail=f"Check the spelling. ({exc})",
            )
            return

        # Compare the new target coordinates with existing targets in the SPICA catalog
        coo_existing_targets = Table(
            [SkyCoord(self.spica_catg["ra"], self.spica_catg["dec"], unit="deg")],
            names=["sc"],
        )
        coo_existing_targets["index"] = np.arange(len(coo_existing_targets))
        tables_joined = join(
            coo_new_target,
            coo_existing_targets,
            join_funcs={"sc": join_skycoord(join_same_target_angle)},
        )

        # If the target is found in the SPICA catalog
        if len(tables_joined) > 0:
            # Check if the target is already selected
            target_main_id = self.spica_catg["target_main_id"][
                tables_joined["index"][0]
            ]

            if (
                target_main_id
                in self.spica_catg["target_main_id"][self.indexList_Targets]
            ):
                messagebox.showwarning(
                    title="Warning",
                    message=f"{name} is already in selection.\n(Simbad Main-ID: {self.name_simbad})",
                )
            else:
                # Confirm the addition of the target
                confirm = messagebox.askyesno(
                    title="Confirmation",
                    message=f"Do you confirm to add {name}?\n(Simbad Main-ID: {self.name_simbad})",
                )

                if confirm:

                    if len(tables_joined["index"]) == 1:
                        j = tables_joined["index"][0]
                        if self.index_AddTarget is None:
                            self.index_AddTarget = [j]
                        elif j not in self.index_AddTarget:
                            self.index_AddTarget.append(j)

                        self.getSelectedTargets()
                        self.plot_radec()

                        if self.popupInfoTargets:
                            self.clear_all(self.my_treeInfoTargets)
                            self.insert_popupInfoTargets()
                            self.tree_frameInfoTargets.title(
                                f"List of selected targets ({len(self.indexList_Targets)} objects)"
                            )

                    else:
                        self.tree_frameAddTarget = Toplevel(self.root)
                        tree_scroll = Scrollbar(self.tree_frameAddTarget)
                        tree_scroll.pack(side=RIGHT, fill=Y)
                        self.tree_frameAddTarget.title(
                            f"Please double-click to select the ProgName for {self.name_simbad}"
                        )

                        # Create the Treeview
                        self.lenAddTarget = len(tables_joined["index"])
                        self.my_treeAddTarget = ttk.Treeview(
                            self.tree_frameAddTarget,
                            yscrollcommand=tree_scroll.set,
                            selectmode="extended",
                            height=self.lenAddTarget,
                        )
                        self.my_treeAddTarget.pack(expand=True, fill="y")

                        # Configure the Scrollbar
                        tree_scroll.config(command=self.my_treeAddTarget.yview)

                        # Define Our Columns
                        self.my_treeAddTarget["columns"] = (
                            "SPICA-DB ID",
                            "Target Main ID",
                            "Progname",
                            "Spica Mode",
                        )

                        # Format Our Columns
                        self.my_treeAddTarget.column("#0", width=0, stretch=NO)
                        self.my_treeAddTarget.column(
                            "SPICA-DB ID", anchor=CENTER, width=100
                        )
                        self.my_treeAddTarget.column(
                            "Target Main ID", anchor=CENTER, width=150
                        )
                        self.my_treeAddTarget.column(
                            "Progname", anchor=CENTER, width=150
                        )
                        self.my_treeAddTarget.column(
                            "Spica Mode", anchor=CENTER, width=150
                        )

                        # Create Headings
                        self.my_treeAddTarget.heading("#0", text="", anchor=W)
                        self.my_treeAddTarget.heading(
                            "SPICA-DB ID", text="SPICA-DB ID", anchor=CENTER
                        )
                        self.my_treeAddTarget.heading(
                            "Target Main ID", text="Target Main ID", anchor=CENTER
                        )
                        self.my_treeAddTarget.heading(
                            "Progname", text="Progname", anchor=CENTER
                        )
                        self.my_treeAddTarget.heading(
                            "Spica Mode", text="Spica Mode", anchor=CENTER
                        )

                        # Create Striped Row Tags
                        self.my_treeAddTarget.tag_configure(
                            "oddrow", background="white"
                        )
                        self.my_treeAddTarget.tag_configure(
                            "evenrow", background="lightblue"
                        )
                        self.my_treeAddTarget.focus()

                        count = 0
                        for k in np.sort(tables_joined["index"]):
                            row = (
                                self.spica_catg["spicadb_id"][k],
                                self.spica_catg["target_main_id"][k],
                                self.spica_catg["progname"][k],
                                self.spica_catg["spica_mode"][k],
                            )

                            if count % 2 == 0:
                                (
                                    self.my_treeAddTarget.insert(
                                        parent="",
                                        index="end",
                                        iid=count,
                                        text="",
                                        values=row,
                                        tags=("evenrow",),
                                    )
                                )
                            else:
                                (
                                    self.my_treeAddTarget.insert(
                                        parent="",
                                        index="end",
                                        iid=count,
                                        text="",
                                        values=row,
                                        tags=("oddrow",),
                                    )
                                )
                            count += 1

                            # self.my_treeAddTarget.insert(parent="", index="end", text="", values=row)

                        self.my_treeAddTarget.bind("<Double-1>", self.OnDoubleClick)

        else:
            # Query the full SPICA database if the target is not found initially
            all_spica_catg = self.dbquery_tap()
            found_in_db = False
            name_lower = name.lower().replace(" ", "")
            for j, starname in enumerate(all_spica_catg["target_main_id"]):
                if name_lower in starname.lower().replace(" ", ""):
                    print("yes", j, starname, name)
                    messagebox.showwarning(
                        title="Warning",
                        message=f"{name_lower} is not observable on {self.date}.",
                    )
                    found_in_db = True
                    break

            if not found_in_db:
                messagebox.showwarning(
                    title="Warning", message=f"{name} is not in SPICA-DB."
                )

    # def open_popupAddTarget2(self):
    #     """

    #     Returns
    #     -------
    #     None.

    #     """
    #     if self.popupAddTarget:
    #         self.tree_frameAddTarget.destroy()
    #         self.popupAddTarget = True
    #     else:
    #         self.popupAddTarget = True

    #     self.tree_frameAddTarget = Toplevel(self.root)
    #     tree_scroll = Scrollbar(self.tree_frameAddTarget)
    #     tree_scroll.pack(side=RIGHT, fill=Y)
    #     self.tree_frameAddTarget.title(f'List of possible fainter stars to include ({len(self.indexList_AddTarget)} objects)')

    #     # Create the Treeview
    #     if len(self.indexList_AddTarget) < 10:
    #         self.lenAddTarget = len(self.indexList_AddTarget)
    #     else:
    #         self.lenAddTarget = 10
    #     self.my_treeAddTarget = ttk.Treeview(self.tree_frameAddTarget, yscrollcommand=tree_scroll.set, selectmode="extended", height=self.lenAddTarget)
    #     self.my_treeAddTarget.pack(expand=True, fill="y")

    #     # Configure the Scrollbar
    #     tree_scroll.config(command=self.my_treeAddTarget.yview)

    #     # Define Our Columns
    #     self.my_treeAddTarget["columns"] = ("SPICA-DB ID", "Target Main ID",
    #                                           "Spec. Type", "Progname", "Spica Mode",
    #                                           "Final Priority", "Completion Rate", "Ra", "Dec", "Diameter", "Vmag", "Hmag")

    #     # Format Our Columns
    #     self.my_treeAddTarget.column("#0", width=0, stretch=NO)
    #     self.my_treeAddTarget.column("SPICA-DB ID", anchor=CENTER, width=100)
    #     self.my_treeAddTarget.column("Target Main ID", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Spec. Type", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Progname", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Spica Mode", anchor=CENTER, width=150)
    #     #self.my_treeInfoTargets.column("Priority PI", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Final Priority", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Completion Rate", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Ra", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Dec", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Diameter", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Vmag", anchor=CENTER, width=150)
    #     self.my_treeAddTarget.column("Hmag", anchor=CENTER, width=150)

    #     # Create Headings
    #     self.my_treeAddTarget.heading("#0", text="", anchor=W)
    #     self.my_treeAddTarget.heading("SPICA-DB ID", text="SPICA-DB ID", anchor=CENTER)
    #     self.my_treeAddTarget.heading("Target Main ID", text="Target Main ID", anchor=CENTER)
    #     self.my_treeAddTarget.heading("Spec. Type", text="Spec. Type", anchor=CENTER)
    #     self.my_treeAddTarget.heading("Progname", text="Progname", anchor=CENTER)
    #     self.my_treeAddTarget.heading("Spica Mode", text="Spica Mode", anchor=CENTER)
    #     #self.my_treeInfoTargets.heading("Priority PI", text="Priority PI", anchor=CENTER)
    #     self.my_treeAddTarget.heading("Final Priority", text="Final Priority", anchor=CENTER)
    #     self.my_treeAddTarget.heading("Completion Rate", text="Completion Rate", anchor=CENTER)
    #     #self.my_treeInfoTargets.heading("Dec", text="Dec", anchor=CENTER)
    #     self.my_treeAddTarget.heading("Ra", text="Ra", anchor=CENTER, command=lambda _col="Ra": \
    #                                    self.treeview_sort_column(self.my_treeAddTarget, _col, False))
    #     self.my_treeAddTarget.heading("Dec", text="Dec", anchor=CENTER, command=lambda _col="Dec": \
    #                                    self.treeview_sort_column(self.my_treeAddTarget, _col, False))
    #     #for col in self.my_treeInfoTargets["columns"]:
    #     #    self.my_treeInfoTargets.heading(col, text=col, anchor=CENTER, command=lambda _col=col: \
    #     #                                    self.treeview_sort_column(self.my_treeInfoTargets, _col, False))
    #     self.my_treeAddTarget.heading("Diameter", text="Diameter", anchor=CENTER, command=lambda _col="Diameter": \
    #                                    self.treeview_sort_column(self.my_treeAddTarget, _col, False))
    #     self.my_treeAddTarget.heading("Vmag", text="Vmag", anchor=CENTER, command=lambda _col="Vmag": \
    #                                    self.treeview_sort_column(self.my_treeAddTarget, _col, False))
    #     self.my_treeAddTarget.heading("Hmag", text="Hmag", anchor=CENTER, command=lambda _col="Hmag": \
    #                                    self.treeview_sort_column(self.my_treeAddTarget, _col, False))
    #     #self.my_treeInfoTargets.heading("Vmag", text="Vmag", anchor=CENTER)

    #     # Insert AddTargets in popup
    #     self.insert_popupAddTarget()

    #     #
    #     self.tree_frameAddTarget.protocol('WM_DELETE_WINDOW', self.closeframeAddTarget)

    def open_popupInfoTargets(self):
        """


        Returns
        -------
        None.

        """
        if self.spica_catg is None:
            messagebox.showinfo(title="Info", message="Please query SPICA-DB first.")
            return

        if self.popupInfoTargets:
            self.tree_frameInfoTargets.destroy()
            self.popupInfoTargets = True
        else:
            self.popupInfoTargets = True

        self.tree_frameInfoTargets = Toplevel(self.root)
        tree_scroll = Scrollbar(self.tree_frameInfoTargets)
        tree_scroll.pack(side=RIGHT, fill=Y)
        self.tree_frameInfoTargets.title(
            f"List of selected targets ({len(self.indexList_Targets)} objects)"
        )

        # Create the Treeview
        if len(self.indexList_Targets) < 10:
            self.lenInfoTargets = len(self.indexList_Targets)
        else:
            self.lenInfoTargets = 10
        self.my_treeInfoTargets = ttk.Treeview(
            self.tree_frameInfoTargets,
            yscrollcommand=tree_scroll.set,
            selectmode="extended",
            height=self.lenInfoTargets,
        )
        self.my_treeInfoTargets.pack(expand=True, fill="y")

        # Configure the Scrollbar
        tree_scroll.config(command=self.my_treeInfoTargets.yview)

        # Define Our Columns
        self.my_treeInfoTargets["columns"] = (
            "SPICA-DB ID",
            "Target Main ID",
            "Spec. Type",
            "Progname",
            "Spica Mode",
            "Final Priority",
            "Completion Rate",
            "Ra",
            "Dec",
            "Diameter",
            "Vmag",
            "Hmag",
        )

        # Format Our Columns
        self.my_treeInfoTargets.column("#0", width=0, stretch=NO)
        self.my_treeInfoTargets.column("SPICA-DB ID", anchor=CENTER, width=100)
        self.my_treeInfoTargets.column("Target Main ID", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Spec. Type", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Progname", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Spica Mode", anchor=CENTER, width=150)
        # self.my_treeInfoTargets.column("Priority PI", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Final Priority", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Completion Rate", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Ra", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Dec", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Diameter", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Vmag", anchor=CENTER, width=150)
        self.my_treeInfoTargets.column("Hmag", anchor=CENTER, width=150)

        # Create Headings
        self.my_treeInfoTargets.heading("#0", text="", anchor=W)
        self.my_treeInfoTargets.heading(
            "SPICA-DB ID", text="SPICA-DB ID", anchor=CENTER
        )
        self.my_treeInfoTargets.heading(
            "Target Main ID", text="Target Main ID", anchor=CENTER
        )
        self.my_treeInfoTargets.heading("Spec. Type", text="Spec. Type", anchor=CENTER)
        self.my_treeInfoTargets.heading("Progname", text="Progname", anchor=CENTER)
        self.my_treeInfoTargets.heading("Spica Mode", text="Spica Mode", anchor=CENTER)
        # self.my_treeInfoTargets.heading("Priority PI", text="Priority PI", anchor=CENTER)
        self.my_treeInfoTargets.heading(
            "Final Priority", text="Final Priority", anchor=CENTER
        )
        self.my_treeInfoTargets.heading(
            "Completion Rate", text="Completion Rate", anchor=CENTER
        )
        # self.my_treeInfoTargets.heading("Dec", text="Dec", anchor=CENTER)
        self.my_treeInfoTargets.heading(
            "Completion Rate",
            text="Completion Rate",
            anchor=CENTER,
            command=lambda _col="Completion Rate": self.treeview_sort_column(
                self.my_treeInfoTargets, _col, False
            ),
        )
        self.my_treeInfoTargets.heading(
            "Ra",
            text="Ra",
            anchor=CENTER,
            command=lambda _col="Ra": self.treeview_sort_column(
                self.my_treeInfoTargets, _col, False
            ),
        )
        self.my_treeInfoTargets.heading(
            "Dec",
            text="Dec",
            anchor=CENTER,
            command=lambda _col="Dec": self.treeview_sort_column(
                self.my_treeInfoTargets, _col, False
            ),
        )
        # for col in self.my_treeInfoTargets["columns"]:
        #    self.my_treeInfoTargets.heading(col, text=col, anchor=CENTER, command=lambda _col=col: \
        #                                    self.treeview_sort_column(self.my_treeInfoTargets, _col, False))
        self.my_treeInfoTargets.heading(
            "Diameter",
            text="Diameter",
            anchor=CENTER,
            command=lambda _col="Diameter": self.treeview_sort_column(
                self.my_treeInfoTargets, _col, False
            ),
        )
        self.my_treeInfoTargets.heading(
            "Vmag",
            text="Vmag",
            anchor=CENTER,
            command=lambda _col="Vmag": self.treeview_sort_column(
                self.my_treeInfoTargets, _col, False
            ),
        )
        self.my_treeInfoTargets.heading(
            "Hmag",
            text="Hmag",
            anchor=CENTER,
            command=lambda _col="Hmag": self.treeview_sort_column(
                self.my_treeInfoTargets, _col, False
            ),
        )
        # self.my_treeInfoTargets.heading("Vmag", text="Vmag", anchor=CENTER)

        # Insert InfoTargets in popup
        self.insert_popupInfoTargets()

        #
        self.tree_frameInfoTargets.protocol(
            "WM_DELETE_WINDOW", self.closeframeInfoTargets
        )

    def closeframeInfoTargets(self):
        self.tree_frameInfoTargets.destroy()
        self.popupInfoTargets = False

    def closeframeAddTarget(self):
        self.tree_frameAddTarget.destroy()
        self.popupAddTarget = False

    # Define a function to clear all the items present in Treeview
    def clear_all(self, tree):
        for item in tree.get_children():
            tree.delete(item)

    def copy(self, event):
        sel = self.my_treeInfoTargets.selection()  # get selected items
        self.root.clipboard_clear()  # clear clipboard
        # copy headers
        headings = [
            self.my_treeInfoTargets.heading("#{}".format(i), "text")
            for i in range(len(self.my_treeInfoTargets.cget("columns")) + 1)
        ]
        self.root.clipboard_append("\t".join(headings) + "\n")
        for item in sel:
            # retrieve the values of the row
            values = [self.my_treeInfoTargets.item(item, "text")]
            values.extend(self.my_treeInfoTargets.item(item, "values"))
            # append the values separated by \t to the clipboard
            self.root.clipboard_append("\t".join(values) + "\n")

    def insert_popupAddTarget(self):
        # Create Striped Row Tags
        self.my_treeAddTarget.tag_configure("oddrow", background="white")
        self.my_treeAddTarget.tag_configure("evenrow", background="lightsalmon1")

        self.my_treeAddTarget.focus()

        count = 0
        targetInfoList = self.spica_catg[self.indexList_AddTarget]
        targetInfoList.sort("vmag")

        for record in targetInfoList:
            model = record["model"]
            model = model.replace("}{", "},{")
            modeltype = json.loads(model)[0]["type"]
            if modeltype == "disk":
                modeldiam = json.loads(model)[0]["diameter"]
            elif modeltype == "elong_disk":
                modeldiam = json.loads(model)[0]["minor_axis_diameter"]

            ra = record["ra"] / 15.0
            if self.ra_sunrise < self.ra_sunset:
                midnight_offset = 12
            else:
                midnight_offset = 0

            if ra < self.ra_sunset / 15:
                ra = ra + midnight_offset
            if ra > self.ra_sunset / 15:
                ra = ra - midnight_offset

            # if  ma.is_masked(record['hmag']):
            #    record['hmag'] = ''
            #    pdb.set_trace()

            if not ma.is_masked(record["priority_pi2"]):
                flag_completion2 = self.update_flag_completion(
                    record["completion_rate"], record["progname2"]
                )
                priority_final2 = self.update_priority_final2(
                    flag_completion2, record["priority_pi2"], record["progname2"]
                )

            if not record["progname2"]:
                row = (
                    record["spicadb_id"],
                    record["target_main_id"],
                    record["spt"],
                    record["progname"],
                    record["spica_mode"],
                    record["priority_final"],
                    str(record["completion_rate"]),
                    round(record["ra"], 2),
                    round(record["dec"], 2),
                    str(round(modeldiam, 2)),
                    str(round(record["vmag"], 2)),
                    str(record["hmag"]),
                )
            else:
                row = (
                    record["spicadb_id"],
                    record["target_main_id"],
                    record["spt"],
                    record["progname"] + " (" + record["progname2"] + ")",
                    record["spica_mode"],
                    str(record["priority_final"]) + " (" + str(priority_final2) + ")",
                    str(record["completion_rate"]),
                    round(record["ra"], 2),
                    round(record["dec"], 2),
                    str(round(modeldiam, 2)),
                    str(round(record["vmag"], 2)),
                    str(record["hmag"]),
                )

            if count % 2 == 0:
                (
                    self.my_treeAddTarget.insert(
                        parent="",
                        index="end",
                        iid=count,
                        text="",
                        values=row,
                        tags=("evenrow",),
                    )
                )
            else:
                (
                    self.my_treeAddTarget.insert(
                        parent="",
                        index="end",
                        iid=count,
                        text="",
                        values=row,
                        tags=("oddrow",),
                    )
                )
            count += 1

        self.my_treeAddTarget.bind("<Control-c>", self.copy)
        self.my_treeAddTarget.bind("<Double-1>", self.OnDoubleClick)

    def insert_popupInfoTargets(self):
        # Create Striped Row Tags
        self.my_treeInfoTargets.tag_configure("oddrow", background="white")
        self.my_treeInfoTargets.tag_configure("evenrow", background="lightblue")

        self.my_treeInfoTargets.focus()

        count = 0
        targetInfoList = self.spica_catg[self.indexList_Targets]
        targetInfoList.sort("vmag")

        for record in targetInfoList:
            model = record["model"]
            model = model.replace("}{", "},{")
            modeltype = json.loads(model)[0]["type"]
            if modeltype == "disk":
                modeldiam = json.loads(model)[0]["diameter"]
            elif modeltype == "elong_disk":
                modeldiam = json.loads(model)[0]["minor_axis_diameter"]

            ra = record["ra"] / 15.0
            if self.ra_sunrise < self.ra_sunset:
                midnight_offset = 12
            else:
                midnight_offset = 0

            if ra < self.ra_sunset / 15:
                ra = ra + midnight_offset
            if ra > self.ra_sunset / 15:
                ra = ra - midnight_offset

            # if  ma.is_masked(record['hmag']):
            #    record['hmag'] = ''
            #    pdb.set_trace()

            if not ma.is_masked(record["priority_pi2"]):
                flag_completion2 = self.update_flag_completion(
                    record["completion_rate"], record["progname2"]
                )
                priority_final2 = self.update_priority_final2(
                    flag_completion2, record["priority_pi2"], record["progname2"]
                )

            if not record["progname2"]:
                row = (
                    record["spicadb_id"],
                    record["target_main_id"],
                    record["spt"],
                    record["progname"],
                    record["spica_mode"],
                    record["priority_final"],
                    str(record["completion_rate"]),
                    round(record["ra"], 2),
                    round(record["dec"], 2),
                    str(round(modeldiam, 2)),
                    str(round(record["vmag"], 2)),
                    str(record["hmag"]),
                )
            else:
                row = (
                    record["spicadb_id"],
                    record["target_main_id"],
                    record["spt"],
                    record["progname"] + " (" + record["progname2"] + ")",
                    record["spica_mode"],
                    str(record["priority_final"]) + " (" + str(priority_final2) + ")",
                    str(record["completion_rate"]),
                    round(record["ra"], 2),
                    round(record["dec"], 2),
                    str(round(modeldiam, 2)),
                    str(round(record["vmag"], 2)),
                    str(record["hmag"]),
                )

            if count % 2 == 0:
                (
                    self.my_treeInfoTargets.insert(
                        parent="",
                        index="end",
                        iid=count,
                        text="",
                        values=row,
                        tags=("evenrow",),
                    )
                )
            else:
                (
                    self.my_treeInfoTargets.insert(
                        parent="",
                        index="end",
                        iid=count,
                        text="",
                        values=row,
                        tags=("oddrow",),
                    )
                )
            count += 1

        self.my_treeInfoTargets.bind("<Control-c>", self.copy)
        self.my_treeInfoTargets.bind("<Double-1>", self.OnDoubleClick)

    def OnDoubleClick(self, event):
        selected_record = self.my_treeAddTarget.item(self.my_treeAddTarget.focus())
        selected_spicadbid = selected_record["values"][0]
        selected_targetmainid = selected_record["values"][1]
        selected_progname = selected_record["values"][2]
        message = (
            f"Do you want to add {selected_targetmainid} from {selected_progname}?"
        )
        confirm = messagebox.askquestion(title="Confirmation", message=message)

        if confirm == "yes":
            j = list(self.spica_catg["spicadb_id"]).index(selected_spicadbid)

            if self.index_AddTarget is None:
                self.index_AddTarget = [j]
            elif j not in self.index_AddTarget:
                self.index_AddTarget.append(j)

            self.getSelectedTargets()
            self.plot_radec()

            if self.popupInfoTargets:
                self.clear_all(self.my_treeInfoTargets)
                self.insert_popupInfoTargets()
                self.tree_frameInfoTargets.title(
                    f"List of selected targets ({len(self.indexList_Targets)} objects)"
                )

        self.tree_frameAddTarget.destroy()

    def retag(self, theTreeToSort):
        tag = "oddrow"
        for iid in theTreeToSort.get_children(""):
            tag = "oddrow" if tag == "evenrow" else "evenrow"
            theTreeToSort.item(iid, tags=(tag,))

    def treeview_sort_column(self, tv, col, reverse):
        if col == "Count":
            items = [(int(tv.set(k, col)), k) for k in tv.get_children("")]
        else:
            items = [(tv.set(k, col), k) for k in tv.get_children("")]
        items.sort(reverse=reverse)

        # rearrange items in sorted positions
        for index, (val, k) in enumerate(items):
            tv.move(k, "", index)

        # retag the treeview items
        self.retag(tv)  # self.my_treeInfoTargets)

        # reverse sort next time
        tv.heading(col, command=lambda: self.treeview_sort_column(tv, col, not reverse))

    def open_popupBestDec(self):
        """


        Returns
        -------
        None.

        """
        if self.spica_catg is None:
            messagebox.showinfo(title="Info", message="Please query SPICA-DB first.")
            return

        if self.popupBestDec:
            self.tree_frameBestDec.destroy()
            self.popupBestDec = True
        else:
            self.popupBestDec = True
        """
        # Query all spica-db catalog content
        if self.ra_sunrise < self.ra_sunset:
            ra_min = self.ra_sunset-180
            ra_max = self.ra_sunrise+180
        else:
            ra_min = self.ra_sunset
            ra_max = self.ra_sunrise
        """
        dec = self.targetListInit["dec"]
        vmag = self.targetListInit["vmag"]
        decSum = []
        decTotal = 0
        for i in range(-30, 90, 5):
            decSum.append(
                np.size(
                    np.where(
                        (dec > i)
                        & (dec <= (i + 5))
                        & (vmag > float(self.strVmagMin.get()))
                        & (vmag < float(self.strVmagMax.get()))
                    )
                )
            )
            decTotal = decTotal + np.size(
                np.where(
                    (dec > i)
                    & (dec <= (i + 5))
                    & (vmag > float(self.strVmagMin.get()))
                    & (vmag < float(self.strVmagMax.get()))
                )
            )

        decLow = np.arange(-30, 90, 5)
        decHigh = np.arange(-25, 95, 5)
        decBestRanges = [list(decLow), list(decHigh), list(decSum)]

        # adqlQuery = ("SELECT round(s.dec / 5, 0)*5, round(s.dec / 5, 0)*5+5, count (*) FROM " + self.catalogSpica + " s where s.vmag < " + self.strVmagMean.get() +
        #             " and s.ra between " + str(ra_min) +" and " + str(ra_max) +
        #             " group by round(s.dec / 5, 0)*5 " +
        #             " order by 3 desc")
        # service = vo.dal.TAPService(oidbTapUrl)
        # results = service.search(adqlQuery)

        # # store new information into data variable
        # data = results.to_table()

        # # Rename colums and sort
        # data.rename_column("mult","dec_low")
        # data.rename_column("sum","dec_high")
        # data.rename_column("count_all","count")
        # data.sort('dec_low')

        self.tree_frameBestDec = Toplevel(self.root)
        self.tree_frameBestDec.title(
            f"# objects per DecRange ({self.strVmagMin.get()} < Vmag < {self.strVmagMax.get()})"
        )
        tree_scroll = Scrollbar(self.tree_frameBestDec)
        tree_scroll.pack(side=RIGHT, fill=Y)

        # Create the Treeview
        self.my_treeBestDec = ttk.Treeview(
            self.tree_frameBestDec,
            yscrollcommand=tree_scroll.set,
            selectmode="extended",
            height=len(decLow),
        )
        self.my_treeBestDec.pack(expand=True, fill="y")

        # Configure the Scrollbar
        tree_scroll.config(command=self.my_treeBestDec.yview)

        # Define Our Columns
        self.my_treeBestDec["columns"] = ("Dec_low", "Dec_high", "Count")

        # Format Our Columns
        self.my_treeBestDec.column("#0", width=0, stretch=NO)
        self.my_treeBestDec.column("Dec_low", anchor=CENTER, width=100)
        self.my_treeBestDec.column("Dec_high", anchor=CENTER, width=150)
        self.my_treeBestDec.column("Count", anchor=CENTER, width=150)

        # Create Headings
        self.my_treeBestDec.heading("#0", text="", anchor=W)
        self.my_treeBestDec.heading(
            "Dec_low", text="Dec. Low (\N{DEGREE SIGN})", anchor=CENTER
        )
        self.my_treeBestDec.heading(
            "Dec_high", text="Dec. High (\N{DEGREE SIGN})", anchor=CENTER
        )
        # self.my_treeBestDec.heading("Count", text=f"Total ({str(decTotal)} objects)", anchor=CENTER)
        self.my_treeBestDec.heading(
            "Count",
            text=f"Total ({str(decTotal)} objects)",
            anchor=CENTER,
            command=lambda _col="Count": self.treeview_sort_column(
                self.my_treeBestDec, _col, False
            ),
        )

        # Create Striped Row Tags
        self.my_treeBestDec.tag_configure("oddrow", background="white")
        self.my_treeBestDec.tag_configure("evenrow", background="lightblue")

        self.my_treeBestDec.focus()

        count = 0
        for i in range(np.shape(decBestRanges)[1]):
            row = (decBestRanges[0][i], decBestRanges[1][i], decBestRanges[2][i])

            if count % 2 == 0:
                (
                    self.my_treeBestDec.insert(
                        parent="",
                        index="end",
                        iid=count,
                        text="",
                        values=row,
                        tags=("evenrow",),
                    )
                )
            else:
                (
                    self.my_treeBestDec.insert(
                        parent="",
                        index="end",
                        iid=count,
                        text="",
                        values=row,
                        tags=("oddrow",),
                    )
                )
            count += 1

    def plotSelectedProgName(self):
        """

        Returns
        -------
        None.

        """
        self.indexProgName = []
        self.indexProgName2 = []
        for j in self.SelectedProgName:
            progname = j.get()
            if progname != str(0):
                # TODO build a filter function to apply on the table atht will split every progname to search for a match
                self.indexProgName += list(
                    filter(
                        lambda x: progname in self.spica_catg["progname"][x].split("-"),
                        range(len(self.spica_catg)),
                    )
                )
                self.indexProgName2 += list(
                    filter(
                        lambda x: progname
                        in self.spica_catg["progname2"][x].split("-"),
                        range(len(self.spica_catg)),
                    )
                )
        """
        self.indexList_Targets = reduce(np.intersect1d, [self.indexProgName, self.indexInstMode, self.indexFinalPriority,
                                                          self.indexDecMin, self.indexDecMax,
                                                          self.indexVmagMin, self.indexVmagMax])
        """

        if self.indexProgName2:
            self.indexProgName = self.indexProgName + self.indexProgName2

        self.getSelectedTargets()
        self.getAddTarget()
        self.plot_radec()

        # Refresh treeview items
        if self.popupInfoTargets:
            self.clear_all(self.my_treeInfoTargets)
            self.insert_popupInfoTargets()
            self.tree_frameInfoTargets.title(
                f"List of selected targets ({len(self.indexList_Targets)} objects)"
            )

        if self.popupAddTarget:
            self.clear_all(self.my_treeAddTarget)
            self.insert_popupAddTarget()
            self.tree_frameAddTarget.title(
                f"List of possible fainter stars to include ({len(self.indexList_AddTarget)} objects)"
            )

    """
    def state_buttonFinalPriority(self, nb_pId, inb_pId):
        if nb_pId != 0:
            self.buttonFinalPriority[inb_pId]['state']='normal'
            self.buttonFinalPriority[inb_pId]['onvalue']=inb_pId+1
        else:
            self.buttonFinalPriority[inb_pId]['state']='normal'
            self.buttonFinalPriority[inb_pId]['onvalue']=0
    """

    def plot_radec(self):  # , indexList, indexList1, indexList2):
        """
        Plots Right Ascension (RA) and Declination (Dec) of observable targets with different priorities and modes.

        Returns
        -------
        None.
        """
        list_targets = self.spica_catg[self.indexList_Targets]
        list_targets = list_targets[
            list(
                filter(
                    lambda x: list_targets["completion_rate"][x] < 1,
                    range(len(list_targets)),
                )
            )
        ]

        self.DecMean = np.median(list_targets["dec"])
        self.VmagMean = np.median(list_targets["vmag"])

        if np.size(self.indexList_Targets) > 0:
            self.strDecMean.set(round(np.median(list_targets["dec"]), 2))
            self.strVmagMean.set(round(np.median(list_targets["vmag"]), 2))
        elif np.size(self.indexList_Targets) == 0:
            self.strDecMean.set(0)
            self.strVmagMean.set(0)

        model = list_targets["model"]
        modeltype = []
        modeldiam = []
        for j in range(len(list_targets)):
            a = model[j].replace("}{", "},{")
            # print(list_targets['spicadb_id'][j], list_targets['target_main_id'][j], list_targets['piname'][j])
            modeltype.append(json.loads(a)[0]["type"])
            if modeltype[j] == "disk":
                modeldiam.append(json.loads(a)[0]["diameter"])
            elif modeltype[j] == "elong_disk":
                modeldiam.append(json.loads(a)[0]["minor_axis_diameter"])
        modeldiam = MaskedColumn(modeldiam)
        list_targets["modeldiam"] = modeldiam

        fig, self.ax = plt.subplots(
            3, figsize=(11, 6), sharex=True
        )  # Figure(figsize = (6.5,4), dpi = 100) sharey=True, sharex=True)
        self.nbtotal_p1 = np.count_nonzero([(list_targets["priority_final"] == 1)])
        self.nbtotal_p2 = np.count_nonzero([(list_targets["priority_final"] == 2)])
        self.nbtotal_p3 = np.count_nonzero([(list_targets["priority_final"] == 3)])
        self.nbtotal_p4 = np.count_nonzero([(list_targets["priority_final"] == 4)])
        self.ax[0].set_title(
            str(np.shape(list_targets)[0])
            + " ("
            + str(self.nbtotal_p1)
            + "/"
            + str(self.nbtotal_p2)
            + "/"
            + str(self.nbtotal_p3)
            + "/"
            + str(self.nbtotal_p4)
            + ")"
            + " observable targets on "
            + self.date
            + " at CHARA"
        )

        print(
            "[INFO] #P1:",
            self.nbtotal_p1,
            ", #P2:",
            self.nbtotal_p2,
            ", #P3:",
            self.nbtotal_p3,
            ", #P4:",
            self.nbtotal_p4,
        )
        # self.state_buttonFinalPriority(nb_p1, 0)
        # self.state_buttonFinalPriority(nb_p2, 1)
        # self.state_buttonFinalPriority(nb_p3, 2)
        """
        if nb_p1 != 0:
            self.buttonFinalPriority[0]['state']='normal'
        else:
            #self.buttonFinalPriority[0]['state']='disabled'
            self.buttonFinalPriority[0]['onvalue']=0
        if nb_p2 != 0:
            self.buttonFinalPriority[1]['state']='normal'
        else:
            #self.buttonFinalPriority[0]['state']='disabled'
            self.buttonFinalPriority[1]['onvalue']=0
        if nb_p3 != 0:
            self.buttonFinalPriority[2]['state']='normal'
        else:
            #self.buttonFinalPriority[0]['state']='disabled'
            self.buttonFinalPriority[2]['onvalue']=0
        """

        modename = ["DIA", "DLD", "IMA", "TMP", "SPI"]
        markername = ["s", "o", "^", "D", "v"]
        colorname = ["red", "blue", "green"]

        # Sort the Inst. Modes by alphabetical order
        uniqMode = sorted(set(list_targets["spica_mode"]))
        if "SPI" in uniqMode:
            iSPI = uniqMode.index("SPI")
            uniqMode.remove("SPI")
            uniqMode.insert(len(uniqMode), "SPI")

        if self.ra_sunrise < self.ra_sunset:
            midnight_offset = 12
        else:
            midnight_offset = 0

        for iMode in range(len(uniqMode)):

            ind_p1 = [
                (list_targets["priority_final"] == 1)
                & (list_targets["spica_mode"] == uniqMode[iMode])
            ]
            if np.any(np.array(ind_p1)):
                ra = list_targets["ra"][tuple(ind_p1)] / 15
                dec = list_targets["dec"][tuple(ind_p1)]
                vmag = list_targets["vmag"][tuple(ind_p1)]
                diam = list_targets["modeldiam"][tuple(ind_p1)]

                ind_ra = list(
                    filter(lambda x: ra[x] < self.ra_sunset / 15, range(len(ra)))
                )
                if ind_ra:
                    ra_morning = ra[ind_ra]
                    dec_morning = dec[ind_ra]
                    vmag_morning = vmag[ind_ra]
                    diam_morning = diam[ind_ra]

                    self.ax[0].scatter(
                        ra_morning + midnight_offset,
                        dec_morning,
                        s=1 / vmag_morning * 100,
                        color="red",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[1].scatter(
                        ra_morning + midnight_offset,
                        vmag_morning,
                        color="red",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[2].scatter(
                        ra_morning + midnight_offset,
                        diam_morning,
                        color="red",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )

                ind_ra = list(
                    filter(lambda x: ra[x] > self.ra_sunset / 15, range(len(ra)))
                )
                if ind_ra:
                    ra_afternoon = ra[ind_ra]
                    dec_afternoon = dec[ind_ra]
                    vmag_afternoon = vmag[ind_ra]
                    diam_afternoon = diam[ind_ra]

                    self.ax[0].scatter(
                        ra_afternoon - midnight_offset,
                        dec_afternoon,
                        s=1 / vmag_afternoon * 100,
                        color="red",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[1].scatter(
                        ra_afternoon - midnight_offset,
                        vmag_afternoon,
                        color="red",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[2].scatter(
                        ra_afternoon - midnight_offset,
                        diam_afternoon,
                        color="red",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )

            ind_p2 = [
                (list_targets["priority_final"] == 2)
                & (list_targets["spica_mode"] == uniqMode[iMode])
            ]
            if np.any(np.array(ind_p2)):
                ra = list_targets["ra"][tuple(ind_p2)] / 15
                dec = list_targets["dec"][tuple(ind_p2)]
                vmag = list_targets["vmag"][tuple(ind_p2)]
                diam = list_targets["modeldiam"][tuple(ind_p2)]

                ind_ra = list(
                    filter(lambda x: ra[x] < self.ra_sunset / 15, range(len(ra)))
                )
                if ind_ra:
                    ra_morning = ra[ind_ra]
                    dec_morning = dec[ind_ra]
                    vmag_morning = vmag[ind_ra]
                    diam_morning = diam[ind_ra]

                    self.ax[0].scatter(
                        ra_morning + midnight_offset,
                        dec_morning,
                        s=1 / vmag_morning * 100,
                        color="blue",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[1].scatter(
                        ra_morning + midnight_offset,
                        vmag_morning,
                        color="blue",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[2].scatter(
                        ra_morning + midnight_offset,
                        diam_morning,
                        color="blue",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )

                ind_ra = list(
                    filter(lambda x: ra[x] > self.ra_sunset / 15, range(len(ra)))
                )
                if ind_ra:
                    ra_afternoon = ra[ind_ra]
                    dec_afternoon = dec[ind_ra]
                    vmag_afternoon = vmag[ind_ra]
                    diam_afternoon = diam[ind_ra]

                    self.ax[0].scatter(
                        ra_afternoon - midnight_offset,
                        dec_afternoon,
                        s=1 / vmag_afternoon * 100,
                        color="blue",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[1].scatter(
                        ra_afternoon - midnight_offset,
                        vmag_afternoon,
                        color="blue",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[2].scatter(
                        ra_afternoon - midnight_offset,
                        diam_afternoon,
                        color="blue",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )

            ind_p3 = [
                (list_targets["priority_final"] == 3)
                & (list_targets["spica_mode"] == uniqMode[iMode])
            ]
            if np.any(np.array(ind_p3)):
                ra = list_targets["ra"][tuple(ind_p3)] / 15
                dec = list_targets["dec"][tuple(ind_p3)]
                vmag = list_targets["vmag"][tuple(ind_p3)]
                diam = list_targets["modeldiam"][tuple(ind_p3)]

                ind_ra = list(
                    filter(lambda x: ra[x] < self.ra_sunset / 15, range(len(ra)))
                )
                if ind_ra:
                    ra_morning = ra[ind_ra]
                    dec_morning = dec[ind_ra]
                    vmag_morning = vmag[ind_ra]
                    diam_morning = diam[ind_ra]

                    self.ax[0].scatter(
                        ra_morning + midnight_offset,
                        dec_morning,
                        s=1 / vmag_morning * 100,
                        color="green",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[1].scatter(
                        ra_morning + midnight_offset,
                        vmag_morning,
                        color="green",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[2].scatter(
                        ra_morning + midnight_offset,
                        diam_morning,
                        color="green",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )

                ind_ra = list(
                    filter(lambda x: ra[x] > self.ra_sunset / 15, range(len(ra)))
                )
                if ind_ra:
                    ra_afternoon = ra[ind_ra]
                    dec_afternoon = dec[ind_ra]
                    vmag_afternoon = vmag[ind_ra]
                    diam_afternoon = diam[ind_ra]

                    self.ax[0].scatter(
                        ra_afternoon - midnight_offset,
                        dec_afternoon,
                        s=1 / vmag_afternoon * 100,
                        color="green",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[1].scatter(
                        ra_afternoon - midnight_offset,
                        vmag_afternoon,
                        color="green",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[2].scatter(
                        ra_afternoon - midnight_offset,
                        diam_afternoon,
                        color="green",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )

            ind_p4 = [
                (list_targets["priority_final"] == 4)
                & (list_targets["spica_mode"] == uniqMode[iMode])
            ]
            if np.any(np.array(ind_p4)):
                ra = list_targets["ra"][tuple(ind_p4)] / 15
                dec = list_targets["dec"][tuple(ind_p4)]
                vmag = list_targets["vmag"][tuple(ind_p4)]
                diam = list_targets["modeldiam"][tuple(ind_p4)]

                ind_ra = list(
                    filter(lambda x: ra[x] < self.ra_sunset / 15, range(len(ra)))
                )
                if ind_ra:
                    ra_morning = ra[ind_ra]
                    dec_morning = dec[ind_ra]
                    vmag_morning = vmag[ind_ra]
                    diam_morning = diam[ind_ra]

                    self.ax[0].scatter(
                        ra_morning + midnight_offset,
                        dec_morning,
                        s=1 / vmag_morning * 100,
                        color="orange",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[1].scatter(
                        ra_morning + midnight_offset,
                        vmag_morning,
                        color="orange",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[2].scatter(
                        ra_morning + midnight_offset,
                        diam_morning,
                        color="orange",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )

                ind_ra = list(
                    filter(lambda x: ra[x] > self.ra_sunset / 15, range(len(ra)))
                )
                if ind_ra:
                    ra_afternoon = ra[ind_ra]
                    dec_afternoon = dec[ind_ra]
                    vmag_afternoon = vmag[ind_ra]
                    diam_afternoon = diam[ind_ra]

                    self.ax[0].scatter(
                        ra_afternoon - midnight_offset,
                        dec_afternoon,
                        s=1 / vmag_afternoon * 100,
                        color="orange",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[1].scatter(
                        ra_afternoon - midnight_offset,
                        vmag_afternoon,
                        color="orange",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )
                    self.ax[2].scatter(
                        ra_afternoon - midnight_offset,
                        diam_afternoon,
                        color="orange",
                        marker=markername[modename.index(uniqMode[iMode])],
                        alpha=0.5,
                        edgecolor="black",
                    )

        self.ax[0].set_xlabel("L.S.T. (hour)")
        self.ax[0].set_ylabel("DEC ($^\circ$)")
        xticks = [i for i in range(0, 26, 2)]
        xtick_labels = [
            (
                "{}".format(t + midnight_offset)
                if t < 12
                else "{}".format(t - midnight_offset)
            )
            for t in xticks
        ]
        self.ax[0].set_xticks(xticks)
        self.ax[0].set_xticklabels(xtick_labels)
        self.ax[0].set_xlim(
            [
                self.ra_sunset / 15 - midnight_offset - 1,
                self.ra_sunrise / 15 + midnight_offset + 1 + max(self.dif) / 15,
            ]
        )
        self.ax[0].set_ylim([-90, 90])
        self.ax[0].set_facecolor("gainsboro")
        self.ax[0].grid(True, linestyle="--", linewidth=0.5)
        self.ax[0].plot(
            [
                self.ra_sunset / 15 - midnight_offset,
                self.ra_sunset / 15 - midnight_offset,
            ],
            [-90, 90],
            "k:",
        )
        self.ax[0].plot(
            [
                self.ra_sunrise / 15 + midnight_offset,
                self.ra_sunrise / 15 + midnight_offset,
            ],
            [-90, 90],
            "k:",
        )

        self.ax[1].set_ylabel("Vmag")
        self.ax[1].set_xticks(xticks)
        self.ax[1].set_xticklabels(xtick_labels)
        self.ax[1].set_xlim(
            [
                self.ra_sunset / 15 - midnight_offset - 1,
                self.ra_sunrise / 15 + midnight_offset + 1 + max(self.dif) / 15,
            ]
        )
        self.ax[1].set_ylim([-3, 14.9])
        self.ax[1].set_facecolor("gainsboro")
        self.ax[1].grid(True, linestyle="--", linewidth=0.5)
        self.ax[1].plot(
            [
                self.ra_sunset / 15 - midnight_offset,
                self.ra_sunset / 15 - midnight_offset,
            ],
            [-3, 14.9],
            "k:",
        )
        self.ax[1].plot(
            [
                self.ra_sunrise / 15 + midnight_offset,
                self.ra_sunrise / 15 + midnight_offset,
            ],
            [-3, 14.9],
            "k:",
        )

        self.ax[2].set_xlabel("L.S.T. (hour)")
        self.ax[2].set_ylabel("Diam (mas)")
        self.ax[2].set_xticks(xticks)
        self.ax[2].set_xticklabels(xtick_labels)
        self.ax[2].set_xlim(
            [
                self.ra_sunset / 15 - midnight_offset - 1,
                self.ra_sunrise / 15 + midnight_offset + 1 + max(self.dif) / 15,
            ]
        )
        self.ax[2].set_ylim([-3, 14.9])
        self.ax[2].set_facecolor("gainsboro")
        self.ax[2].grid(True, linestyle="--", linewidth=0.5)
        self.ax[2].plot(
            [
                self.ra_sunset / 15 - midnight_offset,
                self.ra_sunset / 15 - midnight_offset,
            ],
            [-3, 14.9],
            "k:",
        )
        self.ax[2].plot(
            [
                self.ra_sunrise / 15 + midnight_offset,
                self.ra_sunrise / 15 + midnight_offset,
            ],
            [-3, 14.9],
            "k:",
        )

        # text1 and text2 contain the actual text, created by TextArea
        # text1 and text2 are then packed vertically into a box using VPacker
        text1 = TextArea("Priority_final = 1", textprops=dict(color="red"))
        text2 = TextArea("Priority_final = 2", textprops=dict(color="blue"))
        text3 = TextArea("Priority_final = 3", textprops=dict(color="green"))
        text4 = TextArea("Priority_final = 4", textprops=dict(color="orange"))
        box = VPacker(children=[text1, text2, text3, text4], align="left", pad=2, sep=4)

        # anchored_box creates the text box outside of the plot
        anchored_box = AnchoredOffsetbox(
            loc=3,
            child=box,
            pad=0.2,
            frameon=True,
            bbox_to_anchor=(1.01, 0.7),
            bbox_transform=self.ax[0].transAxes,
            borderpad=0.0,
        )
        anchored_box.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        anchored_box.patch.set_facecolor("wheat")
        self.ax[0].add_artist(anchored_box)

        nb_p1 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 1)
                & (list_targets["spica_mode"] == "DIA")
            ]
        )
        nb_p2 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 2)
                & (list_targets["spica_mode"] == "DIA")
            ]
        )
        nb_p3 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 3)
                & (list_targets["spica_mode"] == "DIA")
            ]
        )
        nb_p4 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 4)
                & (list_targets["spica_mode"] == "DIA")
            ]
        )
        text1 = TextArea(
            "\u25A0 DIA"
            + " ("
            + str(nb_p1)
            + "/"
            + str(nb_p2)
            + "/"
            + str(nb_p3)
            + "/"
            + str(nb_p4)
            + ")"
        )

        nb_p1 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 1)
                & (list_targets["spica_mode"] == "DLD")
            ]
        )
        nb_p2 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 2)
                & (list_targets["spica_mode"] == "DLD")
            ]
        )
        nb_p3 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 3)
                & (list_targets["spica_mode"] == "DLD")
            ]
        )
        nb_p4 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 4)
                & (list_targets["spica_mode"] == "DLD")
            ]
        )
        text2 = TextArea(
            "\u25CF DLD"
            + " ("
            + str(nb_p1)
            + "/"
            + str(nb_p2)
            + "/"
            + str(nb_p3)
            + "/"
            + str(nb_p4)
            + ")"
        )

        nb_p1 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 1)
                & (list_targets["spica_mode"] == "IMA")
            ]
        )
        nb_p2 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 2)
                & (list_targets["spica_mode"] == "IMA")
            ]
        )
        nb_p3 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 3)
                & (list_targets["spica_mode"] == "IMA")
            ]
        )
        nb_p4 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 4)
                & (list_targets["spica_mode"] == "IMA")
            ]
        )
        text3 = TextArea(
            "\u25B2 IMA"
            + " ("
            + str(nb_p1)
            + "/"
            + str(nb_p2)
            + "/"
            + str(nb_p3)
            + "/"
            + str(nb_p4)
            + ")"
        )

        nb_p1 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 1)
                & (list_targets["spica_mode"] == "TMP")
            ]
        )
        nb_p2 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 2)
                & (list_targets["spica_mode"] == "TMP")
            ]
        )
        nb_p3 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 3)
                & (list_targets["spica_mode"] == "TMP")
            ]
        )
        nb_p4 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 4)
                & (list_targets["spica_mode"] == "TMP")
            ]
        )
        text4 = TextArea(
            "\u25C6 TMP"
            + " ("
            + str(nb_p1)
            + "/"
            + str(nb_p2)
            + "/"
            + str(nb_p3)
            + "/"
            + str(nb_p4)
            + ")"
        )

        nb_p1 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 1)
                & (list_targets["spica_mode"] == "SPI")
            ]
        )
        nb_p2 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 2)
                & (list_targets["spica_mode"] == "SPI")
            ]
        )
        nb_p3 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 3)
                & (list_targets["spica_mode"] == "SPI")
            ]
        )
        nb_p4 = np.count_nonzero(
            [
                (list_targets["priority_final"] == 4)
                & (list_targets["spica_mode"] == "SPI")
            ]
        )
        text5 = TextArea(
            "\u25BC SPI"
            + " ("
            + str(nb_p1)
            + "/"
            + str(nb_p2)
            + "/"
            + str(nb_p3)
            + "/"
            + str(nb_p4)
            + ")"
        )

        box = VPacker(
            children=[text1, text2, text3, text4, text5], align="left", pad=2, sep=4
        )

        # anchored_box creates the text box outside of the plot
        anchored_box = AnchoredOffsetbox(
            loc=3,
            child=box,
            pad=0.2,
            frameon=True,
            bbox_to_anchor=(1.01, 0.0),
            bbox_transform=self.ax[0].transAxes,
            borderpad=0.0,
        )
        anchored_box.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        anchored_box.patch.set_facecolor("lightskyblue")
        self.ax[0].add_artist(anchored_box)

        # text1 and text2 contain the actual text, created by TextArea
        # text1 and text2 are then packed vertically into a box using VPacker
        if self.indexList_CalPrim is not None:  # self.calprim_catg:
            nb_calprim = len(self.indexList_CalPrim)  # len(self.calprim_catg)
        else:
            nb_calprim = 0
        if self.indexList_CalSec is not None:  # self.calsec_catg:
            nb_calsec = len(self.indexList_CalSec)  # len(self.calsec_catg)
        else:
            nb_calsec = 0
        text1 = TextArea(
            "CalPRIM: " + str(nb_calprim) + " \u2605",
            textprops=dict(color="dodgerblue"),
        )
        text2 = TextArea(
            "CalSEC: " + str(nb_calsec) + " \u2605", textprops=dict(color="orange")
        )
        box = VPacker(children=[text1, text2], align="left", pad=2, sep=4)

        # anchored_box creates the text box outside of the plot
        anchored_box = AnchoredOffsetbox(
            loc=3,
            child=box,
            pad=0.0,
            frameon=True,
            bbox_to_anchor=(1.01, -0.3),
            bbox_transform=self.ax[0].transAxes,
            borderpad=0.0,
        )
        anchored_box.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        anchored_box.patch.set_facecolor("whitesmoke")
        self.ax[0].add_artist(anchored_box)

        fig.subplots_adjust(hspace=0, right=0.845, top=0.9)

        # Plot primary calibrators if any
        if self.indexList_CalPrim is not None:
            list_calprim = self.calprim_catg[self.indexList_CalPrim]  # indexList1]
            ra = list_calprim["ra"] / 15
            dec = list_calprim["dec"]
            vmag = list_calprim["vmag"]
            diam = list_calprim["ldd"]

            ind_ra = list(filter(lambda x: ra[x] < self.ra_sunset / 15, range(len(ra))))
            if ind_ra:
                ra_morning = ra[ind_ra]
                dec_morning = dec[ind_ra]
                vmag_morning = vmag[ind_ra]
                diam_morning = diam[ind_ra]
                self.ax[0].scatter(
                    ra_morning + midnight_offset,
                    dec_morning,
                    color="dodgerblue",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )
                self.ax[1].scatter(
                    ra_morning + midnight_offset,
                    vmag_morning,
                    color="dodgerblue",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )
                self.ax[2].scatter(
                    ra_morning + midnight_offset,
                    diam_morning,
                    color="dodgerblue",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )

            ind_ra = list(filter(lambda x: ra[x] > self.ra_sunset / 15, range(len(ra))))
            if ind_ra:
                ra_afternoon = ra[ind_ra]
                dec_afternoon = dec[ind_ra]
                vmag_afternoon = vmag[ind_ra]
                diam_afternoon = diam[ind_ra]
                self.ax[0].scatter(
                    ra_afternoon - midnight_offset,
                    dec_afternoon,
                    color="dodgerblue",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )
                self.ax[1].scatter(
                    ra_afternoon - midnight_offset,
                    vmag_afternoon,
                    color="dodgerblue",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )
                self.ax[2].scatter(
                    ra_afternoon - midnight_offset,
                    diam_afternoon,
                    color="dodgerblue",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )

        # Plot secondary calibrators if any
        if self.indexList_CalSec is not None:
            list_calsec = self.calsec_catg[self.indexList_CalSec]
            ra = list_calsec["ra"] / 15
            dec = list_calsec["dec"]
            vmag = list_calsec["Vmag"]
            diam = list_calsec["UDDR"]

            ind_ra = list(filter(lambda x: ra[x] < self.ra_sunset / 15, range(len(ra))))
            if ind_ra:
                ra_morning = ra[ind_ra]
                dec_morning = dec[ind_ra]
                vmag_morning = vmag[ind_ra]
                diam_morning = diam[ind_ra]
                self.ax[0].scatter(
                    ra_morning + midnight_offset,
                    dec_morning,
                    color="orange",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )
                self.ax[1].scatter(
                    ra_morning + midnight_offset,
                    vmag_morning,
                    color="orange",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )
                self.ax[2].scatter(
                    ra_morning + midnight_offset,
                    diam_morning,
                    color="orange",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )

            ind_ra = list(filter(lambda x: ra[x] > self.ra_sunset / 15, range(len(ra))))
            if ind_ra:
                ra_afternoon = ra[ind_ra]
                dec_afternoon = dec[ind_ra]
                vmag_afternoon = vmag[ind_ra]
                diam_afternoon = diam[ind_ra]
                self.ax[0].scatter(
                    ra_afternoon - midnight_offset,
                    dec_afternoon,
                    color="orange",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )
                self.ax[1].scatter(
                    ra_afternoon - midnight_offset,
                    vmag_afternoon,
                    color="orange",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )
                self.ax[2].scatter(
                    ra_afternoon - midnight_offset,
                    diam_afternoon,
                    color="orange",
                    marker="*",
                    alpha=0.5,
                    edgecolor="black",
                    s=150,
                )

        framePlot = LabelFrame(
            self.root, text="Selected stars", font=self.myFontLabelFrame, padx=5
        )
        framePlot.grid(row=8, column=0, columnspan=2, pady=5, ipady=1)
        canvas = FigureCanvasTkAgg(fig, framePlot)
        canvas.draw()
        canvas.get_tk_widget().grid()  # row=0,column=0)

        # cursor = DataCursorPlot.FollowDotCursor(self.ax[0], ra_morning+midnight_offset, dec_morning)

        # #canvas.draw()
        # # Add navigation toolbar
        # toolbarFrame = Frame(framePlot)
        # toolbarFrame.grid()#row=1,column=0)
        # toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

        """
        for i in self.tree_frameInfoTargets.get_children():
            self.tree_frameInfoTargets.delete(i)
        """

        """
        if self.popupInfoTargets:
            self.tree_frameInfoTargets.destroy()
            self.popupInfoTargets = False
            #self.open_popupInfoTargets()
         """
        """
        if self.popupBestDec:
            self.tree_frameBestDec.destroy()
            self.popupBestDec = False
            #self.open_popupBestDec()
        """

    def open_popupProgName(self):
        """


        Returns
        -------
        None.

        """
        top = Toplevel(self.root)
        # self.top.geometry('+850+50')
        top.title("Available workpackages")
        Label(top, text="Workpackages: ", font=self.myFont).grid()

        # Search for distinct prognames (split WPN-WPM to distinct ones)
        wps = [
            p.split("-")
            for p in np.unique(self.spica_catg["progname"][self.indexInstMode]).data
        ]
        self.uniqProgName = np.sort(np.unique(np.concatenate(wps)))
        self.buttonProgName = []
        self.SelectedProgName = []
        for iProgName in self.uniqProgName:
            ProgNameId = StringVar()
            ProgNameId.set("0")
            self.buttonProgName.append(
                Checkbutton(
                    top,
                    text=iProgName,
                    onvalue=iProgName,
                    offvalue=0,
                    indicatoron=1,
                    variable=ProgNameId,
                    selectcolor="red",
                    activebackground="green",
                    command=self.plotSelectedProgName,
                )
            )
            self.SelectedProgName.append(ProgNameId)

        buttonSelectAllProgName = Button(
            top,
            text="SELECT ALL",
            font=self.myFont,
            activebackground="green",
            command=self.select_allProgName,
        )
        buttonSelectAllProgName.grid(column=1, row=0)

        for iProgName in range(len(self.uniqProgName)):
            self.buttonProgName[iProgName].grid(column=0, row=iProgName + 1, sticky=W)

        buttonSubmitProgName = Button(
            top,
            text="SUBMIT",
            font=self.myFont,
            activebackground="green",
            command=self.plotSelectedProgName,
        )
        buttonSubmitProgName.grid(column=0, row=iProgName + 2)

        buttonResetProgName = Button(
            top,
            text="RESET",
            font=self.myFont,
            activebackground="orange",
            command=self.deselect_allProgName,
        )
        buttonResetProgName.grid(column=1, row=iProgName + 2)

        buttonCloseProgName = Button(
            top,
            text="CLOSE",
            font=self.myFont,
            activebackground="red",
            command=top.destroy,
        )
        buttonCloseProgName.grid(column=2, row=iProgName + 2)

    def select_allProgName(self):
        """

        Returns
        -------
        None.

        """
        for j in self.buttonProgName:
            j.select()

    def deselect_allProgName(self):
        self.onQuery()

    def plotSelectedInstMode(self):
        """

        Returns
        -------
        None.

        """
        """
        self.indexInstMode = []
        count= 0

        for j in self.SelectedInstMode:
            if (j.get() != str(0)):
                if count == 0:
                    self.indexInstMode = list(filter(lambda x: self.spica_catg['spica_mode'][x] == j.get(), range(len(self.spica_catg))))
                else:
                    self.indexInstMode = np.concatenate([self.indexInstMode, list(filter(lambda x: self.spica_catg['spica_mode'][x] == j.get(), range(len(self.spica_catg))))])
                count+=1
        self.indexList_Targets = reduce(np.intersect1d, [self.indexProgName, self.indexInstMode, self.indexFinalPriority,
                                                          self.indexDecMin, self.indexDecMax,
                                                          self.indexVmagMin, self.indexVmagMax])
        self.plot_radec()#self.indexList_Targets, self.indexList_CalPrim, self.indexList_CalSec)
        """
        self.indexInstMode = []
        count = 0
        for j in self.SelectedInstMode:
            print("[INFO] Inst. Modes", j.get())
            if count == 0:
                self.indexInstMode = list(
                    filter(
                        lambda x: self.spica_catg["spica_mode"][x] == j.get(),
                        range(len(self.spica_catg)),
                    )
                )
            else:
                self.indexInstMode = np.concatenate(
                    [
                        self.indexInstMode,
                        list(
                            filter(
                                lambda x: self.spica_catg["spica_mode"][x] == j.get(),
                                range(len(self.spica_catg)),
                            )
                        ),
                    ]
                )
            count += 1
        self.indexInstMode = self.indexInstMode.astype("int")
        """
        self.indexList_Targets = reduce(np.intersect1d, [self.indexProgName, self.indexInstMode, self.indexFinalPriority,
                                                          self.indexDecMin, self.indexDecMax,
                                                          self.indexVmagMin, self.indexVmagMax])
        """
        self.getSelectedTargets()
        self.getAddTarget()
        self.plot_radec()

        # Refresh treeview items
        if self.popupInfoTargets:
            self.clear_all(self.my_treeInfoTargets)
            self.insert_popupInfoTargets()
            self.tree_frameInfoTargets.title(
                f"List of selected targets ({len(self.indexList_Targets)} objects)"
            )

        if self.popupAddTarget:
            self.clear_all(self.my_treeAddTarget)
            self.insert_popupAddTarget()
            self.tree_frameAddTarget.title(
                f"List of possible fainter stars to include ({len(self.indexList_AddTarget)} objects)"
            )

    def plotSelectedFinalPriority(self):
        """

        Returns
        -------
        None.

        """
        self.indexFinalPriority = []
        count = 0
        for j in self.SelectedFinalPriority:
            print("[INFO] Priority", j.get())
            if count == 0:
                self.indexFinalPriority = list(
                    filter(
                        lambda x: self.spica_catg["priority_final"][x] == j.get(),
                        range(len(self.spica_catg)),
                    )
                )
            else:
                self.indexFinalPriority = np.concatenate(
                    [
                        self.indexFinalPriority,
                        list(
                            filter(
                                lambda x: self.spica_catg["priority_final"][x]
                                == j.get(),
                                range(len(self.spica_catg)),
                            )
                        ),
                    ]
                )
            count += 1
        self.indexFinalPriority = self.indexFinalPriority.astype("int")
        """
        self.indexList_Targets = reduce(np.intersect1d, [self.indexProgName, self.indexInstMode, self.indexFinalPriority,
                                                          self.indexDecMin, self.indexDecMax,
                                                          self.indexVmagMin, self.indexVmagMax])
        """
        self.getSelectedTargets()
        self.getAddTarget()
        self.plot_radec()

        # Refresh treeview items
        if self.popupInfoTargets:
            self.clear_all(self.my_treeInfoTargets)
            self.insert_popupInfoTargets()
            self.tree_frameInfoTargets.title(
                f"List of selected targets ({len(self.indexList_Targets)} objects)"
            )

        if self.popupAddTarget:
            self.clear_all(self.my_treeAddTarget)
            self.insert_popupAddTarget()
            self.tree_frameAddTarget.title(
                f"List of possible fainter stars to include ({len(self.indexList_AddTarget)} objects)"
            )

    def open_popupInstMode(self):
        """

        Returns
        -------
        None.

        """
        top = Toplevel(self.root)
        # top.geometry('250x750+750+50')
        top.title("Inst. Modes")
        Label(top, text="Modes:", font=self.myFont).grid()

        self.uniqInstMode = sorted(
            set(self.spica_catg["spica_mode"][self.indexProgName])
        )
        buttonInstMode = []
        self.SelectedInstMode = []
        if "SPI" in self.uniqInstMode:
            iSPI = self.uniqInstMode.index("SPI")
            self.uniqInstMode.remove("SPI")
            self.uniqInstMode.insert(len(self.uniqInstMode), "SPI")

        for iInstMode in self.uniqInstMode:
            InstModeId = StringVar()
            InstModeId.set("0")
            buttonInstMode.append(
                Checkbutton(
                    top,
                    text=iInstMode,
                    onvalue=iInstMode,
                    offvalue=0,
                    indicatoron=1,
                    variable=InstModeId,
                    selectcolor="red",
                    activebackground="green",
                    command=self.plotSelectedInstMode,
                )
            )
            self.SelectedInstMode.append(InstModeId)

        for iInstMode in range(len(self.uniqInstMode)):
            buttonInstMode[iInstMode].grid(column=0, row=iInstMode + 1, sticky=W)

        buttonSubmitInstMode = Button(
            top,
            text="SUBMIT",
            font=self.myFont,
            activebackground="green",
            command=self.plotSelectedInstMode,
        )
        buttonSubmitInstMode.grid(column=0, row=iInstMode + 2)

        buttonCloseInstMode = Button(
            top,
            text="CLOSE",
            font=self.myFont,
            activebackground="red",
            command=top.destroy,
        )
        buttonCloseInstMode.grid(column=1, row=iInstMode + 2)

    def entryDateCallback(self, strDate):
        self.date = strDate.get()
        self.labelValDate.config(text=self.date)
        self.onQuery()

    def entryDecMinCallback(self, strDecMin):
        self.indexDecMin = []
        self.decmin = float(strDecMin.get())
        if self.decmin <= -30.0:
            self.decmin = -30.0
        self.indexDecMin = list(
            filter(
                lambda x: self.spica_catg["dec"][x] > self.decmin,
                range(len(self.spica_catg)),
            )
        )
        """
        self.indexList_Targets = reduce(np.intersect1d, [self.indexProgName, self.indexInstMode, self.indexFinalPriority,
                                                          self.indexDecMin, self.indexDecMax,
                                                          self.indexVmagMin, self.indexVmagMax])
        """
        self.getSelectedTargets()
        self.getAddTarget()
        self.plot_radec()

        # Refresh treeview items
        if self.popupInfoTargets:
            self.clear_all(self.my_treeInfoTargets)
            self.insert_popupInfoTargets()
            self.tree_frameInfoTargets.title(
                f"List of selected targets ({len(self.indexList_Targets)} objects)"
            )

        if self.popupAddTarget:
            self.clear_all(self.my_treeAddTarget)
            self.insert_popupAddTarget()
            self.tree_frameAddTarget.title(
                f"List of possible fainter stars to include ({len(self.indexList_AddTarget)} objects)"
            )

    def entryDecMaxCallback(self, strDecMax):
        self.indexDecMax = []
        self.decmax = float(strDecMax.get())
        if self.decmax >= 90.0:
            self.decmax = 90.0
        self.indexDecMax = list(
            filter(
                lambda x: self.spica_catg["dec"][x] < self.decmax,
                range(len(self.spica_catg)),
            )
        )
        """
        self.indexList_Targets = reduce(np.intersect1d, [self.indexProgName, self.indexInstMode, self.indexFinalPriority,
                                                          self.indexDecMin, self.indexDecMax,
                                                          self.indexVmagMin, self.indexVmagMax])
        """
        self.getSelectedTargets()
        self.getAddTarget()
        self.plot_radec()

        # Refresh treeview items
        if self.popupInfoTargets:
            self.clear_all(self.my_treeInfoTargets)
            self.insert_popupInfoTargets()
            self.tree_frameInfoTargets.title(
                f"List of selected targets ({len(self.indexList_Targets)} objects)"
            )

        if self.popupAddTarget:
            self.clear_all(self.my_treeAddTarget)
            self.insert_popupAddTarget()
            self.tree_frameAddTarget.title(
                f"List of possible fainter stars to include ({len(self.indexList_AddTarget)} objects)"
            )

    def entryVmagMinCallback(self, strVmagMin):
        self.indexVmagMin = []
        self.vmagmin = float(strVmagMin.get())
        self.indexVmagMin = list(
            filter(
                lambda x: self.spica_catg["vmag"][x] > self.vmagmin,
                range(len(self.spica_catg)),
            )
        )
        """
        self.indexList_Targets = reduce(np.intersect1d, [self.indexProgName, self.indexInstMode, self.indexFinalPriority,
                                                          self.indexDecMin, self.indexDecMax,
                                                          self.indexVmagMin, self.indexVmagMax])
        """
        self.getSelectedTargets()
        self.getAddTarget()
        self.plot_radec()

        # Refresh treeview items
        if self.popupInfoTargets:
            self.clear_all(self.my_treeInfoTargets)
            self.insert_popupInfoTargets()
            self.tree_frameInfoTargets.title(
                f"List of selected targets ({len(self.indexList_Targets)} objects)"
            )

        if self.popupAddTarget:
            self.clear_all(self.my_treeAddTarget)
            self.insert_popupAddTarget()
            self.tree_frameAddTarget.title(
                f"List of possible fainter stars to include ({len(self.indexList_AddTarget)} objects)"
            )

    def entryVmagMaxCallback(self, strVmagMax):
        self.indexVmagMax = []
        self.indexVmagToAddTarget = []
        self.vmagmax = float(strVmagMax.get())
        self.indexVmagMax = list(
            filter(
                lambda x: self.spica_catg["vmag"][x] < self.vmagmax,
                range(len(self.spica_catg)),
            )
        )
        self.indexVmagToAddTarget = list(
            filter(
                lambda x: self.spica_catg["vmag"][x] >= self.vmagmax
                and self.spica_catg["vmag"][x] < self.vmagmax + 1,
                range(len(self.spica_catg)),
            )
        )

        """
        self.indexList_Targets = reduce(np.intersect1d, [self.indexProgName, self.indexInstMode, self.indexFinalPriority,
                                                          self.indexDecMin, self.indexDecMax,
                                                          self.indexVmagMin, self.indexVmagMax])
        """
        self.getSelectedTargets()
        self.getAddTarget()
        self.plot_radec()

        # Refresh treeview items
        if self.popupInfoTargets:
            self.clear_all(self.my_treeInfoTargets)
            self.insert_popupInfoTargets()
            self.tree_frameInfoTargets.title(
                f"List of selected targets ({len(self.indexList_Targets)} objects)"
            )

        if self.popupAddTarget:
            self.clear_all(self.my_treeAddTarget)
            self.insert_popupAddTarget()
            self.tree_frameAddTarget.title(
                f"List of possible fainter stars to include ({len(self.indexList_AddTarget)} objects)"
            )

    def onLog(self):
        if self.LogOpened:
            self.onCloseLog()
        self.LogOpened = True

        # Load spicadb catalog
        print("\n")
        spicadb_load = self.dbquery_tap()

        nb_completion_ok = [
            x
            for x in range(len(spicadb_load))
            if spicadb_load["completion_rate"][x] >= 1
        ]
        print(f"[INFO] Number of completed SPICA-DB targets: {len(nb_completion_ok)}")

        uniqProgName = list(np.unique(spicadb_load["progname"]))
        # perc_Progselected = np.zeros(len(uniqProgName))
        perc_Progcompleted = np.zeros(len(uniqProgName))

        for i, prog_name in enumerate(uniqProgName):
            nb_total = sum(
                1
                for x in range(len(spicadb_load))
                if spicadb_load["progname"][x] == prog_name
            )
            nb_completed = sum(
                1
                for x in range(len(spicadb_load))
                if spicadb_load["progname"][x] == prog_name
                and spicadb_load["completion_rate"][x] >= 1
            )
            # nb_selected = sum(1 for x in range(len(list_targets)) if list_targets['progname'][x] == prog_name)

            print(
                f"[INFO] Program: {prog_name} -> Total: {nb_total}, Completed: {nb_completed}"
            )

            if nb_total > 0:
                # perc_Progselected[i] = nb_selected / nb_total * 100
                perc_Progcompleted[i] = nb_completed / nb_total * 100

        uniqInstMode = list(np.unique(spicadb_load["spica_mode"]))
        # perc_Modeselected = np.zeros(len(self.uniqInstMode_jmmc))
        perc_Modecompleted = np.zeros(len(uniqInstMode))

        for i, mode in enumerate(uniqInstMode):
            nb_total = sum(
                1
                for x in range(len(spicadb_load))
                if spicadb_load["spica_mode"][x] == mode
            )
            nb_completed = sum(
                1
                for x in range(len(spicadb_load))
                if spicadb_load["spica_mode"][x] == mode
                and spicadb_load["completion_rate"][x] >= 1
            )
            # nb_selected = sum(1 for x in range(len(list_targets)) if list_targets['spica_mode'][x] == mode)

            print(
                f"[INFO] Inst. mode: {mode} -> Total: {nb_total}, Completed: {nb_completed}"
            )

            if nb_total > 0:
                # perc_Modeselected[i] = nb_selected / nb_total * 100
                perc_Modecompleted[i] = nb_completed / nb_total * 100

        print(f"[INFO] Perc. Completed Programs (%): {np.round(perc_Progcompleted,2)}")
        print(
            f"[INFO] Perc. Completed Inst. Modes (%): {np.round(perc_Modecompleted,2)}"
        )

        self.topLog = Toplevel(self.root)
        self.topLog.title("LOG")
        fig2, axs = plt.subplots(2, 1, figsize=(15, 10))  # , sharex=True)

        for iplot, (data, labels) in enumerate(
            [(perc_Progcompleted, uniqProgName), (perc_Modecompleted, uniqInstMode)]
        ):
            axs[iplot].bar(labels, data, label="Completed")
            axs[iplot].bar(labels, 100 - data, bottom=data, label="Not completed")

            # Annotations for 'Completed' bars
            for xpos, ypos, yval in zip(labels, data, data):
                axs[iplot].text(xpos, ypos / 2, f"{yval:.1f}", ha="center", va="center")

            # Annotations for 'Not completed' bars
            for xpos, ypos, yval in zip(labels, data, 100 - data):
                axs[iplot].text(
                    xpos,
                    ypos + (100 - ypos) / 2,
                    f"{yval:.1f}",
                    ha="center",
                    va="center",
                )

            if iplot == 0:
                axs[iplot].set_title("Percentages of completed stars per program")
            else:
                axs[iplot].set_title(
                    "Percentages of completed stars per instrumental mode"
                )

            axs[iplot].set_ylabel("%")
            axs[iplot].set_facecolor("gainsboro")
            axs[iplot].grid(True, linestyle="--", linewidth=0.5)
            axs[iplot].set_ylim([0, 110])
            axs[iplot].legend(bbox_to_anchor=(1.01, 0.5), loc="center left")

        canvas2 = FigureCanvasTkAgg(fig2, self.topLog)
        canvas2.get_tk_widget().grid(row=1, column=0)
        canvas2.draw()

        buttonCloseProgName = Button(
            self.topLog,
            text="CLOSE",
            font=self.myFont,
            fg="white",
            bg="red",
            command=self.onCloseLog,
            cursor="X_cursor",
        )
        buttonCloseProgName.grid(column=0, row=2)

    def onQuit(self):
        if messagebox.askyesno("Exit", "Do you want to quit the SPICA-NSS Tool?"):
            """
            if self.client:
                self.client.disconnect()
            """
            self.root.destroy()
            plt.close("all")
            print("\nBye!")

    def onCloseCalsec(self):
        self.topCalsec.destroy()
        self.CalsecOpened = False

    def onCloseLog(self):
        self.topLog.destroy()
        self.LogOpened = False

    def onAspro(self):
        """

        Returns
        -------
        None.

        """
        """
        #print(datetime.now().strftime('%H:%M:%S >>> '), 'Sending the selected targets to ASPRO...')
        self.topAspro= Toplevel(self.root)
        #top.geometry('250x750+750+50')
        Label(self.topAspro, text= f'You are about to import to ASPRO2:\n', font=self.myFont).grid(column=0,row=0,padx=5,pady=4)#,columnspan=0)
        Label(self.topAspro, text= f'- {len(self.spica_catg[self.indexList_Targets])} science targets,\n', font=self.myFont).grid(column=0,row=1,sticky='W',padx=10)

        if self.calprim_catg:
            Label(self.topAspro, text= f'- {len(self.calprim_catg[self.indexList_CalPrim])} primary calibrators,\n', font=self.myFont).grid(column=0,row=2,sticky='W',padx=10)
        else:
            Label(self.topAspro, text= f'- 0 primary calibrator,\n', font=self.myFont).grid(column=0,row=2,sticky='W',padx=10)
        if self.calsec_catg:
            Label(self.topAspro, text= f'- {len(self.calsec_catg[self.indexList_CalSec])} secondary calibrators.', font=self.myFont).grid(column=0,row=3,sticky='W',padx=10)
        else:
            Label(self.topAspro, text= f'- 0 secondary calibrator.', font=self.myFont).grid(column=0,row=3,sticky='W',padx=10)
        Label(self.topAspro, text= f'\n').grid(column=0,row=4)

        frame = Frame(self.topAspro)
        frame.grid(pady=5)
        Label(frame, text= f'Are you sure?', font=self.myFont).grid(column=0,row=5,columnspan=2)
        Button(frame, text='Yes', font=self.myFont, fg='white', bg='green', cursor='hand1', command=self.import2aspro).grid(column=0,row=6)
        Button(frame, text='No', font=self.myFont, fg='white', bg='red', cursor='hand1', command=self.topAspro.destroy).grid(column=1,row=6)
        """

        if self.spica_catg is None:
            messagebox.showinfo(title="Info", message="Please query SPICA-DB first.")
            return

        self.import2aspro()

    def import2aspro(self):

        targets2aspro = self.spica_catg[self.indexList_Targets]
        # Replace by mainID by HD name if needed
        targets2aspro = self.replacebyHDname(targets2aspro)

        print("[INFO] * Sending to ASPRO2 *")
        print("[INFO] ---------------------")
        print(f"[INFO] {len(targets2aspro)} targets")
        if len(targets2aspro) > 1:
            targets_sort = "transit"
            if targets_sort == "ra":
                targets2aspro = targets2aspro[np.argsort(targets2aspro["ra"])]
            elif targets_sort == "transit":
                coord = SkyCoord(
                    ra=list(targets2aspro["ra"]) * u.deg,
                    dec=list(targets2aspro["dec"]) * u.deg,
                )
                targets = FixedTarget(
                    coord=coord, name=list(targets2aspro["target_main_id"])
                )
                transit_times = self.chara.target_meridian_transit_time(
                    self.date, targets, which="next"
                )
                targets2aspro = targets2aspro[np.argsort(transit_times)]

        calsprim2aspro = None
        calssec2aspro = None
        if self.indexList_CalPrim is not None:
            calsprim2aspro = self.calprim_catg[self.indexList_CalPrim]
            print(f"[INFO] {len(calsprim2aspro)} primary calibrators")
        else:
            print("[INFO] 0 primary calibrator")
        if self.indexList_CalSec is not None:
            calssec2aspro = self.calsec_catg[self.indexList_CalSec]
            print(f"[INFO] {len(calssec2aspro)} secondary calibrators")
        else:
            print("[INFO] 0 secondary calibrator")
        print("")

        self.addDiamModel(targets2aspro, calsprim2aspro, calssec2aspro)
        self.addNssType(targets2aspro, calsprim2aspro, calssec2aspro)
        self.samp_votable(targets2aspro, calsprim2aspro, calssec2aspro)

        # self.topAspro.destroy()

    def addDiamModel(self, targets, calibrators1=None, calibrators2=None):
        """


        Parameters
        ----------
        targets : TYPE
            DESCRIPTION.
        calibrators1 : TYPE, optional
            DESCRIPTION. The default is None.
        calibrators2 : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        model = targets["model"]
        modeltype = []
        modeldiam = []
        for j in range(len(targets)):
            a = model[j].replace("}{", "},{")
            modeltype.append(json.loads(a)[0]["type"])
            if modeltype[j] == "disk":
                modeldiam.append(json.loads(a)[0]["diameter"])
            elif modeltype[j] == "elong_disk":
                modeldiam.append(json.loads(a)[0]["minor_axis_diameter"])
        modeldiam = MaskedColumn(modeldiam)
        targets["diam"] = modeldiam

        if self.indexList_CalPrim is not None:
            diam_calprim = calibrators1["ldd"]
            diam_calprim = MaskedColumn(diam_calprim)
            calibrators1["diam"] = diam_calprim

        if self.indexList_CalSec is not None:
            diam_calsec = calibrators2["UDDR"]
            diam_calsec = MaskedColumn(diam_calsec)
            calibrators2["diam"] = diam_calsec

        return targets, calibrators1, calibrators2

    def addNssType(self, targets, calibrators1=None, calibrators2=None):
        """


        Parameters
        ----------
        targets : TYPE
            DESCRIPTION.
        calibrators1 : TYPE, optional
            DESCRIPTION. The default is None.
        calibrators2 : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        targets["nss_type"] = ["Science"] * len(targets)
        if calibrators1:
            calibrators1["nss_type"] = ["CalPrim"] * len(calibrators1)
        if calibrators2:
            calibrators2["nss_type"] = ["CalSec"] * len(calibrators2)

        return targets, calibrators1, calibrators2

    def onQuery(self):
        """

        Returns
        -------
        None.

        """
        print(f"\n[INFO] Querying the SPICA-DB catalog [{self.catalogSpica}]...")

        # Close already opened widgets
        if self.popupInfoTargets:
            self.tree_frameInfoTargets.destroy()
            self.popupInfoTargets = False

        if self.popupBestDec:
            self.tree_frameBestDec.destroy()
            self.popupBestDec = False

        # Delete previous selection if any
        if self.index_AddTarget is not None:
            self.index_AddTarget = None

        # Create variables to contain the list of indices in each catalog
        self.indexList_Targets = None
        self.indexList_CalPrim = None
        self.indexList_CalSec = None

        # Get the CHARA coordinates
        self.chara = Observer.at_site("CHARA")

        # Get the observable RA domain between sunset and sunrise
        # print(self.iter)
        self.ra_sunset, self.ra_sunrise = self.observable_domain()
        if self.iter == 0:

            # Select only stars with a completion_rate < 1
            tableTargets = self.dbquery_tap()[self.dbquery_tap()["completion_rate"] < 1]

            self.iter += 1
            self.spica_catg_jmmc = tableTargets
            self.uniqProgName_jmmc = sorted(set(tableTargets["progname"]))
            self.uniqInstMode_jmmc = sorted(set(tableTargets["spica_mode"]))
            if "SPI" in self.uniqInstMode_jmmc:
                iSPI = self.uniqInstMode_jmmc.index("SPI")
                self.uniqInstMode_jmmc.remove("SPI")
                self.uniqInstMode_jmmc.insert(len(self.uniqInstMode_jmmc), "SPI")
        else:
            tableTargets = self.spica_catg_jmmc

            # Reset entries for Declination and Magnitude of the Science Objects
            self.entryDecMin = Entry(
                self.FrameObjects,
                textvariable=self.strDecMin,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryDecMin.delete(0, END)
            self.entryDecMin.insert(END, str(self.decmin_default))
            self.entryDecMax = Entry(
                self.FrameObjects,
                textvariable=self.strDecMax,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryDecMax.delete(0, END)
            self.entryDecMax.insert(END, str(self.decmax_default))
            self.entryVmagMin = Entry(
                self.FrameObjects,
                textvariable=self.strVmagMin,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryVmagMin.delete(0, END)
            self.entryVmagMin.insert(END, str(self.vmagmin_default))
            self.entryVmagMax = Entry(
                self.FrameObjects,
                textvariable=self.strVmagMax,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryVmagMax.delete(0, END)
            self.entryVmagMax.insert(END, str(self.vmagmax_default))

            # Reset entries for the Primary Calibrators
            self.entryRaRangePrim = Entry(
                self.FrameCalPrims,
                textvariable=self.strRaRangePrim,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryRaRangePrim.delete(0, END)
            self.entryRaRangePrim.insert(END, str(self.rarangeprim_default))
            self.entryDecRangePrim = Entry(
                self.FrameCalPrims,
                textvariable=self.strDecRangePrim,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryDecRangePrim.delete(0, END)
            self.entryDecRangePrim.insert(END, str(self.decrangeprim_default))
            self.entryVmagRangePrim = Entry(
                self.FrameCalPrims,
                textvariable=self.strVmagRangePrim,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryVmagRangePrim.delete(0, END)
            self.entryVmagRangePrim.insert(END, str(self.vmagrangeprim_default))

            # Reset entries for the Secondary Calibrators
            self.entryRaRangeSec = Entry(
                self.FrameCalSecs,
                textvariable=self.strRaRangeSec,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryRaRangeSec.delete(0, END)
            self.entryRaRangeSec.insert(END, str(self.rarangesec_default))
            self.entryDecRangeSec = Entry(
                self.FrameCalSecs,
                textvariable=self.strDecRangeSec,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryDecRangeSec.delete(0, END)
            self.entryDecRangeSec.insert(END, str(self.decrangesec_default))
            self.entryVmagRangeSec = Entry(
                self.FrameCalSecs,
                textvariable=self.strVmagRangeSec,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryVmagRangeSec.delete(0, END)
            self.entryVmagRangeSec.insert(END, str(self.vmagrangesec_default))
            self.entryLDDChiSec = Entry(
                self.FrameCalSecs,
                textvariable=self.strLDDChiSec,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryLDDChiSec.delete(0, END)
            self.entryLDDChiSec.insert(END, str(self.lddchisec_default))
            self.entryRelErrorSec = Entry(
                self.FrameCalSecs,
                textvariable=self.strRelErrorSec,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryRelErrorSec.delete(0, END)
            self.entryRelErrorSec.insert(END, str(self.relerrorsec_default))
            self.entryMinVisSec = Entry(
                self.FrameCalSecs,
                textvariable=self.strMinVisSec,
                justify="right",
                font=self.myFont,
                cursor="pencil",
                width=5,
            )
            self.entryMinVisSec.delete(0, END)
            self.entryMinVisSec.insert(END, str(self.minvissec_default))
            self.entryMaxBaseline.set(self.maxbaseline_default)

        # nb_completion_ok = list(filter(lambda x: tableTargets['completion_rate'][x] >= 1, range(len(tableTargets))))
        # print(f'[INFO] Number of completed SPICA-DB targets: {len(nb_completion_ok)}')

        # Filtering on RA
        chara_geolat = 34.220920783 * np.pi / 180.0
        dec = tableTargets["dec"] * np.pi / 180.0
        cor = (np.sin(70.0 * np.pi / 180.0) - np.sin(dec) * np.sin(chara_geolat)) / (
            np.cos(dec) * np.cos(chara_geolat)
        )
        self.dif = np.where(cor < 1, np.arccos(cor), 0) * 180.0 / np.pi

        if self.ra_sunset < self.ra_sunrise:
            tableTargets = tableTargets[
                (tableTargets["ra"] + self.dif > self.ra_sunset)
                & (tableTargets["ra"] - self.dif < self.ra_sunrise)
            ]
        else:
            tableTargets = tableTargets[
                (tableTargets["ra"] + self.dif > self.ra_sunset)
                | (tableTargets["ra"] - self.dif < self.ra_sunrise)
            ]
        """
        if self.ra_sunset < self.ra_sunrise:
            tableTargets = tableTargets[(tableTargets['ra'] > self.ra_sunset) & (tableTargets['ra'] < self.ra_sunrise)]
        else:
            tableTargets = tableTargets[(tableTargets['ra'] > self.ra_sunset) | (tableTargets['ra'] < self.ra_sunrise)]
        """

        self.uniqProgName = sorted(set(tableTargets["progname"]))
        self.uniqInstMode = sorted(set(tableTargets["spica_mode"]))
        if "SPI" in self.uniqInstMode:
            iSPI = self.uniqInstMode.index("SPI")
            self.uniqInstMode.remove("SPI")
            self.uniqInstMode.insert(len(self.uniqInstMode), "SPI")

        self.indexProgName = list(range(len(tableTargets)))
        self.indexInstMode = list(range(len(tableTargets)))
        self.indexFinalPriority = list(range(len(tableTargets)))
        self.indexDecMin = list(range(len(tableTargets)))
        self.indexDecMax = list(range(len(tableTargets)))
        self.indexVmagMin = list(range(len(tableTargets)))
        self.indexVmagMax = list(range(len(tableTargets)))
        self.indexVmagToAddTarget = list(range(len(tableTargets)))

        # Calculate the final priority of each star in the database
        self.flag_completion = []
        self.priority_final = []
        for i in range(len(tableTargets)):
            self.flag_completion.append(
                self.update_flag_completion(
                    tableTargets["completion_rate"][i], tableTargets["spica_mode"][i]
                )
            )

            # self.priority_final.append(self.update_priority_final(self.flag_completion[i],
            #                                                     tableTargets['priority_pi'][i]))

            tableTargets["priority_final"][i] = self.update_priority_final2(
                self.flag_completion[i],
                tableTargets["priority_pi"][i],
                tableTargets["progname2"][i],
            )
        for iProgName in list(range(len(self.ProgName))):
            # print(self.ProgName[iProgName])
            if np.any([(tableTargets["progname"] == self.ProgName[iProgName])]):
                self.buttonProgName[iProgName]["state"] = "normal"
                self.buttonProgName[iProgName].select()

        for iInstMode in list(range(len(self.InstMode))):
            # print(self.InstMode[iInstMode])
            if np.any([(tableTargets["spica_mode"] == self.InstMode[iInstMode])]):
                self.buttonInstMode[iInstMode]["state"] = "normal"
                self.buttonInstMode[iInstMode].select()

        for iFinalPriority in list(range(4)):
            if np.any([(tableTargets["priority_final"] == iFinalPriority + 1)]):
                self.buttonFinalPriority[iFinalPriority]["state"] = "normal"
                self.buttonFinalPriority[iFinalPriority].select()

        # Replace by mainID by HD name if needed
        tableTargets = self.replacebyHDname(tableTargets)
        self.spica_catg = tableTargets
        self.indexList_Targets = list(range(len(tableTargets)))
        self.targetListInit = self.spica_catg[self.indexList_Targets]
        self.plot_radec()

    def onFilters(self):
        self.onQuery()
        self.plot_radec()

    def onReset(self):
        """

        Returns
        -------
        None.

        """
        if self.spica_catg is None:
            messagebox.showinfo(title="Info", message="Please query SPICA-DB first.")
            return

        self.iter = 0
        self.onQuery()

        # Reset entries for Declination and Magnitude of the Science Objects
        self.entryDecMin = Entry(
            self.FrameObjects,
            textvariable=self.strDecMin,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryDecMin.delete(0, END)
        self.entryDecMin.insert(END, str(self.decmin_default))
        self.entryDecMax = Entry(
            self.FrameObjects,
            textvariable=self.strDecMax,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryDecMax.delete(0, END)
        self.entryDecMax.insert(END, str(self.decmax_default))
        self.entryVmagMin = Entry(
            self.FrameObjects,
            textvariable=self.strVmagMin,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryVmagMin.delete(0, END)
        self.entryVmagMin.insert(END, str(self.vmagmin_default))
        self.entryVmagMax = Entry(
            self.FrameObjects,
            textvariable=self.strVmagMax,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryVmagMax.delete(0, END)
        self.entryVmagMax.insert(END, str(self.vmagmax_default))

        # Reset entries for the Primary Calibrators
        self.entryRaRangePrim = Entry(
            self.FrameCalPrims,
            textvariable=self.strRaRangePrim,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryRaRangePrim.delete(0, END)
        self.entryRaRangePrim.insert(END, str(self.rarangeprim_default))
        self.entryDecRangePrim = Entry(
            self.FrameCalPrims,
            textvariable=self.strDecRangePrim,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryDecRangePrim.delete(0, END)
        self.entryDecRangePrim.insert(END, str(self.decrangeprim_default))
        self.entryVmagRangePrim = Entry(
            self.FrameCalPrims,
            textvariable=self.strVmagRangePrim,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryVmagRangePrim.delete(0, END)
        self.entryVmagRangePrim.insert(END, str(self.vmagrangeprim_default))

        # Reset entries for the Secondary Calibrators
        self.entryRaRangeSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strRaRangeSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryRaRangeSec.delete(0, END)
        self.entryRaRangeSec.insert(END, str(self.rarangesec_default))
        self.entryDecRangeSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strDecRangeSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryDecRangeSec.delete(0, END)
        self.entryDecRangeSec.insert(END, str(self.decrangesec_default))
        self.entryVmagRangeSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strVmagRangeSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryVmagRangeSec.delete(0, END)
        self.entryVmagRangeSec.insert(END, str(self.vmagrangesec_default))
        self.entryLDDChiSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strLDDChiSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryLDDChiSec.delete(0, END)
        self.entryLDDChiSec.insert(END, str(self.lddchisec_default))
        self.entryRelErrorSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strRelErrorSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryRelErrorSec.delete(0, END)
        self.entryRelErrorSec.insert(END, str(self.relerrorsec_default))
        self.entryMinVisSec = Entry(
            self.FrameCalSecs,
            textvariable=self.strMinVisSec,
            justify="right",
            font=self.myFont,
            cursor="pencil",
            width=5,
        )
        self.entryMinVisSec.delete(0, END)
        self.entryMinVisSec.insert(END, str(self.minvissec_default))
        self.entryMaxBaseline.set(self.maxbaseline_default)

        print("[INFO] nbtotal", self.nbtotal_p1, self.nbtotal_p2, self.nbtotal_p3)
        # To reinitialize the fake DRS
        # for j in self.SelectedFinalPriority:
        #    j.set(0)

        if self.popupInfoTargets:
            self.tree_frameInfoTargets.destroy()
            self.popupInfoTargets = False

        if self.popupBestDec:
            self.tree_frameBestDec.destroy()
            self.popupBestDec = False

    def onDRS(self):
        """

        Returns
        -------
        None.

        """
        print(f"\n[INFO] Updating the catalog after DRS...")

        targets2drs = self.spica_catg[self.indexList_Targets]

        # Fake DRS/QCS
        for k in range(len(targets2drs)):
            targets2drs["completion_rate"][k] = self.update_completion_rate(
                targets2drs["completion_rate"][k], targets2drs["spica_mode"][k]
            )
            targets2drs["ob_refs"][k] = str(int(targets2drs["ob_refs"][k]) + 1)
            targets2drs["qcs_flag"][k] = str(0)

        # Update local spica_db after fake DRS
        for k in range(len(self.indexList_Targets)):
            pos = [(self.spica_catg_jmmc["spicadb_id"] == targets2drs["spicadb_id"][k])]
            self.spica_catg_jmmc["completion_rate"][tuple(pos)] = targets2drs[
                "completion_rate"
            ][k]
            self.spica_catg_jmmc["ob_refs"][tuple(pos)] = targets2drs["ob_refs"][k]
            self.spica_catg_jmmc["qcs_flag"][tuple(pos)] = targets2drs["qcs_flag"][k]

        print("[INFO] Done.")

    def dbquery_tap(self):
        """
        This returns the SPICA-DB queried via TAP protocol. A first selection of
        stars can be made, e.g. in declination

        Parameters
        ----------
        None.

        Returns
        -------
        data : Table
            The SPICA_DB as data variable.
        """
        # tapServerUrl = 'http://tap-preprod.jmmc.fr/vollt/tap/' # url subject to change
        print("tapServerUrl:", tapServerUrl)

        # Query all spica-db catalog content
        # catalog = 'spica_2023_06_07' # let use a snapshot on a specific date for test
        adqlQuery = (
            "SELECT * FROM " + self.catalogSpica
        )  # + " WHERE dec > 20 AND dec < 40" # " + " ORDER BY progname #" #+ " WHERE ra > " + str(ra_sun_set) + " AND ra < " + str(ra_sun_rise)

        service = vo.dal.TAPService(tapServerUrl)
        results = service.search(adqlQuery)

        # store new information into data variable
        data = results.to_table()

        # Sort data in ascending RA
        # data = data[np.argsort(data['ra'])]

        # Select only stars with a completion_rate < 1
        # data = data[data['completion_rate'] < 1]

        # Remove stars without magnitude information
        # data = data[data['vmag'].mask == False]

        # Correct the masked values which can befound in S03 for the completion_rate and the priority_pi
        for i in range(len(data["completion_rate"])):
            if ma.is_masked(data["completion_rate"][i]):
                data["completion_rate"][i] = float(
                    ma.getdata(data["completion_rate"][i])
                )

        for i in range(len(data["priority_pi"])):
            if ma.is_masked(data["priority_pi"][i]):
                data["priority_pi"][i] = ma.getdata(data["priority_pi"][i])

        return data

    def calquery_tap(self):
        """
        This returns the SPICA-DB calibrators queried via TAP protocol. A first selection of
        stars can be made, e.g. in declination

        Parameters
        ----------
        None.

        Returns
        -------
        data : Table
            The SPICA_DB as data variable.
        """
        # tapServerUrl = 'http://tap.jmmc.fr/vollt/tap/' # url subject to change

        # query all calprim content
        calprim = "spica_calprim"
        calAdqlQuery = "SELECT * FROM " + calprim

        service = vo.dal.TAPService(tapServerUrl)
        results = service.search(calAdqlQuery)

        # store new information into data variable
        calibrators = results.to_table()

        return calibrators

    def observable_domain(self):
        """
        This returns the observable RA domain for a given observatory.

        Parameters
        ----------
        dateobs : Class
            The date of the observations.
        observer : Class
            The observer's location

        Returns
        -------
        alpha_sun_set : Float
            The RA (hour angle or alpha) of the sun set.
        alpha_sun_rise : Float
            The RA (hour angle or alpha) of the sun rise.
        """
        # Get the times of the sunset and the sunrise minus/plus one hour (15 degrees)
        sun_set = self.chara.sun_set_time(
            Time(self.date, scale="utc"), which="nearest", horizon=0 * u.deg
        )
        sun_rise = self.chara.sun_rise_time(
            Time(self.date, scale="utc"), which="next", horizon=0 * u.deg
        )
        print("[INFO] Sunset: ISO: {0.iso}, JD: {0.jd}".format(sun_set))
        print("[INFO] Sunrise: ISO: {0.iso}, JD: {0.jd}".format(sun_rise))
        print("")

        # Convert sunset and surise times into local sidereal times
        lmst_sun_set = self.chara.local_sidereal_time(sun_set)
        lmst_sun_rise = self.chara.local_sidereal_time(sun_rise)
        print("[INFO] local_sidereal_time_Sunset:", lmst_sun_set)
        print("[INFO] local_sidereal_time_Sunrise:", lmst_sun_rise)
        print("")

        # Convert the lmst into hour angles (degrees)
        alpha_sun_set = (
            lmst_sun_set.hms.h + lmst_sun_set.hms.m / 60.0 + lmst_sun_set.hms.s / 3600.0
        ) * 15.0
        alpha_sun_rise = (
            lmst_sun_rise.hms.h
            + lmst_sun_rise.hms.m / 60.0
            + lmst_sun_rise.hms.s / 3600.0
        ) * 15.0
        print(
            f"[INFO] On day {self.date}: Observable RA domain [{round(alpha_sun_set,2)}, {round(alpha_sun_rise,2)}] (deg)"
        )
        print("")

        return alpha_sun_set, alpha_sun_rise

    def fixColumnTypes(self, ts):
        """Fixes some datatypes (object -> str) of given table list inplace so we can vstack cal and sci tables."""
        # fix ucds on client side (this does not modify parent's table column's meta)
        ucd4col = {
            # consider target_main_id column as target id (instead of spicadb_id)
            "target_main_id": "meta.id;meta.main",
            "spicadb_id": "meta.code.class",
            ## declare main group (instead of progname)
            #'priority_pi':'meta.code.class;meta.id',
            # do not declare progname as main group  column
            "progname": "meta.code.class",
            "progname2": "meta.code.class",
            #  and fallback some columns to be shown as an extra-info in the target editor
            "piname": "meta.code.class",
            "piname2": "meta.code.class",
            "CALIBRATOR_NAME": "meta.id",
            "SCIENCE_TARGET_NAME": "meta.id",
            # Convert ucd Vizier to ucd Aspro
            "vmag": "phot.mag;em.opt.V",
            "rmag": "phot.mag;em.opt.R",
            "hmag": "phot.mag;em.optIR.H",
            # Add diameter from model for extra infos
            "diam": "meta.code.class",
            # Add NSS_Type for extra infos (NSS_Type = Science, CalPrim, ou CalSec)
            "nss_type": "meta.code.class",
            # "comments":'meta.code.class',
            # NAME # "": "meta.id;meta.main"
            # RA # "": "pos.eq.ra;meta.main"
            # DEC # "": "pos.eq.dec;meta.main"
            # RV # "": "spect.dopplerVeloc.opt"
            # PMRA # "": "pos.pm;pos.eq.ra"
            # PMDEC # "": "pos.pm;pos.eq.dec"
            # PLX # "": "pos.parallax.trig"
            # e_PLX # "": "stat.error;pos.parallax.trig"
            # HD # "": "meta.id"
            # HIP # "": "meta.id"
            # 2MASS # "": "meta.id"
            # OTYPES # "": "src.class"
            # SP_TYPES
            "sptype": "src.spType",
            # FLUX_B"
            "bmag": "phot.mag;em.opt.B",
            # FLUX_V # "": "phot.mag;em.opt.V"
            # FLUX_G
            "gmag": "phot.mag;em.opt.G",
            # FLUX_R # "": "phot.mag;em.opt.R"
            # FLUX_I
            "imag": "phot.mag;em.opt.I",
            # FLUX_J # "": "phot.mag;em.IR.J"
            # FLUX_H
            "hmag": "phot.mag;em.IR.H",
            # FLUX_K # "": "phot.mag;em.IR.K"
            # FLUX_L
            "lmag": "phot.mag;em.IR.3-4um",
            # FLUX_M
            "mmag": "phot.mag;em.IR.4-8um",
            # FLUX_N
            "nmag": "phot.mag;em.IR.8-15um",
            # GROUP # "": "meta.code.class;meta.id"
            # EXTRA_INFORMATION # "": "meta.code.class"
        }
        # print("FIXME remove comments ucd change so we can display them in Aspro2 as soon as calibrator+desc is supported")

        # priority-pi can then be used with free values
        strCols = ["priority_pi"]

        if not isinstance(ts, list):
            ts = [ts]
        for t in ts:
            for c in t.colnames:
                if t[c].dtype == "O" or c in strCols:
                    t[c] = t[c].astype(str)
                if c in ucd4col.keys():
                    t[c].meta["ucd"] = ucd4col[c]

    def normalizeColumnNames(
        self,
        input,
        colNames={
            "ra": "ra",
            "dec": "dec",
            "name": "target_main_id",
            "ldd": "ld_jsdc2",
            "vmag": "vmag",
            "hmag": "hmag",
            "jmag": "jmag",
            "kmag": "kmag",
            "rmag": "rmag",
            "bmag": "bmag",
            "imag": "imag",
            "hmag": "hmag",
            "lmag": "lmag",
            "mmag": "mmag",
            "nmag": "nmag",
            "diam": "diam",
            "sptype": "spt",
            "nss_type": "nss_type",
        },
    ):

        # We could imagine a mapping on top of ucd but because of UCD1+(spica_survey) vs UCD1(VizieR)
        # -> Use a manual mapping
        """Return a new table. Each columns given by colnames keys are copied and renamed to associated values for output table.
        Default colNames values are selected so we can use a calibrator table from spica-calprim or Vizier JSDC2.
        """

        output = Table()
        for c in input.colnames:
            if c in colNames:
                output.add_column(input[c], name=colNames[c])
            elif c.lower() in colNames:
                output.add_column(input[c], name=colNames[c.lower()])
            elif c in colNames.values():
                output.add_column(input[colNames[c]])
        self.fixColumnTypes(output)

        return output

    def samp_votable(self, targets, calibrators1=None, calibrators2=None):
        """
        This sends the filtered list of stars from SPICA-DB to ASPRO2 as a VOtable
        through samp protocol.

        Parameters
        ----------
        targets : Table
            The filtered SPICA-DB as a data variable.

        Returns
        -------
        None.
        """

        self.client = SAMPIntegratedClient(name="SPICA-NSS")

        try:
            self.client.connect()
            tmpname = tempfile.NamedTemporaryFile(delete=False).name

            calPrims = None
            calPrims4sci = None
            calSeconds = None
            calSeconds4sci = None

            # declare some constants
            COLNAME_GRP = "grp"

            # Add primary calibrators if any
            if calibrators1:
                # TODO filter calibrator list given current science data selection
                calPrims = self.normalizeColumnNames(calibrators1)
                calPrims.add_column(
                    "calprim", name=COLNAME_GRP
                )  # set main grouptargets2aspro

                """
                # create a fake star that will gather all cal prim to avoid orphan calibrators
                FAKECALPRIMNAME="lesCALPRIM"
                r = calPrims[0].as_void()
                r["target_main_id"]=FAKECALPRIMNAME
                r["ra"]=r["dec"]=0
                r["ld_jsdc2"]=float("NaN") # TODO fix code so we do not compute a disk for this fake star
                r["vmag"]=r["hmag"]=r["jmag"]=r["kmag"]=r["rmag"]=0
                calPrims.insert_row(0,r)
                """

                # second table to declare some calibrators through Aspro's votable
                calPrims4sci = Table()
                calPrims4sci.add_column(
                    calPrims["target_main_id"], name="CALIBRATOR_NAME"
                )
                # calPrims4sci.add_column(FAKECALPRIMNAME, name="SCIENCE_TARGET_NAME")

                # calibrators1["vmag"].meta["ucd"]='phot.mag;em.opt.V'

            # Add secondary calibrators if any
            if calibrators2:
                calSeconds = self.normalizeColumnNames(calibrators2)
                calSeconds.add_column("calsecond", name=COLNAME_GRP)  # set main group

                # We should not have orphan calibrators since queryJsdc2 perform the association
                """
                """
                """
                # Prepare associations
                FAKECALSECONDNAME="lesCALSECOND"
                # associate all cal to a fake star
                calSeconds4sci=Table()
                calSeconds4sci.add_column(calSeconds["target_main_id"], name="CALIBRATOR_NAME")
                calSeconds4sci.add_column(FAKECALSECONDNAME, name="SCIENCE_TARGET_NAME")


                # and create this fake star that will gather all cal secondary
                r = calSeconds[0].as_void()
                r["target_main_id"]=FAKECALSECONDNAME
                r["ra"]=r["dec"]=0
                r["ld_jsdc2"]=float("NaN") # TODO fix code so we do not compute a disk for this fake star
                r["vmag"]=r["hmag"]=r["jmag"]=r["kmag"]=r["rmag"]=0
                calSeconds.insert_row(0,r)
                """

            if COLNAME_GRP in targets.colnames:
                targets.remove_column(COLNAME_GRP)
            targets.add_column(None, name=COLNAME_GRP)
            targets[COLNAME_GRP].meta["ucd"] = "meta.code.class;meta.id"

            for r in targets:
                r[COLNAME_GRP] = f"priority_pi={r['priority_pi']}"

            # Modify the column 'comments' to include the PI names, the prognames, and the instrumental modes
            for i in range(len(targets["comments"])):
                if not targets["piname2"][i]:
                    targets["comments"][i] = (
                        targets["piname"][i]
                        + ", "
                        + targets["progname"][i]
                        + ", "
                        + targets["spica_mode"][i]
                        + ", "
                        + targets["comments"][i]
                    )
                else:
                    targets["comments"][i] = (
                        targets["piname"][i]
                        + ", "
                        + targets["progname"][i]
                        + ", "
                        + targets["spica_mode"][i]
                        + ", "
                        + targets["comments"][i]
                        + "\n"
                        + targets["piname2"][i]
                        + ", "
                        + targets["progname2"][i]
                        + ", "
                        + targets["comments2"][i]
                    )

            # Add toy model to the tables
            stars2aspro = targets
            self.add_targetmodel(stars2aspro)
            self.fixColumnTypes([stars2aspro])

            if calPrims:
                self.add_toymodel(calPrims)
                self.fixColumnTypes([calPrims])
                stars2aspro = vstack([stars2aspro, calPrims])
            if calSeconds:
                self.add_toymodel(calSeconds)
                self.fixColumnTypes([calSeconds])
                stars2aspro = vstack([stars2aspro, calSeconds])

            # Sort by RA
            if len(stars2aspro) > 1:
                coord = SkyCoord(
                    ra=list(stars2aspro["ra"]) * u.deg,
                    dec=list(stars2aspro["dec"]) * u.deg,
                )
                targets = FixedTarget(
                    coord=coord, name=list(stars2aspro["target_main_id"])
                )
                transit_times = self.chara.target_meridian_transit_time(
                    self.date, targets, which="next"
                )
                stars2aspro = stars2aspro[np.argsort(transit_times)]

            print(stars2aspro)

            # create votable from data table and store it in the temporary file
            vot = astropy.io.votable.from_table(stars2aspro)

            # prepare a table for group's color mapping
            colors = [
                ["calprim", "calsecond", "priority_pi=0", "priority_pi=1"],
                ["#3232ff", "#b2b2ff", "#ff4c4c", "#ffb2b2"],
            ]
            colorsTable = Table(data=colors, names=("GROUP_NAME", "GROUP_COLOR"))

            # Modify
            # add main observation metadata

            year = int(self.date[0:4])
            month = int(self.date[5:7])
            if month < 8:
                semester = "A"
            else:
                semester = "B"
            charaPeriod = f"CHARA {year}{semester}"

            for r in vot.resources:
                r.params.extend(
                    [
                        Param(
                            vot,
                            name="OPERATION",
                            datatype="char",
                            arraysize="*",
                            value="NEW",
                        ),
                        Param(
                            vot,
                            name="INTERFEROMETER",
                            datatype="char",
                            arraysize="*",
                            value="CHARA",
                        ),
                        Param(
                            vot,
                            name="PERIOD",
                            datatype="char",
                            arraysize="*",
                            value=f"{charaPeriod}",
                        ),  # CHARA Future"),
                        Param(
                            vot,
                            name="INSTRUMENT",
                            datatype="char",
                            arraysize="*",
                            value="SPICA_6T",
                        ),
                        Param(
                            vot,
                            name="CONFIGURATIONS",
                            datatype="char",
                            arraysize="*",
                            value="S1 S2 W1 W2 E1 E2",
                        ),
                        Param(vot, name="NIGHT", datatype="boolean", value="true"),
                        Param(
                            vot,
                            name="DATE",
                            datatype="char",
                            arraysize="*",
                            value=self.date,
                        ),
                        Param(vot, name="MIN_ELEVATION", datatype="int", value="30"),
                    ]
                )
            r.tables.append(astropy.io.votable.tree.Table.from_table(vot, colorsTable))

            # store data in temp votable file
            vot.to_xml(tmpname)

            print(
                datetime.now().strftime("%H:%M:%S >>> "), f"'{tmpname}' votable ready"
            )
            # prepare message and broadcast it
            self.message = {
                "samp.mtype": "table.load.votable",
                "samp.params": {"url": Path(tmpname).as_uri()},
            }

            # receiver_ids = client.notify('c1', message)
            # receiver_ids = client.notify_all(message)
            self.connectedAsproClients = []
            for sampClient in self.client.get_registered_clients():
                if self.client.get_metadata(sampClient)["samp.name"] == "Aspro2":
                    print(
                        "[INFO] Opened Aspro2 clients",
                        sampClient,
                        self.client.get_metadata(sampClient)["samp.name"],
                    )
                    self.connectedAsproClients.append(sampClient)

            if self.topSAMP:
                self.topSAMP.destroy()
                self.topSAMP = False

            self.topSAMP = Toplevel(self.root)
            Label(
                self.topSAMP,
                text=f"You are about to import to ASPRO2:\n",
                font=self.myFont,
            ).grid(
                column=0, row=0, padx=5, pady=4
            )  # ,columnspan=0)
            Label(
                self.topSAMP,
                text=f"- {len(self.spica_catg[self.indexList_Targets])} science targets,\n",
                font=self.myFont,
            ).grid(column=0, row=1, sticky="W", padx=10)

            if self.indexList_CalPrim:  # self.calprim_catg:
                Label(
                    self.topSAMP,
                    text=f"- {len(self.calprim_catg[self.indexList_CalPrim])} primary calibrators,\n",
                    font=self.myFont,
                ).grid(column=0, row=2, sticky="W", padx=10)
            else:
                Label(
                    self.topSAMP, text=f"- 0 primary calibrator,\n", font=self.myFont
                ).grid(column=0, row=2, sticky="W", padx=10)
            if self.calsec_catg:
                Label(
                    self.topSAMP,
                    text=f"- {len(self.calsec_catg[self.indexList_CalSec])} secondary calibrators.",
                    font=self.myFont,
                ).grid(column=0, row=3, sticky="W", padx=10)
            else:
                Label(
                    self.topSAMP, text=f"- 0 secondary calibrator.", font=self.myFont
                ).grid(column=0, row=3, sticky="W", padx=10)
            Label(self.topSAMP, text=f"\n").grid(column=0, row=4)

            print(len(self.connectedAsproClients))

            if len(self.connectedAsproClients) == 0:
                Label(
                    self.topSAMP,
                    text="WARNING",
                    font=self.myFont,
                    fg="white",
                    bg="red",
                    relief="solid",
                    borderwidth=2,
                ).grid(column=0, row=5, padx=5, pady=4, columnspan=5)
                Label(
                    self.topSAMP,
                    text="Unable to find a running SAMP Hub. Please launch Aspro2.",
                    font=self.myFont,
                ).grid(column=0, row=6, padx=5, pady=4, columnspan=5)
                Button(
                    self.topSAMP,
                    text="Ok",
                    font=self.myFont,
                    fg="white",
                    bg="red",
                    cursor="hand1",
                    command=self.topSAMP.destroy,
                ).grid(column=0, row=7, columnspan=5)
                self.client.disconnect()
            elif len(self.connectedAsproClients) == 1:
                frame = Frame(self.topSAMP)
                frame.grid(pady=5)
                Label(frame, text="Are you sure?", font=self.myFont).grid(
                    column=0, row=5, columnspan=2
                )
                Button(
                    frame,
                    text="Yes",
                    font=self.myFont,
                    fg="white",
                    bg="green",
                    cursor="hand1",
                    command=self.importClients,
                ).grid(
                    column=0, row=6
                )  # self.import2aspro).grid(column=0,row=9)
                # Button(frame, text='No', font=self.myFont, fg='white', bg='red', cursor='hand1', command=self.topSAMP.destroy).grid(column=1,row=6)
                Button(
                    frame,
                    text="No",
                    font=self.myFont,
                    fg="white",
                    bg="red",
                    cursor="hand1",
                    command=self.cancelClickClients,
                ).grid(column=1, row=6)
            else:
                Label(
                    self.topSAMP,
                    text="WARNING",
                    font=self.myFont,
                    fg="white",
                    bg="red",
                    relief="solid",
                    borderwidth=2,
                ).grid(column=0, row=5, padx=5, pady=4, columnspan=5)
                Label(
                    self.topSAMP,
                    text="Several ASPRO2 clients are opened.",
                    font=self.myFont,
                ).grid(column=0, row=6, padx=5, pady=4, columnspan=5)
                Label(
                    self.topSAMP, text="Choose the ID receiver(s):\n", font=self.myFont
                ).grid(column=0, row=7, padx=5, pady=4, columnspan=5)

                frame = Frame(self.topSAMP)
                frame.grid(pady=5)
                self.SelectedClient = []
                for x, n in enumerate(self.connectedAsproClients):
                    self.SelectedClient.append(StringVar(frame, value=0))
                    l = Checkbutton(frame, text=n, variable=self.SelectedClient[x])
                    l.grid(row=1, column=x + 1, sticky=W)
                Button(
                    frame,
                    text="Send2Aspro",
                    font=self.myFont,
                    fg="white",
                    bg="green",
                    cursor="hand1",
                    command=self.importClients,
                ).grid(column=x + 2, row=3, padx=5)
                Button(
                    frame,
                    text="Cancel",
                    font=self.myFont,
                    fg="white",
                    bg="red",
                    cursor="hand1",
                    command=self.cancelClickClients,
                ).grid(column=x + 3, row=3, padx=5)

        except:
            print("[INFO] Unable to find a running SAMP Hub. Please launch Aspro2.")
            messagebox.showinfo(title="Info", message="Please launch Aspro2 first.")

        # print(datetime.now().strftime('%H:%M:%S >>> '), f"{tmpname} sent by samp")

    def cancelClickClients(self):
        self.topSAMP.destroy()
        self.client.disconnect()

    def importClients(self):
        """


        Returns
        -------
        None.

        """
        if len(self.connectedAsproClients) == 1:
            self.client.notify(self.connectedAsproClients[0], self.message)
        elif len(self.connectedAsproClients) > 1:
            for i in range(len(self.SelectedClient)):
                if self.SelectedClient[i].get() == "1":
                    self.client.notify(self.connectedAsproClients[i], self.message)
        self.topSAMP.destroy()
        self.client.disconnect()

    def update_completion_rate(self, completion_rate: float, spica_mode: str) -> float:
        """
        Updates the completion rate based on the SPICA mode.

        Parameters
        ----------
        completion_rate : float
            The initial completion rate value.
        spica_mode : str
            The SPICA mode which determines how the completion rate is adjusted.

        Returns
        -------
        float
            The updated completion rate after adjustments based on the SPICA mode.
        """
        try:
            # Handle masked completion rate values
            if ma.is_masked(completion_rate):
                completion_rate = 0.0

            # Adjust completion rate based on SPICA mode
            if completion_rate < 1:
                if spica_mode == "DIA":
                    completion_rate += 0.5
                elif spica_mode == "DLD":
                    completion_rate += 0.25
                elif spica_mode in ("IMA", "TMP"):
                    completion_rate += 0.1
                elif spica_mode == "SPI":
                    completion_rate += 0.2

            return completion_rate

        except Exception as e:
            print(f"[ERROR] An unexpected error occurred: {e}")
            return completion_rate

    def update_flag_completion(self, completion_rate, spica_mode):
        """

        Parameters
        ----------
        completion_rate : TYPE
            DESCRIPTION.
        spica_mode : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        flag = []
        # if ma.is_masked(completion_rate):
        #    completion_rate = float(ma.getdata(completion_rate))

        if (spica_mode == "DIA") & (completion_rate > 0):
            flag = 1

        if (spica_mode == "DLD") & (completion_rate >= 0.5):
            flag = 1

        if (spica_mode == "IMA") & (completion_rate >= 0.5):
            flag = 1

        if (spica_mode == "TMP") & (completion_rate >= 0.5):
            flag = 1

        if (spica_mode == "DLD") & (0 < completion_rate < 0.5):
            flag = 2

        if (spica_mode == "IMA") & (0 < completion_rate < 0.5):
            flag = 2

        if (spica_mode == "TMP") & (0 < completion_rate < 0.5):
            flag = 2

        if completion_rate == 0:
            flag = 3

        if spica_mode == "SPI":
            flag = 3

        return flag

    def update_priority_final(self, flag_completion, priority_pi):
        """


        Parameters
        ----------
        flag_completion : TYPE
            DESCRIPTION.
        priority_pi : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        print(flag_completion, priority_pi)

        priority_final = []

        if (flag_completion == 1) & (priority_pi == 0):
            priority_final = 1

        if (flag_completion == 1) & (priority_pi == 1):
            priority_final = 1

        if (flag_completion == 2) & (priority_pi == 0):
            priority_final = 1

        if (flag_completion == 2) & (priority_pi == 1):
            priority_final = 2

        if flag_completion >= 3:
            priority_final = 3

        return priority_final

    # TODO
    def update_priority_final2(self, flag_completion, priority_pi, progname2):
        """


        Parameters
        ----------
        flag_completion : TYPE
            DESCRIPTION.
        priority_pi : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        priority_final = []

        if flag_completion == 1:
            priority_final = 1

        if (flag_completion == 2) & (priority_pi == 0):
            priority_final = 1

        if (flag_completion == 3) & (priority_pi == 0):
            priority_final = 2

        if (priority_pi == 1) & (progname2 == True):
            priority_final = 2

        if (flag_completion == 2) & (priority_pi == 1):
            priority_final = 3

        if (flag_completion == 3) & (priority_pi == 1):
            priority_final = 4

        return priority_final

    def add_toymodel(self, data):
        """

        Parameters
        ----------
        data : TYPE
            DESCRIPTION.

        Returns
        -------
        data : TYPE
            DESCRIPTION.

        """
        """
        # quickshow of relevant content for us
        if True:
            help(Models.gaussian)
            print(Models.gaussian('toto', fwhm=3))
            print([f for f in dir(Models) if not f.startswith("_")])
        """
        # Append dynamically a 'toy' model to our catalog
        # prepare new column
        COL_MODEL = "model"
        if COL_MODEL in data.colnames:
            data.remove_column(COL_MODEL)
        data.add_column(None, index=0, copy=False, name=COL_MODEL)
        data[COL_MODEL].meta["ucd"] = Models.SAMP_UCD_MODEL

        # loop to set a first disk model
        for row in data:
            diam = row["ld_jsdc2"]
            if diam:
                data[COL_MODEL][row.index] = Models.disk("disk2", diameter=diam)

        return data

    def add_targetmodel(self, data):
        """


        Parameters
        ----------
        data : TYPE
            DESCRIPTION.

        Returns
        -------
        data : TYPE
            DESCRIPTION.

        """
        """
        # quickshow of relevant content for us
        if True:
            help(Models.gaussian)
            print(Models.gaussian('toto', fwhm=3))
            print([f for f in dir(Models) if not f.startswith("_")])
        """

        # prepare new column if not present
        COLNAME_XMLMODEL = "model"
        if not COLNAME_XMLMODEL in data.colnames:
            data.add_column(None, index=0, copy=False, name=COLNAME_XMLMODEL)
            data[COLNAME_XMLMODEL].meta["ucd"] = Models.SAMP_UCD_MODEL

        # Let's search for a json model in model columns
        hasJsonModelColumn = COLNAME_XMLMODEL in data.colnames

        # loop to set a first disk model if model is empty
        for row in data:
            diam = row["diam"]  # ["ld_jsdc2"]
            row[COLNAME_XMLMODEL] = row[COLNAME_XMLMODEL].replace("}{", "},{")

            if hasJsonModelColumn:
                model = row[COLNAME_XMLMODEL]
            else:
                model = None

            if model:
                data[COLNAME_XMLMODEL][row.index] = _model(model)
            elif diam:
                data[COLNAME_XMLMODEL][row.index] = Models.disk("disk1", diameter=diam)
                # should be stored as [{type:disk, diameter:123},{type:punct, x=.3, y=-.4, diameter:123}]
            # else : # use other ld columns or a more sophisticated model ....
            else:
                pass

        # prefer str dtype
        data[COLNAME_XMLMODEL] = data[COLNAME_XMLMODEL].astype(str)

        return data

    # *** PRIMARY CALIBRATORS ***
    def query_calprim(self):
        """

        Returns
        -------
        None.

        """
        tableTargets = self.spica_catg[self.indexList_Targets]
        tableCalprim = None

        minRa = min(tableTargets["ra"]) - self.rarangeprim / 60 * 15
        maxRa = max(tableTargets["ra"]) + self.rarangeprim / 60 * 15
        if minRa < 0:
            minRa = 0
        if maxRa > 360:
            maxRa = 360
        minDec = np.median(tableTargets["dec"]) - self.decrangeprim
        maxDec = np.median(tableTargets["dec"]) + self.decrangeprim
        minVmag = np.median(tableTargets["vmag"]) - self.vmagrangeprim
        maxVmag = np.median(tableTargets["vmag"]) + self.vmagrangeprim

        print("[INFO] ---")
        print("[INFO] * Selection of the primary calibrators *")
        print(
            f"[INFO] RA range (deg): {self.rarangeprim/60*15} [{round(minRa,2)}, {round(maxRa,2)}]"
        )
        print(
            f"[INFO] DEC range (deg): {self.decrangeprim} [{round(minDec,2)}, {round(maxDec,2)}]"
        )
        print(
            f"[INFO] Vmag range: {self.vmagrangeprim} [{round(minVmag,2)}, {round(maxVmag,2)}]"
        )

        # Query catalog (TAP) of the primary calibrators
        tableCalprim = self.calquery_tap()

        # Filter primary calibrators on RA between sunset and sunrise
        if self.ra_sunset < self.ra_sunrise:
            tableCalprim = tableCalprim[
                (tableCalprim["ra"] > self.ra_sunset)
                & (tableCalprim["ra"] < self.ra_sunrise)
            ]
        else:
            tableCalprim = tableCalprim[
                (tableCalprim["ra"] > self.ra_sunset)
                | (tableCalprim["ra"] < self.ra_sunrise)
            ]

        # Filter primary calibrators on DEC and VMAG
        tableCalprim = tableCalprim[
            (tableCalprim["dec"] > minDec)
            & (tableCalprim["dec"] < maxDec)
            & (tableCalprim["vmag"] > minVmag)
            & (tableCalprim["vmag"] < maxVmag)
        ]

        if tableCalprim:
            # Remove duplicated targets
            JOIN_SAME_TARGET_ANGLE = 5 * u.arcsec
            coo_calprim = Table(
                [
                    SkyCoord(tableCalprim["ra"], tableCalprim["dec"], unit="deg"),
                    np.arange(len(tableCalprim)),
                ],
                names=["sc", "idx"],
            )
            coo_target = Table(
                [SkyCoord(tableTargets["ra"], tableTargets["dec"], unit="deg")],
                names=["sc"],
            )
            tablesJoined = join(
                coo_calprim,
                coo_target,
                join_funcs={"sc": join_skycoord(JOIN_SAME_TARGET_ANGLE)},
            )
            idxToRemove = np.unique(tablesJoined["idx"])
            if idxToRemove.size > 0:
                print(
                    f"[INFO] Duplicated primary in science targets ({len(tableTargets)}): {len(idxToRemove)}"
                )
                tableCalprim.remove_rows(idxToRemove)
            else:
                print(
                    f"[INFO] No duplicated primary in science targets at {JOIN_SAME_TARGET_ANGLE}."
                )

            # Replace mainID by HD name if needed
            tableCalprim = self.replacebyHDname(tableCalprim)

            self.calprim_catg = tableCalprim
            self.indexList_CalPrim = list(np.arange(np.size(tableCalprim)))
            print(
                f"[INFO] Final number of primary calibrators: {np.size(tableCalprim)}"
            )
        else:
            self.indexList_CalPrim = None
            print(
                "[Warning] There are no primary calibrators left after filtering in Dec and Vmag."
            )
        print("[INFO] ---")
        print("")

        self.plot_radec()

    def replacebyHDname(self, tableObjects):
        inspicadb = 0
        if "target_main_id" in tableObjects.colnames:
            inspicadb = 1
            tableObjects.rename_column("target_main_id", "Name")
        if "name" in tableObjects.colnames:
            tableObjects.rename_column("name", "Name")

        # Loop over all star names

        for j, starname in enumerate(tableObjects["Name"]):
            objectnames = (
                None  # Initialize objectnames to None at the start of each iteration
            )

            # Check if the star name already contains an HD number
            if "HD" not in starname:
                objectnames = Simbad.query_objectids(
                    starname
                )  # Query Simbad for alternate names

            if objectnames:
                if "ID" in objectnames.colnames:
                    objectnames.rename_column("ID", "id")

                # Filter for names containing "HD"
                hdname = [i for i in list(objectnames["id"]) if "HD" in i]
                if hdname:
                    oldname = tableObjects["Name"][j]
                    tableObjects["Name"][j] = hdname[0]
                    print(j, oldname, "->", tableObjects["Name"][j], "(", hdname, ")")
                else:
                    print(f"{j} No HD name found in Simbad for {starname}")
            else:  # Only print this if we actually queried Simbad
                print(f"{j} Already in HD name: {starname}")

        # Rename the column back to 'target_main_id' if it was originally that
        if inspicadb == 1:
            tableObjects.rename_column("Name", "target_main_id")

        return tableObjects

    def entryRaRangePrimCallback(self, strRaRangePrim):
        self.rarangeprim = float(strRaRangePrim.get())
        self.query_calprim()

    def entryDecRangePrimCallback(self, strDecRangePrim):
        self.decrangeprim = float(strDecRangePrim.get())
        self.query_calprim()

    def entryVmagRangePrimCallback(self, strVmagRangePrim):
        self.vmagrangeprim = float(strVmagRangePrim.get())
        self.query_calprim()

    def goCalPrim(self):
        self.query_calprim()

    def delCalPrim(self):
        print("[INFO] Delete primary calibrators")
        del self.indexList_CalPrim
        self.calprim_catg = None
        self.indexList_CalPrim = None
        self.plot_radec()

    def delCalSec(self):
        print("[INFO] Delete secondary calibrators")
        del self.indexList_CalSec
        self.calsec_catg = None
        self.indexList_CalSec = None
        self.plot_radec()

    def check4BeStars(self, tableObjects):
        rowsBestars = []
        for j, starname in enumerate(tableObjects["Name"]):
            result_table = Simbad.query_object(starname)
            result_table.rename_columns(
                result_table.colnames, [col.lower() for col in result_table.colnames]
            )

            if result_table["otype"] == "Be*":
                rowsBestars.append(j)
                print(
                    f"[INFO] {starname} is a Be star and is removed from the calibrators."
                )
        if rowsBestars:
            tableObjects.remove_rows(rowsBestars)

        return tableObjects

    def check4BadCal(self, tableObjects):
        rowsBadcal = []
        adqlQuery = "SELECT * FROM badcal"

        service = vo.dal.TAPService(tapServerUrl)
        results = service.search(adqlQuery)
        tableBadcal = results.to_table()

        JOIN_SAME_TARGET_ANGLE = 5 * u.arcsec
        coo_calsec = Table(
            [
                SkyCoord(tableObjects["ra"], tableObjects["dec"], unit="deg"),
                np.arange(len(tableObjects)),
            ],
            names=["sc", "idx"],
        )
        coo_badcal = Table(
            [SkyCoord(tableBadcal["ra"], tableBadcal["dec"], unit="deg")], names=["sc"]
        )
        tablesJoined = join(
            coo_calsec,
            coo_badcal,
            join_funcs={"sc": join_skycoord(JOIN_SAME_TARGET_ANGLE)},
        )
        idxToRemove = np.unique(tablesJoined["idx"])
        if idxToRemove.size > 0:
            print(f"[INFO] BadCal to remove ({len(tableObjects)}): {len(idxToRemove)}")
            tableObjects.remove_rows(idxToRemove)
        else:
            print(f"[INFO] No BadCal to remove at {JOIN_SAME_TARGET_ANGLE}.")

        return tableObjects

    # *** SECONDARY CALIBRATORS ***
    def query_calsec(self):
        if self.CalsecOpened:
            self.onCloseCalsec()

        tableTargets = self.spica_catg[self.indexList_Targets]
        if self.indexList_CalPrim:
            tableCalprim = self.calprim_catg[self.indexList_CalPrim]
        tableCalsec = None

        minRa = min(tableTargets["ra"]) - self.rarangesec / 60 * 15
        maxRa = max(tableTargets["ra"]) + self.rarangesec / 60 * 15
        if minRa < 0:
            minRa = 0
        if maxRa > 360:
            maxRa = 360
        minDec = np.median(tableTargets["dec"]) - self.decrangesec
        maxDec = np.median(tableTargets["dec"]) + self.decrangesec
        if self.vmagrangesec == 0:
            self.vmagrangesec = 0.1
        minVmag = np.median(tableTargets["vmag"]) - self.vmagrangesec
        maxVmag = np.median(tableTargets["vmag"]) + self.vmagrangesec

        print("[INFO] ---")
        print("[INFO] ** Selection of the secondary calibrators **")
        print(
            f"[INFO] RA range (deg): {self.rarangesec/60*15} [{round(minRa,2)}, {round(maxRa,2)}]"
        )
        print(
            f"[INFO] DEC range (deg): {self.decrangesec} [{round(minDec,2)}, {round(maxDec,2)}]"
        )
        print(
            f"[INFO] Vmag range: {self.vmagrangesec} [{round(minVmag,2)}, {round(maxVmag,2)}]"
        )
        print(f"[INFO] Max. LDD Chi: {self.lddchisec}")
        print(f"[INFO] Max. rel. error (%): {self.relerrorsec}")
        print(f"[INFO] Min. vis2: {self.minvissec}")

        # Query JSDC2 (Vizier) for a list of secondary calibrators based on DEC and Vmag
        # TODO move next constant in a common part ?
        CALSEC_VIZIER_ROW_LIMIT = -1
        JOIN_SAME_TARGET_ANGLE = 5 * u.arcsec

        """tableCalsec = Vizier(catalog = "II/346/jsdc_v2",columns=['*', '_RAJ2000', '_DEJ2000', '+LDD'], column_filters={
                                          'RAJ2000' : f' > {minRa} & < {maxRa}',
                                          'DEJ2000' : f' > {minDec} & < {maxDec}',
                                          'Vmag' : f' > {minVmag} & < {maxVmag}',
                                          'CalFlag' : '0',
                                          'LDDCHI' : f' < {self.lddchisec}',
                                          }, row_limit=CALSEC_VIZIER_ROW_LIMIT).query_constraints()[0]"""

        tableCalsec = Vizier(
            catalog="II/346/jsdc_v2",
            columns=["*", "_RAJ2000", "_DEJ2000", "+LDD"],
            column_filters={
                "DEJ2000": f" > {minDec} & < {maxDec}",
                "Vmag": f"< {maxVmag}",
                "CalFlag": "0",
                "LDDCHI": f" < {self.lddchisec}",
            },
            row_limit=CALSEC_VIZIER_ROW_LIMIT,
        ).query_constraints()[
            0
        ]  # > {minVmag} &

        tableCalsec.rename_column("_DEJ2000", "dec")
        tableCalsec.rename_column("_RAJ2000", "ra")

        """
        if len(tableCalsec) == 0: #>= CALSEC_VIZIER_ROW_LIMIT:
            self.CalsecOpened = True
            self.topCalsec = Toplevel(self.root)
            self.topCalsec.resizable(FALSE, FALSE)
            Label(self.topCalsec, text= 'WARNING', font=self.myFont, bg='red', fg='white').grid(column=0,row=0,padx=5,pady=4)
            Label(self.topCalsec, text= f'Too many secondary calibrators in selection ({CALSEC_VIZIER_ROW_LIMIT})', font=self.myFont).grid(column=0,row=1,padx=5,pady=4)
            Label(self.topCalsec, text= 'Please, refine your selection.', font=self.myFont).grid(column=0,row=2,padx=5,pady=4)
            frame = Frame(self.topCalsec)
            frame.grid(pady=5)
            Button(frame, text='Close', font=self.myFont, fg='white', bg='red', cursor='hand1', command=self.onCloseCalsec).grid(column=1,row=3)
        else:
        """

        # Filter secondary calibrators on RA between sunset and sunrise
        if self.ra_sunset < self.ra_sunrise:
            tableCalsec = tableCalsec[
                (tableCalsec["ra"] > self.ra_sunset)
                & (tableCalsec["ra"] < self.ra_sunrise)
            ]
        else:
            tableCalsec = tableCalsec[
                (tableCalsec["ra"] > self.ra_sunset)
                | (tableCalsec["ra"] < self.ra_sunrise)
            ]

            """
            indCal1 = [(tableCalsec['ra'] > minRa) & (tableCalsec['ra'] < self.ra_sunrise)]
            indCal2 = [(tableCalsec['ra'] > self.ra_sunset)]# & (tableCalsec['ra'] < maxRa)]
            indCal = np.logical_or(indCal1, indCal2)
        """
        # tableCalsec = tableCalsec[np.where(indCal)[1][:]]

        # Calculate the visibility
        wavel = 0.75
        maxbaseline = self.intMaxBaseline.get()  # 330.66
        diam = tableCalsec["UDDR"] * 1.0e-3 * np.pi / (3600.0 * 180.0)
        zstar = np.pi * diam * maxbaseline / (wavel * 1.0e-6)
        vis2 = (2.0 * jv(1, zstar) / zstar) ** 2.0
        print(
            f"[INFO] Selected MaxBaseline {maxbaseline} m for minVis2 {self.minvissec}"
        )

        # Filter in diameter relative error and on minimum visibility
        filterSec = tableCalsec["Name"].mask
        filterSec = filterSec | (
            (tableCalsec["e_LDD"] / tableCalsec["UDDR"] < self.relerrorsec / 100)
            & (vis2 > self.minvissec)
        )

        if np.any(filterSec):
            tableCalsec = tableCalsec[filterSec]

            # Remove duplicated targets
            coo_calsec = Table(
                [
                    SkyCoord(tableCalsec["ra"], tableCalsec["dec"], unit="deg"),
                    np.arange(len(tableCalsec)),
                ],
                names=["sc", "idx"],
            )
            coo_target = Table(
                [SkyCoord(tableTargets["ra"], tableTargets["dec"], unit="deg")],
                names=["sc"],
            )
            tablesJoined = join(
                coo_calsec,
                coo_target,
                join_funcs={"sc": join_skycoord(JOIN_SAME_TARGET_ANGLE)},
            )
            idxToRemove = np.unique(tablesJoined["idx"])
            if idxToRemove.size > 0:
                print(
                    f"[INFO] Duplicated secondary in science targets ({len(tableTargets)}): {len(idxToRemove)} (out of {len(tableCalsec)})"
                )
                tableCalsec.remove_rows(idxToRemove)
                coo_calsec = Table(
                    [
                        SkyCoord(tableCalsec["ra"], tableCalsec["dec"], unit="deg"),
                        np.arange(len(tableCalsec)),
                    ],
                    names=["sc", "idx"],
                )
            else:
                print(
                    f"[INFO] No duplicated secondary in science targets at {JOIN_SAME_TARGET_ANGLE}"
                )

            # Remove duplicated primary calibrators
            if self.indexList_CalPrim:
                coo_calprim = Table(
                    [SkyCoord(tableCalprim["ra"], tableCalprim["dec"], unit="deg")],
                    names=["sc"],
                )
                tablesJoined = join(
                    coo_calsec,
                    coo_calprim,
                    keys="sc",
                    join_funcs={"sc": join_skycoord(JOIN_SAME_TARGET_ANGLE)},
                )
                idxToRemove = np.unique(tablesJoined["idx"])
                if idxToRemove.size > 0:
                    print(
                        f"[INFO] Duplicated secondary in primary calibrators ({len(tableCalprim)}): {len(idxToRemove)} (out of {len(tableCalsec)})"
                    )
                    tableCalsec.remove_rows(idxToRemove)
                else:
                    print(
                        f"[INFO] No duplicated secondary in primary calibrators at {JOIN_SAME_TARGET_ANGLE}"
                    )

            if len(tableCalsec) > 50:
                self.CalsecOpened = True
                self.topCalsec = Toplevel(self.root)
                self.topCalsec.resizable(FALSE, FALSE)
                Label(
                    self.topCalsec,
                    text="WARNING",
                    font=self.myFont,
                    bg="red",
                    fg="white",
                ).grid(column=0, row=0, padx=5, pady=4)
                Label(
                    self.topCalsec,
                    text=f"Too many secondary calibrators in selection ({str(len(tableCalsec))})",
                    font=self.myFont,
                ).grid(column=0, row=1, padx=5, pady=4)
                Label(
                    self.topCalsec,
                    text="Please, refine your selection to 50.",
                    font=self.myFont,
                ).grid(column=0, row=2, padx=5, pady=4)
                frame = Frame(self.topCalsec)
                frame.grid(pady=5)
                Button(
                    frame,
                    text="Close",
                    font=self.myFont,
                    fg="white",
                    bg="red",
                    cursor="hand1",
                    command=self.onCloseCalsec,
                ).grid(column=1, row=3)
            else:
                # Replace by mainID by HD name if needed
                tableCalsec = self.replacebyHDname(tableCalsec)

                # Check for Be stars
                tableCalsec = self.check4BeStars(tableCalsec)

                # Check for BadCalibrators
                tableCalsec = self.check4BadCal(tableCalsec)

                self.calsec_catg = tableCalsec
                self.indexList_CalSec = np.arange(np.size(tableCalsec))
                print(
                    f"[INFO] Final number of secondary calibrators: {np.size(tableCalsec)}"
                )
        else:
            self.indexList_CalSec = None
            print(
                "[Warning] There are no secondary calibrators left after filtering in Relative Error and Visibility."
            )
        # else:
        #    print('[Warning] There are no secondary calibrators in initial Vizier query.')
        print("[INFO] ---")
        print("")
        self.plot_radec()

    def entryRaRangeSecCallback(self, strRaRangeSec):
        self.rarangesec = float(strRaRangeSec.get())
        self.query_calsec()

    def entryDecRangeSecCallback(self, strDecRangeSec):
        self.decrangesec = float(strDecRangeSec.get())
        self.query_calsec()

    def entryVmagRangeSecCallback(self, strVmagRangeSec):
        self.vmagrangesec = float(strVmagRangeSec.get())
        self.query_calsec()

    def entryLDDChiSecCallback(self, strLDDChiSec):
        self.lddchisec = float(strLDDChiSec.get())
        self.query_calsec()

    def entryRelErrorSecCallback(self, strRelErrorSec):
        self.relerrorsec = float(strRelErrorSec.get())
        self.query_calsec()

    def entryMinVisSecCallback(self, strMinVisSec):
        self.minvissec = float(strMinVisSec.get())
        self.query_calsec()

    def entryMinVisSecCallback2(self, strMinVisSec):
        # self.minvissec = float(strMinVisSec.get())
        self.query_calsec()

    def goCalSec(self):
        self.query_calsec()

    def getSelectedTargets(self):
        self.indexList_Targets = reduce(
            np.intersect1d,
            [
                self.indexProgName,
                self.indexInstMode,
                self.indexFinalPriority,
                self.indexDecMin,
                self.indexDecMax,
                self.indexVmagMin,
                self.indexVmagMax,
            ],
        )
        if self.index_AddTarget is not None:
            self.indexList_Targets = np.append(
                self.indexList_Targets, self.index_AddTarget
            )

    def getAddTarget(self):
        self.indexList_AddTarget = reduce(
            np.intersect1d,
            [
                self.indexProgName,
                self.indexInstMode,
                self.indexFinalPriority,
                self.indexDecMin,
                self.indexDecMax,
                self.indexVmagMin,
                self.indexVmagToAddTarget,
            ],
        )

    """
    def getSelectedTargets_with_AddTarget(self, indexAddTarget):
        self.getSelectedTargets()
        self.indexList_Targets = np.append(self.indexList_Targets, indexAddTarget)
        pdb.set_trace()
    """
    # def on_closing(self):
    #     if messagebox.askokcancel("Quit", "Do you want to quit?"):
    #         self.tree_frameInfoTargets.destroy()
    #         self.popupInfoTargets = False

    # Class Variables
    # date="2020-10-10" #date.today()
    CalsecOpened = False
    LogOpened = False
    popupInfoTargets = False
    popupAddTarget = False
    popupBestDec = False
    topSAMP = False

    date = date.today()
    date = date.strftime("%Y-%m-%d")
    # date = Time(date, scale='utc')

    # nstmodeName = ['DIA', 'DLD', 'IMA', 'TMP', 'SPI']
    spica_catg = None
    index_AddTarget = None

    calprim_catg = None
    indexList_CalPrim = None
    calsec_catg = None
    indexList_CalSec = None

    decmin = -30.0
    decmax = 90.0
    decmin_default = decmin
    decmax_default = decmax

    vmagmin = -3.0
    vmagmax = 13.0
    vmagmin_default = vmagmin
    vmagmax_default = vmagmax

    rarangeprim = 60
    decrangeprim = 5
    vmagrangeprim = 2
    rarangeprim_default = rarangeprim
    decrangeprim_default = decrangeprim
    vmagrangeprim_default = vmagrangeprim

    rarangesec = 60
    decrangesec = 2
    vmagrangesec = 2
    lddchisec = 2
    relerrorsec = 10
    minvissec = 0.7
    maxbaseline = 330
    rarangesec_default = rarangesec
    decrangesec_default = decrangesec
    vmagrangesec_default = vmagrangesec
    lddchisec_default = lddchisec
    relerrorsec_default = relerrorsec
    minvissec_default = minvissec
    maxbaseline_default = maxbaseline

    catalogSpica = "spica"

    strDecMean = 0.0


oidbTapUrl = "http://tap.jmmc.fr/vollt/tap/"
tapServerUrl = "http://tap.jmmc.fr/vollt/tap/"
nssVersion = "v2025-09"

if __name__ == "__main__":
    # Start NSS gui
    app = spica_NSS()
