#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on december 2023

@author: dsalabert
"""

import json
import platform
import tempfile
import warnings

from datetime import date, datetime
from functools import reduce
from idlelib.tooltip import Hovertip
from pathlib import Path
from tkinter import *
from tkinter import messagebox, simpledialog, ttk

import astropy.io.votable
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
from astropy.table import MaskedColumn, Table, join, join_skycoord, vstack
from astropy.time import Time
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from scipy.special import jv

if platform.system() == "Darwin":
    from tkmacosx import Button

Simbad.add_votable_fields("otype")
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message=".*partition.*mask.*")
plt.ioff()


class spica_NSS:
    """
    SPICA-NSS graphical observing preparation tool.

    The application provides interactive selection of SPICA
    science targets, calibrator identification, observability
    analysis for CHARA, and interoperability with JMMC Aspro2
    through SAMP and VOtable exchanges.

    Main Capabilities:
        * Query SPICA targets through TAP/ADQL
        * Filter targets by work package, observing mode and priority
        * Search primary and secondary calibrators
        * Visualize sky distributions and observing statistics
        * Export observing lists to Aspro2
    """

    def __init__(self):
        """
        Initialize the SPICA-NSS GUI application.

        Sets up the main Tkinter window, creates all interface frames
        (date selector, action buttons, work-package and instrumental-mode
        filters, priority selectors, science-object filter entries, primary
        and secondary calibrator parameter entries, and the matplotlib canvas),
        configures layout via ``grid``, and starts the Tkinter event loop.

        Class-level default values for declination and magnitude ranges as
        well as calibrator search parameters are used to pre-populate all
        entry widgets.

        Returns:
            None
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
        )
        self.root.eval("::msgcat::mclocale en")  # Message in english

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

        # Workpackages (checkbuttons)
        labelProgName = Label(
            FrameWorkPackages,
            text="ProgNames",
            font=self.myFont,
            fg="black",
            bg="orange",
            width=15,
        )

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
        )

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
        )
        self.entryMaxBaseline.set(self.maxbaseline)
        self.entryMaxBaseline.bind(
            "<ButtonRelease-1>",
            (lambda _: self.entryMinVisSecCallback2(self.entryMaxBaseline)),
        )

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

        FrameAddTarget.rowconfigure(0, weight=1)
        FrameAddTarget.columnconfigure(0, weight=1)
        FrameAddTarget.grid(sticky=E, row=4, column=1, padx=15, pady=1, ipady=1)

        FrameInstModes.grid(sticky=W, row=2, column=0, columnspan=2, padx=15, ipady=1)
        FramePriorities.grid(
            sticky=W, row=3, column=0, columnspan=2, padx=15, pady=1, ipady=1
        )

        self.FrameCalPrims.grid(
            sticky=W, row=5, column=0, columnspan=2, padx=15, ipadx=3, pady=1, ipady=1
        )
        self.FrameCalSecs.grid(
            sticky=W, row=6, column=0, columnspan=2, padx=15, ipadx=3, pady=1, ipady=1
        )
        FrameLog.grid(sticky=W, row=8, column=0, columnspan=2, padx=15, pady=1, ipady=1)

        entryDate.grid(column=0, row=0, padx=5)
        self.labelValDate.grid(column=1, row=0, padx=5)


        buttonQuery.grid(column=0, row=0, padx=5)
        buttonBestDec.grid(column=1, row=0, padx=5)
        buttonInfoTargets.grid(column=2, row=0, padx=5)
        buttonAspro.grid(column=3, row=0, padx=5)
        buttonReset.grid(column=4, row=0, padx=5)
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

        labelRaRangePrim.grid(column=0, row=0, padx=5)
        self.entryRaRangePrim.grid(column=1, row=0)
        labelDecRangePrim.grid(column=2, row=0, padx=5)
        self.entryDecRangePrim.grid(column=3, row=0)
        labelVmagRangePrim.grid(column=4, row=0, padx=5)
        self.entryVmagRangePrim.grid(column=5, row=0)

        Button(
            self.FrameCalPrims,
            text="UNDO",
            font=self.myFont,
            fg="white",
            bg="red",
            command=self.delCalPrim,
            cursor="hand1",
        ).grid(column=6, row=0, padx=5)

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

        Button(
            self.FrameCalSecs,
            text="UNDO",
            font=self.myFont,
            fg="white",
            bg="red",
            command=self.delCalSec,
            cursor="hand1",
        ).grid(column=6, row=0, rowspan=1, padx=5)

        self.root.mainloop()

    def open_popupAddTarget(self):
        """Add an arbitrary science target to the current selection.

        Prompts the user for a target name, resolves it through SIMBAD to
        obtain sky coordinates, and then cross-matches against the SPICA
        catalog.  If a unique match is found the target is added to
        ``self.index_AddTarget`` and the plot is refreshed.  If multiple
        matches are found a Treeview popup lets the user pick the desired
        programme entry.  If the target is not in the observable window for
        the selected date but exists in the full database, a warning is
        shown.  An error dialog is raised when SIMBAD cannot resolve the
        name.

        Returns:
            None
        """
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

    def open_popupInfoTargets(self):
        """
        Open a popup window listing the currently selected science targets.

        Creates a ``Toplevel`` Treeview window displaying key properties of
        every target in ``self.indexList_Targets`` (SPICA-DB ID, main
        identifier, spectral type, programme name, SPICA mode, final
        priority, completion rate, RA, Dec, angular diameter, V magnitude
        and H magnitude).  If the popup is already open it is first
        destroyed and rebuilt so the contents are always up to date.

        Column headings for Completion Rate, Ra, Dec, Diameter, Vmag and
        Hmag are clickable and trigger ascending/descending sort via
        :meth:`treeview_sort_column`.  The window title shows the total
        number of currently selected targets.

        Returns:
            None
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
        self.my_treeInfoTargets.heading(
            "Final Priority", text="Final Priority", anchor=CENTER
        )
        self.my_treeInfoTargets.heading(
            "Completion Rate", text="Completion Rate", anchor=CENTER
        )
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

        # Insert InfoTargets in popup
        self.insert_popupInfoTargets()

        #
        self.tree_frameInfoTargets.protocol(
            "WM_DELETE_WINDOW", self.closeframeInfoTargets
        )

    def closeframeInfoTargets(self):
        """Close the target-info popup and reset its open flag.

        Returns:
            None
        """
        self.tree_frameInfoTargets.destroy()
        self.popupInfoTargets = False

    def closeframeAddTarget(self):
        """Close the add-target popup and reset its open flag.

        Returns:
            None
        """
        self.tree_frameAddTarget.destroy()
        self.popupAddTarget = False

    # Define a function to clear all the items present in Treeview
    def clear_all(self, tree):
        """Remove all rows from a Treeview widget.

        Args:
            tree (ttk.Treeview): The Treeview whose contents should be cleared.

        Returns:
            None
        """
        for item in tree.get_children():
            tree.delete(item)

    def copy(self, event):
        """Copy the selected Treeview rows to the system clipboard.

        Copies column headers followed by the values of every selected row,
        tab-separated, to the Tkinter clipboard so the content can be pasted
        into external applications.  Bound to ``<Control-c>`` on the target
        info Treeview.

        Args:
            event: The Tkinter key-press event that triggered this callback.

        Returns:
            None
        """
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
        """Populate the add-target Treeview with candidate fainter stars.

        Fetches the subset of the SPICA catalog indexed by
        ``self.indexList_AddTarget``, sorts it by V magnitude, and inserts
        each record as a striped row into ``self.my_treeAddTarget``.  For
        targets belonging to a second programme the programme name and final
        priority columns show both values in parentheses.  Also binds
        ``<Control-c>`` and ``<Double-1>`` events on the Treeview.

        Returns:
            None
        """
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
        """Populate the target-info Treeview with the currently selected targets.

        Fetches the subset of the SPICA catalog indexed by
        ``self.indexList_Targets``, sorts it by V magnitude, and inserts each
        record as a striped row into ``self.my_treeInfoTargets``.  For targets
        belonging to a second programme the programme name and final priority
        columns show both values in parentheses.  Also binds ``<Control-c>``
        and ``<Double-1>`` events on the Treeview.

        Returns:
            None
        """
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
        """Handle a double-click on a row in the add-target Treeview.

        Retrieves the SPICA-DB ID, main identifier and programme name of the
        focused row and asks the user for confirmation.  If confirmed, the
        corresponding catalog index is appended to ``self.index_AddTarget``,
        the target selection and sky plot are refreshed, and the add-target
        popup is closed.

        Args:
            event: The Tkinter double-click event that triggered this callback.

        Returns:
            None
        """
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
        """Reapply alternating row tags after a Treeview sort.

        Iterates over all top-level items of the given Treeview and
        alternately assigns the ``"oddrow"`` and ``"evenrow"`` tags so that
        the striped background pattern is preserved after rows are reordered.

        Args:
            theTreeToSort (ttk.Treeview): The Treeview widget to retag.

        Returns:
            None
        """
        tag = "oddrow"
        for iid in theTreeToSort.get_children(""):
            tag = "oddrow" if tag == "evenrow" else "evenrow"
            theTreeToSort.item(iid, tags=(tag,))

    def treeview_sort_column(self, tv, col, reverse):
        """Sort a Treeview by the given column and toggle sort direction.

        Extracts the values of ``col`` for all rows, sorts them (converting
        to ``int`` for the ``"Count"`` column), rearranges the rows
        accordingly, reapplies striped row tags via :meth:`retag`, and
        reconfigures the column heading so the next click reverses the sort.

        Args:
            tv (ttk.Treeview): The Treeview widget to sort.
            col (str): The column identifier to sort by.
            reverse (bool): If ``True``, sort in descending order.

        Returns:
            None
        """
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
        Open a popup window showing the target count per declination band.

        Bins all observable science targets (from ``self.targetListInit``)
        into 5-degree declination strips from −30° to 90°, applying the
        current V-magnitude limits.  The result is displayed in a sortable
        ``Toplevel`` Treeview with three columns: Dec Low, Dec High, and
        Count.  The Count column is clickable for ascending/descending sort.

        If the popup is already open it is first destroyed and rebuilt.

        Returns:
            None
        """
        if self.spica_catg is None:
            messagebox.showinfo(title="Info", message="Please query SPICA-DB first.")
            return

        if self.popupBestDec:
            self.tree_frameBestDec.destroy()
            self.popupBestDec = True
        else:
            self.popupBestDec = True

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
        Update the target selection when the programme-name filter changes.

        Iterates over all checkbutton states in ``self.SelectedProgName`` and
        rebuilds ``self.indexProgName`` as the union of row indices matching
        any enabled programme code (matched against both ``progname`` and
        ``progname2`` columns with hyphen-separated splitting).  The combined
        index list is then propagated through :meth:`getSelectedTargets`,
        :meth:`getAddTarget` and :meth:`plot_radec`, and any open target
        info or add-target popup windows are refreshed.

        Returns:
            None
        """
        self.indexProgName = []
        self.indexProgName2 = []
        for j in self.SelectedProgName:
            progname = j.get()
            if progname != str(0):
                # TODO build a filter function to apply on the table that will split every progname to search for a match
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

    def plot_radec(self):
        """
        Plots Right Ascension (RA) and Declination (Dec) of observable targets with different priorities and modes.

        Returns:
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
            modeltype.append(json.loads(a)[0]["type"])
            if modeltype[j] == "disk":
                modeldiam.append(json.loads(a)[0]["diameter"])
            elif modeltype[j] == "elong_disk":
                modeldiam.append(json.loads(a)[0]["minor_axis_diameter"])
        modeldiam = MaskedColumn(modeldiam)
        list_targets["modeldiam"] = modeldiam

        fig, self.ax = plt.subplots(
            3, figsize=(11, 6), sharex=True
        )
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
        if self.indexList_CalSec is not None:
            nb_calsec = len(self.indexList_CalSec)
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


    def open_popupProgName(self):
        """
        Open a popup window for selecting SPICA work-package programme names.

        Derives the unique programme codes available under the currently
        active instrumental-mode filter and presents them as a set of
        checkbuttons inside a ``Toplevel`` window.  Buttons labelled
        SELECT ALL, SUBMIT, RESET and CLOSE are also provided.

        Returns:
            None
        """
        top = Toplevel(self.root)
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
        Select all programme-name checkbuttons in the current popup.

        Iterates over every button in ``self.buttonProgName`` and calls its
        ``select()`` method, effectively enabling all available work packages.

        Returns:
            None
        """
        for j in self.buttonProgName:
            j.select()

    def deselect_allProgName(self):
        """Reset the programme-name filter by re-running the catalog query.

        Delegates to :meth:`onQuery`, which reloads the catalog from the
        cached data and resets all filters including programme names to their
        default (all selected) state.

        Returns:
            None
        """
        self.onQuery()

    def plotSelectedInstMode(self):
        """
        Update the target selection when the instrumental-mode filter changes.

        Rebuilds ``self.indexInstMode`` as the union of row indices matching
        any enabled SPICA mode (DIA, DLD, IMA, TMP, SPI) selected in
        ``self.SelectedInstMode``.  The combined index is propagated through
        :meth:`getSelectedTargets`, :meth:`getAddTarget` and
        :meth:`plot_radec`, and any open target info or add-target popups are
        refreshed.

        Returns:
            None
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
        Update the target selection when the final-priority filter changes.

        Rebuilds ``self.indexFinalPriority`` as the union of row indices
        matching any enabled priority level (1–4) selected in
        ``self.SelectedFinalPriority``.  The combined index is propagated
        through :meth:`getSelectedTargets`, :meth:`getAddTarget` and
        :meth:`plot_radec`, and any open target info or add-target popups are
        refreshed.

        Returns:
            None
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
        Open a popup window for selecting SPICA instrumental modes.

        Derives the unique instrumental modes available under the current
        programme-name filter and presents them as checkbuttons inside a
        ``Toplevel`` window.  SPI is always placed last.  SUBMIT and CLOSE
        buttons are also provided.

        Returns:
            None
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
        """Handle a change in the observation date entry.

        Updates ``self.date``, refreshes the date label, and triggers a full
        catalog re-query via :meth:`onQuery`.

        Args:
            strDate (tk.Entry): The entry widget holding the new date string.

        Returns:
            None
        """
        self.date = strDate.get()
        self.labelValDate.config(text=self.date)
        self.onQuery()

    def entryDecMinCallback(self, strDecMin):
        """Handle a change in the minimum declination filter entry.

        Updates ``self.decmin`` (clamped to −30°), rebuilds
        ``self.indexDecMin``, and refreshes the target selection and plot.
        Any open target info or add-target popups are also updated.

        Args:
            strDecMin (tk.Entry): Entry widget holding the new minimum
                declination value in degrees.

        Returns:
            None
        """
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
        """Handle a change in the maximum declination filter entry.

        Updates ``self.decmax`` (clamped to 90°), rebuilds
        ``self.indexDecMax``, and refreshes the target selection and plot.
        Any open target info or add-target popups are also updated.

        Args:
            strDecMax (tk.Entry): Entry widget holding the new maximum
                declination value in degrees.

        Returns:
            None
        """
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
        """Handle a change in the minimum V-magnitude filter entry.

        Updates ``self.vmagmin``, rebuilds ``self.indexVmagMin``, and
        refreshes the target selection and plot.  Any open target info or
        add-target popups are also updated.

        Args:
            strVmagMin (tk.Entry): Entry widget holding the new minimum
                V magnitude.

        Returns:
            None
        """
        self.indexVmagMin = []
        self.vmagmin = float(strVmagMin.get())
        self.indexVmagMin = list(
            filter(
                lambda x: self.spica_catg["vmag"][x] > self.vmagmin,
                range(len(self.spica_catg)),
            )
        )

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
        """Handle a change in the maximum V-magnitude filter entry.

        Updates ``self.vmagmax``, rebuilds ``self.indexVmagMax`` for the
        main selection and ``self.indexVmagToAddTarget`` for the one-magnitude
        fainter add-target pool, then refreshes the target selection and plot.
        Any open target info or add-target popups are also updated.

        Args:
            strVmagMax (tk.Entry): Entry widget holding the new maximum
                V magnitude.

        Returns:
            None
        """
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
        """Open the survey statistics (LOG) popup.

        Queries the full SPICA database, computes the completion percentage
        for each programme name and each instrumental mode, and displays the
        results as two stacked bar charts (programme completion and mode
        completion) in a ``Toplevel`` window with annotated bar segments.
        If the LOG window is already open it is first closed and rebuilt.

        Returns:
            None
        """
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

            print(
                f"[INFO] Program: {prog_name} -> Total: {nb_total}, Completed: {nb_completed}"
            )

            if nb_total > 0:
                perc_Progcompleted[i] = nb_completed / nb_total * 100

        uniqInstMode = list(np.unique(spicadb_load["spica_mode"]))
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

            print(
                f"[INFO] Inst. mode: {mode} -> Total: {nb_total}, Completed: {nb_completed}"
            )

            if nb_total > 0:
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
        """Quit the SPICA-NSS application after user confirmation.

        Shows a yes/no dialog; if confirmed, destroys the main Tkinter window,
        closes all matplotlib figures, and prints a farewell message.

        Returns:
            None
        """
        if messagebox.askyesno("Exit", "Do you want to quit the SPICA-NSS Tool?"):
            self.root.destroy()
            plt.close("all")
            print("\nBye!")

    def onCloseCalsec(self):
        """Close the secondary-calibrator warning popup and reset its flag.

        Returns:
            None
        """
        self.topCalsec.destroy()
        self.CalsecOpened = False

    def onCloseLog(self):
        """Close the LOG statistics popup and reset its flag.

        Returns:
            None
        """
        self.topLog.destroy()
        self.LogOpened = False

    def onAspro(self):
        """
        Trigger the export of the current selection to Aspro2 via SAMP.

        Validates that the SPICA catalog has been loaded, then delegates to
        :meth:`import2aspro` which sorts the targets by meridian transit
        time, attaches diameter models, assigns NSS type labels and calls
        :meth:`samp_votable` to broadcast the VOtable to any running Aspro2
        instance.

        Returns:
            None
        """
        if self.spica_catg is None:
            messagebox.showinfo(title="Info", message="Please query SPICA-DB first.")
            return

        self.import2aspro()

    def import2aspro(self):
        """Prepare and dispatch the current selection to Aspro2.

        Sorts the selected science targets by meridian transit time at CHARA,
        gathers any primary and secondary calibrators, attaches angular
        diameter models (:meth:`addDiamModel`) and NSS type labels
        (:meth:`addNssType`), then calls :meth:`samp_votable` to broadcast
        the assembled VOtable to a running Aspro2 instance over SAMP.

        Returns:
            None
        """
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
        Add a unified ``diam`` column to targets and calibrator tables.

        Parses the JSON model stored in the ``model`` column of each science
        target to extract the angular diameter (``diameter`` for uniform
        disks or ``minor_axis_diameter`` for elongated disks) and stores it
        in a new ``diam`` column.  For primary calibrators the ``ldd`` column
        is used, and for secondary calibrators the ``UDDR`` column is used.

        Args:
            targets (astropy.table.Table): Science target table; must contain a ``model`` JSON column.
            calibrators1 (astropy.table.Table, optional): Primary calibrator table; must contain an ``ldd`` column.
            calibrators2 (astropy.table.Table, optional): Secondary calibrator table; must contain a ``UDDR`` column.

        Returns:
            tuple of astropy.table.Table: The three input tables with a ``diam`` column added in-place.
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
        Add an ``nss_type`` classification column to each table.

        Marks each row with a string label indicating its role in the
        observing sequence: ``"Science"`` for science targets,
        ``"CalPrim"`` for primary calibrators, and ``"CalSec"`` for
        secondary calibrators.

        Args:
            targets (astropy.table.Table): Science target table.
            calibrators1 (astropy.table.Table, optional): Primary calibrator table.
            calibrators2 (astropy.table.Table, optional): Secondary calibrator table.

        Returns:
            tuple of astropy.table.Table: The three input tables with the ``nss_type`` column added.
        """
        targets["nss_type"] = ["Science"] * len(targets)
        if calibrators1:
            calibrators1["nss_type"] = ["CalPrim"] * len(calibrators1)
        if calibrators2:
            calibrators2["nss_type"] = ["CalSec"] * len(calibrators2)

        return targets, calibrators1, calibrators2

    def onQuery(self):
        """
        Query the SPICA database and initialise the GUI state.

        On the first call (``self.iter == 0``) the full SPICA catalog is
        fetched via :meth:`dbquery_tap` and targets with a completion rate of
        1 (fully observed) are excluded.  On subsequent calls the previously
        downloaded catalog is reused, avoiding redundant network requests,
        and all filter entry widgets are reset to their default values.

        The observable RA window for the selected date is computed from
        CHARA sunset/sunrise times.  Each target's final priority and
        completion flag are then calculated and programme-name, mode and
        priority checkbuttons are configured accordingly.  Object names are
        optionally replaced by their HD identifiers before the sky plot is
        refreshed via :meth:`plot_radec`.

        Returns:
            None
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

            tableTargets["priority_final"][i] = self.update_priority_final2(
                self.flag_completion[i],
                tableTargets["priority_pi"][i],
                tableTargets["progname2"][i],
            )
        for iProgName in list(range(len(self.ProgName))):
            if np.any([(tableTargets["progname"] == self.ProgName[iProgName])]):
                self.buttonProgName[iProgName]["state"] = "normal"
                self.buttonProgName[iProgName].select()

        for iInstMode in list(range(len(self.InstMode))):
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
        """Re-apply all active filters and refresh the sky plot.

        Delegates to :meth:`onQuery` to recompute the filtered target list,
        then calls :meth:`plot_radec` to update the display.

        Returns:
            None
        """
        self.onQuery()
        self.plot_radec()

    def onReset(self):
        """
        Reset the application to its initial state.

        Forces a fresh query of the SPICA database by resetting
        ``self.iter`` to 0 before calling :meth:`onQuery`.  All entry
        widgets for science-object filters and calibrator search parameters
        are restored to their class-level defaults, and any open target info
        or best-declination popup windows are closed.

        Does nothing if the SPICA catalog has not yet been loaded.

        Returns:
            None
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

        if self.popupInfoTargets:
            self.tree_frameInfoTargets.destroy()
            self.popupInfoTargets = False

        if self.popupBestDec:
            self.tree_frameBestDec.destroy()
            self.popupBestDec = False

    def onDRS(self):
        """
        Simulate a DRS/QCS pipeline reduction on the selected targets.

        Iterates over the currently selected science targets and increments
        their completion rate (via :meth:`update_completion_rate`), their OB
        reference counter and resets their QCS flag to 0.  The local cached
        catalog (``self.spica_catg_jmmc``) is then updated accordingly so
        that subsequent queries reflect the simulated observations without
        hitting the network.

        This method is intended as a development aid to test the evolution of
        target priorities after hypothetical DRS/QCS runs.

        Returns:
            None
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
        Query the SPICA catalog through the TAP service.

        The full content of the SPICA database table defined by
        ``self.catalogSpica`` is retrieved using an ADQL query.
        Masked values in selected columns are converted to standard
        Python numerical values to simplify subsequent filtering and
        priority calculations.

        Returns:
            astropy.table.Table: SPICA catalog containing science targets and associated survey metadata.
        """

        print("tapServerUrl:", tapServerUrl)

        # Query all spica-db catalog content
        adqlQuery = (
            "SELECT * FROM " + self.catalogSpica
        )
        service = vo.dal.TAPService(tapServerUrl)
        results = service.search(adqlQuery)

        # store new information into data variable
        data = results.to_table()

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
        """Query the SPICA primary calibrator catalog through the TAP service.

        Fetches the full content of the ``spica_calprim`` table.  An initial
        declination pre-selection can be applied within the query if needed.

        Returns:
            astropy.table.Table: Primary calibrator catalog as returned by the TAP service.
        """
        # query all calprim content
        calprim = "spica_calprim"
        calAdqlQuery = "SELECT * FROM " + calprim

        service = vo.dal.TAPService(tapServerUrl)
        results = service.search(calAdqlQuery)

        # store new information into data variable
        calibrators = results.to_table()

        return calibrators

    def observable_domain(self):
        """Compute the observable RA domain for CHARA on the selected night.

        Derives the local sidereal time at sunset and sunrise (both referred
        to the CHARA site) and converts them to right-ascension values in
        degrees, defining the RA window within which targets are accessible.

        Returns:
            tuple: Pair ``(alpha_sun_set, alpha_sun_rise)`` where both values
                are floats giving the RA (in degrees) of sunset and sunrise
                respectively.
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
        """
        Convert catalog-specific column names into a standardized
        SPICA/Aspro naming scheme.

        This method allows catalogs originating from different
        sources (SPICA database, JSDC, Vizier, calibrator catalogs)
        to be merged and exported using a consistent set of column
        names.

        Args:
            input (astropy.table.Table): Input catalog.
            colNames (dict, optional): Mapping between source column names and standardized output names.

        Returns:
            astropy.table.Table: Catalog with normalized column names and corrected column types.
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
        Export selected targets and calibrators to Aspro2.

        Science targets, primary calibrators and secondary
        calibrators are converted into a VOtable and transmitted
        through the SAMP protocol to a running Aspro2 instance.

        Args:
            targets (astropy.table.Table): Selected science targets.
            calibrators1 (astropy.table.Table, optional): Primary calibrators.
            calibrators2 (astropy.table.Table, optional): Secondary calibrators.

        Returns:
            None
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
                )

                # second table to declare some calibrators through Aspro's votable
                calPrims4sci = Table()
                calPrims4sci.add_column(
                    calPrims["target_main_id"], name="CALIBRATOR_NAME"
                )

            # Add secondary calibrators if any
            if calibrators2:
                calSeconds = self.normalizeColumnNames(calibrators2)
                calSeconds.add_column("calsecond", name=COLNAME_GRP)  # set main group

                # We should not have orphan calibrators since queryJsdc2 perform the association

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
            )
            Label(
                self.topSAMP,
                text=f"- {len(self.spica_catg[self.indexList_Targets])} science targets,\n",
                font=self.myFont,
            ).grid(column=0, row=1, sticky="W", padx=10)

            if self.indexList_CalPrim:
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
                )
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

    def cancelClickClients(self):
        """Cancel the pending Aspro2 SAMP export.

        Closes the SAMP confirmation popup and disconnects the SAMP client
        without sending any message.

        Returns:
            None
        """
        self.topSAMP.destroy()
        self.client.disconnect()

    def importClients(self):
        """
        Send the prepared VOtable to the selected Aspro2 SAMP client(s).

        If exactly one Aspro2 client is connected the VOtable SAMP message
        is sent directly to that client.  If several clients are open, the
        message is sent only to the subset chosen by the user via the
        checkbox selection in ``self.SelectedClient``.  After broadcasting,
        the SAMP confirmation popup is closed and the client disconnects from
        the SAMP hub.

        Returns:
            None
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

        Args:
            completion_rate (float): The initial completion rate value.
            spica_mode (str): The SPICA mode which determines how the completion rate is adjusted.

        Returns:
            float: The updated completion rate after adjustments based on the SPICA mode.
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
        Convert a completion rate into a survey completion flag.

        The returned flag represents the observing status of a target
        and depends on both the completion rate and the SPICA
        observing mode.

        Args:
            completion_rate (float): Fraction of observations completed for the target.
            spica_mode (str): SPICA observing mode (DIA, DLD, IMA, TMP or SPI).

        Returns:
            int: Completion flag where 1 = sufficiently completed,
                2 = partially completed, 3 = not completed or special mode.
        """

        flag = []
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
        Compute the final scheduling priority.

        The final priority combines the scientific priority assigned
        by the principal investigator with the current completion
        status of the target.

        Args:
            flag_completion (int): Completion category returned by ``update_flag_completion``.
            priority_pi (int): Priority assigned by the principal investigator.

        Returns:
            int: Final priority rank used by the SPICA-NSS target selection process.
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

    def update_priority_final2(self, flag_completion, priority_pi, progname2):
        """
        Compute the final scheduling priority for targets shared between programmes.

        An extended priority formula that accounts for membership in a second
        programme (``progname2``).  The returned value ranges from 1 (highest
        urgency) to 4 (lowest urgency).

        Args:
            flag_completion (int): Completion category as returned by :meth:`update_flag_completion`.
            priority_pi (int): Priority assigned by the principal investigator (0 = high, 1 = low).
            progname2 (str or bool): Secondary programme identifier; a truthy value indicates the target belongs to a second programme.

        Returns:
            int: Final priority rank (1–4).
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
        Attach a simple uniform-disk model to calibrator rows.

        Creates (or replaces) a ``model`` column and sets each row's model
        to an Aspro-compatible uniform-disk representation built from the
        limb-darkened diameter stored in the ``ld_jsdc2`` column.  Rows
        where ``ld_jsdc2`` is falsy are left with a ``None`` model.

        Args:
            data (astropy.table.Table): Calibrator table containing at least a ``ld_jsdc2`` column.

        Returns:
            astropy.table.Table: The input table with the ``model`` column added or updated.
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
        Attach Aspro-compatible geometric models to targets.

        Existing serialized models are converted into Aspro model
        objects. If no model is available, a simple uniform-disk
        model is generated using the target angular diameter.

        Args:
            data (astropy.table.Table): Table containing target properties and optional model definitions.

        Returns:
            astropy.table.Table: Input table augmented with a ``model`` column suitable for export through SAMP to Aspro2.
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
        Select primary calibrator candidates for the currently
        selected science targets.

        Candidate calibrators are extracted from the SPICA primary
        calibrator catalog and filtered according to configurable
        constraints in right ascension, declination and visual
        magnitude relative to the science target sample.

        Results are stored in ``self.calprim_catg`` and
        ``self.indexList_CalPrim``.

        Returns:
            None
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
        """
        Replace object names by their HD identifiers.

        For each object in the input table, SIMBAD is queried for
        alternative identifiers when the current name is not already
        an HD designation. If an HD identifier is found, the object
        name is replaced by the first matching HD entry returned by
        SIMBAD.

        Args:
            tableObjects (astropy.table.Table): Input catalog containing object identifiers. The table must include either a ``target_main_id`` column, a ``name`` column, or a ``Name`` column.

        Returns:
            astropy.table.Table: Catalog with object names replaced by HD identifiers whenever available.

        Note:
            Objects already identified by an HD name are not queried.
            If no HD identifier is found in SIMBAD, the original object
            name is preserved.

            The original column naming convention is restored before
            returning the table.
        """
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
        """Handle a change in the primary-calibrator RA search range entry.

        Updates ``self.rarangeprim`` and triggers :meth:`query_calprim`.

        Args:
            strRaRangePrim (tk.Entry): Entry widget holding the new RA range
                in arcminutes.

        Returns:
            None
        """
        self.rarangeprim = float(strRaRangePrim.get())
        self.query_calprim()

    def entryDecRangePrimCallback(self, strDecRangePrim):
        """Handle a change in the primary-calibrator Dec search range entry.

        Updates ``self.decrangeprim`` and triggers :meth:`query_calprim`.

        Args:
            strDecRangePrim (tk.Entry): Entry widget holding the new
                declination range in degrees.

        Returns:
            None
        """
        self.decrangeprim = float(strDecRangePrim.get())
        self.query_calprim()

    def entryVmagRangePrimCallback(self, strVmagRangePrim):
        """Handle a change in the primary-calibrator V-magnitude range entry.

        Updates ``self.vmagrangeprim`` and triggers :meth:`query_calprim`.

        Args:
            strVmagRangePrim (tk.Entry): Entry widget holding the new
                V-magnitude search range (±).

        Returns:
            None
        """
        self.vmagrangeprim = float(strVmagRangePrim.get())
        self.query_calprim()

    def goCalPrim(self):
        """Manually trigger the primary-calibrator query.

        Delegates directly to :meth:`query_calprim` using the current search
        parameter values.

        Returns:
            None
        """
        self.query_calprim()

    def delCalPrim(self):
        """Remove the current primary-calibrator selection and refresh the plot.

        Clears ``self.calprim_catg`` and ``self.indexList_CalPrim``, then
        calls :meth:`plot_radec` to remove calibrator markers from the sky
        distribution plot.

        Returns:
            None
        """
        print("[INFO] Delete primary calibrators")
        del self.indexList_CalPrim
        self.calprim_catg = None
        self.indexList_CalPrim = None
        self.plot_radec()

    def delCalSec(self):
        """Remove the current secondary-calibrator selection and refresh the plot.

        Clears ``self.calsec_catg`` and ``self.indexList_CalSec``, then
        calls :meth:`plot_radec` to remove secondary calibrator markers from
        the sky distribution plot.

        Returns:
            None
        """
        print("[INFO] Delete secondary calibrators")
        del self.indexList_CalSec
        self.calsec_catg = None
        self.indexList_CalSec = None
        self.plot_radec()

    def check4BeStars(self, tableObjects):
        """Remove Be stars from a calibrator table.

        Queries SIMBAD for the object type of each entry and removes any
        rows classified as ``"Be*"``, which are unsuitable calibrators due to
        their variable circumstellar emission.

        Args:
            tableObjects (astropy.table.Table): Calibrator table containing
                a ``Name`` column.

        Returns:
            astropy.table.Table: Input table with Be-star rows removed.
        """
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
        """Remove known bad calibrators from a calibrator table.

        Queries the ``badcal`` table from the JMMC TAP service and
        cross-matches it against the input table using a 5 arcsec tolerance.
        Any matching rows are removed from the calibrator list.

        Args:
            tableObjects (astropy.table.Table): Calibrator table containing
                ``ra`` and ``dec`` columns in degrees.

        Returns:
            astropy.table.Table: Input table with known bad calibrators
                removed.
        """
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
        """Select secondary calibrator candidates for the current target selection.

        Queries the JSDC2 catalog (Vizier ``II/346/jsdc_v2``) for stars
        within the configurable RA, Dec and V-magnitude windows centred on
        the median position of the selected science targets.  Candidates are
        further filtered by:

        - Observable RA window (between CHARA sunset and sunrise).
        - LDD chi² quality criterion (``self.lddchisec``).
        - Diameter relative error (``self.relerrorsec``).
        - Minimum squared visibility at ``self.maxbaseline`` (``self.minvissec``).

        Duplicates already present in the science-target list or among the
        primary calibrators are removed by sky cross-match.  If more than 50
        candidates survive, a warning popup is shown.  Otherwise the table is
        screened for Be stars and known bad calibrators before being stored
        in ``self.calsec_catg`` and ``self.indexList_CalSec``.

        Returns:
            None
        """
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

        print("[INFO] ---")
        print("")
        self.plot_radec()

    def entryRaRangeSecCallback(self, strRaRangeSec):
        """Handle a change in the secondary-calibrator RA search range entry.

        Updates ``self.rarangesec`` and triggers :meth:`query_calsec`.

        Args:
            strRaRangeSec (tk.Entry): Entry widget holding the new RA range
                in arcminutes.

        Returns:
            None
        """
        self.rarangesec = float(strRaRangeSec.get())
        self.query_calsec()

    def entryDecRangeSecCallback(self, strDecRangeSec):
        """Handle a change in the secondary-calibrator Dec search range entry.

        Updates ``self.decrangesec`` and triggers :meth:`query_calsec`.

        Args:
            strDecRangeSec (tk.Entry): Entry widget holding the new
                declination range in degrees.

        Returns:
            None
        """
        self.decrangesec = float(strDecRangeSec.get())
        self.query_calsec()

    def entryVmagRangeSecCallback(self, strVmagRangeSec):
        """Handle a change in the secondary-calibrator V-magnitude range entry.

        Updates ``self.vmagrangesec`` and triggers :meth:`query_calsec`.

        Args:
            strVmagRangeSec (tk.Entry): Entry widget holding the new
                V-magnitude search range (±).

        Returns:
            None
        """
        self.vmagrangesec = float(strVmagRangeSec.get())
        self.query_calsec()

    def entryLDDChiSecCallback(self, strLDDChiSec):
        """Handle a change in the secondary-calibrator LDD chi² threshold entry.

        Updates ``self.lddchisec`` and triggers :meth:`query_calsec`.

        Args:
            strLDDChiSec (tk.Entry): Entry widget holding the new maximum
                LDD fit chi² value.

        Returns:
            None
        """
        self.lddchisec = float(strLDDChiSec.get())
        self.query_calsec()

    def entryRelErrorSecCallback(self, strRelErrorSec):
        """Handle a change in the secondary-calibrator diameter relative-error entry.

        Updates ``self.relerrorsec`` and triggers :meth:`query_calsec`.

        Args:
            strRelErrorSec (tk.Entry): Entry widget holding the new maximum
                relative diameter error in percent.

        Returns:
            None
        """
        self.relerrorsec = float(strRelErrorSec.get())
        self.query_calsec()

    def entryMinVisSecCallback(self, strMinVisSec):
        """Handle a change in the secondary-calibrator minimum visibility entry.

        Updates ``self.minvissec`` and triggers :meth:`query_calsec`.

        Args:
            strMinVisSec (tk.Entry): Entry widget holding the new minimum
                squared visibility threshold.

        Returns:
            None
        """
        self.minvissec = float(strMinVisSec.get())
        self.query_calsec()

    def entryMinVisSecCallback2(self, strMinVisSec):
        """Alternative handler for the minimum-visibility entry that skips the value update.

        Triggers :meth:`query_calsec` without modifying ``self.minvissec``,
        used when the visibility threshold is updated through a separate
        control (e.g. the max-baseline spinbox) rather than the text entry.

        Args:
            strMinVisSec (tk.Entry): The minimum-visibility entry widget
                (value is not read).

        Returns:
            None
        """
        self.query_calsec()

    def goCalSec(self):
        """Manually trigger the secondary-calibrator query.

        Delegates directly to :meth:`query_calsec` using the current search
        parameter values.

        Returns:
            None
        """
        self.query_calsec()

    def getSelectedTargets(self):
        """Compute ``self.indexList_Targets`` from the active filter indices.

        Takes the intersection of all active filter index arrays
        (programme name, instrumental mode, final priority, Dec min/max,
        Vmag min/max) and stores the result in ``self.indexList_Targets``.
        Any manually added target index (``self.index_AddTarget``) is
        appended afterwards so it is always included regardless of the active
        filters.

        Returns:
            None
        """
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
        """Compute ``self.indexList_AddTarget`` for the fainter-star pool.

        Takes the intersection of the programme-name, instrumental-mode,
        final-priority, Dec min/max and Vmag-min filter indices together with
        ``self.indexVmagToAddTarget`` (targets one magnitude fainter than the
        current Vmag maximum) to build the list of candidates shown in the
        add-target popup.

        Returns:
            None
        """
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

    CalsecOpened = False
    LogOpened = False
    popupInfoTargets = False
    popupAddTarget = False
    popupBestDec = False
    topSAMP = False

    date = date.today()
    date = date.strftime("%Y-%m-%d")

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

tapServerUrl = "http://tap.jmmc.fr/vollt/tap/"
nssVersion = "v2025-09"

if __name__ == "__main__":
    # Start NSS gui
    app = spica_NSS()
