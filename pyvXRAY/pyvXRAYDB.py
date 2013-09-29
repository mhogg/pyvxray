# -*- coding: utf-8 -*-

# Copyright (C) 2013 Michael Hogg

# This file is part of pyvXRAY - See LICENSE.txt for information on usage and redistribution

from abaqusConstants import *
from abaqusGui import *
import os

thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)

class PyvXRAYDB(AFXDataDialog):

    def __init__(self, form):

        # Construct the base class.
        AFXDataDialog.__init__(self, form, 'Create virtual x-rays',self.OK|self.CANCEL, DIALOG_ACTIONS_SEPARATOR)
        self.form = form
        
        okBtn = self.getActionButton(self.ID_CLICKED_OK)
        okBtn.setText('OK')
        
        # Define Tab book
        TabBook = FXTabBook(p=self, tgt=None, sel=0, opts=TABBOOK_NORMAL,x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,pt=DEFAULT_SPACING, pb=DEFAULT_SPACING)

        # Define Regions Tab            
        tabItem   = FXTabItem(p=TabBook, text='Select regions', ic=None, opts=TAB_TOP_NORMAL,x=0, y=0, w=0, h=0, pl=6, pr=6, pt=DEFAULT_PAD, pb=DEFAULT_PAD)
        TabItem_1 = FXVerticalFrame(p=TabBook, opts=FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X,
                                    x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
                                    pt=DEFAULT_SPACING, pb=DEFAULT_SPACING, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        # Define Regions - Bone Region
        GroupBox_1 = FXGroupBox(p=TabItem_1, text='Bone region', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_1 = AFXVerticalAligner(p=GroupBox_1, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=0, pb=0)
        AFXTextField(p=VAligner_1, ncols=0, labelText='Part instance name:      ', tgt=form.bPartNameKw, sel=0, w=300, opts=LAYOUT_FIX_WIDTH)
        AFXTextField(p=VAligner_1, ncols=0, labelText='Set name:', tgt=form.bSetNameKw, sel=0, w=300, opts=LAYOUT_FIX_WIDTH)
        AFXTextField(p=VAligner_1, ncols=0, labelText='Density variable:', tgt=form.BMDfonameKw, sel=0, w=300, opts=LAYOUT_FIX_WIDTH)
        # Define Regions - Implant Region
        GroupBox_2 = FXGroupBox(p=TabItem_1, text='Implant region', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        self.cb1   = FXCheckButton(p=GroupBox_2, text='Show implant on x-rays', tgt=form.showImplantKw, sel=0)
        VAligner_2 = AFXVerticalAligner(p=GroupBox_2, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=0, pb=0)
        self.tf1   = AFXTextField(p=VAligner_2, ncols=0, labelText='Part instance name:       ', tgt=form.iPartNameKw, sel=0, w=300, opts=LAYOUT_FIX_WIDTH)
        self.tf2   = AFXTextField(p=VAligner_2, ncols=0, labelText='Set name:', tgt=form.iSetNameKw, sel=0, w=300, opts=LAYOUT_FIX_WIDTH)
        self.tf3   = AFXTextField(p=VAligner_2, ncols=0, labelText='Density (kg/m^3):', tgt=form.iDensityKw, sel=0, w=300, opts=LAYOUT_FIX_WIDTH)

        # Inputs Tab
        tabItem   = FXTabItem(p=TabBook, text='Inputs', ic=None, opts=TAB_TOP_NORMAL, x=0, y=0, w=0, h=0, pl=6, pr=6, pt=DEFAULT_PAD, pb=DEFAULT_PAD)
        TabItem_2 = FXVerticalFrame(p=TabBook, opts=FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X,
                                    x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
                                    pt=DEFAULT_SPACING, pb=DEFAULT_SPACING, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        GroupBox_3 = FXGroupBox(p=TabItem_2, text='Required inputs', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_3 = AFXVerticalAligner(p=GroupBox_3, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=10, pb=10)
        AFXTextField(p=VAligner_3, ncols=21, labelText='Step list:', tgt=form.stepListKw, sel=0, pt=5, pb=5)
        
        ComboBox_1 = AFXComboBox(p=VAligner_3, ncols=19, nvis=1, text='Coordinate system:', tgt=form.csysNameKw, sel=0, pt=5, pb=5)
        csyses = []
        for csysType,csysNames in form.csysList.items():
            for csysName in csysNames:
                listText = '%s (%s)' % (csysName,csysType)
                csyses.append(listText)
        csyses.sort()
        csyses.insert(0,'GLOBAL') # Add global to the start of the sorted list
        self.form.csysNameKw.setValue(csyses[0])
        for csys in csyses:
            ComboBox_1.appendItem(text=csys)
        ComboBox_1.setMaxVisible(5)
                
        AFXTextField(p=VAligner_3, ncols=21, labelText='Mapping resolution (mm):      ', tgt=form.resGridKw, sel=0, pt=5, pb=5)

        # X-ray settings Tab
        tabItem   = FXTabItem(p=TabBook, text='X-ray settings', ic=None, opts=TAB_TOP_NORMAL, x=0, y=0, w=0, h=0, pl=6, pr=6, pt=DEFAULT_PAD, pb=DEFAULT_PAD)
        TabItem_3 = FXVerticalFrame(p=TabBook, opts=FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X,
                                    x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
                                    pt=DEFAULT_SPACING, pb=DEFAULT_SPACING, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        GroupBox_4 = FXGroupBox(p=TabItem_3, text='Settings', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_4 = AFXVerticalAligner(p=GroupBox_4, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=10, pb=0)
        AFXTextField(p=VAligner_4, ncols=13, labelText='Base name of xray file(s):', tgt=form.imageNameBaseKw, sel=0, pt=5, pb=5)
        AFXTextField(p=VAligner_4, ncols=13, labelText='Approx size of x-ray images (in pixels):  ', tgt=form.preferredXraySizeKw, sel=0, pt=5, pb=5)
        ComboBox_2 = AFXComboBox(p=VAligner_4, ncols=11, nvis=1, text='Image file format', tgt=form.imageFormatKw, sel=0, pt=5, pb=5)
        ComboBox_2.setMaxVisible(5)
        for imageFormat in form.imageFormats:
            ComboBox_2.appendItem(text=imageFormat)
        FXCheckButton(p=GroupBox_4, text='Smooth images', tgt=form.smoothKw, sel=0, pt=10, pb=5)
        FXCheckButton(p=GroupBox_4, text='Manual scaling of images', tgt=form.manualScalingKw, sel=0, pt=10, pb=5)
        
    def processUpdates(self):
        # Update form
        self.form.showImplant = self.cb1.getCheck()
        # Disable implant option if show implant not checked
        tfs = [self.tf1,self.tf2,self.tf3]
        if self.cb1.getCheck():
            for tf in tfs: tf.enable() 
        else: 
            for tf in tfs: tf.disable() 
        return          
