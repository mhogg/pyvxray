# -*- coding: utf-8 -*-

# Copyright (C) 2013 Michael Hogg

# This file is part of pyvXRAY - See LICENSE.txt for information on usage and redistribution

from abaqusConstants import *
from abaqusGui import *
import os

thisPath = os.path.abspath(__file__)
thisDir  = os.path.dirname(thisPath)

class PyvXRAYDB(AFXDataDialog):

    def __init__(self, form):

        # Construct the base class.
        AFXDataDialog.__init__(self, form, 'pyvXRAY - Create virtual x-rays',self.OK|self.CANCEL, DIALOG_ACTIONS_SEPARATOR)
        self.form = form
        
        okBtn = self.getActionButton(self.ID_CLICKED_OK)
        okBtn.setText('OK')
        
        # Define Tab book
        TabBook = FXTabBook(p=self, tgt=None, sel=0, opts=TABBOOK_NORMAL,x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,pt=DEFAULT_SPACING, pb=DEFAULT_SPACING)

        # Define Regions Tab            
        FXTabItem(p=TabBook, text='Select regions', ic=None, opts=TAB_TOP_NORMAL,x=0, y=0, w=0, h=0, pl=6, pr=6, pt=DEFAULT_PAD, pb=DEFAULT_PAD)
        TabItem_1 = FXVerticalFrame(p=TabBook, opts=FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X,
                                    x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
                                    pt=DEFAULT_SPACING, pb=DEFAULT_SPACING, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        
        # Select odb
        GroupBox_1 = FXGroupBox(p=TabItem_1, text='Result file', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_1 = AFXVerticalAligner(p=GroupBox_1, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=0, pb=0)  
        
        ComboBox_1 = AFXComboBox(p=VAligner_1, ncols=35, nvis=1, text='%-27s'%'Odb:', tgt=form.odbNameKw, sel=0, pt=5, pb=5)    
        if len(form.odbList)>0:
            for odbName in form.odbList:
                ComboBox_1.appendItem(odbName)
            form.odbNameKw.setValue(form.odbList[0])
            ComboBox_1.setMaxVisible(10)
            self.odbName = form.odbList[0]
        else: self.odbName=None
        
        # Select bone region 
        GroupBox_2 = FXGroupBox(p=TabItem_1, text='Bone region', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_2 = AFXVerticalAligner(p=GroupBox_2, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=0, pb=0)
        
        self.ComboBox_2 = AFXComboBox(p=VAligner_2, ncols=35, nvis=1, text='Bone set:', tgt=form.bSetNameKw, sel=0, pt=5, pb=5)    
        self.ComboBox_2.setMaxVisible(10)      
        self.populateElementListComboBox() 
        
        self.ComboBox_3 = AFXComboBox(p=VAligner_2, ncols=35, nvis=1, text='%-20s'%'Density variable:', tgt=form.BMDfonameKw, sel=0, pt=5, pb=5)    
        self.ComboBox_3.setMaxVisible(10)      
        self.populateScalarListComboBox() 
        
        # Select implant region
        GroupBox_3 = FXGroupBox(p=TabItem_1, text='Implant region', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)        

        self.cb1   = FXCheckButton(p=GroupBox_3, text='Show implant on x-rays', tgt=form.showImplantKw, sel=0)        
        
        VAligner_3 = AFXVerticalAligner(p=GroupBox_3, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=0, pb=0)   
        self.ComboBox_4 = AFXComboBox(p=VAligner_3, ncols=35, nvis=1, text='Implant set:', tgt=form.iSetNameKw, sel=0, pt=5, pb=5)    
        self.ComboBox_4.setMaxVisible(10)      
        self.populateElementListComboBoxImplant()

        self.tf1 = AFXTextField(p=VAligner_3, ncols=0, labelText='Density (kg/m^3):', tgt=form.iDensityKw, sel=0, w=331, opts=LAYOUT_FIX_WIDTH)
        
        # Inputs Tab
        FXTabItem(p=TabBook, text='Inputs', ic=None, opts=TAB_TOP_NORMAL, x=0, y=0, w=0, h=0, pl=6, pr=6, pt=DEFAULT_PAD, pb=DEFAULT_PAD)
        TabItem_2 = FXVerticalFrame(p=TabBook, opts=FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X,
                                    x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
                                    pt=DEFAULT_SPACING, pb=DEFAULT_SPACING, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        GroupBox_4 = FXGroupBox(p=TabItem_2, text='Required inputs', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_4 = AFXVerticalAligner(p=GroupBox_4, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=10, pb=10)
        
        AFXTextField(p=VAligner_4, ncols=26, labelText='Step list:', tgt=form.stepListKw, sel=0, pt=5, pb=5)
        form.stepListKw.setValue(', '.join(form.stepList))
        
        ComboBox_5 = AFXComboBox(p=VAligner_4, ncols=24, nvis=1, text='Coordinate system:', tgt=form.csysNameKw, sel=0, pt=5, pb=5)
        csyses = []
        for csysType,csysNames in form.csysList.items():
            for csysName in csysNames:
                listText = '%s (%s)' % (csysName,csysType)
                csyses.append(listText)
        csyses.sort()
        csyses.insert(0,'GLOBAL') # Add global to the start of the sorted list
        form.csysNameKw.setValue(csyses[0])
        for csys in csyses:
            ComboBox_5.appendItem(text=csys)
        ComboBox_5.setMaxVisible(5)
                
        AFXTextField(p=VAligner_4, ncols=26, labelText='%-30s'%'Mapping resolution (mm):', tgt=form.resGridKw, sel=0, pt=5, pb=5)

        # X-ray settings Tab
        FXTabItem(p=TabBook, text='X-ray settings', ic=None, opts=TAB_TOP_NORMAL, x=0, y=0, w=0, h=0, pl=6, pr=6, pt=DEFAULT_PAD, pb=DEFAULT_PAD)
        TabItem_3 = FXVerticalFrame(p=TabBook, opts=FRAME_RAISED|FRAME_THICK|LAYOUT_FILL_X,
                                    x=0, y=0, w=0, h=0, pl=DEFAULT_SPACING, pr=DEFAULT_SPACING,
                                    pt=DEFAULT_SPACING, pb=DEFAULT_SPACING, hs=DEFAULT_SPACING, vs=DEFAULT_SPACING)
        GroupBox_5 = FXGroupBox(p=TabItem_3, text='Settings', opts=FRAME_GROOVE|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        VAligner_5 = AFXVerticalAligner(p=GroupBox_5, opts=0, x=0, y=0, w=0, h=0, pl=0, pr=0, pt=10, pb=0)
        
        AFXTextField(p=VAligner_5, ncols=18, labelText='Base name of xray file(s):', tgt=form.imageNameBaseKw, sel=0, pt=5, pb=5)

        AFXTextField(p=VAligner_5, ncols=18, labelText='%-42s'%'Approx size of x-ray images (in pixels):', tgt=form.preferredXraySizeKw, sel=0, pt=5, pb=5)

        ComboBox_6 = AFXComboBox(p=VAligner_5, ncols=16, nvis=1, text='Image file format', tgt=form.imageFormatKw, sel=0, pt=5, pb=5)
        for imageFormat in form.imageFormats:
            ComboBox_6.appendItem(text=imageFormat)
        ComboBox_6.setMaxVisible(5)
        
        FXCheckButton(p=GroupBox_5, text='Smooth images', tgt=form.smoothKw, sel=0, pt=10, pb=5)

        FXCheckButton(p=GroupBox_5, text='Manual scaling of images', tgt=form.manualScalingKw, sel=0, pt=10, pb=5)
        
    def populateElementListComboBox(self):
        """Populate comboBox containing element sets for bone"""
        if len(self.form.elementSets)==0: return
        self.ComboBox_2.clearItems()
        for elementSet in self.form.elementSets:
            self.ComboBox_2.appendItem(elementSet)
        self.form.bSetNameKw.setValue(self.form.elementSets[0]) 
        
    def populateScalarListComboBox(self):
        """Populate comboBox containing scalar fieldoutputs"""
        if len(self.form.scalarList)==0: return
        self.ComboBox_3.clearItems()
        for scalar in self.form.scalarList:
            self.ComboBox_3.appendItem(scalar)
        self.form.BMDfonameKw.setValue(self.form.scalarList[0])   
        
    def populateElementListComboBoxImplant(self):
        """Populate comboBox containing element sets for implant"""
        if len(self.form.elementSets)==0: return
        self.ComboBox_4.clearItems()
        for elementSet in self.form.elementSets:
            self.ComboBox_4.appendItem(elementSet)
        self.form.iSetNameKw.setValue(self.form.elementSets[0]) 

    def processUpdates(self):
        """Update form"""
        # If odb name changes, the re-populate the region list
        if self.form.odbNameKw.getValue() != self.odbName:
            self.odbName = self.form.odbNameKw.getValue()
            self.form.setOdb(self.odbName)
            self.form.getElementSetList()
            self.populateElementListComboBox() 
            self.populateElementListComboBoxImplant()
            self.form.getScalarList()
            self.populateScalarListComboBox()
        # Disable implant option if show implant not checked
        tfs = [self.tf1,self.ComboBox_4]
        if self.cb1.getCheck():
            for tf in tfs: tf.enable() 
        else: 
            for tf in tfs: tf.disable() 
        return          
