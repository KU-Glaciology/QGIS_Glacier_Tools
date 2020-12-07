# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProcessingFeedback,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingOutputRasterLayer,
                       QgsProcessingOutputVectorLayer,
                       QgsProcessingParameterRasterDestination,
                       QgsProcessingParameterVectorDestination,
                       QgsRasterLayer,
                       QgsVectorLayer)
from qgis import processing
from scipy import ndimage
import numpy as np
import rasterio
import rasterio.mask
import time
import gdal
import json


class CrevassePolygons(QgsProcessingAlgorithm):
    """
    This processing tool takes a Digital Elevation Model (DEM) and identifies
    crevasses, or holes of any kind really, using a neighborhood mean approach
    with user-defined parameters.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT = 'INPUT'
    OUTRASTER = 'OUTRASTER'
    OUTVECTOR = 'OUTVECTOR'
    MASK = 'MASK'
    OUTPUT = 'OUTPUT'
    NEIGHBORHOOD = 'NEIGHBORHOOD'
    STAT = 'STAT'
    DEVIATION = 'DEVIATION'
    CIRCULAR = 'CIRCULAR'

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return CrevassePolygons()

    def name(self):
        return 'crevassepolygons'

    def displayName(self):
        return self.tr('Crevasse Polygons')

    def group(self):
        return self.tr('Glacier Tools')

    def groupId(self):
        return 'glaciertools'

    def shortHelpString(self):
        return self.tr("Create polygons that indicate low spots in a DEM")

    def initAlgorithm(self, config=None):
        # Input Digital Elevation Model to analyze
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT,
                self.tr('Input DEM'),
                optional=False
            )
        )
        
        # Optional Polygon Mask
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.MASK,
                self.tr('Polygon Mask'),
                [QgsProcessing.TypeVectorPolygon],
                optional=True
            )
        )
        
        # Size of the neighborhood
        self.addParameter(
            QgsProcessingParameterNumber(
                self.NEIGHBORHOOD,
               self.tr('Neighborhood radius, the resulting neighborhood will be 2r + 1 by 2r + 1'),
               type=QgsProcessingParameterNumber.Integer,
               minValue=0, 
               defaultValue=3, 
               optional=False
            )
        )
        
        # Use Circular Neighborhood?
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CIRCULAR, 
                self.tr('Use Circular Neighborhood'), 
                defaultValue=False
            )
        )
           
        # Aggregate Type, mean or max
        self.addParameter(
            QgsProcessingParameterEnum(
                self.STAT, 
                self.tr('Aggregate Type'),
                ['Mean','Max'],
                defaultValue='Mean',
                optional=False
            )
        )
        
        # Deviation Threshold
        self.addParameter(
            QgsProcessingParameterNumber(
                self.DEVIATION,
                self.tr('Deviation Threshold'),
                type=QgsProcessingParameterNumber.Double,
                minValue=0.0, 
                defaultValue=5.0, 
                optional=False
            )
        )
        
        # Output Raster Location
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTRASTER, 
                self.tr('Output binary crevasse raster layer')
            )
        )
        
        # Output Vector Location
        self.addParameter(
            QgsProcessingParameterVectorDestination(
                self.OUTVECTOR, 
                self.tr('Output crevasse polygon layer')
            )
        )
        
        
    def processAlgorithm(self, parameters, context, feedback):
        
        # Grab our parameters from the user input
        maskSource = self.parameterAsSource(parameters, self.MASK, context)
        circular = self.parameterAsBool(parameters, self.CIRCULAR, context)
        neighborhood = self.parameterAsInt(parameters, self.NEIGHBORHOOD, context)
        stat = self.parameterAsInt(parameters, self.STAT, context)
        deviation = self.parameterAsDouble(parameters, self.DEVIATION, context)
        outraster = self.parameterAsOutputLayer(parameters, self.OUTRASTER, context)
        outvector = self.parameterAsOutputLayer(parameters, self.OUTVECTOR, context)
        dem = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        demProvider = dem.dataProvider()
        demUri = demProvider.dataSourceUri()
        
        # Declare functions to be used below
        def _create_footprint(n, circle):
            
            # Don't create a even number of elements
            if(n%2 == 0):
                n = n + 1
            
            # Create our True mask
            mask = np.ones((n,n))
            
            # If we want a circle, we need to turn some of the mask cells offTheRecordChanged
            if(circle):
                x = np.arange(0, n)
                y = np.arange(0, n)
                c = (n - 1) / 2
                r = neighborhood
                zeros = (x[np.newaxis,:]-c)**2 + (y[:,np.newaxis]-c)**2 > r**2
                mask[zeros] = 0
            
            # Return our footprint to use
            return mask
        
        # Get our processing footprint
        feedback.pushConsoleInfo(str("Getting Footprint"))
        footprint = _create_footprint((2 * neighborhood + 1), circular)
        feedback.pushConsoleInfo(str(footprint))
        
        # Process our DEM
        feedback.pushConsoleInfo(str("Opening DEM"))
        with rasterio.open(demUri) as dataset:
            
            demArray = dataset.read(1)
            
            # Mask out the DEM if a mask layer has been provided
            if(maskSource):
                features = maskSource.getFeatures()
                for current, feature in enumerate(features):
                    geom = feature.geometry()
                    geojson = json.loads(geom.asJson())
                    demArray, out_transform = rasterio.mask.mask(dataset, [geojson], crop=False, indexes=1, nodata=-9999)
            
            feedback.pushConsoleInfo(str(dataset.shape))
            w, h = dataset.shape
            cnt = w * h
            
            i = 0
            def _filter_func(values):
                if(feedback.isCanceled()):
                    return exit()
                nonlocal i
                i = i + 1
                p = (i / cnt) * 100
                feedback.setProgress(p)

                idx = int((len(values) - 1) / 2)
                val = values[idx]
                
                nonlocal stat
                nonlocal deviation
                if(stat == 1):
                    calc = np.max(values)
                else:
                    calc = np.mean(values)
                    
                out = 0
                if((calc - val) >= deviation):
                    out = 1
                return out
                
            feedback.pushConsoleInfo(str("Starting Processing"))
            
            feedback.pushConsoleInfo('outraster=' + str(outraster))

            filtered = ndimage.generic_filter(demArray,
                       function=_filter_func,
                       footprint=footprint,
                       mode='nearest')
                       
            feedback.pushConsoleInfo(str("Completed Filter Processing"))
            feedback.pushConsoleInfo(str("Writing Output Raster Dataset"))
            
            with rasterio.open(
                str(outraster),
                'w',
                driver='GTiff',
                height=filtered.shape[0],
                width=filtered.shape[1],
                count=1,
                dtype=filtered.dtype,
                crs=dataset.crs,
                transform=dataset.transform,
            ) as outrasterFile:
                outrasterFile.write(filtered, 1)
            
            feedback.pushConsoleInfo(str("Converting to polygons"))
            
            # Check for cancelation
            if feedback.isCanceled():
                return {}
                
            polygons = processing.run(
                'grass7:r.to.vect', {
                    '-b': False,
                    '-s': True,
                    '-t': False,
                    '-v': False,
                    '-z': False,
                    'GRASS_OUTPUT_TYPE_PARAMETER': 0,
                    'GRASS_REGION_CELLSIZE_PARAMETER': 0,
                    'GRASS_REGION_PARAMETER': None,
                    'GRASS_VECTOR_DSCO': '',
                    'GRASS_VECTOR_EXPORT_NOCAT': False,
                    'GRASS_VECTOR_LCO': '',
                    'column': 'value',
                    'input': QgsRasterLayer(outraster),
                    'type': 2,
                    'output': QgsProcessing.TEMPORARY_OUTPUT
                }, 
                context=context, 
                feedback=feedback, 
                is_child_algorithm=True
            )
            
            feedback.pushConsoleInfo(str("Removing non-crevasse polygons"))
            
            # Check for cancelation
            if feedback.isCanceled():
                return {}
            
            outExtracted = processing.run(
                "native:extractbyattribute", {
                    'INPUT': polygons['output'],
                    'FIELD': "value",
                    'OPERATOR': 0,
                    'VALUE': '1',
                    'OUTPUT': outvector
                },
                context=context, 
                feedback=feedback, 
                is_child_algorithm=True
            )
            
            feedback.pushConsoleInfo(str("Complete"))
            
            return {self.OUTPUT: outvector }
