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
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterEnum,
                       QgsWkbTypes,
                       QgsFeature,
                       QgsGeometry,
                       QgsPoint,
                       QgsPointXY,
                       QgsEllipse)
from qgis import processing
import numpy as np
import numpy.polynomial.polynomial as poly


class MajorMinorAxes(QgsProcessingAlgorithm):
    """
    This script creates a line feature layer containing a representation
    of the major axes of each feature in a polygon layer by fitting a 
    polynomial best-fit line to the vertices of the polygon in cartesian space.
    """
    
    INPUT = 'INPUT'
    DIM = 'DIM'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return MajorMinorAxes()

    def name(self):
        return 'majorminoraxes'

    def displayName(self):
        return self.tr('Calculate major axes for Polygons')

    def group(self):
        return self.tr('Glacier Tools')

    def groupId(self):
        return 'glaciertools'

    def shortHelpString(self):
        return self.tr("Create major axes for polygons")

    def initAlgorithm(self, config=None):
        # Input Polygon Source
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Input Polygon Layer'),
                [QgsProcessing.TypeVectorPolygon]
            )
        )
        
        # Polynomial Dimension
        self.addParameter(
            QgsProcessingParameterEnum(self.DIM, self.tr("Polynomial Dimension"), ['1','2','3','4'], defaultValue='0')
        )
        
        # Output Line Layer
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Output Layer')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):

        source = self.parameterAsSource(
            parameters,
            self.INPUT,
            context
        )
        
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))
            
        dim = self.parameterAsInt(
            parameters,
            self.DIM,
            context
        )
    
        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            source.fields(),
            QgsWkbTypes.LineString,
            source.sourceCrs()
        )

        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))
        
        # Compute the number of steps to display within the progress bar and
        # get features from source
        total = 100.0 / source.featureCount() if source.featureCount() else 0
        features = source.getFeatures()

        for current, feature in enumerate(features):
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break

            geom = feature.geometry()
            
            vertices = list(geom.vertices())
            x = np.zeros(len(vertices))
            y = np.zeros(len(vertices))
            i = 0
            for vert in vertices:
                x[i]=vert.x()
                y[i]=vert.y()
                i = i + 1
            
            x_new = np.linspace(np.min(x), np.max(x), 100)
            
            # Polynomial fitting
            coefs = poly.polyfit(x, y, dim + 1)
            y_new = poly.polyval(x_new, coefs)
            
            outLine = QgsFeature()
            
            outVertices = []
            for idx in range(len(x_new)):
                outVertices.append(QgsPoint(x_new[idx], y_new[idx]))
            
            outLine.setGeometry(QgsGeometry.fromPolyline(outVertices))
            outLine.setAttributes(feature.attributes())
            
            centroid = geom.centroid()
            centerPt = QgsPoint(centroid.asPoint().x(), centroid.asPoint().y())
                
            sink.addFeature(outLine, QgsFeatureSink.FastInsert)
            
            # Update the progress bar
            feedback.setProgress(int(current * total))

        return {self.OUTPUT: dest_id}
