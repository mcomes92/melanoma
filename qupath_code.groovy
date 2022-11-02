setImageType('BRIGHTFIELD_H_E');
setColorDeconvolutionStains('{"Name" : "H&E default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049", "Stain 2" : "Eosin", "Values 2" : "0.2159 0.8012 0.5581", "Background" : " 255 255 255"}');

def imageData = getCurrentImageData()
def annotations = getAnnotationObjects()

//Keep only the annotation which has the largest area
def max_area=annotations[0].getROI().getArea()
def max_indice=0
for (i=0;i<annotations.size();i++){
    if (annotations[i].getROI().getArea()>=max_area){
        max_area=annotations[i].getROI().getArea()
        max_indice=i
    }
} 
for (i in annotations){
    if (i.getROI().getArea()<max_area){
        removeObject(i,true)
    }
}
// this is the area which will be divided in tiles

/***** Tile creation *****/ 
selectAnnotations()
// setting tiles parameters
runPlugin('qupath.lib.algorithms.TilerPlugin', '{"tileSizeMicrons": 110.0,  "trimToROI": true,  "makeAnnotations": true,  "removeParentAnnotation": false}');

def annotations2 = getAnnotationObjects()

//Keep Tiles and exclude first annotations
def max_area2=annotations2[0].getROI().getArea()
def max_indice2=0
for (i=0;i<annotations2.size();i++){
    if (annotations2[i].getROI().getArea()>=max_area){
        max_area2=annotations2[i].getROI().getArea()
        max_indice2=i
    }
} 
for (i in annotations2){
    if (i.getROI().getArea()==max_area2){
        removeObject(i,true)
    }
}
// saving tiles
def server = getCurrentServer()
// Define output path (here, relative to project)
def name = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())

def pathOutput = "D:\\oncologico\\digital pathology\\progetto melanoma_no_cell_class\\Ti2les\\"+name+"\\"
//Define what you want to export here
tiles = getAnnotationObjects()//getDetectionObjects().findAll{it.getPathClass() != null}


mkdirs(pathOutput)
i=1
for (tile in tiles){
    def requestFull = RegionRequest.createInstance(server.getPath(),1,tile.getROI())
    x = tile.getROI().getCentroidX()
    y = tile.getROI().getCentroidY()
//    tileClass = tile.getPathClass().toString()
//    fileName = pathOutput+"Tile "+"x "+x+" y "+y+" "+tileClass+" "+(i++)+".tif"
    fileName = pathOutput+"Tile "+"x "+x+" y "+y+" "+(i++)+".tif"
    //print describe(requestFull)
    writeImageRegion(server, requestFull, fileName)
    
}


// cell detection on tiles
selectAnnotations()
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImageBrightfield": "Hematoxylin OD",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 0.1,  "maxBackground": 2.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 5.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

// save annotations
// Define output path (here, relative to project)
def name2 = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())+ '.csv'

def path2 = buildFilePath(PROJECT_BASE_DIR, 'annotation results')
mkdirs(path2)
path = buildFilePath(path2, name2)
saveAnnotationMeasurements(path2)
print 'Results exported to ' + path2


