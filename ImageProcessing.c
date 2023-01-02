<<<<<<< HEAD
 
//**************************************************************************
//* WARNING: This file was automatically generated.  Any changes you make  *
//*          to this file will be lost if you generate the file again.     *
//**************************************************************************
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <nivision.h>
#include <nimachinevision.h>
#include <windows.h>

// If you call Machine Vision functions in your script, add NIMachineVision.c to the project.

#define IVA_MAX_BUFFERS 10
#define IVA_STORE_RESULT_NAMES

#define VisionErrChk(Function) {if (!(Function)) {success = 0; goto Error;}}

typedef enum IVA_ResultType_Enum {IVA_NUMERIC, IVA_BOOLEAN, IVA_STRING} IVA_ResultType;

typedef union IVA_ResultValue_Struct    // A result in Vision Assistant can be of type double, BOOL or string.
{
    double numVal;
    BOOL   boolVal;
    char*  strVal;
} IVA_ResultValue;

typedef struct IVA_Result_Struct
{
#if defined (IVA_STORE_RESULT_NAMES)
    char resultName[256];           // Result name
#endif
    IVA_ResultType  type;           // Result type
    IVA_ResultValue resultVal;      // Result value
} IVA_Result;

typedef struct IVA_StepResultsStruct
{
#if defined (IVA_STORE_RESULT_NAMES)
    char stepName[256];             // Step name
#endif
    int         numResults;         // number of results created by the step
    IVA_Result* results;            // array of results
} IVA_StepResults;

typedef struct IVA_Data_Struct
{
    Image* buffers[IVA_MAX_BUFFERS];            // Vision Assistant Image Buffers
    IVA_StepResults* stepResults;              // Array of step results
    int numSteps;                               // Number of steps allocated in the stepResults array
    CoordinateSystem *baseCoordinateSystems;    // Base Coordinate Systems
    CoordinateSystem *MeasurementSystems;       // Measurement Coordinate Systems
    int numCoordSys;                            // Number of coordinate systems
} IVA_Data;



static IVA_Data* IVA_InitData(int numSteps, int numCoordSys);
static int IVA_DisposeData(IVA_Data* ivaData);
static int IVA_Classification_Segmentation(Image* image, Image* imageMask, const ROI* roi, ParticleClassifierPreprocessingOptions2 preprocessingOptions);
static int IVA_Classification_Extract_Particles(Image* image, Image* imageMask, ROI** rois, int numParticles);
static int IVA_DisposeStepResults(IVA_Data* ivaData, int stepIndex);
static int IVA_CLRExtractValue(Image* image);
static int IVA_EdgeTool2(Image* image,
                                  ROI* roi,
                                  int pProcessType,
                                  int pPolarity,
                                  int pKernelSize,
                                  int pWidth,
                                  float pMinThreshold,
                                  int pInterpolationType,
                                  int pColumnProcessingMode,
                                  IVA_Data* ivaData,
                                  int stepIndex);
static int IVA_Classification_Multiple(Image* image, char* fileName, ROI* roi);

int IVA_ProcessImage(Image *image)
{
	int success = 1;
    IVA_Data *ivaData;
    BCGOptions redBCGOptions;
    BCGOptions greenBCGOptions;
    BCGOptions blueBCGOptions;
    float kernel[9] = {-1,-1,-1,-1,10,-1,-1,-1,-1};
    ROI *roi;
    ROI *roi1;

    // Initializes internal data (buffers and array of points for caliper measurements)
    VisionErrChk(ivaData = IVA_InitData(6, 0));

    //-------------------------------------------------------------------//
    //                          Color Brightness                         //
    //-------------------------------------------------------------------//

    redBCGOptions.brightness = 164;
    redBCGOptions.contrast = 65.8000009655952;
    redBCGOptions.gamma = 1.86;

    greenBCGOptions.brightness = 164;
    greenBCGOptions.contrast = 65.8000009655952;
    greenBCGOptions.gamma = 1.86;

    blueBCGOptions.brightness = 164;
    blueBCGOptions.contrast = 65.8000009655952;
    blueBCGOptions.gamma = 1.86;

    // Applies brightness, contrast, and gamma correction to the image
    VisionErrChk(imaqColorBCGTransform(image, image, &redBCGOptions, &greenBCGOptions, &blueBCGOptions, NULL));

	VisionErrChk(IVA_CLRExtractValue(image));

    //-------------------------------------------------------------------//
    //                       Lookup Table: Inverse                       //
    //-------------------------------------------------------------------//

    // Inverts the pixel intensities of the image.
    VisionErrChk(imaqInverse(image, image, NULL));

    //-------------------------------------------------------------------//
    //                        Filters: Convolution                       //
    //-------------------------------------------------------------------//

    // Applies a linear filter to an image by convolving the image with a filtering kernel.
    VisionErrChk(imaqConvolve(image, image, kernel, 3, 3, 0, NULL));

    // Creates a new, empty region of interest.
    VisionErrChk(roi = imaqCreateROI());

    // Creates a new line ROI contour and adds the line to the provided ROI.
    VisionErrChk(imaqAddLineContour(roi, imaqMakePoint(34, 418),
                                         imaqMakePoint(810, 446)));

	VisionErrChk(IVA_EdgeTool2(image, roi, IMAQ_BEST, 
		IMAQ_SEARCH_FOR_ALL_EDGES, 5, 5, 10, IMAQ_BILINEAR_FIXED, 
		IMAQ_AVERAGE_COLUMNS, ivaData, 4));

    // Cleans up resources associated with the object
    imaqDispose(roi);

    // Creates a new, empty region of interest.
    VisionErrChk(roi1 = imaqCreateROI());

    // Creates a new rectangle ROI contour and adds the rectangle to the provided ROI.
    VisionErrChk(imaqAddRectContour(roi1, imaqMakeRect(17, 25, 929, 674)));

	VisionErrChk(IVA_Classification_Multiple(image, 
		"C:\\Users\\HP1\\Downloads\\leaf\\Leaf Classifier.clf", roi1));

    // Cleans up resources associated with the object
    imaqDispose(roi1);

    // Releases the memory allocated in the IVA_Data structure.
    IVA_DisposeData(ivaData);

Error:
	return success;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_CLRExtractValue
//
// Description  : Extracts the value plane from a color image.
//
// Parameters   : image  -  Input image
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_CLRExtractValue(Image* image)
{
    int success = 1;
    Image* plane;


    //-------------------------------------------------------------------//
    //                         Extract Color Plane                       //
    //-------------------------------------------------------------------//

    // Creates an 8 bit image that contains the extracted plane.
    VisionErrChk(plane = imaqCreateImage(IMAQ_IMAGE_U8, 7));

    // Extracts the value plane
    VisionErrChk(imaqExtractColorPlanes(image, IMAQ_HSV, NULL, NULL, plane));

    // Copies the color plane in the main image.
    VisionErrChk(imaqDuplicate(image, plane));

Error:
    imaqDispose(plane);

    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_EdgeTool2
//
// Description  : Finds edges along a path in an image.
//
// Parameters   : image                  -  Input image
//                roi                    -  Region of Interest
//                pProcessType           -  Determines which edges the function looks for. 
//                pPolarity              -  Specifies the polarity of the edges to be found
//                pKernelSize            -  Specifies the size of the edge detection kernel
//                pWidth                 -  Specifies the number of pixels averaged perpendicular
//                                          to the search direction to compute the edge profile
//                                          strength at each point along the search ROI
//                pMinThreshold          -  Specifies the minimum edge strength (gradient magnitude)
//                                          required for a detected edge
//                pInterpolationType     -  Specifies the interpolation method used to locate the edge position
//                pColumnProcessingMode  -  Specifies the method used to find the straight edge
//                ivaData                -  Internal Data structure
//                stepIndex              -  Step index (index at which to store
//                                          the results in the resuts array)
//
// Return Value : None (inline function)
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_EdgeTool2(Image* image,
                                  ROI* roi,
                                  int pProcessType,
                                  int pPolarity,
                                  int pKernelSize,
                                  int pWidth,
                                  float pMinThreshold,
                                  int pInterpolationType,
                                  int pColumnProcessingMode,
                                  IVA_Data* ivaData,
                                  int stepIndex)
{
    int success = 1;
    EdgeOptions2 edgeOptions;
    EdgeReport2* edgeReport = NULL;
    int i;
    int numObjectResults;
    IVA_Result* edgeResults;
    unsigned int visionInfo;


    //-------------------------------------------------------------------//
    //                     Edge Detector - Edge Tool3                    //
    //-------------------------------------------------------------------//

    edgeOptions.polarity = pPolarity;
    edgeOptions.kernelSize = pKernelSize;
    edgeOptions.width = pWidth;
    edgeOptions.minThreshold = pMinThreshold;
    edgeOptions.interpolationType = pInterpolationType;
    edgeOptions.columnProcessingMode = pColumnProcessingMode;

    // Finds edges along a path in the image.
    VisionErrChk(edgeReport = imaqEdgeTool4(image, roi, pProcessType, &edgeOptions, FALSE));

    // ////////////////////////////////////////
    // Store the results in the data structure.
    // ////////////////////////////////////////
    
    // First, delete all the results of this step (from a previous iteration)
    IVA_DisposeStepResults(ivaData, stepIndex);

    // Check if the image is calibrated.
    VisionErrChk(imaqGetVisionInfoTypes(image, &visionInfo));

    // If the image is calibrated, we also need to log the calibrated position (x and y) -> 9 results instead of 6
    numObjectResults = (visionInfo & IMAQ_VISIONINFO_CALIBRATION ? 9 : 6);

    ivaData->stepResults[stepIndex].numResults = edgeReport->numEdges * numObjectResults + 1;
    ivaData->stepResults[stepIndex].results = malloc (sizeof(IVA_Result) * ivaData->stepResults[stepIndex].numResults);
    edgeResults = ivaData->stepResults[stepIndex].results;
    
    if (edgeResults)
    {
        #if defined (IVA_STORE_RESULT_NAMES)
            sprintf(edgeResults->resultName, "# of Edges");
        #endif
        edgeResults->type = IVA_NUMERIC;
        edgeResults->resultVal.numVal = edgeReport->numEdges;
        edgeResults++;
        
        for (i = 0 ; i < edgeReport->numEdges ; i++)
        {
            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.X Position (Pix.)", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].position.x;
            edgeResults++;
            
            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Y Position (Pix.)", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].position.y;
            edgeResults++;
            
            if (visionInfo & IMAQ_VISIONINFO_CALIBRATION)
            {
                #if defined (IVA_STORE_RESULT_NAMES)
                    sprintf(edgeResults->resultName, "Edge %d.X Position (World)", i + 1);
                #endif
                edgeResults->type = IVA_NUMERIC;
                edgeResults->resultVal.numVal = edgeReport->edges[i].calibratedPosition.x;
                edgeResults++;

                #if defined (IVA_STORE_RESULT_NAMES)
                    sprintf(edgeResults->resultName, "Edge %d.Y Position (World)", i + 1);
                #endif
                edgeResults->type = IVA_NUMERIC;
                edgeResults->resultVal.numVal = edgeReport->edges[i].calibratedPosition.y;
                edgeResults++;
            }

            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Distance (Pix.)", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].distance;
            edgeResults++;

            if (visionInfo & IMAQ_VISIONINFO_CALIBRATION)
            {
                #if defined (IVA_STORE_RESULT_NAMES)
                    sprintf(edgeResults->resultName, "Edge %d.Distance (World)", i + 1);
                #endif
                edgeResults->type = IVA_NUMERIC;
                edgeResults->resultVal.numVal = edgeReport->edges[i].calibratedDistance;
                edgeResults++;
            }
            
            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Strength", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].magnitude;
            edgeResults++;

            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Noise", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].noisePeak;
            edgeResults++;
            
            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Rising", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].rising;
            edgeResults++;
        }
    }
    
Error:
    // Disposes temporary structures.
    imaqDispose(edgeReport);

    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_Classification_Multiple
//
// Description  : Classifies all the objects located in the given ROI.
//
// Parameters   : image     -  Input image
//                fileName  -  Character Set File Path
//                roi       -  Region of interest specifying the location of
//                             the object in the image.
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_Classification_Multiple(Image* image, char* fileName, ROI* roi)
{
    int    success = 1;
    ClassifierSession* classifierSession = NULL;    // Classifier Session
    char   description[255];                        // Classifier Description
    Image* binaryImage;                             // Binary Image
    ParticleClassifierPreprocessingOptions2 preprocessingOptions;
    int    numParticles;                            // Number of particles to classify
    int    i;
    ROI**  rois;                                    // ROIs around each particle to classify.
    ClassifierReport** classifierReports = NULL;


    //-------------------------------------------------------------------//
    //                             Classify                              //
    //-------------------------------------------------------------------//

    // Reads the classifier file and properties.
    VisionErrChk(classifierSession = imaqReadClassifierFile(NULL, fileName, IMAQ_CLASSIFIER_READ_ALL, NULL, NULL, description));

    // Create a binary image that will contain the segmented image.
    VisionErrChk(binaryImage = imaqCreateImage(IMAQ_IMAGE_U8, 7));

    // Get the preprocessing options from the classifier session.
    VisionErrChk(imaqGetParticleClassifierOptions2(classifierSession, &preprocessingOptions, NULL));

    // Segments the image.
    VisionErrChk(IVA_Classification_Segmentation(image, binaryImage, roi, preprocessingOptions));

    // Counts the number of particles to classify.
    VisionErrChk(imaqCountParticles(binaryImage, TRUE, &numParticles));

    // Allocates an array of region of interests.
    VisionErrChk(rois = (ROI**)malloc(numParticles * sizeof(ROI*)));

    for (i=0 ; i < numParticles ; i++)
    {
        // Creates a region of interest for each particle.
        VisionErrChk(rois[i] = imaqCreateROI());
    }

    // Get the ROIs of all individual particles.
    VisionErrChk(IVA_Classification_Extract_Particles(image, binaryImage, rois, numParticles));

    // Allocates the classifier reports for all objects in the image.
    VisionErrChk(classifierReports = (ClassifierReport**)malloc(numParticles * sizeof(ClassifierReport*)));

    // Classifies the object located in the given ROIs.
    for (i = 0 ; i < numParticles ; i++)
    {
        VisionErrChk(classifierReports[i] = imaqClassify(image, classifierSession, rois[i], NULL, 0));
    }

Error:
    for (i = 0 ; i < numParticles ; i++)
    {
        imaqDispose(rois[i]);
        imaqDispose(classifierReports[i]);
    }
    free(rois);
    free(classifierReports);
    imaqDispose(binaryImage);
    imaqDispose(classifierSession);

    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_InitData
//
// Description  : Initializes data for buffer management and results.
//
// Parameters   : # of steps
//                # of coordinate systems
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static IVA_Data* IVA_InitData(int numSteps, int numCoordSys)
{
    int success = 1;
    IVA_Data* ivaData = NULL;
    int i;


    // Allocate the data structure.
    VisionErrChk(ivaData = (IVA_Data*)malloc(sizeof (IVA_Data)));

    // Initializes the image pointers to NULL.
    for (i = 0 ; i < IVA_MAX_BUFFERS ; i++)
        ivaData->buffers[i] = NULL;

    // Initializes the steo results array to numSteps elements.
    ivaData->numSteps = numSteps;

    ivaData->stepResults = (IVA_StepResults*)malloc(ivaData->numSteps * sizeof(IVA_StepResults));
    for (i = 0 ; i < numSteps ; i++)
    {
        #if defined (IVA_STORE_RESULT_NAMES)
            sprintf(ivaData->stepResults[i].stepName, "");
        #endif
        ivaData->stepResults[i].numResults = 0;
        ivaData->stepResults[i].results = NULL;
    }

    // Create the coordinate systems
	ivaData->baseCoordinateSystems = NULL;
	ivaData->MeasurementSystems = NULL;
	if (numCoordSys)
	{
		ivaData->baseCoordinateSystems = (CoordinateSystem*)malloc(sizeof(CoordinateSystem) * numCoordSys);
		ivaData->MeasurementSystems = (CoordinateSystem*)malloc(sizeof(CoordinateSystem) * numCoordSys);
	}

    ivaData->numCoordSys = numCoordSys;

Error:
    return ivaData;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_DisposeData
//
// Description  : Releases the memory allocated in the IVA_Data structure
//
// Parameters   : ivaData  -  Internal data structure
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_DisposeData(IVA_Data* ivaData)
{
    int i;


    // Releases the memory allocated for the image buffers.
    for (i = 0 ; i < IVA_MAX_BUFFERS ; i++)
        imaqDispose(ivaData->buffers[i]);

    // Releases the memory allocated for the array of measurements.
    for (i = 0 ; i < ivaData->numSteps ; i++)
        IVA_DisposeStepResults(ivaData, i);

    free(ivaData->stepResults);

    // Dispose of coordinate systems
    if (ivaData->numCoordSys)
    {
        free(ivaData->baseCoordinateSystems);
        free(ivaData->MeasurementSystems);
    }

    free(ivaData);

    return TRUE;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_Classification_Segmentation
//
// Description  : Segments the classification image
//
// Parameters   : image                 - Input Image
//                imageMask             - Segmented image
//                roi                   - Region of Interest
//                preprocessingOptions  - Preprocessing options
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_Classification_Segmentation(Image* image, Image* imageMask, const ROI* roi, ParticleClassifierPreprocessingOptions2 preprocessingOptions)
{
    int              success = 1;
    Rect             boundingBox;
    Rect             expandedBox;
    Rect             reducedBox;
    ROI*             tmpROI;
    int              contourCount;
    int              i;
    int              xRes;
    int              yRes;
    int              expandX;
    int              expandY;
    int              leftOut;
    int              topOut;
    int              rightOut;
    int              bottomOut;
    int              leftExpand;
    int              rightExpand;
    int              leftSlack;
    int              rightSlack;
    int              topExpand;
    int              bottomExpand;
    int              topSlack;
    int              bottomSlack;
    int              useExpandedROI = 0;
    ContourID        roiContourID;
    ContourInfo2*    contourInfo;
    CoordinateSystem baseSystem;
    CoordinateSystem newSystem;
    float            thresholdMin;
    float            thresholdMax;
    ThresholdData*   thresholdData = NULL;
    Image*           tmpImageMask;
    short            lookupTable[256];
    ROIProfile*      roiProfile = NULL;


    // Get the bounding box of the region of interest.
    VisionErrChk(imaqGetROIBoundingBox(roi, &boundingBox));

    // Local Threshold method uses a kernel size, and the ROI must be larger than the kernel size for the function to work
    // Take care of making the ROI to extract larger if needed for the local threshold case
    if (preprocessingOptions.thresholdType == IMAQ_THRESHOLD_LOCAL) 
    {
        // Get the image size
        VisionErrChk(imaqGetImageSize(image, &xRes, &yRes));

        // Take into account clipping of ROI. Just get ROI that's within bounds of image.
        if ( min(xRes, (boundingBox.left + boundingBox.width) ) - max(0, boundingBox.left) < preprocessingOptions.localThresholdOptions.windowWidth || min(yRes, (boundingBox.top + boundingBox.height) ) - max(0, boundingBox.top) < preprocessingOptions.localThresholdOptions.windowHeight)
        {
            // The ROI is smaller than the kernel. Try to expand the kernel in the directions needed to meet the minimum size required by the kernel.
            expandX = (int)(preprocessingOptions.localThresholdOptions.windowWidth / 2);
            expandY = (int)(preprocessingOptions.localThresholdOptions.windowHeight / 2); 
            leftExpand = expandX;
            rightExpand = expandX;
            leftSlack = boundingBox.left;
            rightSlack = xRes - (boundingBox.left + boundingBox.width);

            if (leftExpand > leftSlack) 
            {
                rightExpand += (leftExpand - leftSlack);
            }

            if (rightExpand > rightSlack) 
            {
                leftExpand += (rightExpand - rightSlack);
            }

            leftOut = boundingBox.left - leftExpand;
            if (leftOut < 0) 
            {
                leftOut = 0;
            }
            rightOut = boundingBox.left + boundingBox.width + rightExpand;
            if (rightOut > xRes) 
            {
                rightOut = xRes;
            }

            topExpand = expandY;
            bottomExpand = expandY;
            topSlack = boundingBox.top;
            bottomSlack = yRes - (boundingBox.top + boundingBox.height);
            if (topExpand > topSlack) 
            {
                bottomExpand += (topExpand - topSlack);
            }
            if (bottomExpand > bottomSlack) 
            {
                topExpand += (bottomExpand - bottomSlack);
            }
            topOut = boundingBox.top - topExpand;
            if (topOut < 0) 
            {
                topOut = 0;
            }
            bottomOut = boundingBox.top + boundingBox.height + bottomExpand;
            if (bottomOut > yRes) {
                bottomOut = yRes;
            }
            expandedBox.left = leftOut;
            expandedBox.width = rightOut - leftOut;
            expandedBox.top = topOut;
            expandedBox.height = bottomOut - topOut;

            // Create the reduced Rect so after performing the local threshold, we can reduce the size back to the original ROI dimensions.
            reducedBox.top = max(boundingBox.top - topOut, 0);
            reducedBox.left = max(boundingBox.left - leftOut, 0);
            reducedBox.width = boundingBox.width + min(boundingBox.left, 0);
            reducedBox.height = boundingBox.height + min(boundingBox.top, 0);

            // Set this flag so the image can be reduced after performing the local threshold.
            useExpandedROI = 1;
        }

    }

    // if Expanded Box hasn't been updated, use the boundingBox passed in to extract. 
    if (useExpandedROI == 0)
    {
        expandedBox = boundingBox;
    }
    // Extract the region of interest into the mask image.
    VisionErrChk(imaqScale(imageMask, image, 1, 1, IMAQ_SCALE_LARGER, expandedBox));

    // Create a temporary ROI that will be used to mask the extracted image, to get rid of the pixels
    // outside of the rotated rectangle.
    VisionErrChk(tmpROI = imaqCreateROI());

    // Get the number of countours of the ROI
    VisionErrChk(contourCount = imaqGetContourCount(roi));

    // Copy all the contours of the original ROI in the temp ROI.
    for (i = 0 ; i < contourCount ; i++)
    {
        VisionErrChk(roiContourID = imaqGetContour(roi, i));
        VisionErrChk(imaqCopyContour(tmpROI, roi, roiContourID));
    }

    // Get the region of interest contour.
    VisionErrChk(roiContourID = imaqGetContour(roi, 0));
    VisionErrChk(contourInfo = imaqGetContourInfo2(roi, roiContourID));

    // If the ROI is a rotated rectangle, then compute the new location of the search ROI.
    if ((contourInfo->type == IMAQ_ROTATED_RECT) && (contourInfo->structure.rotatedRect->angle > 0.01))
    {
        baseSystem.origin.x = (boundingBox.left < 0 ? 0 : boundingBox.left);
        baseSystem.origin.y = (boundingBox.top < 0 ? 0 : boundingBox.top);
        baseSystem.angle = 0;
        baseSystem.axisOrientation = IMAQ_DIRECT;

        newSystem.origin = imaqMakePointFloat(0.0, 0.0);
        newSystem.angle = 0;
        newSystem.axisOrientation = IMAQ_DIRECT;

        VisionErrChk(imaqTransformROI2(tmpROI, &baseSystem, &newSystem));
    }

    // Create a temporary image.
    VisionErrChk(tmpImageMask = imaqCreateImage(IMAQ_IMAGE_U8, 7));

    // Determine the threshold range.
    switch (preprocessingOptions.thresholdType)
    {
        case IMAQ_THRESHOLD_MANUAL:
            thresholdMin = preprocessingOptions.manualThresholdRange.minValue;
            thresholdMax = preprocessingOptions.manualThresholdRange.maxValue;
            VisionErrChk(imaqThreshold(imageMask, imageMask, thresholdMin, thresholdMax, TRUE, 1.0));
            break;
        case IMAQ_THRESHOLD_AUTO:
            VisionErrChk(thresholdData = imaqAutoThreshold(tmpImageMask, imageMask, 2, preprocessingOptions.autoThresholdOptions.method));
            if (preprocessingOptions.autoThresholdOptions.particleType == IMAQ_PARTICLE_BRIGHT)
            {
                thresholdMin = (thresholdData->rangeMax > preprocessingOptions.autoThresholdOptions.limits.minValue ?
                    thresholdData->rangeMax : preprocessingOptions.autoThresholdOptions.limits.minValue);
                thresholdMax = 255;
            }
            else
            {
                thresholdMin = 0;
                thresholdMax = (thresholdData->rangeMax < preprocessingOptions.autoThresholdOptions.limits.maxValue ?
                    thresholdData->rangeMax : preprocessingOptions.autoThresholdOptions.limits.maxValue);
            }
            VisionErrChk(imaqThreshold(imageMask, imageMask, thresholdMin, thresholdMax, TRUE, 1.0));
            break;
        case IMAQ_THRESHOLD_LOCAL:
            VisionErrChk(imaqLocalThreshold(imageMask, imageMask, preprocessingOptions.localThresholdOptions.windowWidth, preprocessingOptions.localThresholdOptions.windowHeight, preprocessingOptions.localThresholdOptions.method, preprocessingOptions.localThresholdOptions.deviationWeight, preprocessingOptions.localThresholdOptions.particleType, 1.0));
            break;
        default:
            break;
    }

    // If the expanded ROI was used, reduce it so no particles are found outside requested ROI.
    if (useExpandedROI)
    {
        VisionErrChk(imaqScale(imageMask, imageMask, 1, 1, IMAQ_SCALE_LARGER, reducedBox));
    }

    // Cast the image to 8 bit.
    VisionErrChk(imaqCast(imageMask, imageMask, IMAQ_IMAGE_U8, NULL, 0));

    // Eliminates particles that touch the border of the image.
    if (preprocessingOptions.rejectBorder)
    {
        if ((contourInfo->type == IMAQ_ROTATED_RECT) && (contourInfo->structure.rotatedRect->angle > 0.01))
        {
            // Special case for the rotated rectangle.
            VisionErrChk(imaqLabel2(imageMask, imageMask, TRUE, NULL));
            
            lookupTable[0] = 0;
            for (i=1 ; i < 256 ; i++)
                lookupTable[i] = 1;

            VisionErrChk(roiProfile = imaqROIProfile(imageMask, tmpROI));

            for (i=0 ; i < roiProfile->report.dataCount ; i++)
                lookupTable[0] = roiProfile->report.profileData[i];

            VisionErrChk(imaqLookup(imageMask, imageMask, lookupTable, NULL));
        }
        else
            VisionErrChk(imaqRejectBorder(imageMask, imageMask, TRUE));
    }

    // Remove small particles.
    if (preprocessingOptions.numErosions)
        VisionErrChk(imaqSizeFilter(imageMask, imageMask, TRUE, preprocessingOptions.numErosions, FALSE, NULL));

    // If the rectangle is rotated, mask out the areas of the image that are not in the ROI.
    if ((contourInfo->type == IMAQ_ROTATED_RECT) && (contourInfo->structure.rotatedRect->angle > 0.01))
    {
        // Perform the mask
        VisionErrChk(imaqROIToMask(tmpImageMask, tmpROI, 255, imageMask, NULL));

        VisionErrChk(imaqAnd(imageMask, imageMask, tmpImageMask));
    }

    // Sets the mask offset.
    VisionErrChk(imaqSetMaskOffset(imageMask, imaqMakePoint (max(boundingBox.left,0), max(boundingBox.top,0))));

Error:
    imaqDispose(tmpROI);
    imaqDispose(contourInfo);
    imaqDispose(thresholdData);
    imaqDispose(tmpImageMask);
    imaqDispose(roiProfile);

    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_Classification_Extract_Particles
//
// Description  : Extracts the region of interests of the bounding rectangles of
//                all particles.
//
// Parameters   : image         - Input image
//                imageMask     - Image mask
//                rois          - Array of ROIs
//                numParticles  - Number of particles in the image mask
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_Classification_Extract_Particles(Image* image, Image* imageMask, ROI** rois, int numParticles)
{
    int       success = TRUE;
    ImageInfo imageInfo;
    ImageInfo maskInfo;
    int       i;
    Rect      particleBoundingRect;
    double    measurement;


    // Get the original image information.
    VisionErrChk(imaqGetImageInfo(image, &imageInfo));

    // Get the information of the image mask.
    VisionErrChk(imaqGetImageInfo(imageMask, &maskInfo));

    // Compute the region of interests around each particle.
    for (i=0 ; i < numParticles ; i++)
    {
        VisionErrChk(imaqMeasureParticle(imageMask, i, FALSE, IMAQ_MT_BOUNDING_RECT_TOP, &measurement));
        measurement = measurement  + maskInfo.yOffset - 5;
        particleBoundingRect.top = (measurement < 0 ? 0 : measurement);

        VisionErrChk(imaqMeasureParticle(imageMask, i, FALSE, IMAQ_MT_BOUNDING_RECT_LEFT, &measurement));
        measurement = measurement + maskInfo.xOffset - 5;
        particleBoundingRect.left = (measurement < 0 ? 0 : measurement);

        VisionErrChk(imaqMeasureParticle(imageMask, i, FALSE, IMAQ_MT_BOUNDING_RECT_BOTTOM, &measurement));
        measurement = measurement + maskInfo.yOffset + 5;
        measurement = (measurement > imageInfo.yRes - 1 ? imageInfo.yRes - 1 : measurement);
        particleBoundingRect.height = measurement - particleBoundingRect.top + 1;

        VisionErrChk(imaqMeasureParticle(imageMask, i, FALSE, IMAQ_MT_BOUNDING_RECT_RIGHT, &measurement));
        measurement = measurement + maskInfo.xOffset + 5;
        measurement = (measurement > imageInfo.xRes - 1 ? imageInfo.xRes - 1 : measurement);
        particleBoundingRect.width = measurement - particleBoundingRect.left + 1;

        VisionErrChk(imaqAddRectContour(rois[i], particleBoundingRect));
    }

Error:
    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_DisposeStepResults
//
// Description  : Dispose of the results of a specific step.
//
// Parameters   : ivaData    -  Internal data structure
//                stepIndex  -  step index
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_DisposeStepResults(IVA_Data* ivaData, int stepIndex)
{
    int i;

    
    for (i = 0 ; i < ivaData->stepResults[stepIndex].numResults ; i++)
    {
        if (ivaData->stepResults[stepIndex].results[i].type == IVA_STRING)
            free(ivaData->stepResults[stepIndex].results[i].resultVal.strVal);
    }

    free(ivaData->stepResults[stepIndex].results);

    return TRUE;
}


=======
 
//**************************************************************************
//* WARNING: This file was automatically generated.  Any changes you make  *
//*          to this file will be lost if you generate the file again.     *
//**************************************************************************
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <nivision.h>
#include <nimachinevision.h>
#include <windows.h>

// If you call Machine Vision functions in your script, add NIMachineVision.c to the project.

#define IVA_MAX_BUFFERS 10
#define IVA_STORE_RESULT_NAMES

#define VisionErrChk(Function) {if (!(Function)) {success = 0; goto Error;}}

typedef enum IVA_ResultType_Enum {IVA_NUMERIC, IVA_BOOLEAN, IVA_STRING} IVA_ResultType;

typedef union IVA_ResultValue_Struct    // A result in Vision Assistant can be of type double, BOOL or string.
{
    double numVal;
    BOOL   boolVal;
    char*  strVal;
} IVA_ResultValue;

typedef struct IVA_Result_Struct
{
#if defined (IVA_STORE_RESULT_NAMES)
    char resultName[256];           // Result name
#endif
    IVA_ResultType  type;           // Result type
    IVA_ResultValue resultVal;      // Result value
} IVA_Result;

typedef struct IVA_StepResultsStruct
{
#if defined (IVA_STORE_RESULT_NAMES)
    char stepName[256];             // Step name
#endif
    int         numResults;         // number of results created by the step
    IVA_Result* results;            // array of results
} IVA_StepResults;

typedef struct IVA_Data_Struct
{
    Image* buffers[IVA_MAX_BUFFERS];            // Vision Assistant Image Buffers
    IVA_StepResults* stepResults;              // Array of step results
    int numSteps;                               // Number of steps allocated in the stepResults array
    CoordinateSystem *baseCoordinateSystems;    // Base Coordinate Systems
    CoordinateSystem *MeasurementSystems;       // Measurement Coordinate Systems
    int numCoordSys;                            // Number of coordinate systems
} IVA_Data;



static IVA_Data* IVA_InitData(int numSteps, int numCoordSys);
static int IVA_DisposeData(IVA_Data* ivaData);
static int IVA_Classification_Segmentation(Image* image, Image* imageMask, const ROI* roi, ParticleClassifierPreprocessingOptions2 preprocessingOptions);
static int IVA_Classification_Extract_Particles(Image* image, Image* imageMask, ROI** rois, int numParticles);
static int IVA_DisposeStepResults(IVA_Data* ivaData, int stepIndex);
static int IVA_CLRExtractValue(Image* image);
static int IVA_EdgeTool2(Image* image,
                                  ROI* roi,
                                  int pProcessType,
                                  int pPolarity,
                                  int pKernelSize,
                                  int pWidth,
                                  float pMinThreshold,
                                  int pInterpolationType,
                                  int pColumnProcessingMode,
                                  IVA_Data* ivaData,
                                  int stepIndex);
static int IVA_Classification_Multiple(Image* image, char* fileName, ROI* roi);

int IVA_ProcessImage(Image *image)
{
	int success = 1;
    IVA_Data *ivaData;
    BCGOptions redBCGOptions;
    BCGOptions greenBCGOptions;
    BCGOptions blueBCGOptions;
    float kernel[9] = {-1,-1,-1,-1,10,-1,-1,-1,-1};
    ROI *roi;
    ROI *roi1;

    // Initializes internal data (buffers and array of points for caliper measurements)
    VisionErrChk(ivaData = IVA_InitData(6, 0));

    //-------------------------------------------------------------------//
    //                          Color Brightness                         //
    //-------------------------------------------------------------------//

    redBCGOptions.brightness = 164;
    redBCGOptions.contrast = 65.8000009655952;
    redBCGOptions.gamma = 1.86;

    greenBCGOptions.brightness = 164;
    greenBCGOptions.contrast = 65.8000009655952;
    greenBCGOptions.gamma = 1.86;

    blueBCGOptions.brightness = 164;
    blueBCGOptions.contrast = 65.8000009655952;
    blueBCGOptions.gamma = 1.86;

    // Applies brightness, contrast, and gamma correction to the image
    VisionErrChk(imaqColorBCGTransform(image, image, &redBCGOptions, &greenBCGOptions, &blueBCGOptions, NULL));

	VisionErrChk(IVA_CLRExtractValue(image));

    //-------------------------------------------------------------------//
    //                       Lookup Table: Inverse                       //
    //-------------------------------------------------------------------//

    // Inverts the pixel intensities of the image.
    VisionErrChk(imaqInverse(image, image, NULL));

    //-------------------------------------------------------------------//
    //                        Filters: Convolution                       //
    //-------------------------------------------------------------------//

    // Applies a linear filter to an image by convolving the image with a filtering kernel.
    VisionErrChk(imaqConvolve(image, image, kernel, 3, 3, 0, NULL));

    // Creates a new, empty region of interest.
    VisionErrChk(roi = imaqCreateROI());

    // Creates a new line ROI contour and adds the line to the provided ROI.
    VisionErrChk(imaqAddLineContour(roi, imaqMakePoint(34, 418),
                                         imaqMakePoint(810, 446)));

	VisionErrChk(IVA_EdgeTool2(image, roi, IMAQ_BEST, 
		IMAQ_SEARCH_FOR_ALL_EDGES, 5, 5, 10, IMAQ_BILINEAR_FIXED, 
		IMAQ_AVERAGE_COLUMNS, ivaData, 4));

    // Cleans up resources associated with the object
    imaqDispose(roi);

    // Creates a new, empty region of interest.
    VisionErrChk(roi1 = imaqCreateROI());

    // Creates a new rectangle ROI contour and adds the rectangle to the provided ROI.
    VisionErrChk(imaqAddRectContour(roi1, imaqMakeRect(17, 25, 929, 674)));

	VisionErrChk(IVA_Classification_Multiple(image, 
		"C:\\Users\\HP1\\Downloads\\leaf\\Leaf Classifier.clf", roi1));

    // Cleans up resources associated with the object
    imaqDispose(roi1);

    // Releases the memory allocated in the IVA_Data structure.
    IVA_DisposeData(ivaData);

Error:
	return success;
}

////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_CLRExtractValue
//
// Description  : Extracts the value plane from a color image.
//
// Parameters   : image  -  Input image
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_CLRExtractValue(Image* image)
{
    int success = 1;
    Image* plane;


    //-------------------------------------------------------------------//
    //                         Extract Color Plane                       //
    //-------------------------------------------------------------------//

    // Creates an 8 bit image that contains the extracted plane.
    VisionErrChk(plane = imaqCreateImage(IMAQ_IMAGE_U8, 7));

    // Extracts the value plane
    VisionErrChk(imaqExtractColorPlanes(image, IMAQ_HSV, NULL, NULL, plane));

    // Copies the color plane in the main image.
    VisionErrChk(imaqDuplicate(image, plane));

Error:
    imaqDispose(plane);

    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_EdgeTool2
//
// Description  : Finds edges along a path in an image.
//
// Parameters   : image                  -  Input image
//                roi                    -  Region of Interest
//                pProcessType           -  Determines which edges the function looks for. 
//                pPolarity              -  Specifies the polarity of the edges to be found
//                pKernelSize            -  Specifies the size of the edge detection kernel
//                pWidth                 -  Specifies the number of pixels averaged perpendicular
//                                          to the search direction to compute the edge profile
//                                          strength at each point along the search ROI
//                pMinThreshold          -  Specifies the minimum edge strength (gradient magnitude)
//                                          required for a detected edge
//                pInterpolationType     -  Specifies the interpolation method used to locate the edge position
//                pColumnProcessingMode  -  Specifies the method used to find the straight edge
//                ivaData                -  Internal Data structure
//                stepIndex              -  Step index (index at which to store
//                                          the results in the resuts array)
//
// Return Value : None (inline function)
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_EdgeTool2(Image* image,
                                  ROI* roi,
                                  int pProcessType,
                                  int pPolarity,
                                  int pKernelSize,
                                  int pWidth,
                                  float pMinThreshold,
                                  int pInterpolationType,
                                  int pColumnProcessingMode,
                                  IVA_Data* ivaData,
                                  int stepIndex)
{
    int success = 1;
    EdgeOptions2 edgeOptions;
    EdgeReport2* edgeReport = NULL;
    int i;
    int numObjectResults;
    IVA_Result* edgeResults;
    unsigned int visionInfo;


    //-------------------------------------------------------------------//
    //                     Edge Detector - Edge Tool3                    //
    //-------------------------------------------------------------------//

    edgeOptions.polarity = pPolarity;
    edgeOptions.kernelSize = pKernelSize;
    edgeOptions.width = pWidth;
    edgeOptions.minThreshold = pMinThreshold;
    edgeOptions.interpolationType = pInterpolationType;
    edgeOptions.columnProcessingMode = pColumnProcessingMode;

    // Finds edges along a path in the image.
    VisionErrChk(edgeReport = imaqEdgeTool4(image, roi, pProcessType, &edgeOptions, FALSE));

    // ////////////////////////////////////////
    // Store the results in the data structure.
    // ////////////////////////////////////////
    
    // First, delete all the results of this step (from a previous iteration)
    IVA_DisposeStepResults(ivaData, stepIndex);

    // Check if the image is calibrated.
    VisionErrChk(imaqGetVisionInfoTypes(image, &visionInfo));

    // If the image is calibrated, we also need to log the calibrated position (x and y) -> 9 results instead of 6
    numObjectResults = (visionInfo & IMAQ_VISIONINFO_CALIBRATION ? 9 : 6);

    ivaData->stepResults[stepIndex].numResults = edgeReport->numEdges * numObjectResults + 1;
    ivaData->stepResults[stepIndex].results = malloc (sizeof(IVA_Result) * ivaData->stepResults[stepIndex].numResults);
    edgeResults = ivaData->stepResults[stepIndex].results;
    
    if (edgeResults)
    {
        #if defined (IVA_STORE_RESULT_NAMES)
            sprintf(edgeResults->resultName, "# of Edges");
        #endif
        edgeResults->type = IVA_NUMERIC;
        edgeResults->resultVal.numVal = edgeReport->numEdges;
        edgeResults++;
        
        for (i = 0 ; i < edgeReport->numEdges ; i++)
        {
            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.X Position (Pix.)", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].position.x;
            edgeResults++;
            
            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Y Position (Pix.)", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].position.y;
            edgeResults++;
            
            if (visionInfo & IMAQ_VISIONINFO_CALIBRATION)
            {
                #if defined (IVA_STORE_RESULT_NAMES)
                    sprintf(edgeResults->resultName, "Edge %d.X Position (World)", i + 1);
                #endif
                edgeResults->type = IVA_NUMERIC;
                edgeResults->resultVal.numVal = edgeReport->edges[i].calibratedPosition.x;
                edgeResults++;

                #if defined (IVA_STORE_RESULT_NAMES)
                    sprintf(edgeResults->resultName, "Edge %d.Y Position (World)", i + 1);
                #endif
                edgeResults->type = IVA_NUMERIC;
                edgeResults->resultVal.numVal = edgeReport->edges[i].calibratedPosition.y;
                edgeResults++;
            }

            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Distance (Pix.)", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].distance;
            edgeResults++;

            if (visionInfo & IMAQ_VISIONINFO_CALIBRATION)
            {
                #if defined (IVA_STORE_RESULT_NAMES)
                    sprintf(edgeResults->resultName, "Edge %d.Distance (World)", i + 1);
                #endif
                edgeResults->type = IVA_NUMERIC;
                edgeResults->resultVal.numVal = edgeReport->edges[i].calibratedDistance;
                edgeResults++;
            }
            
            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Strength", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].magnitude;
            edgeResults++;

            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Noise", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].noisePeak;
            edgeResults++;
            
            #if defined (IVA_STORE_RESULT_NAMES)
                sprintf(edgeResults->resultName, "Edge %d.Rising", i + 1);
            #endif
            edgeResults->type = IVA_NUMERIC;
            edgeResults->resultVal.numVal = edgeReport->edges[i].rising;
            edgeResults++;
        }
    }
    
Error:
    // Disposes temporary structures.
    imaqDispose(edgeReport);

    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_Classification_Multiple
//
// Description  : Classifies all the objects located in the given ROI.
//
// Parameters   : image     -  Input image
//                fileName  -  Character Set File Path
//                roi       -  Region of interest specifying the location of
//                             the object in the image.
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_Classification_Multiple(Image* image, char* fileName, ROI* roi)
{
    int    success = 1;
    ClassifierSession* classifierSession = NULL;    // Classifier Session
    char   description[255];                        // Classifier Description
    Image* binaryImage;                             // Binary Image
    ParticleClassifierPreprocessingOptions2 preprocessingOptions;
    int    numParticles;                            // Number of particles to classify
    int    i;
    ROI**  rois;                                    // ROIs around each particle to classify.
    ClassifierReport** classifierReports = NULL;


    //-------------------------------------------------------------------//
    //                             Classify                              //
    //-------------------------------------------------------------------//

    // Reads the classifier file and properties.
    VisionErrChk(classifierSession = imaqReadClassifierFile(NULL, fileName, IMAQ_CLASSIFIER_READ_ALL, NULL, NULL, description));

    // Create a binary image that will contain the segmented image.
    VisionErrChk(binaryImage = imaqCreateImage(IMAQ_IMAGE_U8, 7));

    // Get the preprocessing options from the classifier session.
    VisionErrChk(imaqGetParticleClassifierOptions2(classifierSession, &preprocessingOptions, NULL));

    // Segments the image.
    VisionErrChk(IVA_Classification_Segmentation(image, binaryImage, roi, preprocessingOptions));

    // Counts the number of particles to classify.
    VisionErrChk(imaqCountParticles(binaryImage, TRUE, &numParticles));

    // Allocates an array of region of interests.
    VisionErrChk(rois = (ROI**)malloc(numParticles * sizeof(ROI*)));

    for (i=0 ; i < numParticles ; i++)
    {
        // Creates a region of interest for each particle.
        VisionErrChk(rois[i] = imaqCreateROI());
    }

    // Get the ROIs of all individual particles.
    VisionErrChk(IVA_Classification_Extract_Particles(image, binaryImage, rois, numParticles));

    // Allocates the classifier reports for all objects in the image.
    VisionErrChk(classifierReports = (ClassifierReport**)malloc(numParticles * sizeof(ClassifierReport*)));

    // Classifies the object located in the given ROIs.
    for (i = 0 ; i < numParticles ; i++)
    {
        VisionErrChk(classifierReports[i] = imaqClassify(image, classifierSession, rois[i], NULL, 0));
    }

Error:
    for (i = 0 ; i < numParticles ; i++)
    {
        imaqDispose(rois[i]);
        imaqDispose(classifierReports[i]);
    }
    free(rois);
    free(classifierReports);
    imaqDispose(binaryImage);
    imaqDispose(classifierSession);

    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_InitData
//
// Description  : Initializes data for buffer management and results.
//
// Parameters   : # of steps
//                # of coordinate systems
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static IVA_Data* IVA_InitData(int numSteps, int numCoordSys)
{
    int success = 1;
    IVA_Data* ivaData = NULL;
    int i;


    // Allocate the data structure.
    VisionErrChk(ivaData = (IVA_Data*)malloc(sizeof (IVA_Data)));

    // Initializes the image pointers to NULL.
    for (i = 0 ; i < IVA_MAX_BUFFERS ; i++)
        ivaData->buffers[i] = NULL;

    // Initializes the steo results array to numSteps elements.
    ivaData->numSteps = numSteps;

    ivaData->stepResults = (IVA_StepResults*)malloc(ivaData->numSteps * sizeof(IVA_StepResults));
    for (i = 0 ; i < numSteps ; i++)
    {
        #if defined (IVA_STORE_RESULT_NAMES)
            sprintf(ivaData->stepResults[i].stepName, "");
        #endif
        ivaData->stepResults[i].numResults = 0;
        ivaData->stepResults[i].results = NULL;
    }

    // Create the coordinate systems
	ivaData->baseCoordinateSystems = NULL;
	ivaData->MeasurementSystems = NULL;
	if (numCoordSys)
	{
		ivaData->baseCoordinateSystems = (CoordinateSystem*)malloc(sizeof(CoordinateSystem) * numCoordSys);
		ivaData->MeasurementSystems = (CoordinateSystem*)malloc(sizeof(CoordinateSystem) * numCoordSys);
	}

    ivaData->numCoordSys = numCoordSys;

Error:
    return ivaData;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_DisposeData
//
// Description  : Releases the memory allocated in the IVA_Data structure
//
// Parameters   : ivaData  -  Internal data structure
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_DisposeData(IVA_Data* ivaData)
{
    int i;


    // Releases the memory allocated for the image buffers.
    for (i = 0 ; i < IVA_MAX_BUFFERS ; i++)
        imaqDispose(ivaData->buffers[i]);

    // Releases the memory allocated for the array of measurements.
    for (i = 0 ; i < ivaData->numSteps ; i++)
        IVA_DisposeStepResults(ivaData, i);

    free(ivaData->stepResults);

    // Dispose of coordinate systems
    if (ivaData->numCoordSys)
    {
        free(ivaData->baseCoordinateSystems);
        free(ivaData->MeasurementSystems);
    }

    free(ivaData);

    return TRUE;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_Classification_Segmentation
//
// Description  : Segments the classification image
//
// Parameters   : image                 - Input Image
//                imageMask             - Segmented image
//                roi                   - Region of Interest
//                preprocessingOptions  - Preprocessing options
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_Classification_Segmentation(Image* image, Image* imageMask, const ROI* roi, ParticleClassifierPreprocessingOptions2 preprocessingOptions)
{
    int              success = 1;
    Rect             boundingBox;
    Rect             expandedBox;
    Rect             reducedBox;
    ROI*             tmpROI;
    int              contourCount;
    int              i;
    int              xRes;
    int              yRes;
    int              expandX;
    int              expandY;
    int              leftOut;
    int              topOut;
    int              rightOut;
    int              bottomOut;
    int              leftExpand;
    int              rightExpand;
    int              leftSlack;
    int              rightSlack;
    int              topExpand;
    int              bottomExpand;
    int              topSlack;
    int              bottomSlack;
    int              useExpandedROI = 0;
    ContourID        roiContourID;
    ContourInfo2*    contourInfo;
    CoordinateSystem baseSystem;
    CoordinateSystem newSystem;
    float            thresholdMin;
    float            thresholdMax;
    ThresholdData*   thresholdData = NULL;
    Image*           tmpImageMask;
    short            lookupTable[256];
    ROIProfile*      roiProfile = NULL;


    // Get the bounding box of the region of interest.
    VisionErrChk(imaqGetROIBoundingBox(roi, &boundingBox));

    // Local Threshold method uses a kernel size, and the ROI must be larger than the kernel size for the function to work
    // Take care of making the ROI to extract larger if needed for the local threshold case
    if (preprocessingOptions.thresholdType == IMAQ_THRESHOLD_LOCAL) 
    {
        // Get the image size
        VisionErrChk(imaqGetImageSize(image, &xRes, &yRes));

        // Take into account clipping of ROI. Just get ROI that's within bounds of image.
        if ( min(xRes, (boundingBox.left + boundingBox.width) ) - max(0, boundingBox.left) < preprocessingOptions.localThresholdOptions.windowWidth || min(yRes, (boundingBox.top + boundingBox.height) ) - max(0, boundingBox.top) < preprocessingOptions.localThresholdOptions.windowHeight)
        {
            // The ROI is smaller than the kernel. Try to expand the kernel in the directions needed to meet the minimum size required by the kernel.
            expandX = (int)(preprocessingOptions.localThresholdOptions.windowWidth / 2);
            expandY = (int)(preprocessingOptions.localThresholdOptions.windowHeight / 2); 
            leftExpand = expandX;
            rightExpand = expandX;
            leftSlack = boundingBox.left;
            rightSlack = xRes - (boundingBox.left + boundingBox.width);

            if (leftExpand > leftSlack) 
            {
                rightExpand += (leftExpand - leftSlack);
            }

            if (rightExpand > rightSlack) 
            {
                leftExpand += (rightExpand - rightSlack);
            }

            leftOut = boundingBox.left - leftExpand;
            if (leftOut < 0) 
            {
                leftOut = 0;
            }
            rightOut = boundingBox.left + boundingBox.width + rightExpand;
            if (rightOut > xRes) 
            {
                rightOut = xRes;
            }

            topExpand = expandY;
            bottomExpand = expandY;
            topSlack = boundingBox.top;
            bottomSlack = yRes - (boundingBox.top + boundingBox.height);
            if (topExpand > topSlack) 
            {
                bottomExpand += (topExpand - topSlack);
            }
            if (bottomExpand > bottomSlack) 
            {
                topExpand += (bottomExpand - bottomSlack);
            }
            topOut = boundingBox.top - topExpand;
            if (topOut < 0) 
            {
                topOut = 0;
            }
            bottomOut = boundingBox.top + boundingBox.height + bottomExpand;
            if (bottomOut > yRes) {
                bottomOut = yRes;
            }
            expandedBox.left = leftOut;
            expandedBox.width = rightOut - leftOut;
            expandedBox.top = topOut;
            expandedBox.height = bottomOut - topOut;

            // Create the reduced Rect so after performing the local threshold, we can reduce the size back to the original ROI dimensions.
            reducedBox.top = max(boundingBox.top - topOut, 0);
            reducedBox.left = max(boundingBox.left - leftOut, 0);
            reducedBox.width = boundingBox.width + min(boundingBox.left, 0);
            reducedBox.height = boundingBox.height + min(boundingBox.top, 0);

            // Set this flag so the image can be reduced after performing the local threshold.
            useExpandedROI = 1;
        }

    }

    // if Expanded Box hasn't been updated, use the boundingBox passed in to extract. 
    if (useExpandedROI == 0)
    {
        expandedBox = boundingBox;
    }
    // Extract the region of interest into the mask image.
    VisionErrChk(imaqScale(imageMask, image, 1, 1, IMAQ_SCALE_LARGER, expandedBox));

    // Create a temporary ROI that will be used to mask the extracted image, to get rid of the pixels
    // outside of the rotated rectangle.
    VisionErrChk(tmpROI = imaqCreateROI());

    // Get the number of countours of the ROI
    VisionErrChk(contourCount = imaqGetContourCount(roi));

    // Copy all the contours of the original ROI in the temp ROI.
    for (i = 0 ; i < contourCount ; i++)
    {
        VisionErrChk(roiContourID = imaqGetContour(roi, i));
        VisionErrChk(imaqCopyContour(tmpROI, roi, roiContourID));
    }

    // Get the region of interest contour.
    VisionErrChk(roiContourID = imaqGetContour(roi, 0));
    VisionErrChk(contourInfo = imaqGetContourInfo2(roi, roiContourID));

    // If the ROI is a rotated rectangle, then compute the new location of the search ROI.
    if ((contourInfo->type == IMAQ_ROTATED_RECT) && (contourInfo->structure.rotatedRect->angle > 0.01))
    {
        baseSystem.origin.x = (boundingBox.left < 0 ? 0 : boundingBox.left);
        baseSystem.origin.y = (boundingBox.top < 0 ? 0 : boundingBox.top);
        baseSystem.angle = 0;
        baseSystem.axisOrientation = IMAQ_DIRECT;

        newSystem.origin = imaqMakePointFloat(0.0, 0.0);
        newSystem.angle = 0;
        newSystem.axisOrientation = IMAQ_DIRECT;

        VisionErrChk(imaqTransformROI2(tmpROI, &baseSystem, &newSystem));
    }

    // Create a temporary image.
    VisionErrChk(tmpImageMask = imaqCreateImage(IMAQ_IMAGE_U8, 7));

    // Determine the threshold range.
    switch (preprocessingOptions.thresholdType)
    {
        case IMAQ_THRESHOLD_MANUAL:
            thresholdMin = preprocessingOptions.manualThresholdRange.minValue;
            thresholdMax = preprocessingOptions.manualThresholdRange.maxValue;
            VisionErrChk(imaqThreshold(imageMask, imageMask, thresholdMin, thresholdMax, TRUE, 1.0));
            break;
        case IMAQ_THRESHOLD_AUTO:
            VisionErrChk(thresholdData = imaqAutoThreshold(tmpImageMask, imageMask, 2, preprocessingOptions.autoThresholdOptions.method));
            if (preprocessingOptions.autoThresholdOptions.particleType == IMAQ_PARTICLE_BRIGHT)
            {
                thresholdMin = (thresholdData->rangeMax > preprocessingOptions.autoThresholdOptions.limits.minValue ?
                    thresholdData->rangeMax : preprocessingOptions.autoThresholdOptions.limits.minValue);
                thresholdMax = 255;
            }
            else
            {
                thresholdMin = 0;
                thresholdMax = (thresholdData->rangeMax < preprocessingOptions.autoThresholdOptions.limits.maxValue ?
                    thresholdData->rangeMax : preprocessingOptions.autoThresholdOptions.limits.maxValue);
            }
            VisionErrChk(imaqThreshold(imageMask, imageMask, thresholdMin, thresholdMax, TRUE, 1.0));
            break;
        case IMAQ_THRESHOLD_LOCAL:
            VisionErrChk(imaqLocalThreshold(imageMask, imageMask, preprocessingOptions.localThresholdOptions.windowWidth, preprocessingOptions.localThresholdOptions.windowHeight, preprocessingOptions.localThresholdOptions.method, preprocessingOptions.localThresholdOptions.deviationWeight, preprocessingOptions.localThresholdOptions.particleType, 1.0));
            break;
        default:
            break;
    }

    // If the expanded ROI was used, reduce it so no particles are found outside requested ROI.
    if (useExpandedROI)
    {
        VisionErrChk(imaqScale(imageMask, imageMask, 1, 1, IMAQ_SCALE_LARGER, reducedBox));
    }

    // Cast the image to 8 bit.
    VisionErrChk(imaqCast(imageMask, imageMask, IMAQ_IMAGE_U8, NULL, 0));

    // Eliminates particles that touch the border of the image.
    if (preprocessingOptions.rejectBorder)
    {
        if ((contourInfo->type == IMAQ_ROTATED_RECT) && (contourInfo->structure.rotatedRect->angle > 0.01))
        {
            // Special case for the rotated rectangle.
            VisionErrChk(imaqLabel2(imageMask, imageMask, TRUE, NULL));
            
            lookupTable[0] = 0;
            for (i=1 ; i < 256 ; i++)
                lookupTable[i] = 1;

            VisionErrChk(roiProfile = imaqROIProfile(imageMask, tmpROI));

            for (i=0 ; i < roiProfile->report.dataCount ; i++)
                lookupTable[0] = roiProfile->report.profileData[i];

            VisionErrChk(imaqLookup(imageMask, imageMask, lookupTable, NULL));
        }
        else
            VisionErrChk(imaqRejectBorder(imageMask, imageMask, TRUE));
    }

    // Remove small particles.
    if (preprocessingOptions.numErosions)
        VisionErrChk(imaqSizeFilter(imageMask, imageMask, TRUE, preprocessingOptions.numErosions, FALSE, NULL));

    // If the rectangle is rotated, mask out the areas of the image that are not in the ROI.
    if ((contourInfo->type == IMAQ_ROTATED_RECT) && (contourInfo->structure.rotatedRect->angle > 0.01))
    {
        // Perform the mask
        VisionErrChk(imaqROIToMask(tmpImageMask, tmpROI, 255, imageMask, NULL));

        VisionErrChk(imaqAnd(imageMask, imageMask, tmpImageMask));
    }

    // Sets the mask offset.
    VisionErrChk(imaqSetMaskOffset(imageMask, imaqMakePoint (max(boundingBox.left,0), max(boundingBox.top,0))));

Error:
    imaqDispose(tmpROI);
    imaqDispose(contourInfo);
    imaqDispose(thresholdData);
    imaqDispose(tmpImageMask);
    imaqDispose(roiProfile);

    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_Classification_Extract_Particles
//
// Description  : Extracts the region of interests of the bounding rectangles of
//                all particles.
//
// Parameters   : image         - Input image
//                imageMask     - Image mask
//                rois          - Array of ROIs
//                numParticles  - Number of particles in the image mask
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_Classification_Extract_Particles(Image* image, Image* imageMask, ROI** rois, int numParticles)
{
    int       success = TRUE;
    ImageInfo imageInfo;
    ImageInfo maskInfo;
    int       i;
    Rect      particleBoundingRect;
    double    measurement;


    // Get the original image information.
    VisionErrChk(imaqGetImageInfo(image, &imageInfo));

    // Get the information of the image mask.
    VisionErrChk(imaqGetImageInfo(imageMask, &maskInfo));

    // Compute the region of interests around each particle.
    for (i=0 ; i < numParticles ; i++)
    {
        VisionErrChk(imaqMeasureParticle(imageMask, i, FALSE, IMAQ_MT_BOUNDING_RECT_TOP, &measurement));
        measurement = measurement  + maskInfo.yOffset - 5;
        particleBoundingRect.top = (measurement < 0 ? 0 : measurement);

        VisionErrChk(imaqMeasureParticle(imageMask, i, FALSE, IMAQ_MT_BOUNDING_RECT_LEFT, &measurement));
        measurement = measurement + maskInfo.xOffset - 5;
        particleBoundingRect.left = (measurement < 0 ? 0 : measurement);

        VisionErrChk(imaqMeasureParticle(imageMask, i, FALSE, IMAQ_MT_BOUNDING_RECT_BOTTOM, &measurement));
        measurement = measurement + maskInfo.yOffset + 5;
        measurement = (measurement > imageInfo.yRes - 1 ? imageInfo.yRes - 1 : measurement);
        particleBoundingRect.height = measurement - particleBoundingRect.top + 1;

        VisionErrChk(imaqMeasureParticle(imageMask, i, FALSE, IMAQ_MT_BOUNDING_RECT_RIGHT, &measurement));
        measurement = measurement + maskInfo.xOffset + 5;
        measurement = (measurement > imageInfo.xRes - 1 ? imageInfo.xRes - 1 : measurement);
        particleBoundingRect.width = measurement - particleBoundingRect.left + 1;

        VisionErrChk(imaqAddRectContour(rois[i], particleBoundingRect));
    }

Error:
    return success;
}


////////////////////////////////////////////////////////////////////////////////
//
// Function Name: IVA_DisposeStepResults
//
// Description  : Dispose of the results of a specific step.
//
// Parameters   : ivaData    -  Internal data structure
//                stepIndex  -  step index
//
// Return Value : success
//
////////////////////////////////////////////////////////////////////////////////
static int IVA_DisposeStepResults(IVA_Data* ivaData, int stepIndex)
{
    int i;

    
    for (i = 0 ; i < ivaData->stepResults[stepIndex].numResults ; i++)
    {
        if (ivaData->stepResults[stepIndex].results[i].type == IVA_STRING)
            free(ivaData->stepResults[stepIndex].results[i].resultVal.strVal);
    }

    free(ivaData->stepResults[stepIndex].results);

    return TRUE;
}


>>>>>>> d97195687ee6623f7f62a98bbc5b544072e2f0bf
