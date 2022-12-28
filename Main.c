////////////////////////////////////////////////////////////////////////////////
//
// Author : Vision Assistant
// Purpose: This file implements the algorithm prototyped in Vision Assistant.
//
// WARNING: This file was automatically generated.
//          Any changes you make to this file will be lost if you generate the
//          file again.
//
////////////////////////////////////////////////////////////////////////////////

//==============================================================================
//     Includes
//==============================================================================

#include <stdio.h>
#include <nivision.h>
#include "ImageProcessing.h"


//==============================================================================
//  Defines
//==============================================================================

#define DISPLAY_WINDOW 0


//==============================================================================
//  Main Function
//==============================================================================

int main (int argc, char *argv[])
{
    int success = 1;
    int err = 0;
    char** imagePath;       // Image Path
    int cancelled;
    ImageType imageType;    // Image Type
    Image* image;           // Image


    // IMAQ Vision creates windows in a separate thread
    imaqSetWindowThreadPolicy(IMAQ_SEPARATE_THREAD);

    // Display the Load Image dialog
    imagePath = imaqLoadImagePopup(NULL, "*.*", NULL, "Open Image", FALSE, IMAQ_BUTTON_LOAD, 0, 0, 1, 0, &cancelled, NULL);

    if (!cancelled)
    {
        // Get the type of the image file to create an image of the right type
        imaqGetFileInfo(imagePath[0], NULL, NULL, NULL, NULL, NULL, &imageType);

        // Create an IMAQ Vision image
        image = imaqCreateImage(imageType, 7);

        // Read the image from disk
        imaqReadFile(image, imagePath[0], NULL, NULL);

        // Vision Assistant Algorithm
        success = IVA_ProcessImage(image);
        if (!success)
            err = imaqGetLastError();

        // Display the image
        imaqMoveWindow(DISPLAY_WINDOW, imaqMakePoint(0,0));
        imaqSetWindowPalette(DISPLAY_WINDOW, IMAQ_PALETTE_GRAY, NULL, 0);
        imaqDisplayImage(image, DISPLAY_WINDOW, TRUE);

        // Wait for a key press before exiting
        printf ("Press Enter to exit.\n");
        getchar();

        // Dispose resources
        imaqDispose(image);
    }

    imaqDispose(imagePath);

    return 0;
}
