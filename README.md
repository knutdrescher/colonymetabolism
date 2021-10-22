# colonymetabolism

Overview of functions in this repository: 

deadCellFractionByHeight.m is a Matlab function that calculates the ratio between the live and dead cells for the whole colony and for cells in a shell of thickness 30 µm from the colony surface. 

fluorescenceRatioByHeight.m is a Matlab function that calculates the ratio between mRuby2 and GFP fluorescence levels for the whole colony. This function also calculates the biovolume of the whole colony and the biovolume of the shell of thickness 30 µm from the colony surface. 

getColonyMask.m is a Matlab function that performs the segmentation of the colony. This function is called by the two functions deadCellFractionByHeight.m and fluorescenceRatioByHeight.m

getColonyOutline.m is a Matlab function that corrects substrate tilting of the colony and calculates colony shape and perimeter. 

horizontalAndVerticalIntensityProfile.m is a Matlab function that calculates the fluorescence intensity of the vertical distance to colony surface and horizontal distance to colony surface. 
