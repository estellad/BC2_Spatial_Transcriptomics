library(png)
img_xe <- readPNG(paste0(xe_path, "../Xenium_ome_tif.png"))
img_vis <- readPNG(paste0(vis_path, "spatial/tissue_lowres_image.png"))

# library(jpeg) 
library(mmand) 
library(RNiftyReg) 
# Read images and convert to greyscale 
source <- apply(img_xe, 1:2, mean) 
target <- apply(img_vis, 1:2, mean) 

# Register images
result <- niftyreg(source, target) 

# Calculate morphological gradient 
kernel <- shapeKernel(c(3,3), type="diamond") 
gradient <- dilate(result$image,kernel) - erode(result$image,kernel)

# Display the results 
display(target) 
display(threshold(gradient, method="kmeans"), add=TRUE, col="red")


source_skewed <- applyTransform(forward(result), source)
display(source_skewed)

display(target)
display(source)
length(which(gradient != -Inf))
# 505420
# > dim(gradient)
# [1]  848 1356
# > 848*1356
# [1] 1149888