library(png)
img_xe <- readPNG("~/Desktop/pink_source.png")
img_vis <- readPNG("~/Desktop/blue_target.png")

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
length(which(gradient != -Inf))
# 225770
# > dim(gradient)
# [1] 404 594
# > 404*594
# [1] 239976