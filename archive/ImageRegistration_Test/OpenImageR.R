install.packages("OpenImageR")
library(OpenImageR)
library(tiff)

img_vis_path <- paste0(vis_path, "spatial/tissue_lowres_image.png")

img_xe_path <- paste0(xe_path, "../Xenium_ome_tif.png")

img_xe_path <- "~/Downloads/morphology_focus.ome.tif"

img_vis = OpenImageR::readImage(img_vis_path)
img_xe = tiff:readTIFF(img_xe_path)

print(dim(img_vis))
600 543   3

print(dim(img_xe))
582 803   4

# TODO: does not how to fix the 3rd dimension is not 3 but 4, probably due to matte = TRUE for the screenshot, which is a "color" of transparency
dim(img_xe)[3] <- 3

image_random_path <- "~/Desktop/image.png"
img_random = OpenImageR::readImage(img_v)

img_xe_test <-List_2_Array(img_xe)

r = ncol(img_vis)
c = nrow(img_vis)

trans_mtx <- matrix(
  c(
    8.82797498e-02, -1.91831377e+00,  1.63476055e+04,
    1.84141210e+00,  5.96797885e-02,  4.12499099e+03#,
    # -3.95225478e-07, -4.66405945e-06,  1.03706895e+00
  ),
  nrow = 2, byrow= TRUE
)


res_3d = warpAffine(img = img_vis,
                    M = trans_mtx,
                    R = r,
                    C = c,
                    verbose = TRUE)

imageShow(res_3d, clear_viewer = FALSE)
