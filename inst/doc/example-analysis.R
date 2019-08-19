## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE------------------------------------------------
library(dplyr)
library(eSDM)
library(lwgeom)
library(sf)

source(system.file("eSDM_vignette_helper.R", package = "eSDM"), local = TRUE, echo = FALSE)

## ------------------------------------------------------------------------
# Import, process, and plot Model_B predictions
# model.b <- read.csv("Predictions_Beckeretal2016.csv")
model.b.sf <- readRDS(system.file("extdata/Predictions_Beckeretal2016.rds", package = "eSDM")) %>% 
  eSDM::pts2poly_centroids(0.09 / 2, crs = 4326) %>%
  st_wrap_dateline() %>%
  st_set_agr("constant")

model.b.sf

# Make base map
map.world <- eSDM::gshhg.l.L16

# Other option for making base map
# map.world <- st_geometry(st_as_sf(maps::map('world', plot = FALSE, fill = TRUE)))

## ---- fig.width=7--------------------------------------------------------
plot_sf_3panel(model.b.sf, "pred_bm", main.txt = "Model_B - ", map.base = map.world)

## ------------------------------------------------------------------------
# Import, process, and plot Model_H predictions
# model.h <- read.csv("Predictions_Hazenetal2017.csv")
model.h.sf <- readRDS(system.file("extdata/Predictions_Hazenetal2017.rds", package = "eSDM")) %>% 
  dplyr::select(lon, lat, pred_bm, se) %>%
  eSDM::pts2poly_centroids(0.25 / 2, crs = 4326, agr = "constant")

model.h.sf

## ---- fig.width=7--------------------------------------------------------
plot_sf_3panel(model.h.sf, "pred_bm", main.txt = "Model_H - ", map.base = map.world)

## ---- fig.width=7--------------------------------------------------------
# Import, process, and plot Model_R predictions
# model.r <- st_read("Shapefiles/Predictions_Redfernetal2017.shp")
model.r.sf <- readRDS(system.file("extdata/Predictions_Redfernetal2017.rds", package = "eSDM")) %>% 
  st_make_valid() %>% #
  st_set_agr("constant")

model.r.sf

## ---- fig.width=7--------------------------------------------------------
plot_sf_3panel(model.r.sf, "pred_bm", main.txt = "Model_R - ", map.base = map.world)

## ---- eval=FALSE---------------------------------------------------------
#  # Example code for converting raster to sf object; code not run
#  logo <- raster::raster(system.file("external/rlogo.grd", package="raster"))
#  logo.sf <- as(logo, "SpatialPolygonsDataFrame") %>%
#    sf::st_as_sf()

## ------------------------------------------------------------------------
# Study area polygon
poly.study <- st_read(system.file("extdata/Shapefiles/Study_Area_CCE.shp", package = "eSDM")) %>%
  st_geometry() %>% 
  st_transform(st_crs(model.r.sf))

# Erasing polygon; clip to the buffered study area polygon reduces future computation time
poly.erase <- eSDM::gshhg.l.L16 %>%
  st_transform(st_crs(model.r.sf)) %>%
  lwgeom::st_make_valid() %>%
  st_crop(st_buffer(poly.study, 100000))

# Create the base geometry; st_erase() function defined in eSDM_vignette_helper.R
base.geom <- model.r.sf %>%
  st_geometry() %>%
  st_erase(poly.erase) %>% 
  st_intersection(poly.study) %>%
  st_cast("MULTIPOLYGON")

## ---- fig.width=5, fig.height=7------------------------------------------
# Visualize the base geometry
plot(st_transform(base.geom, 4326), col = NA, border = "black", axes = TRUE)
plot(map.world, add = TRUE, col = "tan", border = NA)
graphics::box()

## ------------------------------------------------------------------------
# Convert SE values to variance
model.b.sf <- model.b.sf %>% 
  mutate(variance = se^2) %>% 
  dplyr::select(pred_bm, se, variance)
model.h.sf <- model.h.sf %>% 
  mutate(variance = se^2) %>% 
  dplyr::select(pred_bm, se, variance)
model.r.sf <- model.r.sf %>% 
  mutate(variance = se^2) %>% 
  dplyr::select(pred_bm, se, variance)

## ------------------------------------------------------------------------
# Perform overlay, and convert overlaid uncertainty values to SEs
over1.sf <- eSDM::overlay_sdm(base.geom, st_transform(model.b.sf, st_crs(base.geom)), c("pred_bm", "variance"), 50) %>% 
  mutate(se = sqrt(variance))
over2.sf <- eSDM::overlay_sdm(base.geom, st_transform(model.h.sf, st_crs(base.geom)), c("pred_bm", "variance"), 50) %>% 
  mutate(se = sqrt(variance))
over3.sf <- eSDM::overlay_sdm(base.geom, model.r.sf, c("pred_bm", "variance"), 50) %>% 
  mutate(se = sqrt(variance))

over3.sfb <- model.r.sf %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(pred_bm, variance) %>% 
  st_sf(geometry = base.geom, agr = "constant") %>% 
  dplyr::mutate(se = sqrt(variance))
all.equal(over3.sf, over3.sfb)
rm(over3.sfb)


## ---- fig.width=7, eval=FALSE--------------------------------------------
#  # Plot overlaid predictions; code not run
#  plot_sf_3panel(over1.sf, "pred_bm", main.txt = "Overlaid Model_B - ", map.base = map.world)
#  plot_sf_3panel(over2.sf, "pred_bm", main.txt = "Overlaid Model_H - ", map.base = map.world)
#  plot_sf_3panel(over3.sf, "pred_bm", main.txt = "Overlaid Model_R - ", map.base = map.world)

## ------------------------------------------------------------------------
# Import and process validation data
# valid.data <- read.csv("eSDM_Validation_data_all.csv", stringsAsFactors = FALSE)
valid.data <- readRDS(system.file("extdata/eSDM_Validation_data_all.rds", package = "eSDM"))%>% 
  arrange(source, lat, lon) %>% 
  mutate(pres_abs = ifelse(pres_abs > 0, 1, 0)) %>% #For demonstration purposes; pres_abs column is already binary
  st_as_sf(coords = c("lon", "lat"), crs = 4326, agr = "constant") %>%
  st_transform(st_crs(base.geom))

# Extract the line transect and home range validation data
valid.data.lt <- valid.data %>% filter(source == "Becker_et_al_2016")
valid.data.hr <- valid.data %>% filter(source == "Irvine_et_al_2014")

# Summarize the number of presence and absence points
valid.data %>% 
  st_set_geometry(NULL) %>% 
  group_by(source) %>%  
  summarize(pres = sum(pres_abs == 1), 
            abs = sum(pres_abs == 0)) %>% 
  knitr::kable(caption = "Validation data summary")

## ---- eval=FALSE---------------------------------------------------------
#  # Calculate evaluation metrics with different validation data sets; code not run
#  names.1 <- c(
#    "Model_B_orig", "Model_H_orig", "Model_R_orig",
#    "Model_B_overlaid", "Model_H_overlaid", "Model_R_overlaid"
#  )
#  
#  eval.lt <- data.frame(do.call(rbind, list(
#    eSDM::evaluation_metrics(model.b.sf, 1, st_transform(valid.data.lt, 4326), "pres_abs"),
#    eSDM::evaluation_metrics(model.h.sf, 1, st_transform(valid.data.lt, 4326), "pres_abs"),
#    eSDM::evaluation_metrics(model.r.sf, 1, valid.data.lt, "pres_abs"),
#    eSDM::evaluation_metrics(over1.sf, 1, valid.data.lt, "pres_abs"),
#    eSDM::evaluation_metrics(over2.sf, 1, valid.data.lt, "pres_abs"),
#    eSDM::evaluation_metrics(over3.sf, 1, valid.data.lt, "pres_abs")
#  ))) %>%
#    mutate(Preds = names.1) %>%
#    dplyr::select(Preds, AUC_LT = X1, TSS_LT = X2)
#  
#  eval.hr <- data.frame(do.call(rbind, list(
#    eSDM::evaluation_metrics(model.b.sf, 1, st_transform(valid.data.hr, 4326), "pres_abs"),
#    eSDM::evaluation_metrics(model.h.sf, 1, st_transform(valid.data.hr, 4326), "pres_abs"),
#    eSDM::evaluation_metrics(model.r.sf, 1, valid.data.hr, "pres_abs"),
#    eSDM::evaluation_metrics(over1.sf, 1, valid.data.hr, "pres_abs"),
#    eSDM::evaluation_metrics(over2.sf, 1, valid.data.hr, "pres_abs"),
#    eSDM::evaluation_metrics(over3.sf, 1, valid.data.hr, "pres_abs")
#  ))) %>%
#    mutate(Preds = names.1) %>%
#    dplyr::select(Preds, AUC_HR = X1, TSS_HR = X2)
#  
#  eval.combo <- data.frame(do.call(rbind, list(
#    eSDM::evaluation_metrics(model.b.sf, 1, st_transform(valid.data, 4326), "pres_abs"),
#    eSDM::evaluation_metrics(model.h.sf, 1, st_transform(valid.data, 4326), "pres_abs"),
#    eSDM::evaluation_metrics(model.r.sf, 1, valid.data, "pres_abs"),
#    eSDM::evaluation_metrics(over1.sf, 1, valid.data, "pres_abs"),
#    eSDM::evaluation_metrics(over2.sf, 1, valid.data, "pres_abs"),
#    eSDM::evaluation_metrics(over3.sf, 1, valid.data, "pres_abs")
#  ))) %>%
#    mutate(Preds = names.1) %>%
#    dplyr::select(Preds, AUC = X1, TSS = X2)

## ------------------------------------------------------------------------
read.csv(system.file("extdata/Table3.csv", package = "eSDM")) %>%
  filter(grepl("Model_", Predictions)) %>% 
  dplyr::select(Predictions, AUC, TSS, `AUC-LT` = AUC.LT, `TSS-LT` = TSS.LT, 
         `AUC-HR` = AUC.HR, `TSS-HR` = TSS.HR) %>% 
  knitr::kable(caption = "Evaluation metrics", digits = 3, align = "lcccccc")

## ------------------------------------------------------------------------
# Rescale predictions
over.sf <- bind_cols(
  over1.sf %>% st_set_geometry(NULL) %>% dplyr::select(pred_bm1 = pred_bm, var1 = variance), 
  over2.sf %>% st_set_geometry(NULL) %>% dplyr::select(pred_bm2 = pred_bm, var2 = variance), 
  over3.sf %>% st_set_geometry(NULL) %>% dplyr::select(pred_bm3 = pred_bm, var3 = variance)
) %>% 
  st_sf(geometry = base.geom, agr = "constant")

over.sf.rescaled <- ensemble_rescale(
  over.sf, c("pred_bm1", "pred_bm2", "pred_bm3"), "abundance", 1648, 
  x.var.idx = c("var1", "var2", "var3")
)

# Check that overlaid predictions predict expected abundance
eSDM::model_abundance(over.sf.rescaled, "pred_bm1")
eSDM::model_abundance(over.sf.rescaled, "pred_bm2")
eSDM::model_abundance(over.sf.rescaled, "pred_bm3")

summary(over.sf.rescaled)

## ------------------------------------------------------------------------
# Calculate ensemble weights
e.weights <- list(
  eSDM::evaluation_metrics(over1.sf, 1, valid.data, "pres_abs"), 
  eSDM::evaluation_metrics(over2.sf, 1, valid.data, "pres_abs"), 
  eSDM::evaluation_metrics(over3.sf, 1, valid.data, "pres_abs")
)

over.df.resc.var <- over.sf.rescaled %>% 
  dplyr::select(var1, var2, var3) %>% 
  st_set_geometry(NULL)

e.weights.unw <- c(1, 1, 1) / 3
e.weights.auc <- sapply(e.weights, function(i) i[1]) / sum(sapply(e.weights, function(i) i[1]))
e.weights.tss <- sapply(e.weights, function(i) i[2]) / sum(sapply(e.weights, function(i) i[2]))
e.weights.var <- data.frame(t(apply(
  1 / over.df.resc.var, 1, function(i) {i / sum(i, na.rm = TRUE)}
)))

e.weights.unw
e.weights.auc
e.weights.tss
head(e.weights.var)

## ------------------------------------------------------------------------
### Create ensembles

# Unweighted; calculate CV because it is used in Fig. 4 plot
ens.sf.unw <- eSDM::ensemble_create(
  over.sf.rescaled, c("pred_bm1", "pred_bm2", "pred_bm3"),  w = e.weights.unw, 
  x.var.idx = NULL
) %>% 
  mutate(SE = sqrt(Var_ens), CV = SE / Pred_ens) %>% 
  dplyr::select(Pred_ens, SE, CV) %>% 
  st_set_agr("constant")

# Weights based on AUC
ens.sf.wauc <- eSDM::ensemble_create(
  over.sf.rescaled, c("pred_bm1", "pred_bm2", "pred_bm3"),  w = e.weights.auc, 
  x.var.idx = NULL
) %>% 
  mutate(SE = sqrt(Var_ens)) %>% 
  dplyr::select(Pred_ens, SE) %>% 
  st_set_agr("constant")

# Weights based on TSS
ens.sf.wtss <- eSDM::ensemble_create(
  over.sf.rescaled, c("pred_bm1", "pred_bm2", "pred_bm3"),  w = e.weights.tss, 
  x.var.idx = NULL
) %>% 
  mutate(SE = sqrt(Var_ens)) %>% 
  dplyr::select(Pred_ens, SE) %>% 
  st_set_agr("constant")

# Weights based on the inverse of the variance
ens.sf.wvar <- eSDM::ensemble_create(
  over.sf.rescaled, c("pred_bm1", "pred_bm2", "pred_bm3"),  w = e.weights.var, 
  x.var.idx = NULL
) %>% 
  mutate(SE = sqrt(Var_ens)) %>% 
  dplyr::select(Pred_ens, SE) %>% 
  st_set_agr("constant")

## ---- eval=FALSE---------------------------------------------------------
#  # Create an ensemble and calculate within-model uncertainty; code not run
#  ens.sf.unw.wmv <- eSDM::ensemble_create(
#    over.sf.rescaled, c("pred_bm1", "pred_bm2", "pred_bm3"),  w = e.weights.unw,
#    x.var.idx = c(var1, var2, var3)
#  ) %>%
#    mutate(SE = sqrt(Var_ens)) %>%
#    dplyr::select(Pred_ens , SE)

## ---- eval=FALSE---------------------------------------------------------
#  # Calculate evaluation metrics for ensembles; code not run
#  names.2 <- c(
#    "Ensemble – unweighted", "Ensemble – AUC-based weights",
#    "Ensemble – TSS-based weights", "Ensemble – variance-based weights"
#  )
#  
#  eval.lt.ens <- data.frame(do.call(rbind, list(
#    eSDM::evaluation_metrics(ens.sf.unw,  "Pred_ens", valid.data.lt, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wauc, "Pred_ens", valid.data.lt, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wtss, "Pred_ens", valid.data.lt, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wvar, "Pred_ens", valid.data.lt, "pres_abs")
#  ))) %>%
#    mutate(Preds = names.2) %>%
#    dplyr::select(Preds, AUC_LT = X1, TSS_LT = X2)
#  
#  eval.hr.ens <- data.frame(do.call(rbind, list(
#    eSDM::evaluation_metrics(ens.sf.unw,  "Pred_ens", valid.data.hr, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wauc, "Pred_ens", valid.data.hr, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wtss, "Pred_ens", valid.data.hr, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wvar, "Pred_ens", valid.data.hr, "pres_abs")
#  ))) %>%
#    mutate(Preds = names.2) %>%
#    dplyr::select(Preds, AUC_HR = X1, TSS_HR = X2)
#  
#  eval.combo.ens <- data.frame(do.call(rbind, list(
#    eSDM::evaluation_metrics(ens.sf.unw,  "Pred_ens", valid.data, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wauc, "Pred_ens", valid.data, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wtss, "Pred_ens", valid.data, "pres_abs"),
#    eSDM::evaluation_metrics(ens.sf.wvar, "Pred_ens", valid.data, "pres_abs")
#  ))) %>%
#    mutate(Preds = names.2) %>%
#    dplyr::select(Preds, AUC = X1, TSS = X2)

## ------------------------------------------------------------------------
read.csv(system.file("extdata/Table3.csv", package = "eSDM")) %>%
  dplyr::select(Predictions, AUC, TSS, `AUC-LT` = AUC.LT, `TSS-LT` = TSS.LT, 
         `AUC-HR` = AUC.HR, `TSS-HR` = TSS.HR) %>% 
  knitr::kable(caption = "Evaluation metrics", digits = 3, align = "lcccccc")

## ---- fig.width=7--------------------------------------------------------
# Simple code to visualize ensemble created with weights based on TSS values
plot_sf_3panel(
  rename(ens.sf.wtss, se = SE), "Pred_ens", main.txt = "Ensemble-TSS - ", 
  map.base = map.world
)

## ---- fig.width=7, fig.height=4, eval=FALSE------------------------------
#  ### Figure 4; code not run
#  library(tmap)
#  
#  # Values passed to tmap_sdm - range of map
#  range.poly <- st_sfc(
#    st_polygon(list(matrix(
#      c(-132, -132, -116, -116, -132, 29.5, 49, 49, 29.5, 29.5), ncol = 2
#    ))),
#    crs = 4326
#  )
#  rpoly.mat <- matrix(st_bbox(range.poly), ncol = 2)
#  
#  # Values passed to tmap_sdm - size of text labels and legend width
#  main.size <- 0.8
#  leg.size  <- 0.55
#  leg.width <- 0.43
#  grid.size <- 0.55
#  
#  # Values passed to tmap_sdm - color scale info
#  blp1 <- tmap_sdm_help(ens.sf.unw, "Pred_ens")
#  blp2 <- tmap_sdm_help(ens.sf.unw, "CV")
#  
#  # Plot of predictions (whales / km^-2)
#  tmap.obj1 <- tmap_sdm(
#    ens.sf.unw, "Pred_ens", blp1, map.world, rpoly.mat,
#    "Unweighted ensemble - predictions",
#    main.size, leg.size, leg.width, grid.size
#  )
#  # Plot of SE values (with same color sceme as predictions)
#  tmap.obj2 <- tmap_sdm(
#    ens.sf.unw, "SE", blp1, map.world, rpoly.mat,
#    "Unweighted ensemble - SE",
#    main.size, leg.size, leg.width, grid.size
#  )
#  # Plot of CV values
#  tmap.obj3 <- tmap_sdm(
#    ens.sf.unw, "CV", blp2, map.world, rpoly.mat,
#    "Unweighted ensemble - CV",
#    main.size, leg.size, leg.width, grid.size
#  )
#  
#  # Generate plot
#  tmap_arrange(
#    list(tmap.obj1, tmap.obj2, tmap.obj3), ncol = 3, asp = NULL, outer.margins = 0.05
#  )

## ---- fig.height=9, fig.width=5.7, eval=FALSE----------------------------
#  ### Figure 5; code not run
#  
#  # Values passed to tmap_sdm - size of text labels and legend width
#  main.size <- 1.1
#  leg.size  <- 0.7
#  leg.width <- 0.6
#  grid.size <- 0.7
#  
#  # Values passed to tmap_sdm - color scale info
#  blp1b <- tmap_sdm_help(ens.sf.wtss, "Pred_ens")
#  blp2b <- tmap_sdm_help_perc(ens.sf.wtss, "Pred_ens")
#  
#  # Plot of predictions (whales / km^-2)
#  tmap.obj1 <- tmap_sdm(
#    ens.sf.wtss, "Pred_ens", blp1, map.world, rpoly.mat, "Ensemble-TSS - Predictions",
#    main.size, leg.size, leg.width, grid.size
#  )
#  # Plot of SE values (with same color sceme as predictions)
#  tmap.obj2 <- tmap_sdm(
#    ens.sf.wtss, "SE", blp1, map.world, rpoly.mat, "Ensemble-TSS - SE",
#    main.size, leg.size, leg.width, grid.size
#  )
#  # Plot of predictions (percentiles)
#  tmap.obj3 <- tmap_sdm(
#    ens.sf.wtss, "Pred_ens", blp2b, map.world, rpoly.mat, "Ensemble-TSS - Predictions",
#    main.size, leg.size, leg.width, grid.size
#  )
#  # Plot of predictions (percentiles) with combined validation data presence points
#  tmap.obj4 <- tmap_sdm(
#    ens.sf.wtss, "Pred_ens", blp2b, map.world, rpoly.mat, "Ensemble-TSS - Predictions",
#    main.size, leg.size, leg.width, grid.size
#  ) +
#    tm_shape(filter(valid.data, pres_abs == 1)) +
#    tm_dots(col = "black", size = 0.04, shape = 19)
#  
#  # Generate plot
#  tmap_arrange(
#    list(tmap.obj1, tmap.obj2, tmap.obj3, tmap.obj4), ncol = 2, nrow = 2,
#    asp = NULL, outer.margins = 0.05
#  )

