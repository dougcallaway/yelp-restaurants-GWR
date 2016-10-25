library(RColorBrewer)
library(dismo)
library(ggmap)
library(rgdal)
library(raster)
library(rasterVis)

## GRAPHICS HELPER FUNCTIONS

# Automated model selection procedure; adapted from Gollini, et. al (2015)
getModelList = function(response, predictors, data, bandwidth, distance.matrix=NULL, kernel='gaussian', adaptive=TRUE) {
  require(GWmodel)
  model.candidates = model.selection.gwr(response, predictors, data, 
                                         bandwidth, adaptive = adaptive, 
                                         kernel = kernel, dMat = distance.matrix
  )
  model.candidates = model.sort.gwr(model.candidates, 
                                    numVars = length(predictors), 
                                    ruler.vector = model.candidates[[2]][, 2]
  )
  model.candidates
}

# Model selection AICc plot; adapted from Gollini, et al (2015)
ModelList.plotAICc = function(modelList, main = 'Model Selection Procedure', ylab = 'AICc', xlab = 'Model Number') {
  plot(modelList[[2]][, 2], pch = 20, lty = 5, main = main, ylab = ylab, xlab = xlab, type = 'b')
}

spatialDataFrame.toGCS = function(spatialDataFrame) {
  spatialDataFrame = spTransform(spatialDataFrame, CRS('+init=epsg:4326'))
  colnames(spatialDataFrame@coords) = c('lon', 'lat')
  spatialDataFrame
}

spatialPixelsDataFrameToRaster = function(SDF, colname, 
                                          project = FALSE, 
                                          proj4string = '+init=epsg:4326') {
  r = raster(SDF, colname)
  if(project) {
    r = projectRaster(r, crs = CRS(proj4string))
  }
  r
}

overlay.levelplot = function(SDF, 
                             colname,
                             main = colname,
                             col=rev(brewer.pal(n = 10, name = 'Spectral')),
                             alpha = 0.5, 
                             at = seq(raster@data@min, 
                                      raster@data@max, 
                                      length.out = 10)) {

  raster = spatialPixelsDataFrameToRaster(SDF, colname, 
                                          project = TRUE,
                                          proj4string = '+init=epsg:4326')
  overlay = levelplot(raster, 
                      margin = FALSE, 
                      contour = TRUE,
                      pretty = TRUE,
                      cuts = 10,
                      par.settings = rasterTheme(region = col), 
                      at = at,
                      alpha.regions = alpha,
                      scales = list(cex = 1.1),
                      main = main,
                      xlab = '',
                      ylab = '')
  print(overlay + basemap.levelplot + overlay)
}

display.brewer.all()
palette  = brewer.pal(12,'Paired')

## STUDY AREA OVERVIEW
data.coordinates.gcs = coordinates(spatialDataFrame.toGCS(restaurants))

# Data Overview
basemap.ggmap = ggmap(get_map(location = rowMeans(bbox(data.coordinates.gcs)), 
                        zoom = 9, 
                        maptype = 'roadmap'))

basemap.ggmap + 
  geom_point(alpha=0.2, color = palette[10], data = as.data.frame(data.coordinates.gcs)) +
  labs(x = 'Longitude', y = 'Latitude') +
  theme(text = element_text(size = 14))

# study area grid
plot(grid)
plot(restaurants, add=TRUE, col=adjustcolor('blue', alpha.f = 0.5))


## GLOBAL MODEL RESULTS
model.view.gwr(var.response, var.predictors, model.list = model.candidates[[1]])
ModelList.plotAICc(model.candidates)


## GWR RESULTS

basemap = gmap(data.coordinates.gcs, type = 'roadmap', zoom = 9, scale = 2, lonlat = TRUE)
basemap.levelplot = levelplot(basemap, maxpixels = ncell(basemap), 
                              col.regions = basemap@legend@colortable, 
                              at = 0:255, 
                              panel = panel.levelplot.raster, 
                              interpolate = TRUE, 
                              colorkey = FALSE, 
                              margin = FALSE)

# full model, no collinearity correction
overlay.levelplot(gwr.full$SDF, 'sqrt_sntmt', col = brewer.pal(n = 9, name = 'YlOrRd'), 
                  main = "Square-Root Sentiment Score\n (with collinearity effects)")

# predictor varialbe boxplots (with collinearity correction)
## Model 2A
opar = par(mar = c(10, 8, 6, 2), las = 2)
boxplot(lcrm.2A$SDF@data[colnames(lcrm.2A$SDF@data) %in% var.predictors], 
        col = palette[2],
        pch = 20,
        outcol = setColor(palette[1], alpha = 0.1),
        main = 'Distribution of Predictor Variable Coefficients\nModel 2A',
        ylab = 'coefficient value',
        xlab = '',
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.5)
abline(0, 0, col = palette[8], lty = 'dashed')

## Model 2B
boxplot(lcrm.2B$SDF@data['sqrt_sntmt'], 
        col = palette[2],
        pch = 20,
        outcol = setColor(palette[1], alpha = 0.1),
        main = 'Distribution of Square-Root-Sentiment Coefficients\nModel 2B',
        ylab = '',
        xlab = 'sqrt_sntmt',
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 1.5)
mtext('coefficient value', 2, line = 5, las = 0, cex = 1.5)

par(opar)

# Model 2A GWR plots, condition number corrected
model = lcrm.2A

## contour plots
overlay.levelplot(model$SDF, 'Intercept', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Intercept: Model 2A")
overlay.levelplot(model$SDF, 'ln_rvw_ct', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Log Review Count")
overlay.levelplot(model$SDF, 'revw_span', col = rev(brewer.pal(n = 9, name = 'Blues')), main = "Review Span")
overlay.levelplot(model$SDF, 'ln_unique', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Log Uniqueness")
overlay.levelplot(model$SDF, 'beer_wine', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Beer & Wine Service")
overlay.levelplot(model$SDF, 'full_bar', col = rev(brewer.pal(n = 9, name = 'Blues')), main = "Full Bar Service")
overlay.levelplot(model$SDF, 'price_rng', col = rev(brewer.pal(n = 9, name = 'Spectral')), main = "Price Range")
overlay.levelplot(model$SDF, 'attire', col = rev(brewer.pal(n = 9, name = 'Spectral')), main = "Attire (dress code)")
overlay.levelplot(model$SDF, 'takeout', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Takeout Service")
overlay.levelplot(model$SDF, 'wait_svc', col = rev(brewer.pal(n = 9, name = 'YlGnBu')), main = "Waiter Service")
overlay.levelplot(model$SDF, 'outdr_seat', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Outdoor Seating")
overlay.levelplot(model$SDF, 'sqrt_CBD', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Square-Root Distance from CBD")
overlay.levelplot(model$SDF, 'sqrt_scott', col = rev(brewer.pal(n = 9, name = 'Spectral')), main = "Square-Root Distance from Scottsdale")
overlay.levelplot(model$SDF, 'sqrt_mwext', col = rev(brewer.pal(n = 9, name = 'Spectral')), main = "Square-Root Distance from Motorway Exit")
overlay.levelplot(model$SDF, 'compr_prox', col = rev(brewer.pal(n = 9, name = 'Spectral')), main = "Competitor Proximity")

## correction locations
overlay.levelplot(model$SDF, 'Local_CN', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Local Condition Number (before correction)")
overlay.levelplot(model$SDF, 'Local_Lambda', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Corrections Applied (lambda)")

# secondary model, condition number corrected
model = lcrm.2B

## contour plots
overlay.levelplot(model$SDF, 'Intercept', col = rev(brewer.pal(n = 9, name = 'Blues')), main = "Intercept: Model 2B")
overlay.levelplot(model$SDF, 'sqrt_sntmt', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Square-Root Sentiment Score")
overlay.levelplot(model$SDF, 'ln_rvw_ct', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Log Review Count")

## correction locations
overlay.levelplot(model$SDF, 'Local_CN', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Sub-Model Local Condition Numbers (before correction)")
overlay.levelplot(model$SDF, 'Local_Lambda', col = brewer.pal(n = 9, name = 'YlOrRd'), main = "Corrections Applied (lambda)")
