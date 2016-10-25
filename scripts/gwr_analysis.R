library(GWmodel)
library(car) # for VIF

# LOAD DATA
restaurants = readShapePoints('../../data/shapefiles/restaurants/restaurants', 
                              proj4string = CRS('+init=epsg:32612')
)
colnames(restaurants@coords) = c('east', 'north')
data.coordinates.proj = coordinates(restaurants)

var.response = "log_stars"
var.predictors = c("sqrt_sntmt", "ln_rvw_ct", "revw_span", "ln_unique", "beer_wine",  "full_bar", 
                   "price_rng",  "attire", "takeout", "wait_svc", "outdr_seat", "sqrt_CBD",   
                   "sqrt_scott", "sqrt_mwext", "compr_prox")


# MAKE EVALUATION GRID

# helper functions
study.boundary.buffer = function(coordinates, buffer) {
  boundary = bbox(coordinates)
  boundary[, 'min'] =  boundary[, 'min'] - buffer
  boundary[, 'max'] =  boundary[, 'max'] + buffer
  boundary
}

study.area = function(coordinates) {
  boundary = bbox(coordinates)
  prod(boundary[, 'max'] - boundary[, 'min'])
}

make.grid = function(data.features, resolution=rep(1, 2), buffer=1) {
  
  coordinates = coordinates(data.features)
  boundary = study.boundary.buffer(coordinates, buffer)
  dimensions = boundary[, 'max'] - boundary[, 'min']
  
  SpatialGrid(GridTopology(boundary[, 'min'],
                           cellsize = resolution,
                           cells.dim = (dimensions / resolution) + 1
  ),
  proj4string = data.features@proj4string
  )  
}

# study area parameters; distances in meters
grid.buffer = 1000
grid.resolution = rep(1000, 2)

grid = make.grid(restaurants, grid.resolution, grid.buffer)


# DISTANCE MATRICES (pre-calculate for time savings)
distance.matrix.grid = gw.dist(dp.locat =  data.coordinates.proj, 
                               rp.locat = coordinates(grid)
)

distance.matrix.self = gw.dist(dp.locat = data.coordinates.proj)


## GLOBAL MODEL

model.kernel = 'exponential'
# find optimal bandwidth (apaptive), assuming full model
bandwidth.gwr = bw.gwr(global.formula, data = restaurants, 
                       dMat = distance.matrix.self, 
                       approach = 'AICc', 
                       kernel = model.kernel, 
                       adaptive = TRUE)
bandwidth.gwr 
  # gaussian: Adaptive bandwidth (number of nearest neighbours): 2099 AICc value: -11728.9
  # bisquare: Adaptive bandwidth (number of nearest neighbours): 6981 AICc value: -11725.74
  # exponential: Adaptive bandwidth (number of nearest neighbours): 687 AICc value: -11743.45 
  # tricube: Adaptive bandwidth (number of nearest neighbours): 6981 AICc value: -11725.05 
  # boxcar: Adaptive bandwidth (number of nearest neighbours): 4416 AICc value: -11725.13

# 6:22 per calibration
model.candidates = getModelList(var.response, 
                                var.predictors, 
                                restaurants, 
                                bandwidth = bandwidth.gwr,
                                kernel = model.kernel,
                                distance.matrix = distance.matrix.self
                                )

# global model summary
global.formula = as.formula(paste(var.response, ' ~ ', paste(var.predictors, collapse= "+")))
global.lm = lm(global.formula, data = restaurants@data)
summary(global.lm)

  # Call:
  #   lm(formula = global.formula, data = restaurants@data)
  # 
  # Residuals:
  #   Min       1Q   Median       3Q      Max 
  # -0.81442 -0.04175  0.00644  0.05060  0.66433 
  # 
  # Coefficients:
  #                Estimate   Std. Error t value Pr(>|t|)    
  #   (Intercept) -9.520e-01  1.791e-02 -53.165  < 2e-16 ***
  #   sqrt_sntmt   1.074e+00  7.993e-03 134.358  < 2e-16 ***
  #   ln_rvw_ct    1.842e-02  1.446e-03  12.741  < 2e-16 ***
  #   revw_span   -1.211e-05  1.423e-06  -8.506  < 2e-16 ***
  #   ln_unique    5.463e-03  1.041e-03   5.250    1.56e-07 ***
  #   beer_wine    4.134e-03  4.248e-03   0.973    0.330545    
  #   full_bar    -1.208e-02  4.226e-03  -2.857    0.004285 ** 
  #   price_rng    2.092e-03  3.058e-03   0.684    0.493932    
  #   attire      -8.286e-03  9.908e-03  -0.836    0.403026    
  #   takeout      2.005e-02  5.798e-03   3.458    0.000548 ***
  #   wait_svc    -2.440e-03  3.555e-03  -0.686    0.492590    
  #   outdr_seat  -3.294e-04  2.767e-03  -0.119    0.905214    
  #   sqrt_CBD     2.464e-05  3.790e-05   0.650    0.515566    
  #   sqrt_scott   2.685e-07  3.682e-05   0.007    0.994182    
  #   sqrt_mwext  -1.827e-05  5.273e-05  -0.346    0.729040    
  #   compr_prox  -5.623e-04  4.910e-03  -0.115    0.908837    
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.11 on 7429 degrees of freedom
  # Multiple R-squared:  0.7988,	Adjusted R-squared:  0.7984 
  # F-statistic:  1966 on 15 and 7429 DF,  p-value: < 2.2e-16 


## GWR ANALYSIS

# Full model GWR (initial analysis)
gwr.full = gwr.basic(global.formula, data = restaurants,
                        regression.points = grid,
                        dMat = distance.matrix.grid,
                        bw = bandwidth.gwr, 
                        kernel = model.kernel,
                        adaptive = TRUE,
                        F123.test = TRUE)
gwr.full

  # ***********************************************************************
  #   *                       Package   GWmodel                             *
  #   ***********************************************************************
  #   Program starts at: 2016-10-16 19:37:43 
  # Call:
  #   gwr.basic(formula = global.formula, data = restaurants, regression.points = grid, 
  #             bw = bandwidth.gwr, kernel = model.kernel, adaptive = TRUE, 
  #             dMat = distance.matrix.grid, F123.test = TRUE)
  # 
  # Dependent (y) variable:  log_stars
  # Independent variables:  sqrt_sntmt ln_rvw_ct revw_span ln_unique beer_wine full_bar price_rng attire takeout wait_svc outdr_seat sqrt_CBD sqrt_scott sqrt_mwext compr_prox
  # Number of data points: 7445
  # ***********************************************************************
  #   *                    Results of Global Regression                     *
  #   ***********************************************************************
  #   
  #   Call:
  #   lm(formula = formula, data = data)
  # 
  # Residuals:
  #   Min       1Q   Median       3Q      Max 
  # -0.81442 -0.04175  0.00644  0.05060  0.66433 
  # 
  # Coefficients:
  #                Estimate   Std. Error t value  Pr(>|t|)    
  #   (Intercept) -9.520e-01  1.791e-02 -53.165  < 2e-16 ***
  #   sqrt_sntmt   1.074e+00  7.993e-03 134.358  < 2e-16 ***
  #   ln_rvw_ct    1.842e-02  1.446e-03  12.741  < 2e-16 ***
  #   revw_span   -1.211e-05  1.423e-06  -8.506  < 2e-16 ***
  #   ln_unique    5.463e-03  1.041e-03   5.250    1.56e-07 ***
  #   beer_wine    4.134e-03  4.248e-03   0.973    0.330545    
  #   full_bar    -1.208e-02  4.226e-03  -2.857    0.004285 ** 
  #   price_rng    2.092e-03  3.058e-03   0.684    0.493932    
  #   attire      -8.286e-03  9.908e-03  -0.836    0.403026    
  #   takeout      2.005e-02  5.798e-03   3.458    0.000548 ***
  #   wait_svc    -2.440e-03  3.555e-03  -0.686    0.492590    
  #   outdr_seat  -3.294e-04  2.767e-03  -0.119    0.905214    
  #   sqrt_CBD     2.464e-05  3.790e-05   0.650    0.515566    
  #   sqrt_scott   2.685e-07  3.682e-05   0.007    0.994182    
  #   sqrt_mwext  -1.827e-05  5.273e-05  -0.346    0.729040    
  #   compr_prox  -5.623e-04  4.910e-03  -0.115    0.908837    
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.11 on 7429 degrees of freedom
  # Multiple R-squared:  0.7988,	Adjusted R-squared:  0.7984 
  # F-statistic:  1966 on 15 and 7429 DF,  p-value: < 2.2e-16
  # 
  # ***Extra Diagnostic information
  # Residual sum of squares: 89.92064
  # Sigma(hat): 0.1099147
  # AIC:  -11717.88
  # AICc:  -11717.8
  # ***********************************************************************
  #   *          Results of Geographically Weighted Regression              *
  #   ***********************************************************************
  #   
  #   *********************Model calibration information*********************
  #   Kernel function: exponential 
  # Adaptive bandwidth: 687 (number of nearest neighbours)
  # Regression points: A seperate set of regression points is used.
  # Distance metric: A distance matrix is specified for this model calibration.
  # 
  # ****************Summary of GWR coefficient estimates:******************
  #   Min.    1st Qu.     Median    3rd Qu.    Max.
  # Intercept  -1.049e+00 -9.579e-01 -9.494e-01 -9.436e-01 -0.9124
  # sqrt_sntmt  1.047e+00  1.072e+00  1.073e+00  1.076e+00  1.1100
  # ln_rvw_ct   1.443e-02  1.813e-02  1.852e-02  1.880e-02  0.0255
  # revw_span  -1.524e-05 -1.226e-05 -1.218e-05 -1.210e-05  0.0000
  # ln_unique  -8.068e-04  4.857e-03  5.772e-03  5.995e-03  0.0088
  # beer_wine  -9.473e-03  2.784e-03  3.515e-03  5.569e-03  0.0153
  # full_bar   -2.024e-02 -1.225e-02 -1.192e-02 -1.143e-02 -0.0055
  # price_rng  -1.270e-03  1.460e-03  2.008e-03  2.548e-03  0.0098
  # attire     -5.141e-02 -1.267e-02 -7.556e-03 -4.387e-03  0.0244
  # takeout     7.462e-03  1.810e-02  1.927e-02  2.086e-02  0.0302
  # wait_svc   -1.127e-02 -2.980e-03 -2.723e-03 -2.187e-03  0.0018
  # outdr_seat -6.332e-03 -1.829e-03  5.967e-05  8.544e-04  0.0033
  # sqrt_CBD   -8.998e-05 -5.804e-06  2.240e-05  5.467e-05  0.0002
  # sqrt_scott -9.585e-05 -8.677e-07  3.661e-06  9.597e-06  0.0002
  # sqrt_mwext -2.327e-04 -7.198e-05 -9.952e-06  3.152e-05  0.0002
  # compr_prox -1.728e-02 -4.034e-03  4.250e-03  1.310e-02  0.0756
  # 
  # ***********************************************************************

# collinearity checks

vif(global.lm)

# VIF > 10 indicates collinearity, however can't detect collinearity with intercept (Wheeler 2010)
  
  # sqrt_sntmt  ln_rvw_ct  revw_span  ln_unique  beer_wine   full_bar  price_rng     attire    takeout 
  # 1.454117   1.888033   1.346923   1.681777   1.248046   2.394264   1.815139   1.127545   1.174418 
  
  # wait_svc outdr_seat   sqrt_CBD sqrt_scott sqrt_mwext compr_prox
  # 1.942982   1.176472   1.766905   1.737359   1.153020   1.032996 

BKW = function(predictorVariables) {
  # calculate overall model condition number
  # adapted from Gollini et. al (2015)
  X = as.matrix(cbind(1, predictorVariables))
  p = dim(X)[2]
  Xscale = sweep(X, 2, sqrt(colSums(X^2)), "/")
  Xsvd = svd(Xscale)$d
  Xsvd[1] / Xsvd[p]
}

# full model
BKW(restaurants@data[, var.predictors])
  # 57.87858 
  # (>30 is high) -> collinearity is a problem in global model

# check BKW vs. model's global condition numbers 
# (these should all match the full-model BKW number above)
lcrm.full = gwr.lcr(global.formula, data = restaurants,
                    regression.points = grid,
                    dMat = distance.matrix.grid,
                    bw = length(restaurants), 
                    kernel = 'boxcar', 
                    adaptive = TRUE)

summary(lcrm.full$SDF$Local_CN)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 57.88   57.88   57.88   57.88   57.88   57.88 


# try removing combinations of 1-2 variables
matrix(outer(var.predictors, var.predictors, function(x, y) {
  apply(matrix(c(x, y), ncol = 2), MARGIN = 1, 
        function(z) BKW(restaurants@data[, var.predictors[!var.predictors %in% z]]))
}), ncol = length(var.predictors), dimnames = list(var.predictors, var.predictors))

  #             sqrt_sntmt ln_rvw_ct revw_span ln_unique beer_wine full_bar price_rng   attire  takeout wait_svc outdr_seat sqrt_CBD sqrt_scott sqrt_mwext compr_prox
  # sqrt_sntmt   24.60564  22.81692  23.56657  24.09543  24.20915 23.67512  20.54148 24.52232 19.91437 23.68784   23.86768 23.32316   22.89601   23.49317   24.47943
  # ln_rvw_ct    22.81692  54.87064  52.60536  47.95379  54.14669 53.72835  51.29058 54.81545 51.67914 53.15017   53.02485 52.11130   51.66238   52.57763   54.58411
  # revw_span    23.56657  52.60536  55.53480  48.57318  54.91809 54.37211  51.92837 55.48120 52.25154 53.79909   53.73247 52.79415   52.41041   53.25951   55.25130
  # ln_unique    24.09543  47.95379  48.57318  51.44959  50.82959 50.38177  47.98648 51.39969 47.73716 49.57919   49.89388 48.90705   48.67362   49.41380   51.17937
  # beer_wine    24.20915  54.14669  54.91809  50.82959  57.25144 56.02234  53.78306 57.19789 54.02841 55.57766   55.53207 54.60910   54.23903   55.06670   56.97737
  # full_bar     23.67512  53.72835  54.37211  50.38177  56.02234 56.76622  52.96248 56.71626 53.46080 55.15282   55.07347 54.05309   53.66125   54.52665   56.48597
  # price_rng    20.54148  51.29058  51.92837  47.98648  53.78306 52.96248  54.27830 54.21344 51.21910 52.54384   52.53910 51.54931   51.25334   52.02587   53.99061
  # attire       24.52232  54.81545  55.48120  51.39969  57.19789 56.71626  54.21344 57.82517 54.53718 56.15987   56.12638 55.18666   54.78426   55.64442   57.55240
  # takeout      19.91437  51.67914  52.25154  47.73716  54.02841 53.46080  51.21910 54.53718 54.64160 52.87929   52.90703 51.91992   51.52754   52.38067   54.35858
  # wait_svc     23.68784  53.15017  53.79909  49.57919  55.57766 55.15282  52.54384 56.15987 52.87929 56.21193   54.46788 53.47215   53.08737   53.95019   55.92931
  # outdr_seat   23.86768  53.02485  53.73247  49.89388  55.53207 55.07347  52.53910 56.12638 52.90703 54.46788   56.18100 53.46905   53.03642   53.93150   55.89699
  # sqrt_CBD     23.32316  52.11130  52.79415  48.90705  54.60910 54.05309  51.54931 55.18666 51.91992 53.47215   53.46905 55.24414   51.74494   52.98403   54.96379
  # sqrt_scott   22.89601  51.66238  52.41041  48.67362  54.23903 53.66125  51.25334 54.78426 51.52754 53.08737   53.03642 51.74494   54.84947   52.60759   54.58006
  # sqrt_mwext   23.49317  52.57763  53.25951  49.41380  55.06670 54.52665  52.02587 55.64442 52.38067 53.95019   53.93150 52.98403   52.60759   55.69945   55.42130
  # compr_prox   24.47943  54.58411  55.25130  51.17937  56.97737 56.48597  53.99061 57.55240 54.35858 55.92931   55.89699 54.96379   54.58006   55.42130   57.60637

# without sqrt_sntmt
var.predictors.noSqrt_sntmt = var.predictors[!var.predictors %in% 'sqrt_sntmt']
BKW(restaurants@data[, var.predictors.noSqrt_sntmt])
  # 24.60564

# with only sqrt_sntmt
BKW(restaurants@data[, c('sqrt_sntmt', 'ln_rvw_ct')])
  # 25.27438


## Basic GWR with sqrt_sntmt removed
formula.2A = as.formula(paste(var.response, ' ~ ', 
                              paste(var.predictors.noSqrt_sntmt, collapse= "+")))

bandwidth.2A = bw.gwr.lcr(formula.2A, data = restaurants, 
                          dMat = distance.matrix.self, 
                          kernel = model.kernel, 
                          adaptive = TRUE)
bandwidth.2A
  # 329


# regression at restaurant points
gwr.self.2A = gwr.basic(formula.2A, data = restaurants,
                      dMat = distance.matrix.self,
                      bw = bandwidth.2A, 
                      kernel = model.kernel,
                      adaptive = TRUE,
                      F123.test = TRUE,
                      cv = FALSE)
gwr.self.2A
  # ***********************************************************************
  #   *                       Package   GWmodel                             *
  #   ***********************************************************************
  #   Program starts at: 2016-10-21 00:57:37 
  # Call:
  #   gwr.basic(formula = localA.formula, data = restaurants, bw = bandwidth.2A, 
  #             kernel = model.kernel, adaptive = TRUE, dMat = distance.matrix.self, 
  #             F123.test = TRUE, cv = FALSE)
  # 
  # Dependent (y) variable:  log_stars
  # Independent variables:  ln_rvw_ct revw_span ln_unique beer_wine full_bar price_rng attire takeout wait_svc outdr_seat sqrt_CBD sqrt_scott sqrt_mwext compr_prox
  # Number of data points: 7445
  # ***********************************************************************
  #   *                    Results of Global Regression                     *
  #   ***********************************************************************
  #   
  #   Call:
  #   lm(formula = formula, data = data)
  # 
  # Residuals:
  #   Min       1Q   Median       3Q      Max 
  # -1.18300 -0.09813  0.02447  0.12596  0.62213 
  # 
  # Coefficients:
  #                Estimate   Std. Error t value    Pr(>|t|)    
  #   (Intercept)  1.116e+00  1.695e-02  65.831  < 2e-16 ***
  #   ln_rvw_ct    5.847e-02  2.620e-03  22.314  < 2e-16 ***
  #   revw_span   -2.294e-05  2.631e-06  -8.719  < 2e-16 ***
  #   ln_unique    6.400e-02  1.750e-03  36.571  < 2e-16 ***
  #   beer_wine    3.756e-02  7.854e-03   4.783    1.76e-06 ***
  #   full_bar    -5.108e-02  7.808e-03  -6.542    6.49e-11 ***
  #   price_rng   -1.427e-02  5.658e-03  -2.522    0.011703 *  
  #   attire       3.039e-02  1.834e-02   1.657    0.097575 .  
  #   takeout      3.847e-02  1.073e-02   3.584    0.000341 ***
  #   wait_svc    -1.744e-02  6.581e-03  -2.650    0.008070 ** 
  #   outdr_seat   3.108e-02  5.105e-03   6.088    1.20e-09 ***
  #   sqrt_CBD     1.049e-04  7.017e-05   1.495    0.134982    
  #   sqrt_scott  -2.664e-04  6.808e-05  -3.913    9.18e-05 ***
  #   sqrt_mwext   2.139e-04  9.760e-05   2.192    0.028417 *  
  #   compr_prox   8.262e-03  9.092e-03   0.909    0.363565    
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.2037 on 7430 degrees of freedom
  # Multiple R-squared:  0.3099,	Adjusted R-squared:  0.3086 
  # F-statistic: 238.3 on 14 and 7430 DF,  p-value: < 2.2e-16
  # 
  # ***Extra Diagnostic information
  # Residual sum of squares: 308.4243
  # Sigma(hat): 0.2035637
  # AIC:  -2543.556
  # AICc:  -2543.483
  # ***********************************************************************
  #   *          Results of Geographically Weighted Regression              *
  #   ***********************************************************************
  #   
  #   *********************Model calibration information*********************
  #   Kernel function: exponential 
  # Adaptive bandwidth: 329 (number of nearest neighbours)
  # Regression points: the same locations as observations are used.
  # Distance metric: A distance matrix is specified for this model calibration.
  # 
  # ****************Summary of GWR coefficient estimates:******************
  #                 Min.    1st Qu.     Median    3rd Qu.    Max.
  # Intercept   9.550e-01  1.056e+00  1.091e+00  1.136e+00  1.2300
  # ln_rvw_ct   4.347e-02  5.563e-02  5.998e-02  6.665e-02  0.0797
  # revw_span  -4.150e-05 -2.899e-05 -2.542e-05 -2.209e-05  0.0000
  # ln_unique   4.325e-02  4.962e-02  5.816e-02  6.689e-02  0.0848
  # beer_wine  -2.807e-03  1.780e-02  3.097e-02  4.031e-02  0.0616
  # full_bar   -8.408e-02 -5.491e-02 -4.912e-02 -4.290e-02 -0.0246
  # price_rng  -4.072e-02 -2.369e-02 -1.176e-02  6.788e-04  0.0230
  # attire     -1.793e-01 -7.995e-03  3.313e-02  5.350e-02  0.1364
  # takeout    -8.755e-03  3.130e-02  3.851e-02  4.958e-02  0.0888
  # wait_svc   -5.371e-02 -2.450e-02 -1.973e-02 -1.079e-02  0.0114
  # outdr_seat -9.221e-03  1.952e-02  2.502e-02  3.122e-02  0.0644
  # sqrt_CBD   -1.406e-04  6.971e-05  1.414e-04  2.483e-04  0.0008
  # sqrt_scott -6.219e-04 -3.058e-04 -2.402e-04 -1.590e-04  0.0004
  # sqrt_mwext -1.285e-03  6.133e-05  3.130e-04  4.891e-04  0.0010
  # compr_prox -1.556e-01 -2.875e-02  1.344e-02  9.134e-02  0.3010
  # ************************Diagnostic information*************************
  #   Number of data points: 7445 
  # Effective number of parameters (2trace(S) - trace(S'S)): 298.1077 
  #                                                   Effective degrees of freedom (n-2trace(S) + trace(S'S)): 7146.892 
  # AICc (GWR book, Fotheringham, et al. 2002, p. 61, eq 2.33): -2680.707 
  # AIC (GWR book, Fotheringham, et al. 2002,GWR p. 96, eq. 4.22): -2876.993 
  # Residual sum of squares: 288.9274 
  # R-square value:  0.3534961 
  # Adjusted R-square value:  0.3265257 
  # ******************F test results of GWR calibration********************
  #   ---F1 test (Leung et al. 2000)
  # F1 statistic  Numerator DF  Denominator DF  Pr(>)
  # 0.97389       7326.12599    7430            0.128
  #   ---F2 test (Leung et al. 2000)
  # F2 statistic Numerator DF Denominator DF    Pr(>)    
  # 1.659      740.354           7430           < 2.2e-16 ***
  #   ---F3 test (Leung et al. 2000)
  #               F3 statistic Numerator DF Denominator DF     Pr(>)    
  #   Intercept       2.95453   1514.58486         7326.1  < 2.2e-16 ***
  #   ln_rvw_ct       2.48505   1832.19625         7326.1  < 2.2e-16 ***
  #   revw_span       0.98474   2463.96661         7326.1    0.6776    
  #   ln_unique       8.60183   1260.88550         7326.1  < 2.2e-16 ***
  #   beer_wine       0.88462    981.69238         7326.1    0.9938    
  #   full_bar        0.52743   1292.87249         7326.1    1.0000    
  #   price_rng       2.06126    926.54008         7326.1  < 2.2e-16 ***
  #   attire          2.61541     37.48637         7326.1    2.937e-07 ***
  #   takeout         0.60427    427.84469         7326.1    1.0000    
  #   wait_svc        0.77077   1144.78756         7326.1    1.0000    
  #   outdr_seat      1.17814   2441.97756         7326.1    2.418e-07 ***
  #   sqrt_CBD        2.09078   2707.82597         7326.1  < 2.2e-16 ***
  #   sqrt_scott      1.97846   2080.92013         7326.1  < 2.2e-16 ***
  #   sqrt_mwext      2.17198   1965.74504         7326.1  < 2.2e-16 ***
  #   compr_prox      2.25748    353.68429         7326.1  < 2.2e-16 ***
  #   ---F4 test (GWR book p92)
  # F4 statistic Numerator DF Denominator DF    Pr(>)   
  # 0.93679   7146.89227           7430 0.002669 **
  #   
  #   ---Significance stars
  # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
  # ***********************************************************************
  #   Program stops at: 2016-10-22 23:19:16

vif(gwr.self.2A$lm)

# VIF > 10 indicates collinearity, however can't detect collinearity with intercept (Wheeler, 2010)

  # ln_rvw_ct  revw_span  ln_unique  beer_wine   full_bar  price_rng     attire    takeout 
  # 1.807815   1.342597   1.386973   1.243766   2.382969   1.812261   1.126593   1.173761 
  
  # wait_svc outdr_seat   sqrt_CBD sqrt_scott sqrt_mwext compr_prox 
  # 1.941067   1.168073   1.766466   1.732309   1.151781   1.032812 

# GWR at grid points (no LCR correction)
gwr.grid.2A = gwr.lcr(formula.2A, data = restaurants,
                      regression.points = grid,
                      dMat = distance.matrix.grid,
                      bw = bandwidth.localA, 
                      kernel = model.kernel,
                      adaptive = TRUE)
gwr.grid.2A
  # ***********************************************************************
  #   *                       Package   GWmodel                             *
  #   ***********************************************************************
  #   Program starts at: 2016-10-23 10:20:28 
  # Call:
  #   gwr.lcr(formula = formula.2A, data = restaurants, regression.points = grid, 
  #           bw = bandwidth.2A, kernel = model.kernel, adaptive = TRUE, 
  #           dMat = distance.matrix.grid)
  # 
  # Dependent (y) variable:  log_stars
  # Independent variables:  ln_rvw_ct revw_span ln_unique beer_wine full_bar price_rng attire takeout wait_svc outdr_seat sqrt_CBD sqrt_scott sqrt_mwext compr_prox
  # ***********************************************************************
  #   *        Results of Ridge Geographically Weighted Regression           *
  #   ***********************************************************************
  #   
  #   *********************Model calibration information*********************
  #   Kernel function: exponential 
  # Adaptive bandwidth: 329 (number of nearest neighbours)
  # Regression points: A seperate set of regression points is used.
  # Distance metric: A distance matrix is specified for this model calibration.
  # Lambda(ridge parameter for gwr ridge model): 0
  # **************Summary of Ridge GWR coefficient estimates:***************
  #                Min.    1st Qu.     Median    3rd Qu.    Max.
  # Intercept   9.611e-01  1.105e+00  1.119e+00  1.134e+00  1.2300
  # ln_rvw_ct   4.366e-02  5.796e-02  5.842e-02  5.903e-02  0.0801
  # revw_span  -4.022e-05 -2.349e-05 -2.319e-05 -2.261e-05  0.0000
  # ln_unique   4.316e-02  6.216e-02  6.445e-02  6.596e-02  0.0848
  # beer_wine  -1.843e-03  3.523e-02  3.717e-02  3.986e-02  0.0619
  # full_bar   -8.414e-02 -5.544e-02 -5.238e-02 -4.622e-02 -0.0247
  # price_rng  -4.079e-02 -1.926e-02 -1.561e-02 -9.963e-03  0.0230
  # attire     -1.658e-01  2.514e-02  3.184e-02  3.426e-02  0.1263
  # takeout    -1.067e-02  3.638e-02  3.826e-02  4.028e-02  0.0897
  # wait_svc   -5.373e-02 -1.982e-02 -1.785e-02 -1.549e-02  0.0113
  # outdr_seat -6.287e-03  2.898e-02  3.140e-02  3.322e-02  0.0631
  # sqrt_CBD   -1.404e-04  6.694e-05  9.423e-05  1.271e-04  0.0008
  # sqrt_scott -6.192e-04 -3.177e-04 -2.736e-04 -2.416e-04  0.0003
  # sqrt_mwext -1.279e-03  1.008e-04  2.068e-04  3.116e-04  0.0010
  # compr_prox -1.555e-01  2.629e-03  7.921e-03  1.779e-02  0.2939
  # ************************Diagnostic information*************************
  #   
  #   ***********************************************************************
  #   Program stops at: 2016-10-23 10:24:53

summary(gwr.grid.2A$SDF$Local_CN)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 24.32   24.79   25.21   31.27   28.63  322.50 


## Basic GWR for sqrt_sntmt + ln_rvw_ct model 

# (model must have at least 2 predictors or error is generated:
# Error in array(STATS, dims[perm]) : 'dims' cannot be of length 0)
formula.2B = as.formula(paste(var.response, ' ~ ', 
                              paste(c('sqrt_sntmt', 'ln_rvw_ct'), collapse= "+")))

bandwidth.2B = bw.gwr.lcr(formula.2B, data = restaurants, 
                          dMat = distance.matrix.self, 
                          kernel = model.kernel, 
                          adaptive = TRUE)
bandwidth.2B
  # 1199

# regression at restaurant points
gwr.self.2B = gwr.basic(formula.2B, data = restaurants,
                        dMat = distance.matrix.self,
                        bw = bandwidth.2B, 
                        kernel = model.kernel,
                        adaptive = TRUE,
                        F123.test = TRUE,
                        cv = FALSE)
gwr.self.2B
  # ***********************************************************************
  #   *                       Package   GWmodel                             *
  #   ***********************************************************************
  #   Program starts at: 2016-10-22 23:34:30 
  # Call:
  #   gwr.basic(formula = formula.2B, data = restaurants, bw = bandwidth.2B, 
  #             kernel = model.kernel, adaptive = TRUE, dMat = distance.matrix.self, 
  #             F123.test = TRUE, cv = FALSE)
  # 
  # Dependent (y) variable:  log_stars
  # Independent variables:  sqrt_sntmt ln_rvw_ct
  # Number of data points: 7445
  # ***********************************************************************
  #   *                    Results of Global Regression                     *
  #   ***********************************************************************
  #   
  #   Call:
  #   lm(formula = formula, data = data)
  # 
  # Residuals:
  #   Min       1Q   Median       3Q      Max 
  # -0.80080 -0.04184  0.00582  0.05241  0.64445 
  # 
  # Coefficients:
  #                Estimate   Std. Error t value  Pr(>|t|)    
  #   (Intercept) -1.000861   0.013172  -75.98   <2e-16 ***
  #   sqrt_sntmt   1.103784   0.007106  155.33   <2e-16 ***
  #   ln_rvw_ct    0.013582   0.001128   12.04   <2e-16 ***
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 0.111 on 7442 degrees of freedom
  # Multiple R-squared:  0.7947,	Adjusted R-squared:  0.7947 
  # F-statistic: 1.441e+04 on 2 and 7442 DF,  p-value: < 2.2e-16
  # 
  # ***Extra Diagnostic information
  # Residual sum of squares: 91.72867
  # Sigma(hat): 0.1110142
  # AIC:  -11595.67
  # AICc:  -11595.67
  # ***********************************************************************
  #   *          Results of Geographically Weighted Regression              *
  #   ***********************************************************************
  #   
  #   *********************Model calibration information*********************
  #   Kernel function: exponential 
  # Adaptive bandwidth: 1199 (number of nearest neighbours)
  # Regression points: the same locations as observations are used.
  # Distance metric: A distance matrix is specified for this model calibration.
  # 
  # ****************Summary of GWR coefficient estimates:******************
  #               Min.      1st Qu. Median  3rd Qu.    Max.
  #   Intercept  -1.01900 -1.00300 -0.98940 -0.97980 -0.9617
  #   sqrt_sntmt  1.08300  1.09300  1.10000  1.10700  1.1170
  #   ln_rvw_ct   0.01018  0.01186  0.01298  0.01389  0.0170
  # ************************Diagnostic information*************************
  #   Number of data points: 7445 
  # Effective number of parameters (2trace(S) - trace(S'S)): 24.04565 
  #    Effective degrees of freedom (n-2trace(S) + trace(S'S)): 7420.954 
  # AICc (GWR book, Fotheringham, et al. 2002, p. 61, eq 2.33): -11612.58 
  # AIC (GWR book, Fotheringham, et al. 2002,GWR p. 96, eq. 4.22): -11629.49 
  # Residual sum of squares: 91.22913 
  # R-square value:  0.7958657 
  # Adjusted R-square value:  0.7952042 
  # ******************F test results of GWR calibration********************
  #   ---F1 test (Leung et al. 2000)
  # F1 statistic Numerator DF   Denominator DF  Pr(>)
  # 0.99737      7436.38928     7442            0.4549
  #   ---F2 test (Leung et al. 2000)
  # F2 statistic  Numerator DF   Denominator DF   Pr(>)    
  # 1.9257        78.4935        7442             1.936e-06 ***
  #   ---F3 test (Leung et al. 2000)
  #                 F3 statistic  Numerator  DF Denominator DF   Pr(>)    
  #   Intercept     1.3717        1294.8838  7436.4              6.371e-15 ***
  #   sqrt_sntmt    1.6515        1522.2366  7436.4              < 2.2e-16 ***
  #   ln_rvw_ct     2.2347        2793.8673  7436.4              < 2.2e-16 ***
  #   ---F4 test (GWR book p92)
  # F4 statistic    Numerator DF    Denominator DF  Pr(>)
  # 0.99455         7420.95435      7442            0.407
  # 
  # ---Significance stars
  # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
  # ***********************************************************************
  #   Program stops at: 2016-10-23 06:19:17 

vif(gwr.self.2B$lm)

# VIF > 10 indicates collinearity, however can't detect collinearity with intercept (Wheeler, 2010)

  # sqrt_sntmt  ln_rvw_ct 
  # 1.128554    1.128554  

gwr.grid.2B = gwr.lcr(formula.2B, data = restaurants,
                      regression.points = grid,
                      dMat = distance.matrix.grid,
                      bw = bandwidth.2B, 
                      kernel = model.kernel,
                      adaptive = TRUE)
gwr.grid.2B
   # ***********************************************************************
   # *                       Package   GWmodel                             *
   # ***********************************************************************
   # Program starts at: 2016-10-17 19:15:52 
   # Call:
   # gwr.lcr(formula = localB.formula, data = restaurants, regression.points = grid, 
   #  bw = bandwidth.2B, kernel = model.kernel, adaptive = TRUE, 
   #  dMat = distance.matrix.grid)
   # 
   # Dependent (y) variable:  log_stars
   # Independent variables:  sqrt_sntmt ln_rvw_ct
   # ***********************************************************************
   # *        Results of Ridge Geographically Weighted Regression           *
   # ***********************************************************************
   # 
   # *********************Model calibration information*********************
   # Kernel function: exponential 
   # Adaptive bandwidth: 1199 (number of nearest neighbours)
   # Regression points: A seperate set of regression points is used.
   # Distance metric: A distance matrix is specified for this model calibration.
   # Lambda(ridge parameter for gwr ridge model): 0
   # **************Summary of Ridge GWR coefficient estimates:***************
   #                Min.  1st Qu.   Median  3rd Qu.    Max.
   # Intercept  -1.01900 -1.00300 -1.00100 -0.99780 -0.9623
   # sqrt_sntmt  1.08300  1.10300  1.10400  1.10500  1.1170
   # ln_rvw_ct   0.01021  0.01330  0.01352  0.01371  0.0170
   # ************************Diagnostic information*************************
   # 
   # ***********************************************************************

summary(gwr.grid.2B$SDF$Local_CN)
  #  Min. 1st Qu.  Median    Mean   3rd Qu.    Max. 
  # 21.48   24.46   24.88   25.13   25.74     30.69 


## LCR GW models to correct high condition numbers
cn.threshold = 30

# Model 2A (full model minus sqrt_sntmt)
lcrm.2A = gwr.lcr(formula.2A, data = restaurants,
                  regression.points = grid,
                  dMat = distance.matrix.grid,
                  bw = bandwidth.2A, 
                  kernel = model.kernel,
                  adaptive = TRUE,
                  lambda.adjust = TRUE,
                  cn.thresh = cn.threshold)

lcrm.2A
   # ***********************************************************************
   # *                       Package   GWmodel                             *
   # ***********************************************************************
   # Program starts at: 2016-10-17 20:31:43 
   # Call:
   # gwr.lcr(formula = localA.formula, data = restaurants, regression.points = grid, 
   #  bw = bandwidth.2A, kernel = model.kernel, lambda.adjust = TRUE, 
   #  cn.thresh = cn.threshold, adaptive = TRUE, dMat = distance.matrix.grid)
   # 
   # Dependent (y) variable:  log_stars
   # Independent variables:  ln_rvw_ct revw_span ln_unique beer_wine full_bar price_rng attire takeout wait_svc outdr_seat sqrt_CBD sqrt_scott sqrt_mwext compr_prox
   # ***********************************************************************
   # *        Results of Ridge Geographically Weighted Regression           *
   # ***********************************************************************
   # 
   # *********************Model calibration information*********************
   # Kernel function: exponential 
   # Adaptive bandwidth: 329 (number of nearest neighbours)
   # Regression points: A seperate set of regression points is used.
   # Distance metric: A distance matrix is specified for this model calibration.
   # Lambda(ridge parameter for gwr ridge model): 0   The threshold of condition number is used for adjusting local Lambda: 30
   # **************Summary of Ridge GWR coefficient estimates:***************
   #                  Min.    1st Qu.     Median    3rd Qu.    Max.
   # Intercept   8.506e-01  1.104e+00  1.119e+00  1.134e+00  1.2220
   # ln_rvw_ct   4.425e-02  5.797e-02  5.842e-02  5.904e-02  0.0804
   # revw_span  -3.986e-05 -2.348e-05 -2.319e-05 -2.260e-05  0.0000
   # ln_unique   4.313e-02  6.215e-02  6.445e-02  6.596e-02  0.0847
   # beer_wine  -2.173e-03  3.522e-02  3.715e-02  3.975e-02  0.0613
   # full_bar   -8.552e-02 -5.546e-02 -5.244e-02 -4.627e-02 -0.0254
   # price_rng  -3.954e-02 -1.924e-02 -1.556e-02 -9.933e-03  0.0245
   # attire     -1.647e-01  2.519e-02  3.185e-02  3.427e-02  0.1254
   # takeout    -7.535e-03  3.651e-02  3.830e-02  4.039e-02  0.0932
   # wait_svc   -5.392e-02 -1.984e-02 -1.788e-02 -1.549e-02  0.0112
   # outdr_seat -6.555e-03  2.898e-02  3.141e-02  3.322e-02  0.0631
   # sqrt_CBD   -9.876e-05  6.740e-05  9.552e-05  1.278e-04  0.0015
   # sqrt_scott -5.992e-04 -3.174e-04 -2.730e-04 -2.405e-04  0.0006
   # sqrt_mwext -1.226e-03  1.015e-04  2.080e-04  3.124e-04  0.0012
   # compr_prox -1.550e-01  2.653e-03  7.932e-03  1.781e-02  0.2937
   # ************************Diagnostic information*************************
   # 
   # ***********************************************************************

# Model 2B (sqrt_sntmt + ln_rvw_ct model)
lcrm.2B = gwr.lcr(formula.2B, data = restaurants,
                  regression.points = grid,
                  dMat = distance.matrix.grid,
                  bw = bandwidth.2B, 
                  kernel = model.kernel,
                  adaptive = TRUE,
                  lambda.adjust = TRUE,
                  cn.thresh = cn.threshold)

lcrm.2B
   # ***********************************************************************
   # *                       Package   GWmodel                             *
   # ***********************************************************************
   # Program starts at: 2016-10-17 20:43:50 
   # Call:
   # gwr.lcr(formula = localB.formula, data = restaurants, regression.points = grid, 
   #  bw = bandwidth.2B, kernel = model.kernel, lambda.adjust = TRUE, 
   #  cn.thresh = cn.threshold, adaptive = TRUE, dMat = distance.matrix.grid)
   # 
   # Dependent (y) variable:  log_stars
   # Independent variables:  sqrt_sntmt ln_rvw_ct
   # ***********************************************************************
   # *        Results of Ridge Geographically Weighted Regression           *
   # ***********************************************************************
   # 
   # *********************Model calibration information*********************
   # Kernel function: exponential 
   # Adaptive bandwidth: 1199 (number of nearest neighbours)
   # Regression points: A seperate set of regression points is used.
   # Distance metric: A distance matrix is specified for this model calibration.
   # Lambda(ridge parameter for gwr ridge model): 0   The threshold of condition number is used for adjusting local Lambda: 30
   # **************Summary of Ridge GWR coefficient estimates:***************
   #                Min.  1st Qu.   Median  3rd Qu.    Max.
   # Intercept  -1.01900 -1.00300 -1.00100 -0.99780 -0.9623
   # sqrt_sntmt  1.08300  1.10300  1.10400  1.10500  1.1170
   # ln_rvw_ct   0.01021  0.01330  0.01352  0.01371  0.0170
   # ************************Diagnostic information*************************
   # 
   # ***********************************************************************