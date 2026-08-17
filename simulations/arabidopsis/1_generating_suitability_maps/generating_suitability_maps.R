#borrowing heavily from https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html#select

library(geodata)
library(usdm)
library(terra)
library(sf)
library(blockCV)
library("ENMeval")
library("pastclim")

#envs <- geodata::worldclim_global(var = "bio", res = 5, path = ".")
#names(envs) <- paste("bio", 1:19, sep = "")

set_data_path(
  path_to_nc = "./datasets/Beyer_High-resolutionTerrestrialClimate_2020/data",
  copy_example = FALSE
)
download_dataset(dataset = "Beyer2020")

envs <- region_slice(
  time_bp = 0,
  bio_variables = get_vars_for_dataset(dataset = "Beyer2020"),
  dataset = "Beyer2020"
)

#envs<-current_env
#names(envs) <- paste("bio", 1:19, sep = "")

#Lasky used isothermality (Bio3), minimum temperature of coldest month (Bio6), temperature annual range (Bio7), 
#mean temperature of wettest quarter (Bio8), mean temperature of the warmest quarter (Bio10), precipitation seasonality (Bio15), 
#precipitation of wettest quarter (Bio16) and precipitation of driest quarter (Bio17). 

#lasky.bios<-c("bio3","bio6","bio7","bio8","bio10","bio15","bio16","bio17") ##dropping  for now

#envs<-envs[[names(envs) %in% lasky.bios]]

# jamiemkass originally dropped these variables
#envs.vif <- usdm::vifstep(envs)
#envs.rem <- envs.vif@excluded
#envs <- envs[[!(names(envs) %in% envs.rem)]]

occs<-read.table("datasets/Yim_ClimateBiogeographyOf_2023/data/OccurrenceDataUsedToFitModels_JBI_20March2024.txt",head=TRUE)
colnames(occs)<-c("longitude","latitude")

occs.cells <- terra::extract(envs[[1]], occs, cellnumbers = TRUE, ID = FALSE)
occs.cellDups <- duplicated(occs.cells[,1])
occs <- occs[!occs.cellDups,]

plot(envs[[1]]) # main="Mean annual temperature")
 # Add points for all the occurrence points onto the raster.
   points(occs)


# We'll now experiment with a different spatial R package called sf (simple features).
# Let's make our occs into a sf object -- as the coordinate reference system (crs) for these 
# points is WGS84, a geographic crs (lat/lon) and the same as our envs rasters, we specify it 
# as the SpatRaster's crs.
occs.sf <- sf::st_as_sf(occs, coords = c("longitude","latitude"), crs = terra::crs(envs))

# Now, we project our point data to an equal-area projection, which converts our 
# degrees to meters, which is ideal for buffering (the next step). 
# We use the typical Eckert IV projection.
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)


# Buffer all occurrences by 500 km, union the polygons together 
# (for visualization), and convert back to a form that the terra package 
# can use. Finally, we reproject the buffers back to WGS84 (lat/lon).
# We choose 500 km here to avoid sampling the Caribbean islands.
occs.buf <- sf::st_buffer(occs.sf, dist = 500000) |> 
  sf::st_union() |> 
  sf::st_sf() |>
  sf::st_transform(crs = terra::crs(envs))
plot(envs[[1]], main = names(envs)[1],xlim=range(occs$longitude)+c(-10,10),ylim=range(occs$lat)+c(-10,10))
points(occs)
# To add sf objects to a plot, use add = TRUE
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)


# Crop environmental rasters to match the study extent
envs.bg <- terra::crop(envs, occs.buf)
# Next, mask the rasters to the shape of the buffers
envs.bg <- terra::mask(envs.bg, occs.buf)
plot(envs.bg[[1]], main = names(envs)[1])
points(occs)
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)


# Randomly sample 10,000 background points from one background extent raster 
# (only one per cell without replacement). 
bg <- terra::spatSample(envs.bg, size = 10000, na.rm = TRUE, 
                        values = FALSE, xy = TRUE) |> as.data.frame()
colnames(bg) <- colnames(occs)

# Notice how we have pretty good coverage (every cell).
plot(envs.bg[[1]])
points(bg, pch = 20, cex = 0.2)

cb1 <- get.checkerboard1(occs, envs.bg, bg, aggregation.factor = 25) ##G added 1 here
#evalplot.grps(pts = occs, pts.grp = cb1$occs.grp, envs = envs.bg) ##this doesnt work

# Plotting the background points shows the checkerboard pattern clearly.
#evalplot.grps(pts = bg, pts.grp = cb1$bg.grp, envs = envs.bg)

e.mx.l <- ENMevaluate(occs = occs, envs = envs, bg = bg, 
                      algorithm = 'maxnet', partitions = 'block', 
                      tune.args = list(fc = "L", rm = 1:2))

tune.args <- list(fc = c("L","LQ","LQH"), rm = 1:5)
e.mx <- ENMevaluate(occs = occs, envs = envs, bg = bg, 
                    algorithm = 'maxnet', partitions = 'block', 
                    tune.args = tune.args)

res <- eval.results(e.mx)
# Select the model with delta AICc equal to 0, or the one with the lowest AICc score.
# In practice, models with delta AICc scores less than 2 are usually considered 
# statistically equivalent.
opt.aicc <- res |> filter(res$delta.AICc == 0)
opt.aicc

##complex model
plot(eval.predictions(e.mx)[['fc.LQH_rm.1']], 
     main = 'LQH_1 prediction', col = terrain.colors(100),xlim=range(occs$longitude)+c(-10,10),ylim=range(occs$lat)+c(-10,10))
points(occs,cex=.1)

aicmod <- which(e.mx@results$AICc == min(e.mx@results$AICc))


times <- seq(0, 120000, 1000)
for (i in 1:length(times)) {
  if (times[i] < 21000) {
    climate <- region_slice(
      time_bp = -1*times[i],
      bio_variables = get_vars_for_dataset(dataset = "Beyer2020"),
      dataset = "Beyer2020"
    )
    suitability <- enm.maxnet@predict(e.mx@models[["fc.LQH_rm.1"]], climate,  list(pred.type ="cloglog", doClamp=FALSE))
    writeRaster(suitability, paste0("./datasets/Beyer_High-resolutionTerrestrialClimate_2020/data/suitability_reconstructions/", as.character(times[i]),".grd"))
  } else if (times[i] %% 2000 == 0) {
    climate <- region_slice(
      time_bp = -1*times[i],
      bio_variables = get_vars_for_dataset(dataset = "Beyer2020"),
      dataset = "Beyer2020"
    )
    suitability <- enm.maxnet@predict(e.mx@models[["fc.LQH_rm.1"]], climate,  list(pred.type ="cloglog", doClamp=FALSE))
    writeRaster(suitability, paste0("./datasets/Beyer_High-resolutionTerrestrialClimate_2020/data/suitability_reconstructions/", as.character(times[i]),".grd"))
  }
}








### past
### https://cran.r-project.org/web/packages/pastclim/vignettes/a0_pastclim_overview.html

set_data_path(".")
download_dataset(dataset = "CHELSA_trace21k_1.0_0.5m_vsi", bio_variables = c("bio03","bio06","bio07","bio08","bio10","bio15","bio16","bio17")) 
##https://rdrr.io/cran/pastclim/man/CHELSA_trace21k_1.0.html

get_vars_for_dataset(dataset = "CHELSA_trace21k_1.0_0.5m_vsi")

climate_20k <- region_slice(
  time_bp = -20000,
  bio_variables =  c("bio03","bio06","bio07","bio08","bio10","bio15","bio16","bio17"),  #c("bio06"),
  dataset = "CHELSA_trace21k_1.0_0.5m_vsi"
)


climate_20k_eu_6 <- region_slice(
  time_bp = -20000,
  bio_variables =  c("bio15"),  #c("bio06"),
  dataset = "CHELSA_trace21k_1.0_0.5m_vsi",
  ext = region_extent$Europe
)


#climate_20k_10f<- terra::aggregate(climate_20k,fact=10). ##down sample

#climate_20k <- region_slice(
#  time_bp = -20000,
#  bio_variables =  c("bio03","bio06","bio07","bio08","bio10","bio15","bio16","bio17"),  #c("bio06"),
#  dataset = "CHELSA_trace21k_1.0_0.5m_vsi"
#)

##https://github.com/jamiemkass/ENMeval/issues/112
ref.vals <- raster::extract(envs.bg, occs.buf, method = 'simple', na.rm = TRUE)

#Clamp variables
v <- raster::extract(envs, bg, method = 'simple', na.rm = TRUE)


tmp[is.na(tmp)]<- -50000

f_envs <- clamp.vars(climate_20k, v, left = NULL, right = NULL, categoricals = NULL)
#f_envs <- clamp.vars(climate_20k, ref.vals, left = NULL, right = NULL, categoricals = NULL)


future <- enm.maxnet@predict(e.mx@models[['fc.LQH_rm.1']],f_envs,  list(pred.type ="cloglog",doClamp=TRUE))

remove.nas<-function(env){
  tmp<-env; 
  #recover()
  tmp[is.na(tmp)]<- -50000
  return(as.raster(tmp))
}

no.na.envs<-lapp(envs, fun=remove.nas)

minmax(tmp)
get_time_bp_steps(dataset = "Beyer2020")
plot(past_20k_env[[1]],xlim=range(occs$longitude)+c(-10,10),ylim=range(occs$lat)+c(-10,10))
plot(ice_masks,add=TRUE)

ice_masks <- get_ice_mask(-20000,"Beyer2020")
times<--c(5e3,10e3,14e3,16e3,18e3,20e3,26e3,30e3,40e3,50e3,60e3,80e3)
ice.times<-lapply(times,get_ice_mask,"Beyer2020")


envs_times<-list()
for(i in 1:length(times)){
  envs_times[[i]] <- region_slice(
    time_bp = times[i],
    bio_variables = get_vars_for_dataset(dataset = "Beyer2020"),
    dataset = "Beyer2020"
  )
}

pdf(file="~/Downloads/past_climate.pdf")
for(i in 1:length(times)){
 plot(envs_times[[i]][[1]], main=paste(times[i]/1000," thousand years ago"),
           xlim=range(occs$longitude)+c(-10,10),ylim=range(occs$lat)+c(-10,10));plot(ice.times[[i]],add=TRUE)
  
}