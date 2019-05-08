## Kimberly Whale Project MaxEnt Occupancy modelling

sapply(c("rJava",
         "dismo",
         "ENMeval",
         "raster",
         "sf",
         "tidyverse",
         "lubridate",
         "adehabitatHR",
         "mapview",
         "leaflet",
         "foreach",
         "doParallel"), 
       require, character.only = TRUE)


## Define common CRS for script
ll<-CRS("+init=epsg:4326")
utm<-CRS("+init=epsg:3577")

##----------------------------------------------------------------
## Define and create model extent polygon
##----------------------------------------------------------------
shelf <- 
  st_read("Data/GIS/NWShelf.shp") %>% 
  as_Spatial() %>%
  spTransform(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

poly_points <- 
  tibble(
    lon = c(118.30353, 118.60977, 121.81641, 122.39319, 126.13129, 126.47461, 122.33826),
    lat = c(-17.49608, -16.99244, -13.12763, -12.79305, -14.04600, -15.39808, -18.96084),
    id = "id") 

poly <-
  poly_points %>%
  st_as_sf(coords=c("lon", "lat"), crs = 4326) %>%
  as_Spatial %>%
  mcp(., percent = 100) %>%
  rgeos::gIntersection(., shelf)

projection(poly)

# mpa1 <- read_sf("Data/GIS/Kimberley MPAs/LCSMP.shp")
mpa2 <- read_sf("Data/GIS/Kimberley MPAs/LHFMP+NLMP.shp")
mpa3 <- read_sf("Data/GIS/Kimberley MPAs/NKMP.shp")
mpa4 <- read_sf("Data/GIS/Kimberley MPAs/YNRBMP.shp")
mpa5 <- st_crop(read_sf("Data/GIS/Kimberley MPAs/capad_marine.shp"), 
                xmin=extent(poly)[1], xmax=extent(poly)[2], 
                ymin=extent(poly)[3], ymax=extent(poly)[4])
  

##----------------------------------------------------------------
## Subset detections from satellite tag data -> 
## Only select one random detection per hour ->
## Remerge it with the other survey data to form occ_master dataset

## Function to select 'n' number of random detections per 'step' unit
fun<-function(x, n = 1, step = "day"){
  f<-function(i){ i[round(runif(n,1,nrow(i)),0),c(1:4)]}
  
  dat <- 
    x %>%
    mutate(datehour = as.character(round_date(Date.Time, unit= step))) %>%
    group_by(datehour) %>%
    do(f(.))
  
  return(dat)
}

occ_argos <- 
  read_csv("Data/2018-12-28_FullData/ArgosData_2008.csv") %>%
  filter(`Loc. quality` %in% 2:3) %>%
  transmute(ID = `Platform ID`,
            Date.Time = lubridate::ymd_hms(`Loc. date`),
            Latitude = Latitude,
            Longitude = Longitude) %>%
  bind_rows(., 
            read_csv("Data/2018-12-28_FullData/2011 Exmouth humpbacks filtered argosdata.csv") %>%
              filter(quality %in% 2:3) %>%
              transmute(ID = ptt,
                        Date.Time = lubridate::ymd_hms(times),
                        Latitude = lat,
                        Longitude = long)) %>%
  group_by(ID) %>%
  do(fun(., n = 1, step="hour")) %>%
  transmute(date = as_date(Date.Time),
            adults = 1,
            calves = NA,
            source = "satTag",
            id = as.character(ID),
            lon = Longitude,
            lat = Latitude,
            year = lubridate::year(date),
            month = lubridate::month(date),
            yearmon = paste(year,month, sep="-"))

occ_master<- 
  read_csv("Data/2018-12-28_FullData/2018-12-28_occ_master.csv") %>%
  filter(!source %in% "argos") %>%
  mutate(year = lubridate::year(date),
         month = lubridate::month(date),
         yearmon = paste(year,month, sep="-")) %>%
  bind_rows(., occ_argos) %>%
  st_as_sf(coords=c("lon","lat"), crs = 4326) %>%
  as_Spatial() %>%
  raster::intersect(., poly)

##----------------------------------------------------------------
## Input polygon for coastline
##----------------------------------------------------------------
map <- 
  st_read("Data/GIS/mainland.shp") %>% as_Spatial() %>%
  spTransform(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
  crop(., 
       extent(occ_master %>% 
                st_as_sf(coords=c("lon","lat"), crs=4326) %>% 
                as_Spatial()) + 10) 

map2 <- 
  st_read("Data/GIS/WACoast.shp") %>% as_Spatial() %>%
  spTransform(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
  crop(., 
       extent(occ_master %>% 
                st_as_sf(coords=c("lon","lat"), crs=4326) %>% 
                as_Spatial()) + 10) 
            
##----------------------------------------------------------------
## Input, format and test predictor variables for collinearity
##----------------------------------------------------------------
# Geographic
along <- raster("Data/Predictor rasters/other/along.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)

# Bathy
bath <- raster("Data/Predictor rasters/other/bathymetry.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
slope <- raster("Data/Predictor rasters/other/slope.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)

# Habitat
# DistToFW <- raster("Data/Predictor rasters/other/DistToFW.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
DistToLand <- raster("Data/Predictor rasters/other/DistToLand.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)

# Turbidity
k490_0607 <- raster("Data/Predictor rasters/Turb/k490_0607.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
k490_0800 <- raster("Data/Predictor rasters/Turb/k490_0800.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
k490_0910 <- raster("Data/Predictor rasters/Turb/k490_0910.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
k490_mean <- mean(k490_0607, k490_0800, k490_0910, na.rm=T); names(k490_mean) <- "turbidity"

# Climatic variables
sal_mean <- raster("Data/Predictor rasters/other/sal_mean.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)

sst_0607 <- raster("Data/Predictor rasters/SST/sst_0607.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
sst_0800 <- raster("Data/Predictor rasters/SST/sst_0800.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
sst_0910 <- raster("Data/Predictor rasters/SST/sst_0910.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
sst_mean <- mean(sst_0607, sst_0800, sst_0910, na.rm=T); names(sst_mean) <- "sst_mean"

chl_0607 <- raster("Data/Predictor rasters/Chlor/chlor_0607.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
chl_0800 <- raster("Data/Predictor rasters/Chlor/chlor_0800.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
chl_0910 <- raster("Data/Predictor rasters/Chlor/chlor_0910.asc", crs = ll) %>% crop(., poly) %>% mask(., poly)
chl_mean <- mean(chl_0607, chl_0800, chl_0910, na.rm=T); names(chl_mean) <- "chlorophyll_a"

env_stack <- stack(along, bath, slope,
                   DistToLand, 
                   sal_mean, sst_mean, chl_mean, k490_mean)

##----------------------------------------------------------------
## Test multicolinearity between predictors
##----------------------------------------------------------------
library(corrplot)
## Helper function for spoke plot
source("R Scripts/spoke.R")

ex<-values(env_stack)
ex2<-raster::extract(env_stack, occ_master)
correl<-cor(ex2, method='pearson', use="complete.obs")

corrplot(
  correl,
  method = "shade",
  addgrid.col = "white",
  order = "FPC",
  type = "lower",
  tl.col = 1,
  diag = F,
  shade.col = NA,
  mar = c(1, 0, 1, 5)
)

spoke(
  pos = correl > 0.7,
  neg = correl < -0.7,
  lwdPos = 2,
  lwdNeg = 2,
  ltyNeg = 1,
  ltyPos = 1,
  colNeg = "red",
  colPos = "blue"
)

##----------------------------------------------------------------
## Format occurrence datasets and subset occurrences with calves
##----------------------------------------------------------------
# Regularise datasets
occ_all <-  
  occ_master %>%
  raster::rasterize(., bath, field = 1) %>%
  raster::rasterToPoints(., spatial = T)

occ_all_cal <- 
  occ_master[!is.na(occ_master$calves),] %>%
  raster::rasterize(., bath, field = 1) %>%
  raster::rasterToPoints(., spatial = T)

occ_junjul <-  
  occ_master[occ_master$month %in% c(6,7),] %>%
  raster::rasterize(., bath, field = 1) %>%
  raster::rasterToPoints(., spatial = T)

occ_junjul_cal <-  
  occ_master[occ_master$month %in% c(6,7) & !is.na(occ_master$calves),] %>%
  raster::rasterize(., bath, field = 1) %>%
  raster::rasterToPoints(., spatial = T)

occ_aug <-  
  occ_master[occ_master$month %in% c(8),] %>%
  raster::rasterize(., bath, field = 1) %>%
  raster::rasterToPoints(., spatial = T)

occ_aug_cal <-  
  occ_master[occ_master$month %in% c(8) & !is.na(occ_master$calves),] %>%
  raster::rasterize(., bath, field = 1) %>%
  raster::rasterToPoints(., spatial = T)

occ_sepoct <-  
  occ_master[occ_master$month %in% c(9:11),] %>%
  raster::rasterize(., bath, field = 1) %>%
  raster::rasterToPoints(., spatial = T)

occ_sepoct_cal <-  
  occ_master[occ_master$month %in% c(9:11) & !is.na(occ_master$calves),] %>%
  raster::rasterize(., bath, field = 1) %>%
  raster::rasterToPoints(., spatial = T)

##############################################################
## Generate Bias grid and absence points
##############################################################
bias <-
  adehabitatHR::kernelUD(occ_all, grid = SpatialPixels(SpatialPoints(coordinates(bath))[!is.na(values(bath)),]))[[1]] %>%
  adehabitatHR::getvolumeUD(.) %>%
  raster(.) %>%
  mask(., poly)

values(bias) <- 100 - values(bias)

biaspts<-
  dismo::randomPoints(mask = bias, n = 10000, p = occ_all, prob=TRUE) %>%
  as_tibble() %>%
  mutate(id = "bg") %>%
  st_as_sf(coords=c("x","y"), crs = 4326) %>%
  as_Spatial()


##############################################################
## MaxEnt Model
##############################################################

### Standard MAXENT with automated settings
### Initial maxent model with all variables
occ<-occ_all
env<-stack(bath,
           slope,
           DistToLand,
           sal_mean, sst_mean, chl_mean, k490_mean)
bg<-biaspts
me<-maxent(x=env, p=occ, a=bg)#, nbg=nrow(occ))
pred<-predict(me, env, progress="text", type="cloglog")
mod.eval<-evaluate(me, p=occ, a=bg, x=env)#; mod.eval@auc
TSS.me<-mod.eval@t[which.max((mod.eval@TPR+mod.eval@TNR)-1)];
kappa.me<-mod.eval@t[which.max(mod.eval@kappa)];
LPT.me<-min(raster::extract(pred, occ), na.rm=T)
# mapview(pred, na.color = NA) +
#   # mapview(pred>LPT.me, na.color = NA) +
#   mapview(pred>kappa.me, na.color = NA) +
#   mapview(pred>TSS.me, na.color = NA) +
#   mapview(occ, color = "white", cex =1) +
#   mapview(bg, color = "red", cex =1)

##----------------------------------------------------------------
## Full dataset
# mme_all <- maxent(x = stack(bath, slope,
#                             DistToLand,
#                             sal_mean, sst_mean, chl_mean, k490_mean),
#                   p = occ_all, a = biaspts)
me_all <-
  ENMevaluate(
    occ = coordinates(occ_all),
    env = stack(bath, slope,
                DistToLand,
                sal_mean, sst_mean, chl_mean, k490_mean),
    bg.coords = coordinates(biaspts),
    method = "randomkfold",
    kfolds = 5,
    algorithm='maxent.jar',
    bin.output = T
  )
attributes(me_all)$env <- stack(bath, slope,
                                DistToLand,
                                sal_mean, sst_mean, chl_mean, k490_mean)

saveRDS(me_all, file="Outputs/1. All data/me_all.rds")

# mme_all_cal <- maxent(x = stack(bath, slope,
#                             DistToLand,
#                             sal_mean, sst_mean, chl_mean, k490_mean),
#                   p = occ_all_cal, a = biaspts)
me_all_cal <-
  ENMevaluate(
    occ = coordinates(occ_all_cal),
    env = stack(bath, slope,
                DistToLand,
                sal_mean, sst_mean, chl_mean, k490_mean),
    bg.coords = coordinates(biaspts),
    method = "randomkfold",
    kfolds = 5,
    algorithm='maxent.jar',
    bin.output = T
  )
attributes(me_all_cal)$env <- stack(bath, slope,
                                DistToLand,
                                sal_mean, sst_mean, chl_mean, k490_mean)

saveRDS(me_all_cal, file="Outputs/1. All data/me_all_cal.rds")

##----------------------------------------------------------------
## June July
# mme_junjul <- maxent(x = stack(bath, slope,
#                             DistToLand,
#                             sal_mean, sst_0607, chl_0607, k490_0607),
#                   p = occ_junjul, a = biaspts)
me_junjul <-
  ENMevaluate(
    occ = coordinates(occ_junjul),
    env = stack(bath, slope,
                DistToLand,
                sal_mean, sst_0607, chl_0607, k490_0607),
    bg.coords = coordinates(biaspts),
    method = "randomkfold",
    kfolds = 5,
    algorithm='maxent.jar',
    bin.output = T
  )
attributes(me_junjul)$env <- stack(bath, slope,
                                    DistToLand,
                                    sal_mean, sst_0607, chl_0607, k490_0607)

saveRDS(me_junjul, file="Outputs/2. June-July/me_junjul.rds")

# mme_junjul_cal <- maxent(x = stack(bath, slope,
#                             DistToLand,
#                             sal_mean, sst_0607, chl_0607, k490_0607),
#                   p = occ_junjul_cal, a = biaspts)
me_junjul_cal <-
  ENMevaluate(
    occ = coordinates(occ_junjul_cal),
    env = stack(bath, slope,
                DistToLand,
                sal_mean, sst_0607, chl_0607, k490_0607),
    bg.coords = coordinates(biaspts),
    method = "randomkfold",
    kfolds = 5,
    algorithm='maxent.jar',
    bin.output = T
  )
attributes(me_junjul_cal)$env <- stack(bath, slope,
                                   DistToLand,
                                   sal_mean, sst_0607, chl_0607, k490_0607)

saveRDS(me_junjul_cal, file="Outputs/2. June-July/me_junjul_cal.rds")

##----------------------------------------------------------------
## August
# mme_aug <- maxent(x = stack(bath, slope,
#                             DistToLand,
#                             sal_mean, sst_0800, chl_0800, k490_0800),
#                   p = occ_aug, a = biaspts)
me_aug <-
  ENMevaluate(
    occ = coordinates(occ_aug),
    env = stack(bath, slope,
                DistToLand,
                sal_mean, sst_0800, chl_0800, k490_0800),
    bg.coords = coordinates(biaspts),
    method = "randomkfold",
    kfolds = 5,
    algorithm='maxent.jar',
    bin.output = T
  )
attributes(me_aug)$env <- stack(bath, slope,
                                DistToLand,
                                sal_mean, sst_0800, chl_0800, k490_0800)

saveRDS(me_aug, file="Outputs/3. August/me_aug.rds")

# mme_aug_cal <- maxent(x = stack(bath, slope,
#                             DistToLand,
#                             sal_mean, sst_0800, chl_0800, k490_0800),
#                   p = occ_aug_cal, a = biaspts)
me_aug_cal <-
  ENMevaluate(
    occ = coordinates(occ_aug_cal),
    env = stack(bath, slope,
                DistToLand,
                sal_mean, sst_0800, chl_0800, k490_0800),
    bg.coords = coordinates(biaspts),
    method = "randomkfold",
    kfolds = 5,
    algorithm='maxent.jar',
    bin.output = T
  )
attributes(me_aug_cal)$env <- stack(bath, slope,
                                DistToLand,
                                sal_mean, sst_0800, chl_0800, k490_0800)

saveRDS(me_aug_cal, file="Outputs/3. August/me_aug_cal.rds")

##----------------------------------------------------------------
## September October
# mme_sepoct <- maxent(x = stack(bath, slope,
#                             DistToLand,
#                             sal_mean, sst_0910, chl_0910, k490_0910),
#                   p = occ_sepoct, a = biaspts)
me_sepoct <-
  ENMevaluate(
    occ = coordinates(occ_sepoct),
    env = stack(bath, slope,
                DistToLand,
                sal_mean, sst_0910, chl_0910, k490_0910),
    bg.coords = coordinates(biaspts),
    method = "randomkfold",
    kfolds = 5,
    algorithm='maxent.jar',
    bin.output = T
  )
attributes(me_sepoct)$env <- stack(bath, slope,
                                   DistToLand,
                                   sal_mean, sst_0910, chl_0910, k490_0910)

saveRDS(me_sepoct, file="Outputs/4. September-October/me_sepoct.rds")

# mme_sepoct_cal <- maxent(x = stack(bath, slope,
#                             DistToLand,
#                             sal_mean, sst_0910, chl_0910, k490_0910),
#                   p = occ_sepoct_cal, a = biaspts)
me_sepoct_cal <-
  ENMevaluate(
    occ = coordinates(occ_sepoct_cal),
    env = stack(bath, slope,
                DistToLand,
                sal_mean, sst_0910, chl_0910, k490_0910),
    bg.coords = coordinates(biaspts),
    method = "randomkfold",
    kfolds = 5,
    algorithm='maxent.jar',
    bin.output = T
  )
attributes(me_sepoct_cal)$env <- stack(bath, slope,
                                   DistToLand,
                                   sal_mean, sst_0910, chl_0910, k490_0910)

saveRDS(me_sepoct_cal, file="Outputs/4. September-October/me_sepoct_cal.rds")

##----------------------------------------------------------------
## Evaluation and post process
##----------------------------------------------------------------
## Helper functions
post <- function(eval, plot=TRUE){
  opt.res<-eval@results[which(eval@results$delta.AICc==0),]
  opt.mod<-eval@models[[which(eval@results$delta.AICc==0)[1]]]
  opt.pred<-eval@predictions[[which(eval@results$delta.AICc==0)[1]]]
  
  if(plot){
    eval.plot(eval@results)
    points(opt.res$rm, opt.res$delta.AICc, cex=2, col=2)
  }
  
  clogplot<-opt.pred;
  values(clogplot)<- 1-(exp(values(opt.pred)*(-(exp(opt.mod@results["Entropy",])))))
  
  varimp<-var.importance(opt.mod)
  
  mod.eval<-evaluate(opt.mod, p=eval@occ.pts[eval@occ.grp!=1,], a=eval@bg.pts, x=eval@env)
  TSS.me<-mod.eval@t[which.max((mod.eval@TPR+mod.eval@TNR)-1)];
  kappa.me<-mod.eval@t[which.max(mod.eval@kappa)];
  LPT.me<-min(raster::extract(clogplot, eval@occ.pts), na.rm=T)
  
  out<- 
    list(
      occ = eval@occ.pts,
      bg = eval@bg.pts,
      env = eval@env,
      opt.result = opt.res,
      opt.model = opt.mod,
      opt.plot = clogplot,
      var.imp = varimp,
      kappa = kappa.me,
      TSS = TSS.me,
      LPT = LPT.me
    )
  
  return(out)
}
resptab<-function(eval, xrange="p"){
  cl<-makeCluster(detectCores())
  registerDoParallel(cl)
  res<-foreach(n=1:max(eval@occ.grp), .combine=function(x, ...){mapply(rbind,x,...,SIMPLIFY=FALSE)}, .multicombine=T) %dopar% {
    require(dismo)
    occ_k<-eval@occ.pts[eval@occ.grp!=n,]
    bg_k<-eval@bg.pts
    env_k = eval@env
    me_k<-maxent(x=env_k, p=occ_k, bg=bg_k)
    out<-list()
    for(l in 1:length(names(env_k))){
      nam<-names(env_k)[l]
      out[[nam]]<-data.frame(kfold=n, response(me_k, range=xrange, rug=T, ylim=c(0,1), las=1, var=nam, xlim=range(values(env_k[[nam]]), na.rm=T)))
    }
    return(out)
  }
  stopCluster(cl)
  return(
    list(
      occ = eval@occ.pts,
      bg = eval@bg.pts,
      env = eval@env,
      res = res)
  )
}
plotresp <- function(resp, var = "all", rug = T, xlim = c("default", "crop")) {
    # if (gg) {
    #   require(tidyverse)
    #   dat <- 
    #     mapply(cbind, resp$res, "Var"= names(resp$res), SIMPLIFY=F) %>%
    #     reduce(rbind)
    #   
    #   ex <- 
    #     raster::extract(resp$env, resp$occ) %>%
    #     as_tibble() %>%
    #     gather(Var, V1)
    #   
    #   ggplot(dat, aes(x=V1, y=p)) + 
    #     geom_smooth(formula = y ~ x) +
    #     geom_rug(data = ex, sides = "b") +
    #     facet_wrap(~Var, scales="free") +
    #     theme_classic()
    #   
    #   if (var %in% "all") {
    #     
    #     
    #   }
    # } else {
      if (var %in% "all") {
        require(fifer)
        # dev.new(noRStudioGD=TRUE)
        auto.layout(n = length(resp$res))
        par(mar = c(4.5, 4.5, 1, 1))
        for (p in 1:length(resp$res)) {
          ex <- raster::extract(resp$env[[names(resp$res)[p]]], resp$occ)
          dat <- resp$res[[names(resp$res)[p]]]
          dd <-
            data.frame(x = seq(range(dat$V1)[1], range(dat$V1)[2], l = 100))
          for (d in 1:max(dat$kfold)) {
            dd[, d + 1] <- dat[dat$kfold %in% d, "p"]
          }
          dd$y_mean <- apply(dd[, -1], 1, mean)
          dd$y_sd <- apply(dd[, c(-1, -(max(dat$kfold) + 1))], 1, sd)
          dd$y_up <- dd$y_mean + dd$y_sd
          dd$y_lo <- dd$y_mean - dd$y_sd
          if (xlim %in% "default") {
            xl = range(dd$x, na.rm = T)
          }
          if (xlim %in% "crop") {
            xl = range(ex, na.rm = T)
          }
          plot(
            y_mean ~ x,
            dd,
            type = "l",
            las = 1,
            ylim = range(c(min(dd$y_lo), max(dd$y_up), 0, 1)),
            xlim = xl,
            xlab = names(resp$res)[p],
            ylab = bquote(italic("f") ~ "(" ~ .(names(resp$res)[p]) ~ ")")
          )
          with(dd, polygon(
            c(x, rev(x)),
            c(y_lo, rev(y_up)),
            border = NA,
            col = adjustcolor(8, alpha.f = 0.5)
          ))
          lines(y_up ~ x, dd, lty = 3)
          lines(y_lo ~ x, dd, lty = 3)
          lines(y_mean ~ x, dd, col = 2)
          if (rug) {
            rug(ex, col = "blue")
          }
        }
      } else {
        # dev.new(noRStudioGD=TRUE)
        par(mfrow = c(1, 1), mar = c(4.5, 4.5, 1, 1))
        ex <- raster::extract(resp$env[[var]], resp$occ)
        dat <- resp$res[[var]]
        dd <-
          data.frame(x = seq(range(dat$V1)[1], range(dat$V1)[2], l = 100))
        for (d in 1:max(dat$kfold)) {
          dd[, d + 1] <- dat[dat$kfold %in% d, "p"]
        }
        dd$y_mean <- apply(dd[, -1], 1, mean)
        dd$y_sd <- apply(dd[, c(-1, -(max(dat$kfold) + 1))], 1, sd)
        dd$y_up <- dd$y_mean + dd$y_sd
        dd$y_lo <- dd$y_mean - dd$y_sd
        if (xlim == "default") {
          xl = range(dd$x, na.rm = T)
        }
        if (xlim == "crop") {
          xl = range(ex, na.rm = T)
        }
        plot(
          y_mean ~ x,
          dd,
          type = "l",
          las = 1,
          ylim = range(c(min(dd$y_lo), max(dd$y_up), 0, 1)),
          xlim = xl,
          xlab = var,
          ylab = bquote(italic("f") ~ "(" ~ .(var) ~ ")")
        )
        with(dd, polygon(
          c(x, rev(x)),
          c(y_lo, rev(y_up)),
          border = NA,
          col = adjustcolor(8, alpha.f = 0.5)
        ))
        lines(y_up ~ x, dd, lty = 3)
        lines(y_lo ~ x, dd, lty = 3)
        lines(y_mean ~ x, dd, col = 2)
        if (rug) {
          rug(ex, col = "blue")
        #}
      }
    }
  }
##----------------------------------------------------------------

resp_all<-resptab(me_all)
resp_all_cal<-resptab(me_all_cal)
resp_junjul<-resptab(me_junjul)
resp_junjul_cal<-resptab(me_junjul_cal)
resp_aug<-resptab(me_aug)
resp_aug_cal<-resptab(me_aug_cal)
resp_sepoct<-resptab(me_sepoct)
resp_sepoct_cal<-resptab(me_sepoct_cal)

# pdf("Outputs/resp_sepoct.pdf", width=6, height=7); plotresp(resp_sepoct, xlim="crop"); dev.off()

# me_sepoct_cal<-readRDS(file.choose())

post_all<- post(me_all)
post_all_cal<- post(me_all_cal)
post_junjul<- post(me_junjul)
post_junjul_cal<- post(me_junjul_cal)
post_aug<- post(me_aug)
post_aug_cal<- post(me_aug_cal)
post_sepoct<- post(me_sepoct)
post_sepoct_cal<- post(me_sepoct_cal)

ff<-function(postlist) {
  res <- data.frame(model = c("all","all_cal","junjul","junjul_cal","aug","aug_cal","sepoct","sepoct_cal"),
                    feature = NA,
                    rm = NA,
                    auc = NA,
                    kappa = NA,
                    TSS = NA,
                    LPT = NA)
  for(i in 1:length(postlist)){
    res[i,"feature"] <- postlist[[i]]$opt.result$features[1]
    res[i,"rm"] <- postlist[[i]]$opt.result$rm[1]
    res[i,"auc"] <- postlist[[i]]$opt.result$train.AUC[1]
    res[i,"kappa"] <- postlist[[i]]$kappa
    res[i,"TSS"] <- postlist[[i]]$TSS
    res[i,"LPT"] <- postlist[[i]]$LPT
  }
  return(res)
}

tab <- ff(list(post_all, post_all_cal,
               post_junjul, post_junjul_cal,
               post_aug, post_aug_cal,
               post_sepoct, post_sepoct_cal))
write_csv(tab, "Outputs/Eval_tab.csv")

## 
# ggplot(post_all$var.imp) +
#   geom_col(aes(y=percent.contribution, x=fct_reorder(variable, percent.contribution))) +
#   coord_flip()


##----------------------------------------------------------------
## Save model outputs as raster files (.asc)
# writeRaster(post_sepoct$opt.plot, "Outputs/5. Map outputs/clogplot_sepoct.asc")
# writeRaster(post_sepoct_cal$opt.plot, "Outputs/5. Map outputs/clogplot_sepoct_cal.asc")
plotmaps<- function(post, mpa=F, scales = T, title=NULL, mpa.col="white", threshold=NULL){
  par(mar=c(7,4,1,1))
  plot(as(poly,"SpatialLines"), bg="aliceblue", axes=T, las=1, col=NA)
  
  if(is.null(threshold)){
    plot(
      post$opt.plot,
      add = T,
      col = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Spectral")))(200),
      horizontal = TRUE,
      zlim = c(0, 1)
    )
  } else {
    plot(
      post$opt.plot > post[[threshold]],
      add = T,
      col = c("navajowhite1","firebrick"),
      horizontal = TRUE,
      legend=F,
      zlim = c(0, 1)
    )
  }
  
  if(mpa){
    plot(as(as_Spatial(mpa2),"SpatialLines"), add=T, col=mpa.col, lwd=1)
    plot(as(as_Spatial(mpa3),"SpatialLines"), add=T, col=mpa.col, lwd=1)
    plot(as(as_Spatial(mpa4),"SpatialLines"), add=T, col=mpa.col, lwd=1)
    plot(as(as_Spatial(mpa5),"SpatialLines"), add=T, col=mpa.col, lwd=1)
  }
  plot(map2, add=T, col=8, border = grey(0), lwd=0.1)
  plot(as(poly,"SpatialLines"), bg="aliceblue", axes=T, add=T, lwd=1)
  scalebar(100, type="bar", lonlat=T, label=c("0","50","100"), below="km")
  if(!is.null(threshold)){
    legend(
      'bottomright',
      fill = c("navajowhite1", "firebrick"),
      border = c("navajowhite1", "firebrick"),
      legend = c("Unsuitable Habitat", "Suitable Habitat"),
      bty = "n"
    )
  }
  if(!is.null(title)){
    legend('topleft', pch=NA, col=NA, legend=title, bty="n")
  }
  box()
}

png("Outputs/5. Map outputs/Maps/maps_aug.png", width=12, height=10, units="in", res=300)
par(mfrow=c(2,2))
plotmaps(post_aug, mpa=T,
         title = c(expression(paste(bold("(a)"), " August (Adults)"))))
plotmaps(post_aug, mpa=T, threshold="TSS",
         title = c(expression(paste(bold("(b)"), " August (Adults)"))))

plotmaps(post_aug_cal, mpa=T,
         title = c(expression(paste(bold("(c)"), " August (Mothers and calves)"))))
plotmaps(post_aug_cal, mpa=T, threshold="TSS",
         title = c(expression(paste(bold("(d)"), " August (Mothers and calves)"))))
dev.off()


##----------------------------------------------------------------
## Variable contribution plot
vartab<-post_all$var.imp[c(1:2)]; colnames(vartab)<-c("variable","all")
vartab$variable <- c("Proximity\nTo Land", "Depth", "Chlorophyll a", "Salinity", "Slope", "SST", "Turbidity")

post_junjul$var.imp[c(1:2)]$variable <- c("Proximity\nTo Land", "Depth", "Chlorophyll a", "Turbidity", "Salinity", "Slope", "SST")
post_aug$var.imp[c(1:2)]$variable <- c("Proximity\nTo Land", "Depth", "Chlorophyll a", "Turbidity", "Salinity", "Slope", "SST")
post_sepoct$var.imp[c(1:2)]$variable <- c("Proximity\nTo Land", "Depth", "Chlorophyll a", "Turbidity", "Salinity", "Slope", "SST")

vartab1<-merge(vartab, post_junjul$var.imp[c(1:2)]); colnames(vartab1)<-c("variable","all","JunJul")
vartab2<-merge(vartab1, post_aug$var.imp[c(1:2)]); colnames(vartab2)<-c("variable","all","JunJul","Aug")
vartab3<-merge(vartab2, post_sepoct$var.imp[c(1:2)]); colnames(vartab3)<-c("variable","all","JunJul","Aug","SepOct")
vartab4<-vartab3[order(rowMeans(vartab3[-1])),]
tab<-t(vartab4[c(5,4,3,2)]); colnames(tab)<-vartab4$variable

png("VarImp_all.png",width=5, height=5, units='in', res=500); 
par(mar=c(5,7,1,1))
barplot(
  as.matrix(tab),
  beside = T,
  horiz = T,
  las = 1,
  xlim = c(0, 60),
  ylim = c(0.5, 35.5),
  xlab = expression(bold("Variable contribution (%)")),
  border = "white",
  col = RColorBrewer::brewer.pal(4, "Set2"),
  font = 2,
  yaxs = "i"
)
legend(
  "bottomright",
  fill = rev(RColorBrewer::brewer.pal(4, "Set2")),
  border = rep(NA, 5),
  legend = c(
    "All months combined",
    "June & July",
    "August",
    "September & October"
  ),
  bty = "n",
  cex = 1
)
axis(2, at=seq(0.5, 35.5, by=5), tcl=-0.3, labels=F); axis(2, at=seq(0.5, 35.5, by=5), tcl=0.3, labels=F)
box(bty="l")
dev.off()

##----------------------------------------------------------------
## Leaflet output

m <-
mapview(((post_all$opt.plot > post_all$TSS) + (post_all_cal$opt.plot > post_all_cal$TSS) *2), na.color=NA, layer="All", map.type="Esri.WorldImagery", 
        col.regions= colorRampPalette(c("moccasin","firebrick","forestgreen", "forestgreen")), label=F, homebutton = F) +
  # mapview( post_all$occ %>%  st_as_sf(coords = c("LON", "LAT"), crs = 4326),
  #   color = "white", cex = 0.5, layer = "All", alpha=0.2, alpha.regions=0, homebutton = F) +
  mapview(((post_junjul$opt.plot > post_junjul$TSS) + (post_junjul_cal$opt.plot > post_junjul_cal$TSS) *2), na.color=NA, layer="June-July",
          col.regions= colorRampPalette(c("moccasin","firebrick","forestgreen", "forestgreen")), label=F, homebutton = F) +
  # mapview( post_junjul$occ %>%  st_as_sf(coords = c("LON", "LAT"), crs = 4326),
  #   color = "white", cex = 0.5, layer = "June-July", alpha=0.2, alpha.regions=0, homebutton = F) +
  mapview(((post_aug$opt.plot > post_aug$TSS) + (post_aug_cal$opt.plot > post_aug_cal$TSS) *2), na.color=NA, layer="August",
          col.regions= colorRampPalette(c("moccasin","firebrick","forestgreen", "forestgreen")), label=F, homebutton = F) +
  # mapview( post_aug$occ %>% st_as_sf(coords = c("LON", "LAT"), crs = 4326),
  #   color = "white", cex = 0.5, layer = "August", alpha=0.2, alpha.regions=0, homebutton = F) +
  mapview(((post_sepoct$opt.plot > post_sepoct$TSS) + (post_sepoct_cal$opt.plot > post_sepoct_cal$TSS) *2), na.color=NA, layer="September-October",
          col.regions= colorRampPalette(c("moccasin","firebrick","forestgreen", "forestgreen")), label=F, homebutton = F) +
  # mapview( post_sepoct$occ %>%  st_as_sf(coords = c("LON", "LAT"), crs = 4326),
  #   color = "white", cex = 0.5, layer = "September-October", alpha=0.2, alpha.regions=0, homebutton = F) +
  mapview(poly, color="white", lwd=1, alpha.regions=0, homebutton=F, label=NA, popup = NA, query.type="click") +
  mapview(mpa2, color="black", lwd=1, alpha.regions=0, homebutton=F, layer="Marine Protected Areas", label=NA, popup = NA, query.type="click") +
  mapview(mpa3, color="black", lwd=1, alpha.regions=0, homebutton=F, layer="Marine Protected Areas", label=NA, popup = NA, query.type="click") +
  mapview(mpa4, color="black", lwd=1, alpha.regions=0, homebutton=F, layer="Marine Protected Areas", label=NA, popup = NA, query.type="click") +
  mapview(mpa5, color="black", lwd=1, alpha.regions=0, homebutton=F, layer="Marine Protected Areas", label=NA, popup = NA, query.type="click")

mm<-  
  m@map %>% 
  addLayersControl(baseGroups = c("All", "June-July", "August", "September-October"),
                   overlayGroups = c("Marine Protected Areas"),
                   options = layersControlOptions(collapsed = FALSE)) %>%
  addMiniMap(tiles = providers$OpenStreetMap, toggleDisplay = TRUE,
             position = "bottomleft") %>%
  addMeasure("topleft", primaryLengthUnit = "meters", primaryAreaUnit = "sqmeters") %>%
  addLegend(position = "bottomright", colors = c("moccasin","firebrick","forestgreen"), 
            labels = c("Unsuitable habitat", "Suitable habitat", "Calving habitats"), opacity = 1) %>%
  hideGroup(c("All", "June-July", "August", "September-October", "Marine Protected Areas"))
  
htmlwidgets::saveWidget(mm, file="Kimberley_HumpbackWhale_Models.html")

#-----------------------------------------



