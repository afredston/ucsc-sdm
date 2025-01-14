library(terra) # be sure that the raster package is not loaded, since a lot of terra functions are the same names 
library(geodata)
library(here)
library(mecofun) # devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")
library(RColorBrewer)
library(lattice)
library(randomForest)
#########
# Ring Ouzel: Data preparation 
#########

avi_dat <- read.table('data/Data_SwissBreedingBirds.csv', header=T, sep=',')

summary(avi_dat)

avi_cols <- c('Turdus_torquatus', 'bio_5', 'bio_2', 'bio_14', 'std', 'rad', 'blockCV_tile')
# 5 = max temp of warmest month 
# 2 = mean diurnal range 
# 14 = precipitation of driest month

avi_df <- data.frame(avi_dat)[avi_cols]

summary(avi_df)

bio_curr <- worldclim_country("Switzerland", var="bio", path="data", download=T)
bio_fut_all <- cmip6_tile(model="CNRM-ESM2-1", ssp="370", time="2061-2080", var = "tmax", res = 10, path = here("data"), lon = 8, lat = 46.5)
# there are 12 layers, I'm assuming those are months 
bio_fut <- bio_fut_all[["wc2.1_30s_tmax_07"]]
plot(bio_fut)
names(bio_fut) <- "bio_5"

bio_curr_file <- "data/CHE_wc2.1_30s_bio.tif"
bio_curr <- subset(rast(bio_curr_file), c(2, 5, 14))
plot(bio_curr)
names(bio_curr)

if(!file.exists(here("data","CH_mask.tif"))){
bg <- rast('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif') # spatial mask of Switzerland in Swiss coordinates 
writeRaster(bg, filename=here("data","CH_mask.tif"))
} else {
  bg <- rast(here("data","CH_mask.tif"))
}

ch_ext <- c(5, 11, 45, 48)
bio_curr <- crop(bio_curr, ch_ext)
bio_fut <- crop(bio_fut, ch_ext)
plot(bio_curr)
plot(bio_fut)
# Re-project to Swiss coordinate system and clip to Swiss political boundary
bio_curr <- project(bio_curr, bg)
bio_curr <- resample(bio_curr, bg)
bio_curr <- mask(bio_curr, bg)
plot(bio_curr)
names(bio_curr) <-c("bio_2","bio_5","bio_14")

bio_fut <- project(bio_fut, bg)
bio_fut <- resample(bio_fut, bg)
bio_fut <- mask(bio_fut, bg)
plot(bio_fut)
##########
# 2.3.1 Ring Ouzel: model fitting
##########

m_glm <- glm( Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2), family='binomial', data=avi_df)
summary(m_glm)

# Names of our variables:
pred <- c('bio_2', 'bio_5', 'bio_14')

# We want three panels next to each other:
par(mfrow=c(1,3)) 

# Plot the partial responses
partial_response(m_glm, predictors = avi_df[,pred])


# We prepare the response surface by making a dummy data set where two predictor variables range from their minimum to maximum value, and the remaining predictor is kept constant at its mean:
xyz <- data.frame(expand.grid(seq(min(avi_df[,pred[1]]),max(avi_df[,pred[1]]),length=50), seq(min(avi_df[,pred[2]]),max(avi_df[,pred[2]]),length=50)), mean(avi_df[,pred[3]]))
names(xyz) <- pred

# Make predictions
xyz$z <- predict(m_glm, xyz, type='response')
summary(xyz)

# Make a colour scale
cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)

# plot 3D-surface
wireframe(z ~ bio_2 + bio_5, data = xyz, zlab = list("Occurrence prob.", rot=90),
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1), 
          main='GLM', xlab='bio_2', ylab='bio_5', screen=list(z = 120, x = -70, y = 3))

# Plot inflated response curves:
par(mfrow=c(1,3)) 
inflated_response(m_glm, predictors = avi_df[,pred], method = "stat6", lwd = 3, main='GLM') 

m_rf <- randomForest( x=avi_df[,2:4], y=avi_df[,1], ntree=1000, nodesize=10, importance =T)

importance(m_rf,type=1)

varImpPlot(m_rf)

head(getTree(m_rf,1,T))

par(mfrow=c(1,3)) 
partial_response(m_rf, predictors = avi_df[,pred], main='Random Forest')

xyz$z <- predict(m_rf, xyz)   # Note that we created the xyz data.frame in the GLM example above
wireframe(z ~ bio_2 + bio_5, data = xyz, zlab = list("Occurrence prob.", rot=90),
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1), 
          main='RF', xlab='bio_2', ylab='bio_5', screen=list(z = 120, x = -70, y = 3))


par(mfrow=c(1,3)) 
inflated_response(m_rf, predictors = avi_df[,pred], method = "stat6", lwd = 3, main='RF') 

#########
# 2.4.1 Ring Ouzel: Model assessment
#########
# Make cross-validated predictions for GLM:
crosspred_glm <- mecofun::crossvalSDM(m_glm, traindat= avi_df[!is.na(avi_df$blockCV_tile),], colname_pred=pred, colname_species = "Turdus_torquatus", kfold= avi_df[!is.na(avi_df$blockCV_tile),'blockCV_tile'])

# Make cross-validated predictions for RF:
crosspred_rf <- mecofun::crossvalSDM(m_rf, traindat= avi_df[!is.na(avi_df$blockCV_tile),], colname_pred=pred, colname_species = "Turdus_torquatus", kfold= avi_df[!is.na(avi_df$blockCV_tile),'blockCV_tile'])

# Look at correlation between GLM and RF predictions:
plot(crosspred_glm, crosspred_rf, pch=19, col='grey35')

(eval_glm <- mecofun::evalSDM(observation = avi_df[!is.na(avi_df$blockCV_tile),1], predictions = crosspred_glm))

(eval_rf <- mecofun::evalSDM(observation = avi_df[!is.na(avi_df$blockCV_tile),1], predictions = crosspred_rf))

# Derive median predictions:
crosspred_ens <- apply(data.frame(crosspred_glm, crosspred_rf),1,median)

# Evaluate ensemble predictions
(eval_ens <- mecofun::evalSDM(observation = avi_df[!is.na(avi_df$blockCV_tile),1], predictions = crosspred_ens))


##########
# 2.5.1 Ring Ouzel: Predictions
##########

# Make predictions to current climate:
bio_curr_df <- terra::as.data.frame(bio_curr, xy = TRUE)
colnames(bio_curr_df)

bio_curr_df$pred_glm <- mecofun::predictSDM(m_glm, bio_curr_df)
bio_curr_df$pred_rf <- mecofun::predictSDM(m_rf, bio_curr_df)
bio_curr_df$pred_ens <- apply(bio_curr_df[,-c(1:5)],1,median)

# Make binary predictions:
bio_curr_df$bin_glm <- ifelse(bio_curr_df$pred_glm > eval_glm$thresh, 1, 0)
bio_curr_df$bin_rf <- ifelse(bio_curr_df$pred_rf > eval_rf$thresh, 1, 0)
bio_curr_df$bin_ens <- ifelse(bio_curr_df$pred_ens > eval_ens$thresh, 1, 0)

# Make raster stack of predictions:
r_pred_curr <- rast(bio_curr_df[,-c(3:5)])
plot(r_pred_curr)


avi_df$bio_2 <- NULL
avi_df$bio_14 <- NULL

m_glm <- glm( Turdus_torquatus ~ bio_5 + I(bio_5^2), family='binomial', data=avi_df)
m_rf <- randomForest( x=avi_df[,2], y=avi_df[,1], ntree=1000, nodesize=10, importance =T)

bio_fut_df$pred_glm <- mecofun::predictSDM(m_glm, bio_fut_df)
#bio_fut_df$pred_rf <- mecofun::predictSDM(m_rf, bio_fut_df)
#bio_fut_df$pred_ens <- apply(bio_fut_df[,-c(1:5)],1,median)

bio_fut_rast <- rast(bio_fut_df)
plot(bio_fut_rast)
